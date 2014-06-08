package s2

import (
	"code.google.com/p/gos2/r1"
	"code.google.com/p/gos2/s1"
	"math"
	"sort"
)

type CellEdge struct {
	cell CellID
	edge int
}

type CellEdgeMultimap map[CellID][]int

func (m CellEdgeMultimap) Insert(cid CellID, i int) {
	if _, ok := m[cid]; !ok {
		m[cid] = []int{}
	}
	m[cid] = append(m[cid], i)
}

// Indexing structure to efficiently compute intersections.
type LoopIndex struct {
	loop           *Loop
	index_computed bool
	query_count    int
	min_level_used int
	mapping        CellEdgeMultimap
}

func NewLoopIndex(loop *Loop) LoopIndex {
	return LoopIndex{
		loop:           loop,
		index_computed: false,
		query_count:    0,
		min_level_used: MaxCellLevel,
		mapping:        CellEdgeMultimap{},
	}
}

func (idx *LoopIndex) Reset() {
	idx.min_level_used = MaxCellLevel
	idx.index_computed = false
	idx.query_count = 0
	idx.mapping = CellEdgeMultimap{}
}

func (idx *LoopIndex) IncrementQueryCount() {
	idx.query_count++
}

func (idx *LoopIndex) edge_from(i int) *Point {
	return &idx.loop.vertices[i]
}

func (idx *LoopIndex) edge_to(i int) *Point {
	return &idx.loop.vertices[i+1]
}

func (idx *LoopIndex) num_edges() int {
	return len(idx.loop.vertices)
}

func (idx *LoopIndex) PredictAdditionalCalls(n int) {
	if !idx.index_computed {
		if idx.num_edges() > 100 && (idx.query_count+n) > 30 {
			idx.ComputeIndex()
		}
	}
}

func (idx *LoopIndex) ComputeIndex() {
	for i := 0; i < idx.num_edges(); i++ {
		cover, level := idx.GetCovering(*idx.edge_from(i), *idx.edge_to(i), true)
		idx.min_level_used = int(math.Min(float64(idx.min_level_used), float64(level)))
		for _, cid := range cover {
			idx.mapping.Insert(cid, i)
		}
	}
	idx.index_computed = true
}

func LenientCrossing(a, b, c, d Point) bool {
	maxDetError := 1e-14
	acb := a.PointCross(c).Dot(b.Vector)
	bda := b.PointCross(d).Dot(a.Vector)
	if math.Abs(acb) < maxDetError || math.Abs(bda) < maxDetError {
		return true
	}
	if acb*bda < 0 {
		return false
	}
	cbd := c.PointCross(b).Dot(d.Vector)
	dac := d.PointCross(a).Dot(c.Vector)
	if math.Abs(cbd) < maxDetError || math.Abs(dac) < maxDetError {
		return true
	}
	return acb*cbd >= 0 && acb*dac >= 0
}

func (idx *LoopIndex) EdgeIntersectsCellBoundary(a, b Point, cell Cell) bool {
	vertices := make([]Point, 4)
	for i := 0; i < 4; i++ {
		vertices[i] = cell.Vertex(i)
	}
	for i := 0; i < 4; i++ {
		from := vertices[i]
		to := vertices[(i+1)%4]
		if LenientCrossing(a, b, from, to) {
			return true
		}
	}
	return false
}

// Appends to "candidate_crossings" all edge references which may cross the
// given edge. This is done by covering the edge and then finding all references
// of edges whose coverings overlap this covering. Parent cells are checked
// level by level. Child cells are checked all at once by taking advantage of
// the natural ordering of CellIDs.
func (idx *LoopIndex) FindCandidateCrossings(a, b Point, candidate_crossings *[]int) {
	cover, _ := idx.GetCovering(a, b, false)
	idx.EdgesInParentCells(cover, idx.mapping, idx.min_level_used, candidate_crossings)
	idx.EdgesInChildrenCells(a, b, &cover, idx.mapping, candidate_crossings)
}

func (idx *LoopIndex) EdgesInParentCells(
	cover []CellID,
	mapping CellEdgeMultimap,
	min_level_used int,
	candidate_crossings *[]int) {
	// Find all parent cells of covering cells.
	parent_cells := map[CellID]bool{}
	for _, cid := range cover {
		for parentLevel := cid.Level() - 1; parentLevel >= min_level_used; parentLevel-- {
			if _, ok := parent_cells[cid.Parent(parentLevel)]; !ok {
				parent_cells[cid.Parent(parentLevel)] = true
			} else {
				break
			}
		}
	}

	// Put parent cell edge references into result.
	for pi, _ := range parent_cells {
		sort.IntSlice(mapping[pi]).Sort()
		indexes := mapping[pi]
		for _, index := range indexes {
			*candidate_crossings = append(*candidate_crossings, index)
		}
	}
}

func (idx *LoopIndex) EdgesInChildrenCells(
	a, b Point,
	cover *[]CellID,
	mapping CellEdgeMultimap,
	candidate_crossings *[]int) {
	num_cells := 0

	// Put all the edge references of (covering cells + descendant cells)
	// into result. This relies on the natural ordering of CellIDs.
	keys := []CellID{}
	for k, _ := range mapping {
		keys = append(keys, k)
	}

	sort.Sort(byID(keys))

	rewind := false

	for len(*cover) > 0 {
		back := len(*cover) - 1
		cell := (*cover)[back]
		*cover = append((*cover)[:back], (*cover)[back+1:]...)
		num_cells++
		start := sort.Search(len(keys), func(i int) bool { return keys[i] >= cell.RangeMin() })
		end := sort.Search(len(keys), func(i int) bool { return keys[i] >= cell.RangeMax() })
		if start < len(keys) && end < len(keys) && keys[start] == cell.RangeMin() && keys[end] == cell.RangeMax() {
			num_edges := 0
			for cid := start; cid <= end; cid++ {
				indexes := mapping[keys[cid]]
				for _, index := range indexes {
					*candidate_crossings = append(*candidate_crossings, index)
					num_edges++
					if num_edges == 16 && !cell.IsLeaf() {
						rewind = true
						break
					}
				}
			}

			// If there were too many to insert, uninsert and recurse.
			if rewind {
				*candidate_crossings = (*candidate_crossings)[:num_edges]
				sort.IntSlice(mapping[cell]).Sort()
				for _, id := range mapping[cell] {
					*candidate_crossings = append(*candidate_crossings, id)
				}
				// Recurse on the children -- hopefully some will be empty.
				if mapping[cell][0] != mapping[cell][keys[start]] ||
					mapping[cell][len(mapping[cell])-1] != mapping[cell][keys[end]] {
					children := cell.Children()
					for _, child := range children {
						c := CellFromCellID(child)
						if idx.EdgeIntersectsCellBoundary(a, b, c) {
							*cover = append(*cover, child)
						}
					}
				}
			}
		}
	}
}

// Returns the smallest cell containing all four points, or Sentinel if they
// are not all on the same face. The points don't need to be normalized.
func containingCell4(pa, pb, pc, pd Point) CellID {
	a := cellIDFromPoint(pa)
	b := cellIDFromPoint(pb)
	c := cellIDFromPoint(pc)
	d := cellIDFromPoint(pd)
	if a.Face() != b.Face() || a.Face() != c.Face() || a.Face() != d.Face() {
		return Sentinel()
	}

	for a != b || a != c || a != d {
		a = a.immediateParent()
		b = b.immediateParent()
		c = c.immediateParent()
		d = d.immediateParent()
	}
	return a
}

// Returns the smallest cell containing both points, or Sentinel if they
// are not all on the same face. The points don't need to be normalized.
func containingCell2(pa, pb Point) CellID {
	a := cellIDFromPoint(pa)
	b := cellIDFromPoint(pb)
	if a.Face() != b.Face() {
		return Sentinel()
	}
	for a != b {
		a = a.immediateParent()
		b = b.immediateParent()
	}
	return a
}

func (idx *LoopIndex) GetCovering(a, b Point, thicken_edge bool) ([]CellID, int) {
	covering := []CellID{}
	// Thicken the edge in all directions by roughly 1% of the edge length
	// when thicken_edge is true.
	thickening := 0.01

	// Selects the ideal S2 level at which to cover the edge, this will be
	// the level whose S2 cells have a width roughly commensurate to the
	// length of the edge. We multiply the edge length by 2*thickening to
	// guarantee the thickening is honored when doing the covering-by-cap
	// trick.
	edge_length := float64(a.Distance(b))
	ideal_level := MinWidth.MaxLevel(edge_length * (1 + 2*thickening))
	var containing_cell CellID
	if !thicken_edge {
		containing_cell = containingCell2(a, b)
	} else {
		if ideal_level == MaxCellLevel {
			// If the edge is tiny, instabilities are more likely,
			// so we want to limit the number of operations.
			// We pretend we are in a cell much larger so as to
			// trigger the 'needs covering' case, so we won't try to
			// thicken the edge.
			containing_cell = CellID(0xFFF0).Parent(3)
		} else {
			pq := b.Sub(a.Vector).Mul(thickening)
			ortho := pq.Cross(a.Vector).Normalize().Mul(edge_length * thickening)
			p := a.Sub(pq)
			q := b.Add(pq)
			containing_cell = containingCell4(Point{p.Sub(ortho)}, Point{p.Add(ortho)},
				Point{q.Sub(ortho)}, Point{q.Add(ortho)})
		}
	}

	// Best case: edge is fully contained in a cell that's not too big.
	if containing_cell != Sentinel() && containing_cell.Level() >= ideal_level-2 {
		covering = append(covering, containing_cell)
		return covering, containing_cell.Level()
	}

	if ideal_level == 0 {
		// Edge is very long, maybe even longer than a face width, so
		// the trick below doesn't work. For now, we will add the whole
		// S2 sphere.
		for cid := CellIDBegin(0); cid != CellIDEnd(0); cid = cid.next() {
			covering = append(covering, cid)
		}
		return covering, 0
	}

	// Cover the edge by a cap centered on the edge midpoint, then cover
	// the cap by four big-enough cells around the cell vertex closest to
	// the cap center.
	middle := Point{a.Add(b.Vector).Mul(0.5).Normalize()}
	actual_level := int(math.Min(float64(ideal_level), MaxCellLevel-1))
	cellIDFromPoint(middle).AppendVertexNeighbors(actual_level, &covering)
	return covering, actual_level
}

// An iterator on loops/edges that may cross a query edge (a,b).
type LoopIndexIterator struct {
	current                     int
	num_edges                   int
	is_brute_force              bool
	index                       *LoopIndex
	candidates                  []int
	current_index_in_candidates int
}

func NewLoopIndexIterator(idx *LoopIndex) *LoopIndexIterator {
	return &LoopIndexIterator{current: 0, is_brute_force: false, index: idx}
}

func (it *LoopIndexIterator) GetCandidates(a, b Point) {
	it.index.PredictAdditionalCalls(1)
	it.is_brute_force = !it.index.index_computed
	if it.is_brute_force {
		it.index.IncrementQueryCount()
		it.current = 0
		it.num_edges = it.index.num_edges()
	} else {
		it.candidates = []int{}
		it.index.FindCandidateCrossings(a, b, &it.candidates)
		it.current_index_in_candidates = 0
		if len(it.candidates) != 0 {
			it.current = it.candidates[0]
		}
	}
}

// Index of the current loop in the iteration.
func (it *LoopIndexIterator) Index() int {
	return it.current
}

// Iterate to the next available candidate.
func (it *LoopIndexIterator) Next() {
	if it.is_brute_force {
		it.current++
	} else {
		it.current_index_in_candidates++
		if it.current_index_in_candidates < len(it.candidates) {
			it.current = it.candidates[it.current_index_in_candidates]
		}
	}
}

// True if there are no more candidates.
func (it *LoopIndexIterator) Done() bool {
	if it.is_brute_force {
		return it.current >= it.num_edges
	} else {
		return it.current_index_in_candidates >= len(it.candidates)
	}
}

// A Loop represents a simple spherical polygon. It consists of a single
// chain of vertices where the first vertex is implicitly connected to the
// last. All loops are defined to have a CCW orientation, i.e. the interior
// of the polygon is on the left side of the edges. This implies that a
// clockwise loop enclosing a small area is interpreted to be a CCW loop
// enclosing a very large area.
//
// Loops are not allowed to have any duplicate vertices (whether adjacent or
// not), and non-adjacent edges are not allowed to intersect. Loops must have
// at least 3 vertices. Although these restrictions are not enforced in
// optimized code, you may get unexpected results if they are violated.
//
// Point containment is defined such that if the sphere is subdivided into
// faces (loops), every point is contained by exactly one face. This implies
// that loops do not necessarily contain all (or any) of their vertices.
type Loop struct {
	vertices      []Point
	bound         Rect
	origin_inside bool
	depth         int
	index         LoopIndex

	// Map for speeding up FindVertex: We will compute a map from vertex to
	// index in the vertex array as soon as there has been enough calls.
	num_find_vertex_calls int
	vertex_to_index       map[Point]int
}

func NewLoopFromPath(vertices []Point) *Loop {
	loop := &Loop{
		vertices: make([]Point, len(vertices)),
		bound:    FullRect(),
		depth:    0,
		num_find_vertex_calls: 0,
	}
	loop.ResetMutableFields()
	copy(loop.vertices, vertices)
	loop.index = NewLoopIndex(loop)
	loop.InitOrigin()
	loop.InitBound()
	return loop
}

func NewLoopFromCell(cell Cell) *Loop {
	loop := &Loop{
		vertices: make([]Point, 4),
		bound:    cell.RectBound(),
		depth:    0,
		num_find_vertex_calls: 0,
	}
	for i := 0; i < 4; i++ {
		loop.vertices[i] = cell.Vertex(i)
	}
	loop.index = NewLoopIndex(loop)
	loop.InitOrigin()
	loop.InitBound()
	return loop
}

func (l Loop) IsHole() bool { return (l.depth & 1) != 0 }
func (l Loop) Sign() int {
	if l.IsHole() {
		return -1
	}
	return 1
}

func (l *Loop) Clone() *Loop {
	loop := &Loop{
		vertices:      make([]Point, len(l.vertices)),
		bound:         l.bound,
		origin_inside: l.origin_inside,
		depth:         l.depth,
		num_find_vertex_calls: 0,
	}
	copy(loop.vertices, l.vertices)
	loop.index = NewLoopIndex(loop)
	return loop
}

func (l Loop) Bound() Rect {
	return l.bound
}

func (l Loop) CapBound() Cap {
	return l.bound.CapBound()
}

func (l *Loop) FindVertex(p Point) int {
	l.num_find_vertex_calls++
	if len(l.vertices) < 10 || l.num_find_vertex_calls < 20 {
		// Exhaustive search
		for i := 1; i <= len(l.vertices); i++ {
			if *l.vertex(i) == p {
				return i
			}
		}
		return -1
	}
	if len(l.vertex_to_index) == 0 { // We haven't computed it yet.
		for i := len(l.vertices); i > 0; i-- {
			l.vertex_to_index[*l.vertex(i)] = i
		}
	}

	if idx, ok := l.vertex_to_index[p]; ok {
		return idx
	}
	return -1
}

func (l *Loop) Invert() {
	l.ResetMutableFields()
	// Reverse vertices.
	for i, j := 0, len(l.vertices)-1; i < j; i, j = i+1, j-1 {
		l.vertices[i], l.vertices[j] = l.vertices[j], l.vertices[i]
	}
	l.origin_inside = !l.origin_inside
	if l.bound.Lo().Lat > -(math.Pi/2) && l.bound.Hi().Lat < math.Pi/2 {
		// The complement of this loop contains both poles.
		l.bound = FullRect()
	} else {
		l.InitBound()
	}
}

func (l *Loop) ResetMutableFields() {
	l.index.Reset()
	l.num_find_vertex_calls = 0
	l.vertex_to_index = map[Point]int{}
}

func (l *Loop) InitOrigin() {
	// The bounding box does not need to be correct before calling this
	// function, but it must at least contain vertex(1) since we need to
	// do a Contains() test on this point below.
	//l.bound.Contains(LatLngFromPoint(*l.vertex(1)))

	// To ensure that every point is contained in exactly on face of a
	// subdivision of the sphere, all containment tests are done by counting
	// the edge crossings starting at a fixed point on there sphere
	// (Origin()). We need to know whether this point is inside or outside
	// of the loop. We do this by first guessing that it is outside, and
	// then seeing whether we get the correct containment result for
	// vertex 1. If the result is incorrect, the origin must be inside the
	// loop.
	//
	// A loop with consecutive vertices A,B,C contains vertex B iff the
	// fixed vector R = Ortho(B) is on the left side of the wedge ABC.
	// The test below is written so that B is inside if C=R but not if A=R.
	l.origin_inside = false // Initialize before calling Contains()
	v1_inside := OrderedCCW(Point{l.vertex(1).Ortho()}, *l.vertex(0), *l.vertex(2), *l.vertex(1))
	if v1_inside != l.Contains(*l.vertex(1)) {
		l.origin_inside = true
	}
}

func (l *Loop) InitBound() {
	// The bounding rectangle of a loop is not necessarily the same as the
	// bounding rectangle of its vertices. First, the loop may wrap entirely
	// around the sphere (e.g. a loop that defines two revolutions of a
	// candy-cane stripe). Second, the loop may include one or both poles.
	// Note that a small clockwise loop near the equator contains both
	// poles.

	bounder := NewRectBounder()
	for i := 0; i <= len(l.vertices); i++ {
		bounder.AddPoint(l.vertex(i))
	}

	b := bounder.Bound()
	// Note that we need to initialize l.bound with a temporary value since
	// Contains() does a bounding rectangle check before doing anything
	// else.
	l.bound = FullRect()
	if l.Contains(PointFromCoords(0, 0, 1)) {
		b = Rect{
			Lat: r1.Interval{b.Lat.Lo, math.Pi / 2},
			Lng: s1.FullInterval(),
		}
	}
	// If a loop contains the south pole, then either it wraps entirely
	// around the sphere (full longitude range), or it also contains the
	// north pole in which case b.Lng.IsFull() due to the test above.
	// Either way, we only need to do the south pole containment test if
	// b.Lng.IsFull()
	if b.Lng.IsFull() && l.Contains(PointFromCoords(0, 0, -1)) {
		b.Lat.Lo = -math.Pi / 2
	}
	l.bound = b
}

func (l Loop) vertex(i int) *Point {
	j := i - len(l.vertices)
	if j >= 0 {
		return &l.vertices[j]
	}
	return &l.vertices[i]
}

func (l Loop) IsNormalized() bool {
	// Optimization: if the longitude span is less than 180 degrees, then
	// the loop covers less than half the sphere and is therefore
	// normalized.
	if l.bound.Lng.Length() < math.Pi {
		return true
	}
	// We allow some error so that hemispheres are always considered
	// normalized.
	return l.TurningAngle() >= -1e-14
}

func (l *Loop) Normalize() {
	if !l.IsNormalized() {
		l.Invert()
	}
}

// Return (first, dir) such that first..first+n*dir are valid indices.
func (l Loop) CanonicalFirstVertex() (first, dir int) {
	first = 0
	n := len(l.vertices)
	for i := 1; i < n; i++ {
		if l.vertex(i).LessThan(l.vertex(first).Vector) {
			first = i
		}
	}
	if l.vertex(first + 1).LessThan(l.vertex(first + n - 1).Vector) {
		dir = 1
		// 0 <= first <= n-1, so (first+n*dir) <= 2*n-1.
	} else {
		dir = -1
		// n <= first <= 2*n-1, so (first+n*dir) >= 0.
	}
	return
}

func (l Loop) TurningAngle() float64 {
	// Don't crash even if the loop is not well-defined.
	if len(l.vertices) < 3 {
		return 0
	}
	// To ensure that we get the same result when the loop vertex order is
	// rotated, and that we get the same result with the opposite sign when
	// the vertices are reversed, we need to be careful to add up the
	// individual turn angles in a consistent order. In general, adding up
	// a set of numbers in a different order can change the sum due to
	// rounding errors.
	n := len(l.vertices)
	i, dir := l.CanonicalFirstVertex()
	angle := TurnAngle(*l.vertex((i + n - dir) % n), *l.vertex(i), *l.vertex((i + dir) % n))
	for n = n - 1; n > 0; n-- {
		i += dir
		angle += TurnAngle(*l.vertex(i - dir), *l.vertex(i), *l.vertex(i + dir))
	}
	return float64(dir) * angle
}

func (l Loop) ContainsCell(cell Cell) bool {
	if !l.bound.ContainsPoint(cell.Center()) {
		return false
	}
	cell_loop := NewLoopFromCell(cell)
	return l.ContainsLoop(cell_loop)
}

func (a Loop) ContainsLoop(b *Loop) bool {
	// For this loop A to contain the given loop B, all of the following
	// must be true:
	//
	//  (1) There are no edge crossings between A and B except at vertices.
	//
	//  (2) At every vertex that is shared between A and B, the local edge
	//      ordering implies that A contains B.
	//
	//  (3) If there are no shared vertices, then A must contain a vertex
	//      of B and B must not contain a vertex of A. (An arbitrary vertex
	//      may be chosen in each case.)
	//
	// The second part of (3) is necessary to detect the case of two loops
	// whose union is the entire sphere, i.e. two loops that contain each
	// other's boundaries but not each other's interiors.
	if !a.bound.ContainsRect(b.bound) {
		return false
	}

	// Unless there are shared vertices, we need to check whether A
	// contains a vertex of B. Since shared vertices are rare, it is more
	// efficient to do this test up front as a quick rejection test.
	if !a.Contains(*b.vertex(0)) && a.FindVertex(*b.vertex(0)) < 0 {
		return false
	}

	// Now check whether there are any edge crossings, and also check the
	// loop relationship at any shared vertices.
	var wedge ContainsWedgeProcessor
	if a.AreBoundariesCrossing(b, &wedge) || wedge.doesnt_contain {
		return false
	}

	// At this point we know that the boundaries of A and B do not
	// intersect, and that A contains a vertex of B. However, we still
	// need to check for the case mentioned above, where (A union B) is
	// the entire sphere. Normally this check is very cheap due to the
	// bounding box precondition.
	if a.bound.Union(b.bound).IsFull() {
		if b.Contains(*a.vertex(0)) && b.FindVertex(*a.vertex(0)) < 0 {
			return false
		}
	}
	return true
}

func (a Loop) ContainsNested(b *Loop) bool {
	if !a.bound.ContainsRect(b.bound) {
		return false
	}
	// We are given that A and B do not share any edges, and that either
	// one loop contains the other or they do not intersect.
	m := a.FindVertex(*b.vertex(1))
	if m < 0 {
		// Since b.vertex(1) is not shared, we can check whether A
		// contains it.
		return a.Contains(*b.vertex(1))
	}
	// Check whether the edge order around b.vertex(1) is compatible with
	// A containing B.
	return WedgeContains(*a.vertex(m - 1), *a.vertex(m), *a.vertex(m + 1),
		*b.vertex(0), *b.vertex(2))
}

func (l Loop) Contains(p Point) bool {
	if !l.bound.Contains(LatLngFromPoint(p)) {
		return false
	}
	inside := l.origin_inside
	origin := OriginPoint()
	crosser := NewEdgeCrosser(&origin, &p, l.vertex(0))
	vlen := len(l.vertices)
	if vlen < 2000 {
		for i := 1; i <= vlen; i++ {
			inside = inside != crosser.EdgeOrVertexCrossing(l.vertex(i))
		}
		return inside
	}

	it := NewLoopIndexIterator(&l.index)
	prev_index := -2
	for it.GetCandidates(origin, p); !it.Done(); it.Next() {
		ai := it.Index()
		if prev_index != ai-1 {
			crosser.RestartAt(l.vertex(ai))
		}
		prev_index = ai
		inside = inside != crosser.EdgeOrVertexCrossing(l.vertex(ai+1))
	}
	return inside
}

func (l *Loop) MayIntersect(cell Cell) bool {
	if !l.bound.Intersects(cell.RectBound()) {
		return false
	}
	return NewLoopFromCell(cell).Intersects(l)
}

type WedgeProcessor interface {
	ProcessWedge(a0, ab1, a2, b0, b2 Point) bool
}

// WedgeProcessor to be used to check if loop A intersects loop B.
// Intersects() then returns true when A and B have at least one pair
// of associated wedges that intersect.
type IntersectsWedgeProcessor struct {
	intersects bool
}

func (p *IntersectsWedgeProcessor) ProcessWedge(a0, ab1, a2, b0, b2 Point) bool {
	p.intersects = WedgeIntersects(a0, ab1, a2, b0, b2)
	return p.intersects
}

// WedgeProcessor to be used to check if loop A contains loop B.
// DoesntContain() then returns true if there is a wedge of B not contained
// in the associated wedge of A (and hence loop B is not contained in loop A).
type ContainsWedgeProcessor struct {
	doesnt_contain bool
}

func (p *ContainsWedgeProcessor) ProcessWedge(a0, ab1, a2, b0, b2 Point) bool {
	p.doesnt_contain = !WedgeContains(a0, ab1, a2, b0, b2)
	return p.doesnt_contain
}

// WedgeProcessor to be used to check if the interior of loop A contains the
// interior of loop B, or their boundaries cross each other (therefore they
// have a proper intersection). CrossesOrMayContain() then returns -1 if A
// crossed B, 0 if it is not possible for A to contain B, and 1 otherwise.
type ContainsOrCrossesProcessor struct {
	// True if any crossing on the boundary is discovered.
	has_boundary_crossing bool
	// True if A (B) has a strictly superwedge on a pair of wedges that
	// share a common center point.
	a_has_strictly_super_wedge bool
	b_has_strictly_super_wedge bool
	// True if there is a pair of disjoint wedges with a common center
	// point.
	has_disjoint_wedge bool
}

func (p *ContainsOrCrossesProcessor) CrossesOrMayContain() int {
	if p.has_boundary_crossing {
		return -1
	}
	if p.has_disjoint_wedge || p.b_has_strictly_super_wedge {
		return 0
	}
	return 1
}

func (p *ContainsOrCrossesProcessor) ProcessWedge(a0, ab1, a2, b0, b2 Point) bool {
	wedgeRelation := GetWedgeRelation(a0, ab1, a2, b0, b2)
	if wedgeRelation == WEDGE_PROPERLY_OVERLAPS {
		p.has_boundary_crossing = true
		return true
	}
	p.a_has_strictly_super_wedge = wedgeRelation == WEDGE_PROPERLY_CONTAINS
	p.b_has_strictly_super_wedge = wedgeRelation == WEDGE_IS_PROPERLY_CONTAINED
	if p.a_has_strictly_super_wedge && p.b_has_strictly_super_wedge {
		p.has_boundary_crossing = true
		return true
	}
	p.has_disjoint_wedge = wedgeRelation == WEDGE_IS_DISJOINT
	return false
}

// This method checks all edges of this loop (A) for intersection against
// all edges of B. If there is any shared vertex, the wedges centered at this
// vertex are sent to wedge_processor.
//
// Returns true only when the edges intersect in the sense of RobustCrossing,
// returns false immediately when the wedge_processor returns true: this means
// the wedge processor knows the value of the property that the caller wants
// to compute, and no further inspection is needed. For instance, if the
// property is "loops intersect", then a wedge intersection is all it takes
// to return true.
//
// See Intersects().
func (a *Loop) AreBoundariesCrossing(b *Loop, wedge_processor WedgeProcessor) bool {
	a.index.PredictAdditionalCalls(len(b.vertices))
	it := NewLoopIndexIterator(&a.index)
	for j := 0; j < len(b.vertices); j++ {
		crosser := NewEdgeCrosser(b.vertex(j), b.vertex(j+1), b.vertex(0))
		prev_index := -2
		for it.GetCandidates(*b.vertex(j), *b.vertex(j + 1)); !it.Done(); it.Next() {
			ai := it.Index()
			if prev_index != ai-1 {
				crosser.RestartAt(a.vertex(ai))
			}
			prev_index = ai
			crossing := crosser.RobustCrossing(a.vertex(ai + 1))
			if crossing < 0 {
				continue
			}
			if crossing > 0 {
				return true
			}
			// We only need to check each shared vertex once, so we
			// only consider the case where
			// a.vertex(ai+1) == b.vertex(j+1).
			if *a.vertex(ai + 1) == *b.vertex(j + 1) &&
				wedge_processor.ProcessWedge(*a.vertex(ai), *a.vertex(ai + 1), *a.vertex(ai + 2),
					*b.vertex(j), *b.vertex(j + 1)) {
				return false
			}
		}
	}
	return false
}

func (a *Loop) Intersects(b *Loop) bool {
	// a.Intersects(b) if and only if !a.Complement().Contains(b).
	// This code is similar to Contains(), but is optimized for the case
	// where both loops enclose less than half of the sphere.

	// The largest of the two loops should be edgeindex'd.
	if len(b.vertices) > len(a.vertices) {
		return b.Intersects(a)
	}

	if !a.bound.Intersects(b.bound) {
		return false
	}

	// Unless there are shared vertices, we need to check whether A
	// contains a vertex of B. Since shared vertices are rare, it is more
	// efficient to do this test up front as a quick acceptance test.
	if a.Contains(*b.vertex(0)) && a.FindVertex(*b.vertex(0)) < 0 {
		return true
	}

	// Now check whether there are any edge crossings, and also check the
	// loop relationship at any shared vertices.
	var wedge_processor IntersectsWedgeProcessor
	if a.AreBoundariesCrossing(b, &wedge_processor) || wedge_processor.intersects {
		return true
	}

	// We know that A does not contain a vertex of B, and that there are
	// no edge crossings. Therefore the only way that A can intersect B is
	// if B entirely contains A. We can check this by testing whether B
	// contains an arbitrary non-shared vertex of A. Note that this check
	// is usually cheap because of the bounding box precondition.
	if b.bound.ContainsRect(a.bound) {
		if b.Contains(*a.vertex(0)) && b.FindVertex(*a.vertex(0)) < 0 {
			return true
		}
	}
	return false
}

func (a *Loop) ContainsOrCrosses(b *Loop) int {
	// There can be containment or crossing only if the bounds intersect.
	if a.bound.Intersects(b.bound) {
		return 0
	}
	// Now check whether there are any edge crossings, and also check
	// the loop relationship at any shared vertices. Note that unlike
	// Contains() or Intersects(), we can't do a point containment test
	// as a shortcut because we need to detect whether there are any
	// edge crossings.
	var wedgeProcessor ContainsOrCrossesProcessor
	if a.AreBoundariesCrossing(b, &wedgeProcessor) {
		return -1
	}
	res := wedgeProcessor.CrossesOrMayContain()
	if res <= 0 {
		return res
	}

	// At this point we know that the boundaries do not intersect, and we
	// are given that (A union B) is a proper subset of the sphere.
	// Furthermore, either A contains B or there are no shared vertices
	// (due to the check above). So now we just need to distinguish the
	// case where A contains B from the case where B contains A or the
	// two loops are disjoint.
	if !a.bound.ContainsRect(b.bound) {
		return 0
	}
	if !a.Contains(*b.vertex(0)) && a.FindVertex(*b.vertex(0)) < 0 {
		return 0
	}
	return 1
}

func (a *Loop) BoundaryApproxEquals(b *Loop, maxError float64) bool {
	numVerts := len(a.vertices)
	if numVerts != len(b.vertices) {
		return false
	}
	for offset := 0; offset < numVerts; offset++ {
		if a.vertex(offset).ApproxEqualWithin(*b.vertex(0), maxError) {
			success := true
			for i := 0; i < numVerts; i++ {
				if !a.vertex(i+offset).ApproxEqualWithin(*b.vertex(i), maxError) {
					success = false
					break
				}
			}
			if success {
				return true
			}
		}
	}
	return false
}

func (a *Loop) BoundaryNear(b *Loop, maxError float64) bool {
	for offset := 0; offset < len(a.vertices); offset++ {
		if MatchBoundaries(a, b, offset, maxError) {
			return true
		}
	}
	return false
}

func MatchBoundaries(a, b *Loop, offset int, maxError float64) bool {
	pending := []IntPair{}
	done := map[IntPair]bool{}
	pending = append(pending, IntPair{0, 0})
	alen := len(a.vertices)
	blen := len(b.vertices)
	for len(pending) > 0 {
		back := len(pending) - 1
		i := pending[back].first
		j := pending[back].second
		pending = pending[:back]
		if i == alen && j == blen {
			return true
		}
		done[IntPair{i, j}] = true

		io := i + offset
		if io >= alen {
			io -= alen
		}

		if i < alen {
			if _, ok := done[IntPair{i + 1, j}]; !ok {
				if a.vertex(io+1).DistanceToEdge(
					*b.vertex(j),
					*b.vertex(j + 1)).Radians() <= maxError {
					pending = append(pending, IntPair{i + 1, j})
				}
			}
		}
		if j < blen {
			if _, ok := done[IntPair{i, j + 1}]; !ok {
				if b.vertex(j+1).DistanceToEdge(
					*a.vertex(io),
					*a.vertex(io + 1)).Radians() <= maxError {
					pending = append(pending, IntPair{i, j + 1})
				}
			}
		}
	}
	return false
}
