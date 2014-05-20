package s2

import (
	"code.google.com/p/gos2/s1"
	"math"
	"sort"
)

type Edge struct {
	v0, v1 Point
}

type ByEdgeFirst []Edge

func (a ByEdgeFirst) Len() int { return len(a) }
func (a ByEdgeFirst) Swap(i, j int) {
	a[i].v0, a[j].v0 = a[j].v0, a[i].v0
	a[i].v1, a[j].v1 = a[j].v1, a[i].v1
}
func (a ByEdgeFirst) Less(i, j int) bool {
	return a[i].v0.LessThan(a[j].v0.Vector)
}

type PolygonBuilderOptions struct {
	xor_edges            bool
	undirected_edges     bool
	validate             bool
	vertex_merge_radius  s1.Angle
	edge_splice_fraction float64
}

func DIRECTED_XOR() PolygonBuilderOptions {
	return PolygonBuilderOptions{
		xor_edges:            true,
		undirected_edges:     false,
		validate:             false,
		vertex_merge_radius:  0,
		edge_splice_fraction: 0.866,
	}
}

type PointIndex struct {
	map_          map[CellID][]Point
	vertex_radius float64
	edge_fraction float64
	level         int
}

func NewPointIndex(vertex_radius, edge_fraction float64) PointIndex {
	max_level := MinWidth.MaxLevel(2 * vertex_radius)
	return PointIndex{
		map_:          map[CellID][]Point{Sentinel(): []Point{Point{}}},
		vertex_radius: vertex_radius,
		edge_fraction: edge_fraction,
		level:         int(math.Min(float64(max_level), MaxCellLevel-1)),
	}
}

func (idx PointIndex) FindNearbyPoint(v0, v1 Point, out *Point) bool {
	// Return a point whose distance from the edge (v0, v1) is less than
	// vertex_radius, and which is not equal to v0 or v1. The current
	// implementation returns the closest such point.
	//
	// Strategy: we compute a very cheap covering by approximating the edge
	// as two spherical caps, one around each endpoint, and then computing a
	// 4-cell covering using each one. We could improve the quality of the
	// covering by using some intermediate points along the edge as well.
	length := float64(v0.Distance(v1))
	//	normal := v0.PointCross(v1)
	level := int(math.Min(float64(idx.level), float64(MinWidth.MaxLevel(length))))
	ids := []CellID{}
	cellIDFromPoint(v0).AppendVertexNeighbors(level, &ids)
	cellIDFromPoint(v1).AppendVertexNeighbors(level, &ids)

	// Sort the cell ids so that we can skip duplicates in the loop below.
	sort.Sort(byID(ids))

	// TODO: map_ won't work as is. You need to sort the mapped values
	// before using them.

	//	best_dist := 2 * idx.vertex_radius
	for i := len(ids) - 1; i >= 0; i-- {
		// Skip duplicates
		if i > 0 && ids[i-1] == ids[i] {
			continue
		}
		//		max_id := ids[i].RangeMax()
		for _, v := range idx.map_[ids[i]] {
			if v == v0 || v == v1 {
				continue
			}
		}
	}
	return true
}

func (idx *PointIndex) Insert(p Point) {
	ids := []CellID{}
	cellIDFromPoint(p).AppendVertexNeighbors(idx.level, &ids)
	for i := len(ids) - 1; i >= 0; i-- {
		idx.map_[ids[i]] = append(idx.map_[ids[i]], p)
	}
}

func (idx *PointIndex) Erase(p Point) {
	ids := []CellID{}
	cellIDFromPoint(p).AppendVertexNeighbors(idx.level, &ids)
	for i := len(ids) - 1; i >= 0; i-- {
		for k, v := range idx.map_[ids[i]] {
			if v == p {
				idx.map_[ids[i]] = append(idx.map_[ids[i]][:k], idx.map_[ids[i]][k+1:]...)
				break
			}
		}
	}
}

func (idx PointIndex) QueryCap(axis Point) []Point {
	out := []Point{}
	id := cellIDFromPoint(axis).Parent(idx.level)
	for _, v := range idx.map_[id] {
		if float64(axis.Distance(v)) < idx.vertex_radius {
			out = append(out, v)
		}
	}
	return out
}

type VertexSet []Point
type EdgeSet map[Point]*VertexSet

type ByVertexSet VertexSet

func (a ByVertexSet) Len() int           { return len(a) }
func (a ByVertexSet) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a ByVertexSet) Less(i, j int) bool { return a[i].LessThan(a[j].Vector) }

func (vs VertexSet) Empty() bool {
	return len(vs) == 0
}

func (vs *VertexSet) Insert(p Point) {
	*vs = append(*vs, p)
	sort.Sort(ByVertexSet(*vs))
}

func (vs *VertexSet) Erase(idx int) {
	if idx >= 0 && idx < len(*vs) {
		*vs = append((*vs)[:idx], (*vs)[idx+1:]...)
	}
}

// Returns an index to the found point, or -1
func (vs VertexSet) Find(p Point) int {
	idx := sort.Search(len(vs), func(i int) bool {
		return vs[i].GTE(p.Vector)
	})
	if idx < len(vs) && vs[idx] == p {
		return idx
	}
	return -1
}

type PolygonBuilder struct {
	options             PolygonBuilderOptions
	edges               EdgeSet
	starting_vertices   []Point
	vertex_merge_radius s1.Angle
}

func NewPolygonBuilder(options PolygonBuilderOptions) PolygonBuilder {
	return PolygonBuilder{
		options:           options,
		edges:             EdgeSet{},
		starting_vertices: []Point{},
	}
}

func (b PolygonBuilder) HasEdge(v0, v1 Point) bool {
	candidates, ok := b.edges[v0]
	if ok {
		if candidates.Find(v1) != -1 {
			return true
		}
	}
	return false
}

func (b *PolygonBuilder) AddEdge(v0, v1 Point) bool {
	if v0 == v1 {
		return false
	}
	if b.options.xor_edges && b.HasEdge(v1, v0) {
		b.EraseEdge(v1, v0)
		return false
	}

	if _, ok := b.edges[v0]; !ok {
		b.edges[v0] = &VertexSet{}
		b.starting_vertices = append(b.starting_vertices, v0)
	}
	b.edges[v0].Insert(v1)

	if b.options.undirected_edges {
		if _, ok := b.edges[v1]; !ok {
			b.starting_vertices = append(b.starting_vertices, v1)
		}
		b.edges[v1] = &VertexSet{}
		b.edges[v1].Insert(v0)
	}
	return true
}

func (b *PolygonBuilder) EraseEdge(v0, v1 Point) {
	// Note that there may be more than one copy of an edge if we are not
	// XORing them, so a VertexSet is a multiset.
	vset := b.edges[v0]
	vset.Erase(vset.Find(v1))
	if vset.Empty() {
		delete(b.edges, v0)
	}
	if b.options.undirected_edges {
		vset = b.edges[v1]
		vset.Erase(vset.Find(v0))
		if vset.Empty() {
			delete(b.edges, v1)
		}
	}
}

type MergeMap map[Point]Point

// The overall strategy is to start from each vertex and grow a maximal
// cluster of mergeable vertices. In graph theoretic terms, we find the
// connected components of the undirected graph whose edges connect pairs of
// vertices that are separated by at most vertex_merge_radius.
//
// We then choose a single representative vertex for each cluster, and
// update "merge_map" appropriately. We choose an arbitrary existing vertex
// rather than computing the centroid of all the vertices to avoid creating new
// vertex pairs that need to be merged. (We guarantee that all vertex pairs
// are separated by at least the merge radius in the output.)
func (b *PolygonBuilder) BuildMergeMap(index PointIndex) MergeMap {
	// First, we build the set of all the distinct vertices in the input.
	// We need to include the source and destination of every edge.
	vertices := map[Point]bool{}
	merge_map := MergeMap{}

	for v0, v := range b.edges {
		vertices[v0] = true
		for _, v1 := range *v {
			vertices[v1] = true
		}
	}

	// build a spatial index containing all the distinct vertices
	for p := range vertices {
		index.Insert(p)
	}

	// next we loop through all the vertices and attempt to grow a maximal
	// mergeable group starting from each vertex
	frontier := []Point{}
	for p := range vertices {
		// skip any vertices that have already been merged with
		// another vertex
		if _, ok := merge_map[p]; ok {
			continue
		}
		// grow a maximal mergeable component starting from "p", the
		// canonical representation of the mergeable group
		frontier = append(frontier, p)
		for len(frontier) > 0 {
			i := len(frontier) - 1
			mergeable := index.QueryCap(frontier[i])
			frontier = append(frontier[:i], frontier[i+1:]...)
			for j := len(mergeable) - 1; j >= 0; j-- {
				v1 := mergeable[j]
				if v1 != p {
					// Erase from the index any vertices that will be merged.
					// This ensures that we won't try to merge the same vertex twice.
					index.Erase(v1)
					frontier = append(frontier, v1)
					merge_map[v1] = p
				}
			}
		}
	}
	return merge_map
}

func (b *PolygonBuilder) MoveVertices(merge_map MergeMap) {
	if len(merge_map) == 0 {
		return
	}
	edges := []Edge{}

	for v0, vset := range b.edges {
		for _, v1 := range *vset {
			_, ok0 := merge_map[v0]
			_, ok1 := merge_map[v1]
			if ok0 || ok1 {
				if !b.options.undirected_edges || v0.LessThan(v1.Vector) {
					edges = append(edges, Edge{v0, v1})
				}
			}
		}
	}

	// Now erase all the old edges and add all the new edges. This will
	// automatically take care of any XORing that needs to be done, because
	// EraseEdge also erases the sibling of undirected edges.
	for _, e := range edges {
		v0 := e.v0
		v1 := e.v1
		b.EraseEdge(v0, v1)
		if new0, ok := merge_map[v0]; ok {
			v0 = new0
		}
		if new1, ok := merge_map[v1]; ok {
			v1 = new1
		}
		b.AddEdge(v0, v1)
	}
}

func (b *PolygonBuilder) SpliceEdges(index PointIndex) {
	// We keep a stack of unprocessed edges. Initially all edges are
	// pushed onto the stack.
	edges := []Edge{}

	for v0, vset := range b.edges {
		for _, v1 := range *vset {
			if !b.options.undirected_edges || v0.LessThan(v1.Vector) {
				edges = append(edges, Edge{v0, v1})
			}
		}
	}

	// For each edge, we check whether there are any nearby vertices that
	// should be spliced into it. If there are, we choose one such vertex,
	// split the edge into two pieces, and iterate on each piece.
	for len(edges) > 0 {
		i := len(edges) - 1
		e := edges[i]
		edges = append(edges[:i], edges[i+1:]...)

		// If we are XORing, edges may be erased before we get to them.
		if b.options.xor_edges && !b.HasEdge(e.v0, e.v1) {
			continue
		}

		var vmid Point
		if !index.FindNearbyPoint(e.v0, e.v1, &vmid) {
			continue
		}

		b.EraseEdge(e.v0, e.v1)
		if b.AddEdge(e.v0, vmid) {
			edges = append(edges, Edge{e.v0, vmid})
		}
		if b.AddEdge(vmid, e.v1) {
			edges = append(edges, Edge{vmid, e.v1})
		}
	}
}

func (b *PolygonBuilder) AssembleLoops(loops *[]*Loop, unused_edges *[]Edge) bool {
	if b.options.vertex_merge_radius.Radians() > 0 {
		index := NewPointIndex(
			b.options.vertex_merge_radius.Radians(),
			b.options.edge_splice_fraction)
		mergemap := b.BuildMergeMap(index)
		b.MoveVertices(mergemap)
		if b.options.edge_splice_fraction > 0 {
			b.SpliceEdges(index)
		}
	}

	var dummy_unused_edges []Edge
	if unused_edges == nil {
		dummy_unused_edges = []Edge{}
		unused_edges = &dummy_unused_edges
	}

	// We repeatedly choose an edge and attempt to assemble a loop
	// starting from that edge. (This is always possible unless the
	// input includes extra edges that are not part of any loop.) To
	// maintain a consistent scanning order over b.edges between
	// different machine architectures (e.g. 'clovertown' vs 'opteron'),
	// we follow the order they were added to the builder.
	var loop *Loop

	for i := 0; i < len(b.starting_vertices); {
		v0 := b.starting_vertices[i]
		candidates, ok := b.edges[v0]
		if !ok {
			i++
			continue
		}

		v1 := (*candidates)[0]
		loop = b.AssembleLoop(v0, v1, unused_edges)
		if loop != nil {
			*loops = append(*loops, loop)
			b.EraseLoop(loop.vertices, 0, len(loop.vertices))
		}
	}
	return len(*unused_edges) == 0
}

func (b *PolygonBuilder) EraseLoop(vertices []Point, start, end int) {
	i := end - 1
	for j := 0; j < end; j++ {
		b.EraseEdge(vertices[i], vertices[j])
		i = j
	}
}

func (b *PolygonBuilder) AssembleLoop(v0, v1 Point, unused_edges *[]Edge) *Loop {
	path := []Point{}        // The path so far.
	index := map[Point]int{} // Maps a vertex to its index in "path"
	path = append(path, v0)
	path = append(path, v1)
	index[v1] = 1

	for len(path) >= 2 {
		// Note that "v0" and "v1" become invalid if "path" is modified.
		v0 = path[len(path)-2]
		v1 = path[len(path)-1]
		var v2 Point
		v2_found := false
		if vset, ok := b.edges[v1]; ok {
			for _, v := range *vset {
				// We prefer the leftmost outgoing edge,
				// ignoring any reverse edges.
				if v == v0 {
					continue
				}
				if !v2_found || OrderedCCW(v0, v2, v, v1) {
					v2 = v
				}
				v2_found = true
			}
		}
		if !v2_found {
			// We've hit a dead end. Remove this edge and backtrack.
			*unused_edges = append(*unused_edges, Edge{v0, v1})
			b.EraseEdge(v0, v1)
			delete(index, v1)
			i := len(path) - 1
			path = append(path[:i], path[i+1:]...)
		} else {
			if _, ok := index[v2]; !ok {
				// This is the first time we've visited this
				// vertex.
				index[v2] = len(path)
				path = append(path, v2)
			} else {
				// We've completed a loop. Throw away any
				// initial vertices that are not part of the
				// loop.
				path = path[index[v2]:]
				loop := NewLoopFromPath(path)
				if b.options.undirected_edges && !loop.IsNormalized() {
					return b.AssembleLoop(path[1], path[0], unused_edges)
				}
				return loop
			}
		}
	}
	return nil
}
