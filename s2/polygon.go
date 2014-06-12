package s2

import (
	"log"
)

type LoopMap map[*Loop][]*Loop

type Polygon struct {
	loops       []*Loop
	bound       Rect
	ownsLoops   bool
	hasHoles    bool
	numVertices int
}

func ContainsChild(a, b *Loop, loopMap LoopMap) bool {
	if a == b {
		return true
	}
	children := loopMap[a]
	for _, child := range children {
		if ContainsChild(child, b, loopMap) {
			return true
		}
	}
	return false
}

func (p *Polygon) Init(loops *[]*Loop) {
	p.loops = make([]*Loop, len(*loops))
	if !AreLoopsValid(*loops) {
		log.Println("invalid loops")
	}
	copy(p.loops, *loops)
	for _, loop := range p.loops {
		p.numVertices += len(loop.vertices)
	}

	loopMap := LoopMap{}
	for _, loop := range p.loops {
		p.InsertLoop(loop, nil, loopMap)
	}

	// Reorder the loops in depth-first order.
	p.loops = []*Loop{}
	p.InitLoop(nil, -1, loopMap)

	for i := 0; i < len(p.loops); i++ {
		for j := 0; j < len(p.loops); j++ {
			if i == j {
				continue
			}
			if ContainsChild(p.loops[i], p.loops[j], loopMap) != p.loops[i].ContainsNested(p.loops[j]) {
				log.Println("FAILED EQ TEST")
			}
		}
	}

	p.hasHoles = false
	p.bound = EmptyRect()

	for _, loop := range p.loops {
		if loop.Sign() < 0 {
			p.hasHoles = true
		} else {
			p.bound = p.bound.Union(loop.bound)
		}
	}
}

func NewPolygonFromLoops(loops *[]*Loop) *Polygon {
	p := &Polygon{
		bound:     EmptyRect(),
		ownsLoops: true,
		hasHoles:  false,
	}
	p.Init(loops)
	return p
}

func NewPolygonFromLoop(loop *Loop) *Polygon {
	p := &Polygon{
		loops:       []*Loop{},
		bound:       loop.Bound(),
		ownsLoops:   false,
		hasHoles:    false,
		numVertices: len(loop.vertices),
	}
	p.loops = append(p.loops, loop)
	return p
}

func (p Polygon) CapBound() Cap {
	return p.bound.CapBound()
}

func (a Polygon) ContainsOrCrosses(b *Loop) int {
	inside := false
	for _, loop := range a.loops {
		result := loop.ContainsOrCrosses(b)
		if result < 0 {
			return -1
		}
		if result > 0 {
			inside = inside != true
		}
	}
	if inside {
		return 1
	}
	return 0
}

func (a Polygon) AnyLoopContains(b *Loop) bool {
	for _, loop := range a.loops {
		if loop.ContainsLoop(b) {
			return true
		}
	}
	return false
}

func (a *Polygon) ContainsAllShells(b *Polygon) bool {
	for _, loop := range b.loops {
		if loop.Sign() < 0 {
			continue
		}
		if a.ContainsOrCrosses(loop) <= 0 {
			return false
		}
	}
	return true
}

func (a *Polygon) ExcludesAllHoles(b *Polygon) bool {
	for _, loop := range b.loops {
		if loop.Sign() > 0 {
			continue
		}
		if a.ContainsOrCrosses(loop) != 0 {
			return false
		}
	}
	return true
}

func (a *Polygon) ContainsPolygon(b *Polygon) bool {
	if len(a.loops) == 1 && len(b.loops) == 1 {
		return a.loops[0].ContainsLoop(b.loops[0])
	}
	if !a.bound.ContainsRect(b.bound) {
		if !a.bound.Lng.Union(b.bound.Lng).IsFull() {
			return false
		}
	}
	if !a.hasHoles && !b.hasHoles {
		for _, loop := range b.loops {
			if !a.AnyLoopContains(loop) {
				return false
			}
		}
		return true
	}
	return a.ContainsAllShells(b) && b.ExcludesAllHoles(a)
}

func (a Polygon) ContainsPoint(p Point) bool {
	if len(a.loops) == 1 {
		return a.loops[0].Contains(p)
	}
	if !a.bound.ContainsPoint(p) {
		return false
	}
	inside := false
	for _, loop := range a.loops {
		inside = inside != loop.Contains(p)
		if inside && !a.hasHoles { // Shells are disjoint.
			break
		}
	}
	return inside
}

func (p Polygon) ContainsCell(cell Cell) bool {
	if len(p.loops) == 1 {
		return p.loops[0].ContainsCell(cell)
	}
	if !p.bound.ContainsPoint(cell.Center()) {
		return false
	}
	loop := NewLoopFromCell(cell)
	poly := NewPolygonFromLoop(loop)
	return p.ContainsPolygon(poly)
}

func (p Polygon) MayIntersect(cell Cell) bool {
	if len(p.loops) == 1 {
		return p.loops[0].MayIntersect(cell)
	}
	if !p.bound.Intersects(cell.RectBound()) {
		return false
	}
	loop := NewLoopFromCell(cell)
	poly := NewPolygonFromLoop(loop)
	return p.Intersects(poly)
}

func (a Polygon) Intersects(b *Polygon) bool {
	if len(a.loops) == 1 && len(b.loops) == 1 {
		return a.loops[0].Intersects(b.loops[0])
	}
	if !a.bound.Intersects(b.bound) {
		return false
	}
	if !a.hasHoles && !b.hasHoles {
		for _, l1 := range a.loops {
			for _, l2 := range b.loops {
				if l1.Intersects(l2) {
					return true
				}
			}
		}
	}
	return false
}

func (p *Polygon) Release(loops *[]*Loop) {
	if loops != nil {
		*loops = append(*loops, p.loops...)
	}
	p.loops = []*Loop{}
	p.bound = EmptyRect()
	p.hasHoles = false
}

func (p *Polygon) InitLoop(loop *Loop, depth int, loopMap LoopMap) {
	if loop != nil {
		loop.depth = depth
		p.loops = append(p.loops, loop)
	}
	for _, child := range loopMap[loop] {
		p.InitLoop(child, depth+1, loopMap)
	}
}

func (p *Polygon) InsertLoop(newLoop, parent *Loop, loopMap LoopMap) {
	for _, child := range loopMap[parent] {
		if child.ContainsNested(newLoop) {
			p.InsertLoop(newLoop, child, loopMap)
			return
		}
	}
	// No loop may contain the complement of another loop. (Handling this
	// case is significantly more complicated).
	//
	// Some of the children of the parent loop may now be children of the
	// new loop.
	for i := 0; i < len(loopMap[parent]); {
		child := loopMap[parent][i]
		if newLoop.ContainsNested(child) {
			loopMap[newLoop] = append(loopMap[newLoop], child)
			//			copy(loopMap[parent][i:], loopMap[parent][i+1:])
			//			loopMap[parent][len(loopMap[parent])-1] = nil
			//			loopMap[parent] = loopMap[parent][:len(loopMap[parent])-1]
			loopMap[parent] = append(loopMap[parent][:i], loopMap[parent][i+1:]...)
		} else {
			i++
		}
	}
	loopMap[parent] = append(loopMap[parent], newLoop)
}

type PointPair struct {
	first, second Point
}

type IntPair struct {
	first, second int
}

func AreLoopsValid(loops []*Loop) bool {
	// If a loop contains an edge AB, then no other loop may contain
	// AB or BA.
	if len(loops) > 1 {
		edges := map[PointPair]IntPair{}
		for i, loop := range loops {
			for j := 0; j < len(loop.vertices); j++ {
				key := PointPair{*loop.vertex(j), *loop.vertex(j + 1)}
				if _, ok := edges[key]; !ok {
					edges[key] = IntPair{i, j}
					continue
				}
				return false
			}
		}
	}

	// Verify that no loop covers more than half of the sphere, and that
	// no two loops cross.
	for i, loop := range loops {
		if !loop.IsNormalized() {
			return false
		}
		for j := i + 1; j < len(loops); j++ {
			// This test not only checks for edge crossings, it
			// also detects cases where the two boundaries cross
			// at a shared vertex.
			if loop.ContainsOrCrosses(loops[j]) < 0 {
				return false
			}
		}
	}
	return true
}
