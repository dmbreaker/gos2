package s2

import "fmt"

type LoopMap map[*Loop][]*Loop

type Polygon struct {
	loops       []*Loop
	bound       Rect
	ownsLoops   bool
	hasHoles    bool
	numVertices int
}

func PolygonFromLoops(loops *[]*Loop) Polygon {
	p := Polygon{
		loops:     make([]*Loop, 0, len(*loops)),
		bound:     EmptyRect(),
		ownsLoops: true,
		hasHoles:  false,
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

	for _, loop := range p.loops {
		if loop.Sign() < 0 {
			p.hasHoles = true
		} else {
			p.bound = p.bound.Union(loop.bound)
		}
	}
	return p
}

func (p *Polygon) Release(loops *[]*Loop) {
	if loops != nil {
		fmt.Println("Release", p.loops)
		copy(*loops, p.loops)
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
	children := loopMap[loop]
	for _, child := range children {
		p.InitLoop(child, depth+1, loopMap)
	}
}

func (p *Polygon) InsertLoop(newLoop, parent *Loop, loopMap LoopMap) {
	children := loopMap[parent]
	fmt.Println("children", children)
	for _, child := range children {
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
	newChildren := loopMap[newLoop]
	for i := 0; i < len(children); {
		child := children[i]
		if newLoop.ContainsNested(child) {
			newChildren = append(newChildren, child)
			copy(children[i:], children[i+1:])
			children[len(children)-1] = nil
			children = children[:len(children)-1]
		} else {
			i++
		}
	}
	children = append(children, newLoop)
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
