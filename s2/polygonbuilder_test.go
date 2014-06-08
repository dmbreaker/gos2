package s2

import (
	"code.google.com/p/gos2/r3"
	"code.google.com/p/gos2/s1"
	"fmt"
	"math"
	"math/rand"
	"testing"
)

func TestAddEdge(t *testing.T) {
	pb := NewPolygonBuilder(DIRECTED_XOR())
	p1 := PointFromLatLng(LatLngFromDegrees(47, -122))
	p2 := PointFromLatLng(LatLngFromDegrees(48, -122))
	if ok := pb.AddEdge(p1, p2); !ok {
		t.Errorf("Failed to add edge")
	}
}

func TestExistingDirectedEdge(t *testing.T) {
	pb := NewPolygonBuilder(DIRECTED_XOR())
	p1 := PointFromLatLng(LatLngFromDegrees(47, -122))
	p2 := PointFromLatLng(LatLngFromDegrees(48, -122))
	pb.AddEdge(p1, p2)
	pb.AddEdge(p2, p1) // overwrite edge (p1, p2) with (p2, p1)
	if pb.HasEdge(p1, p2) {
		t.Errorf("Undirected edge found in directed graph")
	}
}

func TestHasEdge(t *testing.T) {
	pb := NewPolygonBuilder(DIRECTED_XOR())
	p1 := PointFromLatLng(LatLngFromDegrees(47, -122))
	p2 := PointFromLatLng(LatLngFromDegrees(48, -122))
	pb.AddEdge(p1, p2)
	if !pb.HasEdge(p1, p2) {
		t.Errorf("Edge not found")
	}
}

func TestSingleLoop(t *testing.T) {
	rect := []float64{
		-172.08984375,
		73.1758971742261,
		-21.4453125,
		-25.641526373065755,
	}
	vertices := []Point{
		PointFromLatLng(LatLngFromDegrees(rect[1], rect[0])),
		PointFromLatLng(LatLngFromDegrees(rect[3], rect[0])),
		PointFromLatLng(LatLngFromDegrees(rect[3], rect[2])),
		PointFromLatLng(LatLngFromDegrees(rect[1], rect[2])),
	}
	pb := NewPolygonBuilder(DIRECTED_XOR())
	for i := 0; i < 4; i++ {
		pb.AddEdge(vertices[i], vertices[(i+1)%len(vertices)])
	}

	loops := []*Loop{}
	edges := []Edge{}
	if !pb.AssembleLoops(&loops, &edges) {
		t.Errorf("unused edges")
	}

	if len(loops) != 1 {
		t.Errorf("len(%v) == %v, want 1", loops, len(loops))
	}
}

type Chain struct {
	str    string
	closed bool
}

type TestCase struct {
	undirectedEdges    int
	xorEdges           int
	canSplit           bool
	minMerge, maxMerge float64
	minVertexAngle     float64
	chainsIn           []Chain
	loopsOut           []string
	numUnusedEdges     int
}

func TestAssembleLoops(t *testing.T) {
	/*
		tests := []TestCase{
			{0, 0, true, 0.0, 10.0, 90.0, []Chain{Chain{"", false}}, []string{""}, 0},

			{0, 0, true, 0.0, 4.0, 15.0, []Chain{
				Chain{"0:0, 0:10, 10:5", true},
				Chain{"0:0, 5:5", false},
				Chain{"10:5, 20:7, 30:10, 40:15, 50:3, 60:-20", false},
			}, []string{
				"0:0, 0:10, 10:5",
			}, 6},

			{0, -1, true, 0.0, 4.0, 90.0, []Chain{
				Chain{"0:0, 0:5, 5:5, 5:0", true},
				Chain{"0:5, 0:10, 5:10, 5:5", true},
				Chain{"5:0, 5:5, 10:5, 10:0", true},
				Chain{"5:5, 5:10, 10:10, 10:5", true},
				Chain{"0:10, 0:15, 0:20", false},
				Chain{"20:0, 15:0, 10:0", false},
			}, []string{
				"0:0, 0:5, 5:5, 5:0",
				"0:5, 0:10, 5:10, 5:5",
				"5:0, 5:5, 10:5, 10:0",
				"5:5, 5:10, 10:10, 10:5",
			}, 4},
		}
		for _, test := range tests {
			if !runTestCase(t, test) {
				t.Errorf("TestBuilder(t, %v) == false", test)
			}
		}
	*/
}

func runTestCase(t *testing.T, test TestCase) bool {
	for iter := 0; iter < 1; iter++ {
		options := PolygonBuilderOptions{edge_splice_fraction: .866}
		options.undirected_edges = evalTristate(test.undirectedEdges)
		options.xor_edges = evalTristate(test.xorEdges)
		minMerge := s1.Angle(test.minMerge * (math.Pi / 180)).Radians()
		maxMerge := s1.Angle(test.maxMerge * (math.Pi / 180)).Radians()
		minSin := math.Sin(s1.Angle(test.minVertexAngle * (math.Pi / 180)).Radians())

		// Half of the time we allow edges to be split into smaller
		// pieces (up to 5 levels, i.e. up to 32 pieces).
		maxSplits := max(0, rand.Intn(10)-4)
		if !test.canSplit {
			maxSplits = 0
		}

		// We choose randomly among two different values for the edge
		// fraction, just to exercise that code.
		edgeFraction := options.edge_splice_fraction
		var vertexMerge, maxPerturb float64
		if minSin < edgeFraction && oneIn(2) {
			edgeFraction = minSin
		}
		if maxSplits == 0 && oneIn(2) {
			// Turn off edge splicing completely.
			edgeFraction = 0
			vertexMerge = minMerge + smallFraction()*(maxMerge-minMerge)
			maxPerturb = 0.5 * math.Min(vertexMerge-minMerge, maxMerge-vertexMerge)
		} else {
			// Splice edges. These bounds also assume that edges
			// may be split.
			//
			// If edges are actually split, need to bump up the
			// minimum merge radius to ensure that split edges
			// in opposite directions are unified. Otherwise
			// there will be tiny degenerate loops created.
			if maxSplits > 0 {
				minMerge += 1e-15
			}
			minMerge /= edgeFraction
			maxMerge *= minSin
			if maxMerge < minMerge {
				t.Errorf("%v < %v", maxMerge, minMerge)
			}

			vertexMerge = minMerge + smallFraction()*(maxMerge-minMerge)
			maxPerturb = 0.5 * math.Min(edgeFraction*(vertexMerge-minMerge), maxMerge-vertexMerge)
		}

		// We can perturb by any amount up to the maximum, but choosing
		// a lower maximum decreases the error bounds when checking the
		// output.
		maxPerturb *= smallFraction()

		// This is the minimum length of a split edge to prevent
		// unexpected merging and/or splicing.
		minEdge := minMerge + (vertexMerge+2*maxPerturb)/minSin
		options.vertex_merge_radius = s1.Angle(vertexMerge)
		options.edge_splice_fraction = edgeFraction
		options.validate = true
		builder := NewPolygonBuilder(options)

		// On each iteration we randomly rotate the test case around
		// the sphere. This causes the PolygonBuilder to choose
		// different first edges when trying to build loops.
		x, y, z := randomFrame()
		m := r3.MatrixFromCols(x.Vector, y.Vector, z.Vector)
		for _, chain := range test.chainsIn {
			addChain(chain, m, maxSplits, maxPerturb, minEdge, &builder)
		}

		var loops []*Loop
		var unusedEdges []Edge
		if test.xorEdges < 0 {
			builder.AssembleLoops(&loops, &unusedEdges)
		} else {
			var polygon Polygon
			builder.AssemblePolygon(&polygon, &unusedEdges)
			polygon.Release(&loops)
		}

		expected := []*Loop{}
		for _, str := range test.loopsOut {
			if str != "" {
				vertices := getVertices(str, m)
				expected = append(expected, NewLoopFromPath(vertices))
			}
		}

		// We assume that the vertex locations in the expected output
		// polygon are separated from the corresponding vertex
		// locations in the input edges by at most half of the
		// minimum merge radius. Essentially this means that the
		// expected output vertices should be near the centroid of the
		// various input vertices.
		//
		// If any edges were split, we need to allow a bit more error
		// due to inaccuracies in the interpolated positions.
		// Similarly, if any vertices were perturbed, we need to bump
		// up the error to allow for numerical errors in the actual
		// perturbation.
		maxError := 0.5*minMerge + maxPerturb
		if maxSplits > 0 || maxPerturb > 0 {
			maxError += 1e-15
		}

		fmt.Println("loops", loops)

		ok0 := findMissingLoops(loops, expected, m, maxSplits, maxError, "Actual")
		ok1 := findMissingLoops(expected, loops, m, maxSplits, maxError, "Expected")
		ok2 := unexpectedUnusedEdgeCount(len(unusedEdges), test.numUnusedEdges, maxSplits)
		if ok0 || ok1 || ok2 {
			// We found a problem.
			dumpUnusedEdges(unusedEdges, m, test.numUnusedEdges)
			fmt.Printf(`During iteration %d:
  undirected: %v
  xor: %v
  maxSplits: %d
  maxPerturb: %.6g
  vertexMergeRadius: %.6g
  edgeSpliceFraction: %.6g
  minEdge: %.6g
  maxError: %.6g

`, iter, options.undirected_edges, options.xor_edges, maxSplits,
				s1.Angle(maxPerturb).Degrees(),
				options.vertex_merge_radius.Degrees(),
				options.edge_splice_fraction,
				s1.Angle(minEdge).Degrees(), s1.Angle(maxError).Degrees())
			return false
		}
	}
	return true
}

func evalTristate(state int) bool {
	if state > 0 {
		return true
	}
	if state < 0 {
		return false
	}
	return oneIn(2)
}

func smallFraction() float64 {
	r := rand.Float64()
	u := rand.Float64()
	if r < 0.3 {
		return 0.0
	}
	if r < 0.6 {
		return u
	}
	return math.Pow(1e-10, u)
}

func randomFrame() (x, y, z Point) {
	x = randomPoint()
	y = Point{x.Cross(randomPoint().Vector).Normalize()}
	z = Point{x.Cross(y.Vector).Normalize()}
	return
}

func samplePoint(s2cap Cap) Point {
	// We consider the cap axis to be the "z" axis. We choose other
	// axes to complete the coordinate frame.
	m := FrameFromPoint(s2cap.center)
	// The surface area of a spherical cap is directly proportional to its
	// height. First we choose a random height, and then we choose a random
	// point along the circle at that height.
	h := rand.Float64() * s2cap.height
	theta := 2 * math.Pi * rand.Float64()
	r := math.Sqrt(h * (2 - h)) // radius of circle
	// The result should already be very close to unit-length, but we might
	// as well make it as accurate as possible.
	return PointFromFrame(m, PointFromCoords(math.Cos(theta)*r, math.Sin(theta)*r, 1-h))
}

func Perturb(x Point, maxPerturb float64) Point {
	// Perturb "x" randomly within the radius of maxPerturb.
	if maxPerturb == 0 {
		return x
	}
	return samplePoint(CapFromCenterHeight(Point{x.Normalize()}, maxPerturb*(math.Pi/180)))
}

func addChain(chain Chain, m r3.Matrix, maxSplits int,
	maxPerturb, minEdge float64, builder *PolygonBuilder) {
	// Transform the given edge chain to the frame (x, y, z), optionally
	// split each edge into pieces and/or perturb the vertices up to the
	// given radius, and add them to the builder.
	vertices := getVertices(chain.str, m)
	if chain.closed {
		vertices = append(vertices, vertices[0])
	}
	for i := 1; i < len(vertices); i++ {
		addEdge(vertices[i-1], vertices[i], maxSplits, maxPerturb, minEdge, builder)
	}
}

func addEdge(v0, v1 Point, maxSplits int, maxPerturb, minEdge float64, builder *PolygonBuilder) {
	length := float64(v0.Angle(v1.Vector))
	if maxSplits > 0 && oneIn(2) && length >= 2*minEdge {
		// Choose an interpolation parameter such that the length of
		// each piece is at least minEdge.
		f := minEdge / length
		t := f + (1-2*f)*rand.Float64()

		// Now add the two sub-edges recursively.
		vmid := EdgeInterpolate(t, v0, v1)
		addEdge(v0, vmid, maxSplits-1, maxPerturb, minEdge, builder)
		addEdge(vmid, v1, maxSplits-1, maxPerturb, minEdge, builder)
	} else {
		builder.AddEdge(Perturb(v0, maxPerturb), Perturb(v1, maxPerturb))
	}
}

func getVertices(s string, m r3.Matrix) (vertices []Point) {
	line := makePolyline(s)
	for _, v := range line.vertices {
		vertices = append(vertices, Point{m.MulVector(v.Vector).Normalize()})
	}
	return
}

func makePolyline(s string) *Polyline {
	vertices := parsePoints(s)
	return PolylineFromPoints(vertices)
}

func findMissingLoops(actual, expected []*Loop,
	m r3.Matrix, maxSplits int, maxError float64, label string) bool {
	// Dump any loops from "actual" that are not present in "expected".
	found := false
	fmt.Println(label, actual, expected)
	for i, loop := range actual {
		if findLoop(loop, expected, maxSplits, maxError) {
			continue
		}
		fmt.Printf("%s loop %d:\n", label, i)
		for _, v := range loop.vertices {
			ll := LatLngFromPoint(Point{m.Transpose().MulVector(v.Vector)})
			fmt.Printf("  [%.6f, %.6f]\n", ll.Lat.Degrees(), ll.Lng.Degrees())
		}
		found = true
	}
	return found
}

func unexpectedUnusedEdgeCount(numActual, numExpected, maxSplits int) bool {
	if maxSplits == 0 {
		return numActual != numExpected
	} else {
		return (numActual > 0) != (numExpected > 0)
	}
}

func dumpUnusedEdges(unusedEdges []Edge, m r3.Matrix, numExpected int) {
	if len(unusedEdges) == numExpected {
		return
	}
	fmt.Printf("Wrong number of unused edges (%d expected, %d actual):\n",
		numExpected, len(unusedEdges))
	for _, edge := range unusedEdges {
		p0 := LatLngFromPoint(Point{m.Transpose().MulVector(edge.v0.Vector)})
		p1 := LatLngFromPoint(Point{m.Transpose().MulVector(edge.v1.Vector)})
		fmt.Printf("  [%.6f, %.6f] -> [%.6f, %.6f]\n",
			p0.Lat.Degrees(), p0.Lng.Degrees(),
			p1.Lat.Degrees(), p1.Lng.Degrees())
	}
}

func findLoop(loop *Loop, candidates []*Loop, maxSplits int, maxError float64) bool {
	for _, c := range candidates {
		if maxSplits == 0 {
			if loop.BoundaryApproxEquals(c, maxError) {
				return true
			}
		} else {
			if loop.BoundaryNear(c, maxError) {
				return true
			}
		}
	}
	return false
}
