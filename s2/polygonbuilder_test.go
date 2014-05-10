package s2

import (
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
