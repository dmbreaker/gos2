package s2

import (
	"math"
	"testing"
)

// returns normalized point
func pc(x, y, z float64) Point {
	return PointFromCoords(x, y, z)
}

func TestEquality(t *testing.T) {
	tests := []struct {
		a    Point
		b    Point
		want bool
	}{
		{pc(0, 0, 0), pc(0, 0, 0), true},
		{pc(.5, 0, 0), pc(.5, 0, 0), true},
	}
	for _, test := range tests {
		got := test.a == test.b
		if got != test.want {
			t.Errorf("%v == %v = %v, want %v", test.a, test.b, got, test.want)
		}
	}
}

func TestWedges(t *testing.T) {
	tests := []struct {
		a0            Point
		ab1           Point
		a2            Point
		b0            Point
		b2            Point
		contains      bool
		intersects    bool
		wedgeRelation int
	}{
		// Intersection in one wedge.
		{pc(-1, 0, 10), pc(0, 0, 1), pc(1, 2, 10),
			pc(0, 1, 10), pc(1, -2, 10),
			false, true, WEDGE_PROPERLY_OVERLAPS},

		// Intersection in two wedges.
		{pc(-1, -1, 10), pc(0, 0, 1), pc(1, -1, 10),
			pc(1, 0, 10), pc(-1, 1, 10),
			false, true, WEDGE_PROPERLY_OVERLAPS},

		// Normal containment.
		{pc(-1, -1, 10), pc(0, 0, 1), pc(1, -1, 10),
			pc(-1, 0, 10), pc(1, 0, 10),
			true, true, WEDGE_PROPERLY_CONTAINS},

		// Containment with equality on one side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(-1, -1, 10),
			pc(2, 1, 10), pc(1, -5, 10),
			true, true, WEDGE_PROPERLY_CONTAINS},

		// Containment with equality on the other side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(-1, -1, 10),
			pc(1, -2, 10), pc(-1, -1, 10),
			true, true, WEDGE_PROPERLY_CONTAINS},

		// Containment with equality on both sides.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(4, -5, 10),
			pc(-2, 3, 10), pc(4, -5, 10),
			true, true, WEDGE_EQUALS},

		// Disjoint with equality on one side.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(4, -5, 10),
			pc(4, -5, 10), pc(-2, -3, 10),
			false, false, WEDGE_IS_DISJOINT},

		// Disjoint with equality on the other side.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(0, 5, 10),
			pc(4, -5, 10), pc(-2, 3, 10),
			false, false, WEDGE_IS_DISJOINT},

		// Disjoint with equality on both sides.
		{pc(-2, 3, 10), pc(0, 0, 1), pc(4, -5, 10),
			pc(4, -5, 10), pc(-2, 3, 10),
			false, false, WEDGE_IS_DISJOINT},

		// B contains A with equality on one side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(1, -5, 10),
			pc(2, 1, 10), pc(-1, -1, 10),
			false, true, WEDGE_IS_PROPERLY_CONTAINED},

		// B contains A with equality on the other side.
		{pc(2, 1, 10), pc(0, 0, 1), pc(1, -5, 10),
			pc(-2, 1, 10), pc(1, -5, 10),
			false, true, WEDGE_IS_PROPERLY_CONTAINED},
	}
	for _, test := range tests {
		got := WedgeContains(test.a0, test.ab1, test.a2, test.b0, test.b2)
		if got != test.contains {
			t.Errorf("WedgeContains(): got %v, want %v", got, test.contains)
		}
		got = WedgeIntersects(test.a0, test.ab1, test.a2, test.b0, test.b2)
		if got != test.intersects {
			t.Errorf("WedgeIntersects(): got %v, want %v", got, test.intersects)
		}
		relation := GetWedgeRelation(test.a0, test.ab1, test.a2, test.b0, test.b2)
		if relation != test.wedgeRelation {
			t.Errorf("GetWedgeRelation(): got %v, want %v", relation, test.wedgeRelation)
		}
	}
}

// Given a point X and an edge AB, check that the distance from X to AB is
// "distance_radians" and the closest point on AB is "expected_closest"
func TestDistanceToEdge(t *testing.T) {
	tests := []struct {
		x                Point
		a                Point
		b                Point
		distance_radians float64
		expected_closest Point
	}{
		{pc(1, 0, 0), pc(1, 0, 0), pc(0, 1, 0), 0, pc(1, 0, 0)},
		{pc(0, 1, 0), pc(1, 0, 0), pc(0, 1, 0), 0, pc(0, 1, 0)},
		{pc(1, 3, 0), pc(1, 0, 0), pc(0, 1, 0), 0, pc(1, 3, 0)},
		{pc(0, 0, 1), pc(1, 0, 0), pc(0, 1, 0), math.Pi / 2, pc(1, 0, 0)},
		{pc(0, 0, -1), pc(1, 0, 0), pc(0, 1, 0), math.Pi / 2, pc(1, 0, 0)},
		{pc(-1, -1, 0), pc(1, 0, 0), pc(0, 1, 0), 0.75 * math.Pi, pc(0, 0, 0)},

		{pc(0, 1, 0), pc(1, 0, 0), pc(1, 1, 0), math.Pi / 4, pc(1, 1, 0)},
		{pc(0, -1, 0), pc(1, 0, 0), pc(1, 1, 0), math.Pi / 2, pc(1, 0, 0)},

		{pc(0, -1, 0), pc(1, 0, 0), pc(-1, 1, 0), math.Pi / 2, pc(1, 0, 0)},
		{pc(-1, -1, 0), pc(1, 0, 0), pc(-1, 1, 0), math.Pi / 2, pc(-1, 1, 0)},

		{pc(1, 1, 1), pc(1, 0, 0), pc(0, 1, 0), math.Asin(math.Sqrt(1. / 3)), pc(1, 1, 0)},
		{pc(1, 1, -1), pc(1, 0, 0), pc(0, 1, 0), math.Asin(math.Sqrt(1. / 3)), pc(1, 1, 0)},

		{pc(-1, 0, 0), pc(1, 1, 0), pc(1, 1, 0), 0.75 * math.Pi, pc(1, 1, 0)},
		{pc(0, 0, -1), pc(1, 1, 0), pc(1, 1, 0), math.Pi / 2, pc(1, 1, 0)},
		{pc(-1, 0, 0), pc(1, 0, 0), pc(1, 0, 0), math.Pi, pc(1, 0, 0)},
	}
	for _, test := range tests {
		got := test.x.DistanceToEdge(test.a, test.b).Radians()
		if math.Abs(got-test.distance_radians) > 1e-14 {
			t.Errorf("%v.DistanceToEdge(%v, %v) = %v, want %v",
				test.x, test.a, test.b, got, test.distance_radians)
		}

		closest := test.x.ClosestPoint(test.a, test.b)
		if test.expected_closest.ApproxEqual(pc(0, 0, 0)) {
			if !closest.ApproxEqual(test.a) && !closest.ApproxEqual(test.b) {
				t.Errorf("%v != %v || %v != %v", closest, test.a, closest, test.b)
			}
		} else {
			if !closest.ApproxEqual(test.expected_closest) {
				t.Errorf("%v != %v", closest, test.expected_closest)
			}
		}
	}
}
