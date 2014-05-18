package s2

import (
	"code.google.com/p/gos2/r3"
	"math"
	"testing"
)

func TestOriginPoint(t *testing.T) {
	if math.Abs(OriginPoint().Norm()-1) > 1e-16 {
		t.Errorf("Origin point norm = %v, want 1", OriginPoint().Norm())
	}
}

func TestPointCross(t *testing.T) {
	tests := []struct {
		p1x, p1y, p1z, p2x, p2y, p2z float64
	}{
		{1, 0, 0, 1, 0, 0},
		{1, 0, 0, 0, 1, 0},
		{0, 1, 0, 1, 0, 0},
		{1, 2, 3, -4, 5, -6},
	}
	for _, test := range tests {
		p1 := PointFromCoords(test.p1x, test.p1y, test.p1z)
		p2 := PointFromCoords(test.p2x, test.p2y, test.p2z)
		result := p1.PointCross(p2)
		if !float64Eq(result.Norm(), 1) {
			t.Errorf("|%v ⨯ %v| = %v, want 1", p1, p2, result.Norm())
		}
		if x := result.Dot(p1.Vector); !float64Eq(x, 0) {
			t.Errorf("|(%v ⨯ %v) · %v| = %v, want 0", p1, p2, p1, x)
		}
		if x := result.Dot(p2.Vector); !float64Eq(x, 0) {
			t.Errorf("|(%v ⨯ %v) · %v| = %v, want 0", p1, p2, p2, x)
		}
	}
}

func TestCCW(t *testing.T) {
	tests := []struct {
		p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z float64
		want                                        bool
	}{
		{1, 0, 0, 0, 1, 0, 0, 0, 1, true},
		{0, 1, 0, 0, 0, 1, 1, 0, 0, true},
		{0, 0, 1, 1, 0, 0, 0, 1, 0, true},
		{1, 1, 0, 0, 1, 1, 1, 0, 1, true},
		{-3, -1, 4, 2, -1, -3, 1, -2, 0, true},

		// All degenerate cases of CCW().  Let M_1, M_2, ... be the sequence of
		// submatrices whose determinant sign is tested by that function.  Then the
		// i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
		//
		//    det(M) = 0
		//    det(M_j) = 0 for j < i
		//    det(M_i) != 0
		//    A < B < C in lexicographic order.
		// det(M_1) = b0*c1 - b1*c0
		{-3, -1, 0, -2, 1, 0, 1, -2, 0, false},
		// det(M_2) = b2*c0 - b0*c2
		{-6, 3, 3, -4, 2, -1, -2, 1, 4, false},
		// det(M_3) = b1*c2 - b2*c1
		{0, -1, -1, 0, 1, -2, 0, 2, 1, false},
		// From this point onward, B or C must be zero, or B is proportional to C.
		// det(M_4) = c0*a1 - c1*a0
		{-1, 2, 7, 2, 1, -4, 4, 2, -8, false},
		// det(M_5) = c0
		{-4, -2, 7, 2, 1, -4, 4, 2, -8, false},
		// det(M_6) = -c1
		{0, -5, 7, 0, -4, 8, 0, -2, 4, false},
		// det(M_7) = c2*a0 - c0*a2
		{-5, -2, 7, 0, 0, -2, 0, 0, -1, false},
		// det(M_8) = c2
		{0, -2, 7, 0, 0, 1, 0, 0, 2, false},
	}

	for _, test := range tests {
		p1 := PointFromCoords(test.p1x, test.p1y, test.p1z)
		p2 := PointFromCoords(test.p2x, test.p2y, test.p2z)
		p3 := PointFromCoords(test.p3x, test.p3y, test.p3z)
		result := CCW(p1, p2, p3)
		if result != test.want {
			t.Errorf("CCW(%v, %v, %v) = %v, want %v", p1, p2, p3, result, test.want)
		}
		if test.want {
			// For these cases we can test the reversibility condition
			result = CCW(p3, p2, p1)
			if result == test.want {
				t.Errorf("CCW(%v, %v, %v) = %v, want %v", p3, p2, p1, result, !test.want)
			}
		}
	}
}

func TestPointDistance(t *testing.T) {
	tests := []struct {
		x1, y1, z1 float64
		x2, y2, z2 float64
		want       float64 // radians
	}{
		{1, 0, 0, 1, 0, 0, 0},
		{1, 0, 0, 0, 1, 0, math.Pi / 2},
		{1, 0, 0, 0, 1, 1, math.Pi / 2},
		{1, 0, 0, -1, 0, 0, math.Pi},
		{1, 2, 3, 2, 3, -1, 1.2055891055045298},
	}
	for _, test := range tests {
		p1 := PointFromCoords(test.x1, test.y1, test.z1)
		p2 := PointFromCoords(test.x2, test.y2, test.z2)
		if a := p1.Distance(p2).Radians(); !float64Eq(a, test.want) {
			t.Errorf("%v.Distance(%v) = %v, want %v", p1, p2, a, test.want)
		}
		if a := p2.Distance(p1).Radians(); !float64Eq(a, test.want) {
			t.Errorf("%v.Distance(%v) = %v, want %v", p2, p1, a, test.want)
		}
	}
}

func TestApproxEqual(t *testing.T) {
	epsilon := 1e-14
	tests := []struct {
		x1, y1, z1 float64
		x2, y2, z2 float64
		want       bool
	}{
		{1, 0, 0, 1, 0, 0, true},
		{1, 0, 0, 0, 1, 0, false},
		{1, 0, 0, 0, 1, 1, false},
		{1, 0, 0, -1, 0, 0, false},
		{1, 2, 3, 2, 3, -1, false},
		{1, 0, 0, 1 * (1 + epsilon), 0, 0, true},
		{1, 0, 0, 1 * (1 - epsilon), 0, 0, true},
		{1, 0, 0, 1 + epsilon, 0, 0, true},
		{1, 0, 0, 1 - epsilon, 0, 0, true},
		{1, 0, 0, 1, epsilon, 0, true},
		{1, 0, 0, 1, epsilon, epsilon, false},
		{1, epsilon, 0, 1, -epsilon, epsilon, false},
	}
	for _, test := range tests {
		p1 := PointFromCoords(test.x1, test.y1, test.z1)
		p2 := PointFromCoords(test.x2, test.y2, test.z2)
		if got := p1.ApproxEqual(p2); got != test.want {
			t.Errorf("%v.ApproxEqual(%v), got %v want %v", p1, p2, got, test.want)
		}
	}
}

func TestColinearPoints(t *testing.T) {
	// The following points happen to be *exactly colinear* along a line
	// that is approximately tangent to the surface of the unit sphere.
	// In fact, "c" is the exact midpoint of the line segment "ab". All of
	// these points are close enough to unit length to satisfy IsUnitLength.
	a := Point{r3.Vector{0.72571927877036835, 0.46058825605889098, 0.51106749730504852}}
	b := Point{r3.Vector{0.7257192746638208, 0.46058826573818168, 0.51106749441312738}}
	c := Point{r3.Vector{0.72571927671709457, 0.46058826089853633, 0.51106749585908795}}
	c_sub_a := c.Sub(a.Vector)
	b_sub_c := b.Sub(c.Vector)
	if c_sub_a != b_sub_c {
		t.Errorf("%v != %v", c_sub_a, b_sub_c)
	}

	if RobustCCW(a, b, c) == 0 {
		t.Errorf("%v == %v", RobustCCW(a, b, c), 0)
	}

	if RobustCCW(a, b, c) != RobustCCW(b, c, a) {
		t.Errorf("%v != %v", RobustCCW(a, b, c), RobustCCW(b, c, a))
	}

	if RobustCCW(a, b, c) != -RobustCCW(c, b, a) {
		t.Errorf("%v != %v", RobustCCW(a, b, c), -RobustCCW(c, b, a))
	}

	// The points "x1" and "x2" are exactly proportional, i.e. they both
	// lie on a common line through the origin. Both points are considered
	// to be normalized, and in fact they both satisfy (x == x.Normalize()).
	// Therefore the triangle (x1, x2, -x1) consists of three distinct
	// points that all lie on a common line through the origin.
	x1 := Point{r3.Vector{0.99999999999999989, 1.4901161193847655e-08, 0}}
	x2 := Point{r3.Vector{1, 1.4901161193847656e-08, 0}}
	neg_x1 := Point{x1.Neg()}
	if x1.Vector != x1.Normalize() {
		t.Errorf("%v != %v", x1, x1.Normalize())
	}
	if x2.Vector != x2.Normalize() {
		t.Errorf("%v != %v", x2, x2.Normalize())
	}

	if RobustCCW(x1, x2, neg_x1) == 0 {
		t.Errorf("%v == %v", RobustCCW(x1, x2, neg_x1), 0)
	}

	if RobustCCW(x1, x2, neg_x1) != RobustCCW(x2, neg_x1, x1) {
		t.Errorf("%v != %v", RobustCCW(x1, x2, neg_x1), RobustCCW(x2, neg_x1, x1))
	}

	if RobustCCW(x1, x2, neg_x1) != -RobustCCW(neg_x1, x2, x1) {
		t.Errorf("%v != %v", RobustCCW(x1, x2, neg_x1), -RobustCCW(neg_x1, x2, x1))
	}

	// Here are two more points that are distinct, exactly proportional,
	// and that satisfy (x == x.Normalize()).
	x3 := PointFromCoords(1, 1, 1)
	x4 := Point{x3.Mul(0.99999999999999989)}
	neg_x3 := Point{x3.Neg()}
	if x3.Vector != x3.Normalize() {
		t.Errorf("%v != %v", x3, x3.Normalize())
	}
	if x4.Vector != x4.Normalize() {
		t.Errorf("%v != %v", x4, x4.Normalize())
	}
	if RobustCCW(x3, x4, neg_x3) == 0 {
		t.Errorf("%v == %v", RobustCCW(x3, x4, neg_x3), 0)
	}
}
