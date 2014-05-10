package s2

import (
	"code.google.com/p/gos2/r1"
	"code.google.com/p/gos2/s1"
	"math"
	"strconv"
	"strings"
	"testing"
)

func makeloop(s string) *Loop {
	points := strings.Split(s, ",")
	path := []Point{}
	for _, p := range points {
		p = strings.Trim(p, " ")
		degs := strings.Split(p, ":")
		lat, _ := strconv.ParseFloat(degs[0], 64)
		lng, _ := strconv.ParseFloat(degs[1], 64)
		ll := LatLngFromDegrees(lat, lng)
		path = append(path, PointFromLatLng(ll))
	}
	return NewLoopFromPath(path)
}

var (
	// The northern hemisphere, defined using two pairs of antipodal points.
	north_hemi = makeloop("0:-180, 0:-90, 0:0, 0:90")

	// The northern hemisphere, defined using three points 120 degrees apart.
	north_hemi3 = makeloop("0:-180, 0:-60, 0:60")

	// The southern hemisphere, defined using two pairs of antipodal points.
	south_hemi = makeloop("0:90, 0:0, 0:-90, 0:-180")

	// The western hemisphere, defined using two pairs of antipodal points.
	west_hemi = makeloop("0:-180, -90:0, 0:0, 90:0")

	// The eastern hemisphere, defined using two pairs of antipodal points.
	east_hemi = makeloop("90:0, 0:0, -90:0, 0:-180")

	// The "near" hemisphere, defined using two pairs of antipodal points.
	near_hemi = makeloop("0:-90, -90:0, 0:90, 90:0")

	// The "far" hemisphere, defined using two pairs of antipodal points.
	far_hemi = makeloop("90:0, 0:90, -90:0, 0:-90")

	// A spiral stripe that slightly over-wraps the equator.
	candy_cane = makeloop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")

	// A small clockwise loop in the northern & eastern hemispheres.
	small_ne_cw = makeloop("35:20, 45:20, 40:25")

	// Loop around the north pole at 80 degrees.
	arctic_80 = makeloop("80:-150, 80:-30, 80:90")

	// Loop around the south pole at 80 degrees.
	antarctic_80 = makeloop("-80:120, -80:0, -80:-120")

	// A completely degenerate triangle along the equator that RobustCCW()
	// considers to be CCW.
	line_triangle = makeloop("0:1, 0:3, 0:2")

	// A nearly-degenerate CCW chevron near the equator with very long sides
	// (about 80 degrees).  Its area is less than 1e-640, which is too small
	// to represent in double precision.
	skinny_chevron = makeloop("0:0, -1e-320:80, 0:1e-320, 1e-320:80")

	// A diamond-shaped loop around the point 0:180.
	loop_a = makeloop("0:178, -1:180, 0:-179, 1:-180")

	// Another diamond-shaped loop around the point 0:180.
	loop_b = makeloop("0:179, -1:180, 0:-178, 1:-180")

	// The intersection of A and B.
	a_intersect_b = makeloop("0:179, -1:180, 0:-179, 1:-180")

	// The union of A and B.
	a_union_b = makeloop("0:178, -1:180, 0:-178, 1:-180")

	// A minus B (concave).
	a_minus_b = makeloop("0:178, -1:180, 0:179, 1:-180")

	// B minus A (concave).
	b_minus_a = makeloop("0:-179, -1:180, 0:-178, 1:-180")

	// A shape gotten from a by adding one triangle to one edge, and
	// subtracting another triangle on an opposite edge.
	loop_c = makeloop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180")

	// A shape gotten from a by adding one triangle to one edge, and
	// adding another triangle on an opposite edge.
	loop_d = makeloop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180")
)

func TestGetRectBound(t *testing.T) {
	if !candy_cane.bound.Lng.IsFull() {
		t.Errorf("%v.IsFull() == false", candy_cane.bound.Lng)
	}
	deg := candy_cane.bound.Lo().Lat.Degrees()
	if deg >= -20 {
		t.Errorf("%v >= -20", deg)
	}
	deg = candy_cane.bound.Hi().Lat.Degrees()
	if deg <= 10 {
		t.Errorf("%v <= 10", deg)
	}

	if !small_ne_cw.bound.IsFull() {
		t.Errorf("%v.IsFull() == false", small_ne_cw.bound)
	}

	var p1, p2 LatLng
	var rect Rect

	p1 = LatLngFromDegrees(80, -180)
	p2 = LatLngFromDegrees(90, 180)
	rect = Rect{
		Lat: r1.Interval{p1.Lat.Radians(), p2.Lat.Radians()},
		Lng: s1.Interval{p1.Lng.Radians(), p2.Lng.Radians()},
	}

	if !arctic_80.bound.Equal(rect) {
		t.Errorf("%v.Equal(%v) == false", arctic_80.bound, rect)
	}

	p1 = LatLngFromDegrees(-90, -180)
	p2 = LatLngFromDegrees(-80, 180)
	rect = Rect{
		Lat: r1.Interval{p1.Lat.Radians(), p2.Lat.Radians()},
		Lng: s1.Interval{p1.Lng.Radians(), p2.Lng.Radians()},
	}

	if !antarctic_80.bound.Equal(rect) {
		t.Errorf("%v.Equal(%v) == false", antarctic_80.bound, rect)
	}

	// Create a loop that contains the complement of the "arctic_80" loop.
	arctic_80_inv := arctic_80.Clone()
	arctic_80_inv.Invert()
	// The highest altitude of each edge is attained at its midpoint
	mid := arctic_80_inv.vertex(0).Add(arctic_80_inv.vertex(1).Vector).Mul(0.5)
	want := arctic_80_inv.bound.Hi().Lat.Radians()
	got := LatLngFromPoint(Point{mid}).Lat.Radians()
	if math.Abs(got-want) > 1e-14 {
		t.Errorf("%v != %v", want, got)
	}

	if !south_hemi.bound.Lng.IsFull() {
		t.Errorf("%v.IsFull() == false", south_hemi.bound.Lng)
	}

	i := r1.Interval{-math.Pi / 2, 0}
	if !south_hemi.bound.Lat.Equal(i) {
		t.Errorf("%v.Equal(%v) == false", south_hemi.bound.Lat, i)
	}
}
