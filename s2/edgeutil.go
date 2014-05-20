package s2

import (
	"code.google.com/p/gos2/s1"
	"math"
)

func WedgeIntersects(a0, ab1, a2, b0, b2 Point) bool {
	// For A not to intersect B (where each loop interior is defined to be
	// its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
	// note that it's important to write these conditions as negatives
	// (!OrderedCCW(a,b,c,o) rather than OrderedCCW(c,b,a,o)) to get
	// correct results when two vertices are the same.
	return !OrderedCCW(a0, b2, b0, ab1) && OrderedCCW(b0, a2, a0, ab1)
}

func WedgeContains(a0, ab1, a2, b0, b2 Point) bool {
	// For A to contain B (where each loop interior is defined to be its
	// left side), the CCW edge order around ab1 must be a2 b2 b0 a0. We
	// split this test into two parts that test three vertices each.
	return OrderedCCW(a2, b2, b0, ab1) && OrderedCCW(b0, a0, a2, ab1)
}

// This is named GetDistance() in the C++ API.
func (x Point) DistanceToEdgeWithNormal(a, b, a_cross_b Point) s1.Angle {
	// There are three cases. If X is located in the spherical wedge
	// defined by A, B, and the axis A x B, then the closest point is on
	// the segment AB. Otherwise the closest point is either A or B; the
	// dividing line between these two cases is the great circle passing
	// through (A x B) and the midpoint of AB.
	if CCW(a_cross_b, a, x) && CCW(x, b, a_cross_b) {
		// The closest point to X lies on the segment AB. We compute
		// the distance to the corresponding great circle. The result
		// is accurate for small distances but not necessarily for
		// large distances (approaching Pi/2).
		//
		// TODO: sanity check a != b
		sin_dist := math.Abs(x.Vector.Dot(a_cross_b.Vector)) / a_cross_b.Vector.Norm()
		return s1.Angle(s1.Angle(math.Asin(math.Min(1.0, sin_dist))).Radians())
	}
	// Otherwise, the closest point is either A or B. The cheapest method is
	// just to compute the minimum of the two linear (as opposed to spherical)
	// distances and convert the result to an angle. Again, this method is
	// accurate for small but not large distances (approaching Pi).
	xa := x.Vector.Sub(a.Vector).Norm2()
	xb := x.Vector.Sub(b.Vector).Norm2()
	linear_dist2 := math.Min(xa, xb)
	return s1.Angle(2 * math.Asin(math.Min(1.0, 0.5*math.Sqrt(linear_dist2))))
}

// This is named GetDistance() in the C++ API.
func (x Point) DistanceToEdge(a, b Point) s1.Angle {
	return x.DistanceToEdgeWithNormal(a, b, a.PointCross(b))
}

func (x Point) ClosestPointWithNormal(a, b, a_cross_b Point) Point {
	// Find the closest point to X along the great circle through AB.
	dx := x.Vector.Dot(a_cross_b.Vector) / a_cross_b.Vector.Norm2()
	p := Point{x.Vector.Sub(a_cross_b.Vector.Mul(dx))}

	// If this point is on the edge AB, then it's the closest point.
	if CCW(a_cross_b, a, p) && CCW(p, b, a_cross_b) {
		return p
	}

	// Otherwise, the closest point is either A or B.
	if x.Vector.Sub(a.Vector).Norm2() <= x.Vector.Sub(b.Vector).Norm2() {
		return a
	}
	return b
}

func (x Point) ClosestPoint(a, b Point) Point {
	return x.ClosestPointWithNormal(a, b, a.PointCross(b))
}

type EdgeCrosser struct {
	a         *Point
	b         *Point
	a_cross_b Point

	// The fields below are updated for each vertex in the chain
	c   *Point // Previous vertex in the vertex chain.
	acb int    // The orientation of the triangle ACB.
}

func NewEdgeCrosser(a, b, c *Point) EdgeCrosser {
	ec := EdgeCrosser{
		a:         a,
		b:         b,
		a_cross_b: Point{a.Cross(b.Vector)},
	}
	ec.RestartAt(c)
	return ec
}

func (e *EdgeCrosser) RestartAt(c *Point) {
	e.c = c
	e.acb = -RobustCCW2(*e.a, *e.b, *e.c, e.a_cross_b)
}

func (e *EdgeCrosser) EdgeOrVertexCrossing(d *Point) bool {
	c := e.c
	crossing := e.RobustCrossing(d)
	if crossing < 0 {
		return false
	}
	if crossing > 0 {
		return true
	}
	return VertexCrossing(*e.a, *e.b, *c, *d)
}

func (e *EdgeCrosser) RobustCrossing(d *Point) int {
	// For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC
	// must all be oriented the same way (CW or CCW). We keep the
	// orientation of ACB as part of our state. When each new point D
	// arrives, we compute the orientation of BDA and check whether it
	// matches ACB. This checks whether the points C and D are on opposite
	// sides of the great circle through AB.

	// Recall that RobustCCW is invariant with respect to rotating its
	// arguments, i.e. ABC has the same orientation as BDA.
	bda := RobustCCW2(*e.a, *e.b, *d, e.a_cross_b)
	var result int
	if bda == -e.acb && bda != 0 {
		result = -1 // Most common case -- triangles have opposite orientations.
	} else if (bda & e.acb) == 0 {
		result = 0 // At least one value is zero -- two vertices are identical.
	} else { // Slow path.
		result = e.RobustCrossingInternal(d)
	}
	// Now save the current vertex D as the next vertex C, and also save the
	// orientation of the new triangle ACB (which is opposite to the current
	// triangle BDA).
	e.c = d
	e.acb = -bda
	return result
}

func (e *EdgeCrosser) RobustCrossingInternal(d *Point) int {
	// ACB and BDA have the appropriate orientations, so now we check the
	// triangles CBD and DAC.
	c_cross_d := e.c.Cross(d.Vector)
	cbd := -RobustCCW2(*e.c, *d, *e.b, Point{c_cross_d})
	if cbd != e.acb {
		return -1
	}
	dac := RobustCCW2(*e.c, *d, *e.a, Point{c_cross_d})
	if dac == e.acb {
		return 1
	}
	return -1
}

func VertexCrossing(a, b, c, d Point) bool {
	// If A == B or C == D there is no intersection. We need to check this
	// case first in case 3 or more input points are identical.
	if a == b || c == d {
		return false
	}

	// If any other pair of vertices is equal, there is a crossing iff
	// OrderedCCW indicates that the edge AB is further CCW around the
	// shared vertex O (either A or B) than the edge CD, starting from an
	// arbitrary fixed reference point.
	if a == d {
		return OrderedCCW(Point{a.Ortho()}, c, b, a)
	}
	if b == c {
		return OrderedCCW(Point{b.Ortho()}, d, a, b)
	}
	if a == c {
		return OrderedCCW(Point{a.Ortho()}, d, b, a)
	}
	if b == d {
		return OrderedCCW(Point{b.Ortho()}, c, a, b)
	}
	return false
}

// The class computes a bounding rectangle that contains all edges
// defined by a vertex chain v0,v1,v2,... All vertices must be unit length.
// Note that the bounding rectangle of an edge can be larger than the
// bounding rectangle of its endpoints, e.g. consider an edge that passes
// through the north pole.
type RectBounder struct {
	a      *Point // The previous vertex in the chain.
	latlng LatLng // The corresponding lat-lng.
	bound  Rect   // The current bounding rectangle.
}

func NewRectBounder() RectBounder {
	return RectBounder{bound: EmptyRect()}
}

func (r RectBounder) Bound() Rect {
	return r.bound
}

func (r *RectBounder) AddPoint(b *Point) {
	ll := LatLngFromPoint(*b)
	if r.bound.IsEmpty() {
		r.bound = r.bound.AddPoint(ll)
	} else {
		// We can't just call bound.AddPoint(ll) here, since we need to
		// ensure that all the longitudes between "a" and "b" are
		// included.
		r.bound = r.bound.Union(RectFromPointPair(r.latlng, ll))

		// Check whether the min/max latitude occurs in the edge
		// interior. We find the normal to the plane containing AB,
		// and then a vector "dir" in this plane that also passes
		// through the equator. We use RobustCrossProd to ensure that
		// the edge normal is accurate even when the two points are
		// very close together.
		a_cross_b := r.a.PointCross(*b)
		dir := a_cross_b.Cross(PointFromCoords(0, 0, 1).Vector)
		da := dir.Dot(r.a.Vector)
		db := dir.Dot(b.Vector)
		if da*db < 0 {
			// min/max latitude occurs in the edge interior.
			abslat := math.Acos(math.Abs(a_cross_b.Z / a_cross_b.Norm()))
			if da < 0 {
				// It's possible that abslat < r.Lat.Lo due to
				// numerical errors.
				r.bound.Lat.Hi = math.Max(abslat, r.bound.Lat.Hi)
			} else {
				r.bound.Lat.Lo = math.Min(-abslat, r.bound.Lat.Lo)
			}

			// If the edge comes very close to the north or south
			// pole then we may not be certain which side of the
			// pole it is on. We handle this by expanding the
			// longitude bounds if the maximum latitude is
			// approximately Pi/2.
			if abslat >= math.Pi/2-1e-15 {
				r.bound.Lng = s1.FullInterval()
			}
		}
	}
	r.a = b
	r.latlng = ll
}