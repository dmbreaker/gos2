package s2

import (
	"code.google.com/p/gos2/r1"
	"code.google.com/p/gos2/s1"
	"math"
)

const (
	maxError = 1.0 / (1 << 51)
)

var (
	poleMinLat = math.Asin(math.Sqrt(1. / 3))
)

type Cell struct {
	id          CellID
	face        int8
	level       int8
	orientation int8
	uv          [2][2]float64
}

func CellFromCellID(id CellID) Cell {
	face, i, j, orientation := id.faceIJOrientation()
	cell := Cell{
		id:          id,
		face:        int8(face),
		level:       int8(id.Level()),
		orientation: int8(orientation),
	}
	ij := [2]int{i, j}
	cellSize := sizeIJ(int(cell.level))
	for d := 0; d < 2; d++ {
		lo := ij[d] & -cellSize
		hi := lo + cellSize
		cell.uv[d][0] = stToUV((1.0 / float64(maxSize)) * float64(lo))
		cell.uv[d][1] = stToUV((1.0 / float64(maxSize)) * float64(hi))
	}
	return cell
}

func (c Cell) CapBound() Cap {
	// Use th cell center in (u,v)-space as the cap axis. This vector
	// is very close to Center() and faster to compute. Neither one of
	// these vectors yields the bounding cap with minimal surface area,
	// but they are both pretty close.
	//
	// It's possible to show that the two vertices that are furthest from
	// the (u,v)-origin never determine the maximum cap size (this is a
	// possible future optimization).
	u := 0.5 * (c.uv[0][0] + c.uv[0][1])
	v := 0.5 * (c.uv[1][0] + c.uv[1][1])
	s2cap := CapFromCenterHeight(Point{faceUVToXYZ(int(c.face), u, v).Normalize()}, 0)
	for k := 0; k < 4; k++ {
		s2cap.AddPoint(c.Vertex(k))
	}
	return s2cap
}

func (c Cell) Subdivide(children *[4]Cell) bool {
	if c.id.IsLeaf() {
		return false
	}

	// Compute the cell midpoint in uv-space.
	u, v := c.id.centerUV()
	// Create four children with the appropriate bounds.
	id := c.id.childBegin()
	for pos := 0; pos < 4; pos, id = pos+1, id.next() {
		child := &children[pos]
		child.face = c.face
		child.level = c.level + 1
		child.orientation = c.orientation ^ int8(posToOrientation[pos])
		child.id = id
		// We want to split the cell in half in "u" and "v". To decide
		// which side to set equal to the midpoint value, we look at
		// cell's (i,j) position within its parent. The index for "i"
		// is in bit 1 of ij.
		ij := posToIJ[c.orientation][pos]
		i := ij >> 1
		j := ij & 1
		child.uv[0][i] = c.uv[0][i]
		child.uv[0][1-i] = u
		child.uv[1][j] = c.uv[1][j]
		child.uv[1][1-j] = v
	}
	return true
}

func AverageArea(level int) float64 {
	return AvgArea.Value(level)
}

func (c Cell) Latitude(i, j int) float64 {
	p := Point{faceUVToXYZ(int(c.face), c.uv[0][i], c.uv[1][j])}
	return latitude(p).Radians()
}

func (c Cell) Longitude(i, j int) float64 {
	p := Point{faceUVToXYZ(int(c.face), c.uv[0][i], c.uv[1][j])}
	return longitude(p).Radians()
}

func (c Cell) VertexRaw(k int) Point {
	// Vertices are returned in the order SW, SE, NE, NW.
	return Point{faceUVToXYZ(int(c.face), c.uv[0][(k>>1)^(k&1)], c.uv[1][k>>1])}
}

func (c Cell) EdgeRaw(k int) Point {
	face := int(c.face)
	switch k {
	case 0:
		return Point{vNorm(face, c.uv[1][0])} // South
	case 1:
		return Point{uNorm(face, c.uv[0][1])} // East
	case 2:
		return Point{vNorm(face, c.uv[1][1]).Neg()} // North
	default:
		return Point{uNorm(face, c.uv[0][0]).Neg()} // West
	}
}

func (c Cell) Vertex(k int) Point {
	return Point{c.VertexRaw(k).Normalize()}
}

func (c Cell) CenterRaw() Point {
	return Point{c.id.rawPoint()}
}

func (c Cell) Center() Point {
	return Point{c.CenterRaw().Normalize()}
}

func (c Cell) MayIntersect(cell Cell) bool {
	return c.id.Intersects(cell.id)
}

func (c Cell) ContainsCell(cell Cell) bool {
	return c.id.Contains(cell.id)
}

func (c Cell) ContainsPoint(p Point) bool {
	// We can't just call xyzToFaceUV, because for points that lie on the
	// boundary between two faces (i.e. u or v is +1/-1) we need to return
	// true for both adjacent cells.
	var u, v float64
	if !faceXYZToUV(int(c.face), p, &u, &v) {
		return false
	}
	return u >= c.uv[0][0] && u <= c.uv[0][1] &&
		v >= c.uv[1][0] && v <= c.uv[1][1]
}

func (c Cell) RectBound() Rect {
	if c.level > 0 {
		// Except for cells at level 0, the latitude and longitude
		// extremes are attained at the vertices. Furthermore, the
		// latitude range is determined by one pair of diagonally
		// opposite vertices and the longitude range is determined by
		// the other pair.
		//
		// We first determine which corner (i,j) of the cell has the
		// largest absolute latitude. To maximize latitude, we want to
		// find the point in the cell that has the largest absolute
		// z-coordinate and the smallest absolute x- and y-coordinates.
		// To do this we look at each coordinate (u and v), and
		// determine whether we want to minimize or maximize that
		// coordinate based on the axis direction and the cell's (u,v)
		// quadrant.
		u := c.uv[0][0] + c.uv[0][1]
		v := c.uv[1][0] + c.uv[1][1]
		i, j := ijFromFaceZ(c.face, u, v)

		// We grow the bounds slightly to make sure that the bounding
		// rectangle also contains the normalized versions of the
		// vertices. Note that the maximum result magnitude is Pi, with
		// a floating-point exponent of 1. Therefore adding or
		// subtracting 2**-51 will always change the result.
		lat := r1.IntervalFromPointPair(c.Latitude(i, j), c.Latitude(1-i, 1-j))
		lat = lat.Expanded(maxError).Intersection(validRectLatRange)
		if lat.Lo == -M_PI_2 || lat.Hi == M_PI_2 {
			return Rect{lat, s1.FullInterval()}
		}

		lng := s1.IntervalFromPointPair(c.Longitude(i, 1-j), c.Longitude(1-i, j))
		return Rect{lat, lng.Expanded(maxError)}
	}

	// The 4 cells around the equator extend to +/-45 degrees latitude at
	// the midpoints of their top and bottom edges. The two cells covering
	// the poles extend down to +/-35.26 degrees at their vertices.
	//
	// The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
	switch c.face {
	case 0:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{-M_PI_4, M_PI_4},
		}
	case 1:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{M_PI_4, 3 * M_PI_4},
		}
	case 2:
		return Rect{
			r1.Interval{poleMinLat, M_PI_2},
			s1.Interval{-math.Pi, math.Pi},
		}
	case 3:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{3 * M_PI_4, -3 * M_PI_4},
		}
	case 4:
		return Rect{
			r1.Interval{-M_PI_4, M_PI_4},
			s1.Interval{-3 * M_PI_4, -M_PI_4},
		}
	default:
		return Rect{
			r1.Interval{-M_PI_2, -poleMinLat},
			s1.Interval{-math.Pi, math.Pi},
		}
	}
}

func ijFromFaceZ(face int8, u, v float64) (i, j int) {
	switch uAxis(int(face)).Z {
	case 0:
		i = bool2int(u < 0)
	default:
		i = bool2int(u > 0)
	}
	switch vAxis(int(face)).Z {
	case 0:
		j = bool2int(v < 0)
	default:
		j = bool2int(v > 0)
	}
	return
}

func bool2int(b bool) int {
	if b {
		return 1
	}
	return 0
}