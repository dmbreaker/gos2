package s2

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
	cellSize := sizeIJ(id.Level())
	for d := 0; d < 2; d++ {
		lo := ij[d] & -cellSize
		hi := lo + cellSize
		cell.uv[d][0] = stToUV((1.0 / float64(maxSize)) * float64(lo))
		cell.uv[d][1] = stToUV((1.0 / float64(maxSize)) * float64(hi))
	}
	return cell
}

func (c Cell) VertexRaw(k int) Point {
	return Point{faceUVToXYZ(int(c.face), c.uv[0][(k>>1)^(k&1)], c.uv[1][k>>1])}
}

func (c Cell) Vertex(k int) Point {
	return Point{c.VertexRaw(k).Normalize()}
}
