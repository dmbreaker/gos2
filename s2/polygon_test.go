package s2

import (
	"math/rand"
	"strings"
	"testing"
)

func makepolygon(s string) *Polygon {
	loopStrs := strings.Split(s, ";")
	loops := []*Loop{}
	for _, str := range loopStrs {
		if str != "" {
			loop := makeloop(str)
			loop.Normalize()
			loops = append(loops, loop)
		}
	}
	return NewPolygonFromLoops(&loops)
}

var (
	// A set of nested loops around the point 0:0 (lat:lng).
	// Every vertex of kNear0 is a vertex of kNear1.
	kNear0    = "-1:0, 0:1, 1:0, 0:-1;"
	kNear1    = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;"
	kNear2    = "5:-2, -2:5, -1:-2;"
	kNear3    = "6:-3, -3:6, -2:-2;"
	kNearHemi = "0:-90, -90:0, 0:90, 90:0;"

	// A set of nested loops around the point 0:180 (lat:lng).
	// Every vertex of kFar0 and kFar2 belongs to kFar1, and all
	// the loops except kFar2 are non-convex.
	kFar0    = "0:179, 1:180, 0:-179, 2:-180;"
	kFar1    = "0:179, -1:179, 1:180, -1:-179, 0:-179, 3:-178, 2:-180, 3:178;"
	kFar2    = "-1:-179, -1:179, 3:178, 3:-178;" // opposite direction
	kFar3    = "-3:-178, -2:179, -3:178, 4:177, 4:-177;"
	kFarHemi = "0:-90, 60:90, -60:90;"

	// A set of nested loops around the point -90:0 (lat:lng).
	kSouthPoint = "-89.9999:0.001"
	kSouth0a    = "-90:0, -89.99:0, -89.99:0.01;"
	kSouth0b    = "-90:0, -89.99:0.02, -89.99:0.03;"
	kSouth0c    = "-90:0, -89.99:0.04, -89.99:0.05;"
	kSouth1     = "-90:0, -89.9:-0.1, -89.9:0.1;"
	kSouth2     = "-90:0, -89.8:-0.2, -89.8:0.2;"
	kSouthHemi  = "0:-180, 0:60, 0:-60;"

	// Two different loops that surround all the Near and Far loops except
	// for the hemispheres.
	kNearFar1 = "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, 1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;"
	kNearFar2 = "-8:-4, 8:-4, 2:15, 2:170, 8:-175, -8:-175, -2:170, -2:15;"
	// Loops that result from intersection of other loops.
	kFarHSouthH = "0:-180, 0:90, -60:90, 0:-90;"
	near_H3210  = makepolygon(kNear0 + kNear2 + kNear3 + kNearHemi + kNear1)

	far_H3210          = makepolygon(kFar2 + kFarHemi + kFar0 + kFar1 + kFar3)
	south_0ab          = makepolygon(kSouth0a + kSouth0b)
	south_210b         = makepolygon(kSouth2 + kSouth0b + kSouth1)
	south_H20abc       = makepolygon(kSouth2 + kSouth0b + kSouthHemi + kSouth0a + kSouth0c)
	nf1_n10_f2_s10abc  = makepolygon(kSouth0c + kFar2 + kNear1 + kNearFar1 + kNear0 + kSouth1 + kSouth0b + kSouth0a)
	nf2_n2_f210_s210ab = makepolygon(kFar2 + kSouth0a + kFar1 + kSouth1 + kFar0 + kSouth0b + kNearFar2 + kSouth2 + kNear2)
	far_H              = makepolygon(kFarHemi)
	south_H            = makepolygon(kSouthHemi)
	far_H_south_H      = makepolygon(kFarHSouthH)
)

func TestSplitting(t *testing.T) {
	SplitAndAssemble(t, near_H3210)
	SplitAndAssemble(t, far_H3210)
	SplitAndAssemble(t, south_0ab)
	SplitAndAssemble(t, south_210b)
	SplitAndAssemble(t, south_H20abc)
	SplitAndAssemble(t, nf1_n10_f2_s10abc)
	SplitAndAssemble(t, nf2_n2_f210_s210ab)
	SplitAndAssemble(t, far_H)
	SplitAndAssemble(t, south_H)
	SplitAndAssemble(t, far_H_south_H)
}

func SplitAndAssemble(t *testing.T, polygon *Polygon) {
	builder := NewPolygonBuilder(DIRECTED_XOR())
	builder.AddPolygon(polygon)
	var expected Polygon
	if !builder.AssemblePolygon(&expected, nil) {
		// TODO: use Fatalf()
		t.Errorf("%v.AssemblePolygon() failed", builder)
	}
	for iter := 0; iter < 10; iter++ {
		coverer := NewRegionCoverer()
		diameter := 2 * polygon.CapBound().Radius().Radians()
		minLevel := MaxWidth.MinLevel(diameter)
		level := minLevel + rand.Intn(4)
		coverer.SetMinLevel(minLevel)
		coverer.SetMaxLevel(level)
		coverer.SetMaxCells(500)

		cells := coverer.Covering(*polygon)
		var covering CellUnion
		covering.Init(cells)
		CheckCompleteCovering(t, *polygon, covering, false, CellID(0))
	}
}

func CheckContains(t *testing.T, astr, bstr string) {
	a := makepolygon(astr)
	b := makepolygon(bstr)
	if !a.ContainsPolygon(b) {
		t.Errorf("!%v.Contains(%v)", a, b)
	}
}

func TestPolygonInit(t *testing.T) {
	CheckContains(t, kNear1, kNear0)
	CheckContains(t, kNear2, kNear1)
	CheckContains(t, kNear3, kNear2)
	CheckContains(t, kNearHemi, kNear3)
}
