package s2

import (
	"math"
)

type LengthMetric struct {
	Deriv float64
	dim   int
}

func NewLengthMetric(deriv float64) LengthMetric {
	return LengthMetric{
		Deriv: deriv,
		dim:   1,
	}
}

// Return the length on the unit sphere for cells at the given level.
func (m LengthMetric) Value(level int) float64 {
	return math.Ldexp(m.Deriv, -m.dim*level)
}

// Return the level at which the length is approximately the given value.
// The return value is always a valid level.
func (m LengthMetric) ClosestLevel(value float64) int {
	return m.MinLevel(math.Sqrt2 * value)
}

// Return the minimum level such that the length is at most the given value,
// or maxLevel if there is no such level. The return value is always a valid
// level.
func (m LengthMetric) MinLevel(value float64) int {
	if value <= 0 {
		return MaxCellLevel
	}
	// This code is equivalent to computing a floating-point "level"
	// value and rounding up. Frexp() returns a fraction in the range
	// [0.5, 1) and the corresponding exponent.
	_, level := math.Frexp(value / m.Deriv)
	level = int(math.Max(0, math.Min(MaxCellLevel, float64(-((level-1)>>uint(m.dim-1))))))
	// TODO: sanity check level
	return level
}

// Return the maximum level such that the length is at least the given value,
// or zero if there is no such level. The return value is always a valid level.
func (m LengthMetric) MaxLevel(value float64) int {
	if value <= 0 {
		return MaxCellLevel
	}
	_, level := math.Frexp(m.Deriv / value)
	level = int(math.Max(0, math.Min(MaxCellLevel, float64((level-1)>>uint(m.dim-1)))))
	// TODO: sanity check level
	return level
}

var (
	MinAngleSpan LengthMetric
	MaxAngleSpan LengthMetric
	AvgAngleSpan LengthMetric
	MinWidth     LengthMetric
	MaxWidth     LengthMetric
	AvgWidth     LengthMetric
	MinEdge      LengthMetric
	MaxEdge      LengthMetric
	AvgEdge      LengthMetric
	MinDiag      LengthMetric
	MaxDiag      LengthMetric
	AvgDiag      LengthMetric
)

// All of the values below were obtained by a combination of hand analysis and
// Mathematica.
// TODO: implement tangent and quadratic projections.
func init() {
	MinAngleSpan = NewLengthMetric(1)
	MaxAngleSpan = NewLengthMetric(2)
	AvgAngleSpan = NewLengthMetric(math.Pi / 2)      // 1.571 (true for all projections)
	MinWidth = NewLengthMetric(math.Sqrt(2. / 3))    // 0.816
	MaxWidth = NewLengthMetric(MaxAngleSpan.Deriv)   // (true for all projections)
	AvgWidth = NewLengthMetric(1.411459345844456965) // 1.411
	MinEdge = NewLengthMetric(2 * math.Sqrt(2) / 3)  // 0.943
	MaxEdge = NewLengthMetric(MaxAngleSpan.Deriv)    // (true for all projections)
	AvgEdge = NewLengthMetric(1.440034192955603643)  // 1.440
	MinDiag = NewLengthMetric(2 * math.Sqrt(2) / 3)  // 0.943
	MaxDiag = NewLengthMetric(2 * math.Sqrt(2))      // 2.828
	AvgDiag = NewLengthMetric(2.031817866418812674)  // 2.032
}
