package s2

import (
	"math"
)

type Metric struct {
	Deriv float64
	dim   int
}

type LengthMetric struct {
	Metric
}

type AreaMetric struct {
	Metric
}

func NewMetric(deriv float64, dim int) Metric    { return Metric{deriv, dim} }
func NewLengthMetric(deriv float64) LengthMetric { return LengthMetric{Metric{deriv, 1}} }
func NewAreaMetric(deriv float64) AreaMetric     { return AreaMetric{Metric{deriv, 2}} }

// Return the length on the unit sphere for cells at the given level.
func (m Metric) Value(level int) float64 {
	return math.Ldexp(m.Deriv, -m.dim*level)
}

// Return the level at which the length is approximately the given value.
// The return value is always a valid level.
func (m Metric) ClosestLevel(value float64) int {
	var val float64
	if m.dim == 1 {
		val = math.Sqrt2 * value
	} else {
		val = 2 * value
	}
	return m.MinLevel(val)
}

// Return the minimum level such that the length is at most the given value,
// or maxLevel if there is no such level. The return value is always a valid
// level.
func (m Metric) MinLevel(value float64) int {
	if value <= 0 {
		return MaxCellLevel
	}
	// This code is equivalent to computing a floating-point "level"
	// value and rounding up. Frexp() returns a fraction in the range
	// [0.5, 1) and the corresponding exponent.
	_, level := math.Frexp(value / m.Deriv)
	level = max(0, min(MaxCellLevel, -((level-1)>>uint(m.dim-1))))
	return level
}

// Return the maximum level such that the length is at least the given value,
// or zero if there is no such level. The return value is always a valid level.
func (m Metric) MaxLevel(value float64) int {
	if value <= 0 {
		return MaxCellLevel
	}
	_, level := math.Frexp(m.Deriv / value)
	level = max(0, min(MaxCellLevel, (level-1)>>uint(m.dim-1)))
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
	AvgArea      AreaMetric
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
	MinEdge = NewLengthMetric(2 * math.Sqrt2 / 3)    // 0.943
	MaxEdge = NewLengthMetric(MaxAngleSpan.Deriv)    // (true for all projections)
	AvgEdge = NewLengthMetric(1.440034192955603643)  // 1.440
	MinDiag = NewLengthMetric(2 * math.Sqrt2 / 3)    // 0.943
	MaxDiag = NewLengthMetric(2 * math.Sqrt2)        // 2.828
	AvgDiag = NewLengthMetric(2.031817866418812674)  // 2.032
	AvgArea = NewAreaMetric(4 * math.Pi / 6)         // 2.094 (true for all projections)
}
