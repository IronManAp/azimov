package tree

import (
	"math"

	"gonum.org/v1/gonum/blas"
	"gonum.org/v1/gonum/mat"
)

type VectorClusterable struct {
	weight float64
	stats  mat.VecDense
	sumsq  float64
}

func NewVectorClusterable(vector mat.VecDense, weights float64) *VectorClusterable {
	var vecTmp mat.VecDense

	vecTmp.CopyVec(&vector)
	vecTmp.ScaleVec(weights, &vector)

	return &VectorClusterable{
		weight: weights,
		stats:  vecTmp,
		sumsq:  blas.Float64Level1.Ddot(nil, vector.RawVector().N, vector.RawVector().Data, 1, vector.RawVector().Data, 1) * weights,
	}
}

func (vec *VectorClusterable) Objf() float64 {
	var direct_sumsq float64

	if vec.weight > math.SmallestNonzeroFloat64 {
		direct_sumsq = blas.Float64Level1.Ddot(nil, vec.stats.Len(), vec.stats.RawVector().Data, 1, vec.stats.RawVector().Data, 1) / vec.weight
	} else {
		direct_sumsq = 0.0
	}

	ans := -(vec.sumsq - direct_sumsq)

	if ans > 0.0 {
		if ans > 1 {
			panic("Positive obective function encountered")
		}

		ans = 0.0
	}

	return ans
}

func (vec *VectorClusterable) Scale(f float64) {
	if f < 0.0 {
		panic("f doit etre > 0.0")
	}

	vec.weight *= f
	vec.stats.ScaleVec(f, &vec.stats)
	vec.sumsq *= f
}

func (vec *VectorClusterable) Copy() *VectorClusterable {
	return &VectorClusterable{
		weight: vec.weight,
		stats:  vec.stats,
		sumsq:  vec.sumsq,
	}
}

func (vec *VectorClusterable) Add(other VectorClusterable) {
	vec.weight += other.weight
	vec.stats.AddVec(&vec.stats, &other.stats)
	vec.sumsq += other.sumsq
}

func (vec *VectorClusterable) Sub(other VectorClusterable) {
	vec.weight -= other.weight
	vec.stats.SubVec(&vec.stats, &other.stats)
	vec.sumsq -= other.sumsq

	if vec.weight < 0.0 {
		if vec.weight < -0.1 && vec.weight < -0.0001*math.Abs(other.weight) {
			panic("Weight negative")
		}

		vec.weight = 0.0
	}

	if vec.weight == 0.0 {
		vec.sumsq = 0.0
		vec.stats.Zero()
	}
}
