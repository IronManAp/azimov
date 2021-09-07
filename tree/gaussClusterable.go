package tree

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

const (
	LOG_2PI float64 = 1.8378770664093454835606594728112
)

type GaussClusterable struct {
	count     float64
	stats     mat.Dense
	var_floor float64
}

func NewGaussClusterable(_dim int, _var_floor float64) *GaussClusterable {
	v := make([]float64, 2*_dim)

	for i := 0; i < 2*_dim; i++ {
		v[i] = 0.0
	}

	return &GaussClusterable{
		count:     0.0,
		stats:     *mat.NewDense(2, _dim, v),
		var_floor: _var_floor,
	}
}

func (gauss *GaussClusterable) Objf() float64 {
	if gauss.count <= 0.0 {
		if gauss.count < -0.1 {
			panic("GaussClusterable::Objf(), count is negative")
		}
	} else {
		dim := gauss.stats.RawMatrix().Cols
		data := make([]float64, dim)

		for i := 0; i < dim; i++ {
			data[i] = 0
		}

		vars := *mat.NewVecDense(dim, data)
		objf_per_frame := 0.0

		for d := 0; d < dim; d++ {
			mean := (gauss.stats.At(0, d) / gauss.count)
			Var := gauss.stats.At(0, d)/gauss.count - mean*mean
			floored_var := math.Max(Var, gauss.var_floor)

			vars.SetVec(d, floored_var)
			objf_per_frame += -0.5 * Var / floored_var
		}

		objf_per_frame += -0.5 * (vars.SumLog() + LOG_2PI*float64(dim))

		if math.IsNaN(objf_per_frame) {
			panic("GaussClusterable::Objf(), objf is NaN")
		} else {
			return objf_per_frame * gauss.count
		}
	}

	return 0.0
}

func (gauss *GaussClusterable) Add(gauss2 GaussClusterable) {
	gauss.count += gauss2.count
	gauss.stats.Add(&gauss.stats, &gauss2.stats)
}

func (gauss *GaussClusterable) Sub(gauss2 GaussClusterable) {
	gauss.count -= gauss2.count
	gauss.stats.Sub(&gauss.stats, &gauss2.stats)
}

func (gauss *GaussClusterable) Copy() *GaussClusterable {
	newGauss := NewGaussClusterable(0, 0.0)
	newGauss.Add(*gauss)

	return newGauss
}

func (gauss *GaussClusterable) Scale(f float64) {
	if f < 0.0 {
		panic("f doit etre >= 0.0")
	}

	gauss.count *= f
	gauss.stats.Scale(f, &gauss.stats)
}

func (gauss *GaussClusterable) X_stats() mat.VecDense {
	return *mat.NewVecDense(gauss.stats.RawMatrix().Cols, gauss.stats.RawRowView(0))
}

func (gauss *GaussClusterable) X2_stats() mat.VecDense {
	return *mat.NewVecDense(gauss.stats.RawMatrix().Cols, gauss.stats.RawRowView(1))
}

func (gauss *GaussClusterable) Count() float64 {
	return gauss.count
}
