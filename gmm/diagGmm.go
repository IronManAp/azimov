package gmm

import (
	"math"

	"azimov.org/azimov/math"
	"azimov.org/azimov/tree"
	"gonum.org/v1/gonum/mat"
)

type DiagGmm struct {
	gConsts       mat.VecDense
	valid_gConsts bool
	weights       mat.VecDense
	inv_vars      mat.Dense
	means_invvars mat.Dense
}

func NewDiagGmm(gc tree.GaussClusterable, var_floor float64) *DiagGmm {
	x := gc.X_stats()
	x2 := gc.X2_stats()
	count := gc.Count()

	if count <= 0.0 {
		panic("Count doit etre > 0.0")
	}

	var dg *DiagGmm

	dg.Resize(1, x.RawVector().N)

	x.ScaleVec(1.0/count, &x)
	x2.ScaleVec(1.0/count, &x2)
	x2.AddScaledVec(&x2, -1.0, &x)
	x2.ApplyFloor(var_floor)
	x2.InvertElement()

	if x2.Min() < 0 {
		panic("x2 a des valeurs negatives.")
	}

	mean := mat.NewDense(1, x.RawVector().N, nil)
	mean.SetRow(0, x.RawVector().Data)

	inv_var := mat.NewDense(1, x.RawVector().N, nil)
	inv_var.SetRow(0, x2.RawVector().Data)
	dg.SetInvVarsAndMean(*inv_var, *mean)

	weights := mat.NewVecDense(1, nil)
	weights.SetVec(0, 1.0)

	dg.SetWeights(*weights)

	return dg
}

func (dgmm *DiagGmm) Resize(nmix, dim int) {
	if nmix <= 0 && dim <= 0 {
		panic("nmix && dim doivent etre > 0")
	}

	if dgmm.gConsts.RawVector().N != nmix {
		dgmm.gConsts.ReuseAsVec(nmix)
	}

	if dgmm.weights.RawVector().Inc != nmix {
		dgmm.weights.ReuseAsVec(nmix)
	}

	if dgmm.inv_vars.RawMatrix().Rows != nmix || dgmm.inv_vars.RawMatrix().Cols != nmix {
		dgmm.inv_vars.ReuseAs(nmix, dim)
		dgmm.inv_vars.SetData(1.0)
	}

	if dgmm.means_invvars.RawMatrix().Rows != nmix || dgmm.means_invvars.RawMatrix().Cols != dim {
		dgmm.means_invvars.ReuseAs(nmix, dim)
	}

	dgmm.valid_gConsts = false
}

func (dgmm *DiagGmm) SetInvVarsAndMean(invvars, means mat.Dense) {
	dgmm.inv_vars.Copy(&invvars)
	new_means_invvars := mat.NewDense(means.RawMatrix().Rows, means.RawMatrix().Cols, means.RawMatrix().Data)
	new_means_invvars.MulElem(new_means_invvars, &invvars)
	dgmm.means_invvars.Copy(new_means_invvars)
	dgmm.valid_gConsts = false
}

func (dgmm *DiagGmm) SetWeights(weights mat.VecDense) {
	dgmm.weights.CopyVec(&weights)
	dgmm.valid_gConsts = false
}

func (dgmm *DiagGmm) ComputeGconsts() int {
	num_mix := dgmm.NumGauss()
	dim := dgmm.Dim()
	offset := -0.5 * tree.LOG_2PI * float64(dim)
	num_bad := 0

	if num_mix != dgmm.gConsts.Len() {
		dgmm.gConsts.ReuseAsVec(num_mix)
	}

	for mix := 0; mix < num_mix; mix++ {
		gc := math.Log(dgmm.weights.At(mix, 0)) + offset
		for d := 0; d < dim; d++ {
			gc += 0.5*math.Log(dgmm.inv_vars.At(mix, d)) - 0.5*dgmm.means_invvars.At(mix, d)*dgmm.means_invvars.At(mix, d)/dgmm.inv_vars.At(mix, d)
		}

		if math.IsNaN(gc) {
			panic("Un element n'est pas un nombre dans gconst computation")
		}

		if math.IsInf(gc, 1) {
			num_bad++

			if gc > 0 {
				gc = -gc
			}
		}

		dgmm.gConsts.RawVector().Data[mix] = gc
	}

	dgmm.valid_gConsts = true
	return num_bad
}

func (dgmm *DiagGmm) NumGauss() int {
	return dgmm.weights.Len()
}

func (dgmm *DiagGmm) Dim() int {
	return dgmm.means_invvars.RawMatrix().Cols
}

func (dgmm *DiagGmm) Generate(output *mat.VecDense) {
	tot := dgmm.weights.Sum()
	r := tot * mathSR.RandUniforme() * 0.99999
	i := 0
	sum := 0.0

	for sum+dgmm.weights.RawVector().Data[i] < r {
		sum += dgmm.weights.RawVector().Data[i]
		i++
	}

	var inv_var mat.Dense

}
