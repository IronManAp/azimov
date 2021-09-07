package gmm

import (
	"gonum.org/v1/gonum/mat"
)

type DiagGmmNormal struct {
	weights mat.Matrix
	means   mat.Matrix
	vars    mat.Matrix
}

func (gmmN *DiagGmmNormal) CopyToDiagGmm(diagGmm DiagGmm) {

}
