package tree

type ScalarClusterable struct {
	x     float64
	x2    float64
	count float64
}

func (scalar *ScalarClusterable) Mean() float64 {
	if scalar.count != 0.0 {
		return scalar.x / scalar.count
	}

	return 0.0
}

func NewScalarClusterable(x float64) *ScalarClusterable {
	return &ScalarClusterable{
		x:     x,
		x2:    x * x,
		count: 1,
	}
}

func (scalar *ScalarClusterable) SetZero() {
	scalar.count = 0.0
	scalar.x = 0.0
	scalar.x2 = 0.0
}

func (scalar *ScalarClusterable) Objf() float64 {
	if scalar.count == 0.0 {
		return 0.0
	}

	if scalar.count < 0 {
		panic("count doit etre > 0")
	}

	return -(scalar.x2 - scalar.x*scalar.x/scalar.count)
}

func (scalar *ScalarClusterable) Add(other ScalarClusterable) {
	scalar.x += other.x
	scalar.x2 += other.x2
	scalar.count += other.count
}

func (scalar *ScalarClusterable) Sub(other ScalarClusterable) {
	scalar.x -= other.x
	scalar.x2 -= other.x2
	scalar.count -= other.count
}

func (scalar *ScalarClusterable) Copy() *ScalarClusterable {
	newScalar := NewScalarClusterable(0)
	newScalar.Add(*scalar)

	return newScalar
}
