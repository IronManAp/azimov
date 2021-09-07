package mathSR

import "math/rand"

func RandUniforme() float64 {
	return (rand.Float64()*32767.0 + 1.0) / (32767.0 + 2.0)
}
