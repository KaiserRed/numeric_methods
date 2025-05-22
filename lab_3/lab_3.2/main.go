package main

import (
	"fmt"
	"math"
)

type Spline struct {
	xi []float64
	fi []float64
	ai []float64
	bi []float64
	ci []float64
	di []float64
}

func main() {
	xi := []float64{0.1, 0.5, 0.9, 1.3, 1.7}
	fi := []float64{10.0, 2.0, 1.1111, 0.76923, 0.58824}
	xStar := 0.8

	spline := BuildNaturalSpline(xi, fi)

	splineValue := spline.Evaluate(xStar)
	trueValue := 1 / xStar
	error := math.Abs(splineValue - trueValue)

	PrintSplineInfo(spline, xi)
	fmt.Printf("\nТочка интерполяции X* = %.1f\n", xStar)
	fmt.Printf("Значение сплайна: %.5f\n", splineValue)
	fmt.Printf("Точное значение: %.5f\n", trueValue)
	fmt.Printf("Погрешность: %.5f\n", error)
	fmt.Println("Проверки")
	if !spline.CheckNodeValues() {
		fmt.Println("Ошибка: сплайн не проходит через узлы")
	}
	fmt.Println("Значения в узлах: Пройдена")

	if !spline.CheckBoundaryConditions() {
		fmt.Println("Ошибка: нарушены граничные условия")
	}
	fmt.Println("Гранияный условий: Пройдена")

}

func BuildNaturalSpline(xi, fi []float64) Spline {
	n := len(xi)
	if n < 3 {
		panic("Для построения сплайна нужно минимум 3 точки")
	}

	h := calculateSteps(xi)
	A, b := buildTridiagonalSystem(xi, fi, h)
	ci := solveTridiagonal(A, b)
	ci = append([]float64{0}, ci...)
	ci = append(ci, 0)

	ai, bi, di := calculateSplineCoefficients(xi, fi, h, ci)

	return Spline{xi, fi, ai, bi, ci, di}
}

func calculateSteps(xi []float64) []float64 {
	h := make([]float64, len(xi)-1)
	for i := 0; i < len(h); i++ {
		h[i] = xi[i+1] - xi[i]
	}
	return h
}

func buildTridiagonalSystem(xi, fi, h []float64) ([][]float64, []float64) {
	n := len(xi)
	A := make([][]float64, n-2)
	for i := range A {
		A[i] = make([]float64, n-2)
	}

	A[0][0] = 2 * (h[0] + h[1])
	A[0][1] = h[1]
	for i := 1; i < n-3; i++ {
		A[i][i-1] = h[i]
		A[i][i] = 2 * (h[i] + h[i+1])
		A[i][i+1] = h[i+1]
	}
	if n > 3 {
		A[n-3][n-4] = h[n-3]
		A[n-3][n-3] = 2 * (h[n-3] + h[n-2])
	}

	b := make([]float64, n-2)
	for i := 0; i < n-2; i++ {
		b[i] = 3 * ((fi[i+2]-fi[i+1])/h[i+1] - (fi[i+1]-fi[i])/h[i])
	}

	return A, b
}

func solveTridiagonal(A [][]float64, b []float64) []float64 {
	n := len(b)
	if n == 0 {
		return nil
	}

	alpha := make([]float64, n)
	beta := make([]float64, n)

	alpha[0] = -A[0][1] / A[0][0]
	beta[0] = b[0] / A[0][0]

	for i := 1; i < n-1; i++ {
		denominator := A[i][i] + A[i][i-1]*alpha[i-1]
		alpha[i] = -A[i][i+1] / denominator
		beta[i] = (b[i] - A[i][i-1]*beta[i-1]) / denominator
	}

	x := make([]float64, n)
	x[n-1] = (b[n-1] - A[n-1][n-2]*beta[n-2]) /
		(A[n-1][n-1] + A[n-1][n-2]*alpha[n-2])

	for i := n - 2; i >= 0; i-- {
		x[i] = alpha[i]*x[i+1] + beta[i]
	}

	return x
}

func calculateSplineCoefficients(xi, fi, h, ci []float64) ([]float64, []float64, []float64) {
	n := len(xi)
	ai := make([]float64, n-1)
	bi := make([]float64, n-1)
	di := make([]float64, n-1)

	for i := 0; i < n-1; i++ {
		ai[i] = fi[i]
		if i < n-2 {
			bi[i] = (fi[i+1]-fi[i])/h[i] - h[i]*(2*ci[i]+ci[i+1])/3
			di[i] = (ci[i+1] - ci[i]) / (3 * h[i])
		} else {
			bi[i] = (fi[i+1]-fi[i])/h[i] - h[i]*2*ci[i]/3
			di[i] = -ci[i] / (3 * h[i])
		}
	}

	return ai, bi, di
}

func (s *Spline) Evaluate(x float64) float64 {
	interval := s.findInterval(x)
	dx := x - s.xi[interval]
	return s.ai[interval] + s.bi[interval]*dx + s.ci[interval]*dx*dx + s.di[interval]*dx*dx*dx
}

func (s *Spline) findInterval(x float64) int {
	for i := 0; i < len(s.xi)-1; i++ {
		if x >= s.xi[i] && x <= s.xi[i+1] {
			return i
		}
	}
	return len(s.xi) - 2
}

func PrintSplineInfo(s Spline, xi []float64) {
	fmt.Println("Построение кубического сплайна с естественными граничными условиями")
	fmt.Println("Узлы интерполяции:")
	fmt.Println("i\tx_i\t\tf_i")
	for i := range xi {
		fmt.Printf("%d\t%.4f\t%.5f\n", i, xi[i], s.fi[i])
	}

	fmt.Println("\nКоэффициенты сплайна:")
	fmt.Println("i\t[x_i-1, x_i]\ta_i\t\tb_i\t\tc_i\t\td_i")
	for i := 1; i < len(s.ai)+1; i++ {
		if i-1 < len(s.ai) {
			fmt.Printf("%d  [%.1f,%.1f]  %.5f  %.5f  %.5f  %.5f\n",
				i, xi[i-1], xi[i], s.ai[i-1], s.bi[i-1], s.ci[i-1], s.di[i-1])
		}
	}
}

func (s *Spline) CheckNodeValues() bool {
	for i, x := range s.xi {
		// Вычисляем значение сплайна в узле x_i
		val := s.Evaluate(x)
		// Сравниваем с заданным значением f_i
		if math.Abs(val-s.fi[i]) > 1e-9 {
			fmt.Printf("Ошибка в узле %d: x=%.2f, сплайн=%.6f, f_i=%.6f\n",
				i, x, val, s.fi[i])
			return false
		}
	}
	return true
}

func (s *Spline) CheckBoundaryConditions() bool {
	n := len(s.xi)
	// Вторая производная в x_0 должна быть 0
	d2_0 := 2*s.ci[0] + 6*s.di[0]*(s.xi[0]-s.xi[0])
	// Вторая производная в x_n должна быть 0
	d2_n := 2*s.ci[n-2] + 6*s.di[n-2]*(s.xi[n-1]-s.xi[n-2])

	if math.Abs(d2_0) > 1e-9 || math.Abs(d2_n) > 1e-9 {
		fmt.Printf("Нарушены граничные условия: d2_0=%.6f, d2_n=%.6f\n",
			d2_0, d2_n)
		return false
	}
	return true
}
