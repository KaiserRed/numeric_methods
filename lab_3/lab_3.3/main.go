package main

import (
	"fmt"
	"math"
)

func main() {
	xi := []float64{0.1, 0.5, 0.9, 1.3, 1.7, 2.1}
	yi := []float64{10.0, 2.0, 1.1111, 0.76923, 0.58824, 0.47618}
	n := len(xi)

	fmt.Println("Исходные данные:")
	fmt.Println("i\tx_i\ty_i")
	for i := 0; i < n; i++ {
		fmt.Printf("%d\t%.1f\t%.5f\n", i, xi[i], yi[i])
	}
	fmt.Println()

	fmt.Println("==============================================")
	fmt.Println("Линейная аппроксимация (1-я степень):")
	fmt.Println("Ищем многочлен P1(x) = a + b*x")
	linearCoeffs, linearErr := leastSquares(xi, yi, 1)
	fmt.Printf("\nПриближающий многочлен: P1(x) = %.4f + %.4f*x\n", linearCoeffs[0], linearCoeffs[1])
	fmt.Printf("Сумма квадратов ошибок: %.6f\n\n", linearErr)

	fmt.Println("==============================================")
	fmt.Println("Квадратичная аппроксимация (2-я степень):")
	fmt.Println("Ищем многочлен P2(x) = a + b*x + c*x^2")
	quadraticCoeffs, quadraticErr := leastSquares(xi, yi, 2)
	fmt.Printf("\nПриближающий многочлен: P2(x) = %.4f + %.4f*x + %.4f*x^2\n",
		quadraticCoeffs[0], quadraticCoeffs[1], quadraticCoeffs[2])
	fmt.Printf("Сумма квадратов ошибок: %.6f\n", quadraticErr)

	fmt.Println("\n==============================================")
	fmt.Println("Таблица значений аппроксимирующих функций:")
	fmt.Println("x\t   y\t    P1(x)\t   P2(x)\t  ΔP1\t  ΔP2")

	for i := 0; i < n; i++ {
		x := xi[i]
		y := yi[i]
		p1 := linearCoeffs[0] + linearCoeffs[1]*x
		p2 := quadraticCoeffs[0] + quadraticCoeffs[1]*x + quadraticCoeffs[2]*x*x
		delta1 := math.Abs(y - p1)
		delta2 := math.Abs(y - p2)
		fmt.Printf("%.2f  %.5f  %.5f  %.5f  %.5f  %.5f\n", x, y, p1, p2, delta1, delta2)
	}
}

func leastSquares(x, y []float64, degree int) ([]float64, float64) {
	n := len(x)
	m := degree + 1

	A := make([][]float64, m)
	for i := range A {
		A[i] = make([]float64, m)
	}
	b := make([]float64, m)

	for k := 0; k < m; k++ {
		for j := 0; j < m; j++ {
			for i := 0; i < n; i++ {
				A[k][j] += math.Pow(x[i], float64(k+j))
			}
		}
		for i := 0; i < n; i++ {
			b[k] += y[i] * math.Pow(x[i], float64(k))
		}
	}

	fmt.Println("\nНормальная система уравнений:")
	switch degree {
	case 1:
		fmt.Printf("%.2fa + %.2fb = %.5f\n", A[0][0], A[0][1], b[0])
		fmt.Printf("%.2fa + %.2fb = %.5f\n", A[1][0], A[1][1], b[1])
	case 2:
		fmt.Printf("%.2fa + %.2fb + %.2fc = %.5f\n", A[0][0], A[0][1], A[0][2], b[0])
		fmt.Printf("%.2fa + %.2fb + %.2fc = %.5f\n", A[1][0], A[1][1], A[1][2], b[1])
		fmt.Printf("%.2fa + %.2fb + %.2fc = %.5f\n", A[2][0], A[2][1], A[2][2], b[2])
	}

	coeffs := gauss(A, b)

	fmt.Println("\nРешение системы:")
	switch degree {
	case 1:
		fmt.Printf("a = %.6f\n", coeffs[0])
		fmt.Printf("b = %.6f\n", coeffs[1])
	case 2:
		fmt.Printf("a = %.6f\n", coeffs[0])
		fmt.Printf("b = %.6f\n", coeffs[1])
		fmt.Printf("c = %.6f\n", coeffs[2])
	}

	err := 0.0
	for i := 0; i < n; i++ {
		pred := 0.0
		for j := 0; j < m; j++ {
			pred += coeffs[j] * math.Pow(x[i], float64(j))
		}
		err += math.Pow(y[i]-pred, 2)
	}

	return coeffs, err
}

func gauss(A [][]float64, b []float64) []float64 {
	n := len(b)
	x := make([]float64, n)

	for k := 0; k < n-1; k++ {
		maxRow := k
		maxVal := math.Abs(A[k][k])
		for i := k + 1; i < n; i++ {
			if math.Abs(A[i][k]) > maxVal {
				maxVal = math.Abs(A[i][k])
				maxRow = i
			}
		}

		if maxRow != k {
			A[k], A[maxRow] = A[maxRow], A[k]
			b[k], b[maxRow] = b[maxRow], b[k]
		}

		for i := k + 1; i < n; i++ {
			factor := A[i][k] / A[k][k]
			for j := k; j < n; j++ {
				A[i][j] -= factor * A[k][j]
			}
			b[i] -= factor * b[k]
		}
	}

	// Обратный ход
	for i := n - 1; i >= 0; i-- {
		sum := 0.0
		for j := i + 1; j < n; j++ {
			sum += A[i][j] * x[j]
		}
		x[i] = (b[i] - sum) / A[i][i]
	}

	return x
}
