package main

import (
	"fmt"
	"math"
	"strings"
)

func main() {
	xiA := []float64{0.1, 0.5, 0.9, 1.3}
	yiA := make([]float64, len(xiA))
	for i := range xiA {
		yiA[i] = 1 / xiA[i]
	}

	xiB := []float64{0.1, 0.5, 1.1, 1.3}
	yiB := make([]float64, len(xiB))
	for i := range xiB {
		yiB[i] = 1 / xiB[i]
	}

	xStar := 0.8
	trueValue := 1 / xStar

	fmt.Println("\n==============================================")
	fmt.Println("ВАРИАНТ A) Узлы: 0.1, 0.5, 0.9, 1.3")
	fmt.Printf("Точка интерполяции X* = %.1f\n", xStar)
	processInterpolation(xiA, yiA, xStar, trueValue)

	fmt.Println("\n==============================================")
	fmt.Println("ВАРИАНТ Б) Узлы: 0.1, 0.5, 1.1, 1.3")
	fmt.Printf("Точка интерполяции X* = %.1f\n", xStar)
	processInterpolation(xiB, yiB, xStar, trueValue)
}

func processInterpolation(xi, yi []float64, xStar, trueValue float64) {
	lagrangeResult := lagrangeInterpolation(xStar, xi, yi)
	lagrangeError := math.Abs(lagrangeResult - trueValue)

	printLagrangeTable(xi, yi, xStar)
	fmt.Println("\nМетод Лагранжа")
	fmt.Println(formatLagrangePolynomial(xi, yi))
	fmt.Printf("L(%.1f) = %.6f\n", xStar, lagrangeResult)
	fmt.Printf("Точное значение f(%.1f) = %.6f\n", xStar, trueValue)
	fmt.Printf("Абсолютная погрешность: %.6f\n", lagrangeError)

	newtonResult, dividedDiffs := newtonInterpolation(xStar, xi, yi)
	newtonError := math.Abs(newtonResult - trueValue)

	printDividedDifferences(xi, yi, dividedDiffs)
	fmt.Println("\nМетод Ньютона:")
	fmt.Println(formatNewtonPolynomial(xi, dividedDiffs))
	fmt.Printf("P(%.1f) = %.6f\n", xStar, newtonResult)
	fmt.Printf("Точное значение f(%.1f) = %.6f\n", xStar, trueValue)
	fmt.Printf("Абсолютная погрешность: %.6f\n", newtonError)
}

func lagrangeInterpolation(x float64, xi, yi []float64) float64 {
	var result float64
	n := len(xi)
	for i := 0; i < n; i++ {
		term := yi[i]
		for j := 0; j < n; j++ {
			if j != i {
				term *= (x - xi[j]) / (xi[i] - xi[j])
			}
		}
		result += term
	}
	return result
}

func newtonInterpolation(x float64, xi, yi []float64) (float64, [][]float64) {
	n := len(xi)

	f := make([][]float64, n)
	for i := range f {
		f[i] = make([]float64, n)
		f[i][0] = yi[i]
	}

	for j := 1; j < n; j++ {
		for i := 0; i < n-j; i++ {
			f[i][j] = (f[i+1][j-1] - f[i][j-1]) / (xi[i+j] - xi[i])
		}
	}

	result := f[0][0]
	product := 1.0
	for j := 1; j < n; j++ {
		product *= (x - xi[j-1])
		result += product * f[0][j]
	}

	return result, f
}

func printLagrangeTable(xi, yi []float64, xStar float64) {
	fmt.Println("\nМетод Лагранжа:")
	fmt.Println(" i    x_i      f_i    ω₄'(x_i)  f_i/ω₄'(x_i)  X*-x_i ")

	// Вычисляем ω₄'(x_i)
	omegaPrimes := make([]float64, len(xi))
	for i := range xi {
		prod := 1.0
		for j := range xi {
			if j != i {
				prod *= (xi[i] - xi[j])
			}
		}
		omegaPrimes[i] = prod
	}

	for i := range xi {
		fmt.Printf(" %d  %6.4f  %6.4f  %8.4f  %12.6f  %6.4f \n",
			i, xi[i], yi[i], omegaPrimes[i], yi[i]/omegaPrimes[i], xStar-xi[i])
	}
}

func printDividedDifferences(xi, yi []float64, f [][]float64) {
	fmt.Println("\nМетод Ньютона (таблица разделённых разностей):")
	fmt.Println(" i    x_i      f[x_i]    f[x_i,x_i+1]  f[x_i,x_i+1,x_i+2]  f[x0,x1,x2,x3] ")

	n := len(xi)
	for i := 0; i < n; i++ {
		fmt.Printf(" %d  %6.4f  %10.6f ", i, xi[i], f[i][0])

		if i < n-1 {
			fmt.Printf(" %12.6f ", f[i][1])
		} else {
			fmt.Printf("              ")
		}

		if i < n-2 {
			fmt.Printf(" %17.6f ", f[i][2])
		} else {
			fmt.Printf("                   ")
		}

		if i == 0 && n >= 4 {
			fmt.Printf(" %14.6f \n", f[0][3])
		} else {
			fmt.Printf("               \n")
		}
	}
}

func formatLagrangePolynomial(xi, yi []float64) string {
	var buf strings.Builder
	buf.WriteString("L(x) = ")
	for i := range xi {
		if i > 0 {
			buf.WriteString(" + ")
		}
		buf.WriteString(fmt.Sprintf("%.4f", yi[i]))
		buf.WriteString(" * [")
		for j := range xi {
			if j != i {
				buf.WriteString(fmt.Sprintf("(x - %.1f)/%.4f", xi[j], xi[i]-xi[j]))
				if j < len(xi)-1 && (j != len(xi)-2 || i != len(xi)-1) {
					buf.WriteString(" * ")
				}
			}
		}
		buf.WriteString("]")
	}
	return buf.String()
}

func formatNewtonPolynomial(xi []float64, diffs [][]float64) string {
	var buf strings.Builder
	buf.WriteString("P(x) = ")
	buf.WriteString(fmt.Sprintf("%.6f", diffs[0][0])) // f[x0]

	for j := 1; j < len(diffs[0]); j++ {
		buf.WriteString(" + ")
		buf.WriteString(fmt.Sprintf("%.6f", diffs[0][j])) // f[x0,x1,...,xj]
		for k := 0; k < j; k++ {
			buf.WriteString(fmt.Sprintf("*(x - %.4f)", xi[k]))
		}
	}

	return buf.String()
}
