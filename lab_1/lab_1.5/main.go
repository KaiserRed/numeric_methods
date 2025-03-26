package main

import (
	"bufio"
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"os"
	"strconv"
	"strings"
)

type Matrix [][]float64

func main() {
	A, epsilon, err := readInput("input.txt")
	if err != nil {
		fmt.Printf("Ошибка чтения: %v\n", err)
		return
	}

	if len(A) != len(A[0]) {
		fmt.Println("Ошибка: матрица должна быть квадратной")
		return
	}

	eigenvalues, iterations := qrAlgorithm(A, epsilon)

	verificationErrors := verifyEigenvalues(A, eigenvalues)

	if err := writeResults("output.txt", A, eigenvalues, iterations, epsilon, verificationErrors); err != nil {
		fmt.Printf("Ошибка записи: %v\n", err)
		return
	}

	fmt.Println("QR-алгоритм успешно завершен. Результаты в output.txt")
}

func readInput(filename string) (Matrix, float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var matrix Matrix
	var epsilon float64
	lineNum := 0

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		if lineNum == 0 {
			epsilon, err = strconv.ParseFloat(line, 64)
			if err != nil {
				return nil, 0, fmt.Errorf("ошибка чтения точности: %v", err)
			}
		} else {
			fields := strings.Fields(line)
			row := make([]float64, len(fields))
			for i, f := range fields {
				val, err := strconv.ParseFloat(f, 64)
				if err != nil {
					return nil, 0, fmt.Errorf("ошибка в строке %d: %v", lineNum, err)
				}
				row[i] = val
			}
			matrix = append(matrix, row)
		}
		lineNum++
	}

	if lineNum < 2 {
		return nil, 0, fmt.Errorf("файл должен содержать минимум 2 строки (точность и матрица)")
	}

	return matrix, epsilon, nil
}

func qrAlgorithm(A Matrix, epsilon float64) ([]complex128, int) {
	//n := len(A)
	currentA := copyMatrix(A)
	iterations := 0

	for {
		Q, R := qrDecomposition(currentA)

		// Обновление матрицы: A = R*Q
		currentA = multiplyMatrices(R, Q)
		iterations++

		if isUpperTriangular(currentA, epsilon) {
			break
		}

		if iterations > 1000 {
			break
		}
	}

	eigenvalues := extractEigenvalues(currentA, epsilon)
	return eigenvalues, iterations
}

func qrDecomposition(A Matrix) (Q, R Matrix) {
	n := len(A)
	Q = make(Matrix, n)
	R = make(Matrix, n)
	for i := range Q {
		Q[i] = make([]float64, n)
		R[i] = make([]float64, n)
	}

	for j := 0; j < n; j++ {
		v := make([]float64, n)
		for i := 0; i < n; i++ {
			v[i] = A[i][j]
		}

		for k := 0; k < j; k++ {
			dot := 0.0
			for i := 0; i < n; i++ {
				dot += v[i] * Q[i][k]
			}
			R[k][j] = dot

			for i := 0; i < n; i++ {
				v[i] -= dot * Q[i][k]
			}
		}

		norm := 0.0
		for i := 0; i < n; i++ {
			norm += v[i] * v[i]
		}
		norm = math.Sqrt(norm)

		for i := 0; i < n; i++ {
			Q[i][j] = v[i] / norm
		}
		R[j][j] = norm
	}

	return Q, R
}

func isUpperTriangular(A Matrix, epsilon float64) bool {
	n := len(A)
	for i := 1; i < n; i++ {
		for j := 0; j < i; j++ {
			if math.Abs(A[i][j]) > epsilon {
				return false
			}
		}
	}
	return true
}

func extractEigenvalues(A Matrix, epsilon float64) []complex128 {
	n := len(A)
	eigenvalues := make([]complex128, n)

	for i := 0; i < n; {
		if i < n-1 && math.Abs(A[i+1][i]) > epsilon {
			// Обработка блока 2x2 для комплексных собственных значений
			a := A[i][i]
			b := A[i][i+1]
			c := A[i+1][i]
			d := A[i+1][i+1]

			trace := a + d
			det := a*d - b*c
			discriminant := trace*trace - 4*det

			if discriminant < 0 {
				realPart := trace / 2
				imagPart := math.Sqrt(-discriminant) / 2
				eigenvalues[i] = complex(realPart, imagPart)
				eigenvalues[i+1] = complex(realPart, -imagPart)
				i += 2
			} else {
				eigenvalues[i] = complex((trace+math.Sqrt(discriminant))/2, 0)
				eigenvalues[i+1] = complex((trace-math.Sqrt(discriminant))/2, 0)
				i += 2
			}
		} else {
			eigenvalues[i] = complex(A[i][i], 0)
			i++
		}
	}

	return eigenvalues
}

func verifyEigenvalues(A Matrix, eigenvalues []complex128) []float64 {
	n := len(A)
	errors := make([]float64, n)

	for k := 0; k < n; k++ {
		maxError := 0.0
		v := make([]complex128, n)
		for i := range v {
			v[i] = complex(rand.Float64(), 0)
		}

		// Вычисляем A*v
		Av := make([]complex128, n)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				Av[i] += complex(A[i][j], 0) * v[j]
			}
		}

		// Вычисляем λ*v
		lambdaV := make([]complex128, n)
		for i := range lambdaV {
			lambdaV[i] = eigenvalues[k] * v[i]
		}

		// Вычисляем максимальную относительную ошибку
		for i := 0; i < n; i++ {
			error := cmplx.Abs(Av[i]-lambdaV[i]) / cmplx.Abs(lambdaV[i])
			if error > maxError {
				maxError = error
			}
		}
		errors[k] = maxError
	}

	return errors
}

func multiplyMatrices(A, B Matrix) Matrix {
	n := len(A)
	result := make(Matrix, n)
	for i := range result {
		result[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				result[i][j] += A[i][k] * B[k][j]
			}
		}
	}
	return result
}

func copyMatrix(A Matrix) Matrix {
	n := len(A)
	copy := make(Matrix, n)
	for i := range copy {
		copy[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			copy[i][j] = A[i][j]
		}
	}
	return copy
}

func writeResults(filename string, A Matrix, eigenvalues []complex128, iterations int, epsilon float64, errors []float64) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	writer.WriteString("Исходная матрица:\n")
	for _, row := range A {
		for _, val := range row {
			writer.WriteString(fmt.Sprintf("%10.6f ", val))
		}
		writer.WriteString("\n")
	}

	writer.WriteString(fmt.Sprintf("\nТочность вычислений: %.0e\n", epsilon))
	writer.WriteString(fmt.Sprintf("Количество итераций: %d\n", iterations))

	writer.WriteString("\nСобственные значения:\n")
	for i, val := range eigenvalues {
		if imag(val) == 0 {
			writer.WriteString(fmt.Sprintf("λ%d = %.6f", i+1, real(val)))
		} else {
			writer.WriteString(fmt.Sprintf("λ%d = %.6f + %.6fi", i+1, real(val), imag(val)))
		}
		writer.WriteString(fmt.Sprintf(" (ошибка: %.3e)\n", errors[i]))
	}

	return writer.Flush()
}
