package main

import (
	"bufio"
	"fmt"
	"math"
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

	if !isSymmetric(A) {
		fmt.Println("Ошибка: матрица не симметрическая!")
		return
	}

	eigenvalues, eigenvectors, iterErrors, iterations := jacobiRotation(A, epsilon)

	if err := writeResults("output.txt", A, eigenvalues, eigenvectors, iterErrors, iterations, epsilon); err != nil {
		fmt.Printf("Ошибка записи: %v\n", err)
		return
	}

	fmt.Println("Вычисления успешно завершены. Результаты в output.txt")
}

func readInput(filename string) (A [][]float64, epsilon float64, err error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	lines := make([]string, 0)

	for scanner.Scan() {
		line := scanner.Text()
		if line != "" {
			lines = append(lines, line)
		}
	}

	if len(lines) < 2 {
		return nil, 0, fmt.Errorf("файл должен содержать минимум 2 строки (матрица и точность)")
	}

	A = make([][]float64, len(lines)-1)
	for i := 0; i < len(lines)-1; i++ {
		parts := strings.Fields(lines[i])
		A[i] = make([]float64, len(parts))
		for j, p := range parts {
			val, err := strconv.ParseFloat(p, 64)
			if err != nil {
				return nil, 0, fmt.Errorf("ошибка в строке %d матрицы: %v", i+1, err)
			}
			A[i][j] = val
		}
	}

	epsilon, err = strconv.ParseFloat(strings.TrimSpace(lines[len(lines)-1]), 64)
	if err != nil {
		return nil, 0, fmt.Errorf("ошибка чтения точности: %v", err)
	}

	n := len(A)
	for _, row := range A {
		if len(row) != n {
			return nil, 0, fmt.Errorf("матрица должна быть квадратной")
		}
	}

	return A, epsilon, nil
}

func isSymmetric(A Matrix) bool {
	n := len(A)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if math.Abs(A[i][j]-A[j][i]) > 1e-6 {
				return false
			}
		}
	}
	return true
}

func jacobiRotation(A Matrix, epsilon float64) ([]float64, Matrix, []float64, int) {
	n := len(A)
	V := make(Matrix, n)
	for i := range V {
		V[i] = make([]float64, n)
		V[i][i] = 1.0
	}

	iterations := 0
	maxOffDiag := maxOffDiagonal(A)
	errors := []float64{maxOffDiag}

	for maxOffDiag > epsilon && iterations < 1000 {
		p, q := findMaxOffDiagonal(A)
		phi := 0.5 * math.Atan2(2*A[p][q], A[q][q]-A[p][p])
		c := math.Cos(phi)
		s := math.Sin(phi)

		newA := make(Matrix, n)
		for i := range newA {
			newA[i] = make([]float64, n)
			copy(newA[i], A[i])
		}

		for i := 0; i < n; i++ {
			if i != p && i != q {
				newA[i][p] = c*A[i][p] - s*A[i][q]
				newA[p][i] = newA[i][p]
				newA[i][q] = s*A[i][p] + c*A[i][q]
				newA[q][i] = newA[i][q]
			}
		}

		newA[p][p] = c*c*A[p][p] - 2*c*s*A[p][q] + s*s*A[q][q]
		newA[q][q] = s*s*A[p][p] + 2*c*s*A[p][q] + c*c*A[q][q]
		newA[p][q] = 0
		newA[q][p] = 0

		A = newA

		for i := 0; i < n; i++ {
			vip := V[i][p]
			viq := V[i][q]
			V[i][p] = c*vip - s*viq
			V[i][q] = s*vip + c*viq
		}

		maxOffDiag = maxOffDiagonal(A)
		errors = append(errors, maxOffDiag)
		iterations++
	}

	eigenvalues := make([]float64, n)
	for i := 0; i < n; i++ {
		eigenvalues[i] = A[i][i]
	}

	return eigenvalues, V, errors, iterations
}

func maxOffDiagonal(A Matrix) float64 {
	max := 0.0
	n := len(A)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if val := math.Abs(A[i][j]); val > max {
				max = val
			}
		}
	}
	return max
}

func findMaxOffDiagonal(A Matrix) (int, int) {
	max := 0.0
	p, q := 0, 0
	n := len(A)
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if val := math.Abs(A[i][j]); val > max {
				max = val
				p, q = i, j
			}
		}
	}
	return p, q
}

func verifyEigen(A Matrix, eigenvalues []float64, eigenvectors Matrix) []float64 {
	n := len(A)
	errors := make([]float64, n)

	for k := 0; k < n; k++ {
		maxError := 0.0
		for i := 0; i < n; i++ {
			sum := 0.0
			for j := 0; j < n; j++ {
				sum += A[i][j] * eigenvectors[j][k]
			}
			error := math.Abs(sum - eigenvalues[k]*eigenvectors[i][k])
			if error > maxError {
				maxError = error
			}
		}
		errors[k] = maxError
	}

	return errors
}

func writeResults(filename string, A Matrix, eigenvalues []float64, eigenvectors Matrix, iterErrors []float64, iterations int, epsilon float64) error {
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

	writer.WriteString("\nСобственные значения:\n")
	for i, val := range eigenvalues {
		writer.WriteString(fmt.Sprintf("λ%d = %.6f\n", i+1, val))
	}

	writer.WriteString("\nСобственные векторы (по столбцам):\n")
	for i := 0; i < len(eigenvectors); i++ {
		for j := 0; j < len(eigenvectors); j++ {
			writer.WriteString(fmt.Sprintf("%10.6f ", eigenvectors[i][j]))
		}
		writer.WriteString("\n")
	}

	verificationErrors := verifyEigen(A, eigenvalues, eigenvectors)
	writer.WriteString("\nПроверка точности (A*v - λ*v):\n")
	for i, err := range verificationErrors {
		writer.WriteString(fmt.Sprintf("Для λ%d: %.3e\n", i+1, err))
	}

	writer.WriteString("\nЗависимость погрешности от итераций:\n")
	for i, err := range iterErrors {
		writer.WriteString(fmt.Sprintf("%4d: %.3e\n", i, err))
	}

	writer.WriteString(fmt.Sprintf("\nВсего итераций: %d\n", iterations))

	return writer.Flush()
}
