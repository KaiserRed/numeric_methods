package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

const maxIterations = 1000

func simpleIteration(A [][]float64, b []float64, epsilon float64) ([]float64, int, error) {
	n := len(b)
	x := make([]float64, n)
	xNew := make([]float64, n)
	iterations := 0

	for iterations < maxIterations {
		iterations++
		for i := 0; i < n; i++ {
			sum := 0.0
			for j := 0; j < n; j++ {
				if i != j {
					sum += A[i][j] * x[j]
				}
			}
			xNew[i] = (b[i] - sum) / A[i][i]
		}

		diff := 0.0
		for i := 0; i < n; i++ {
			diff += math.Pow(xNew[i]-x[i], 2)
		}
		diff = math.Sqrt(diff)

		if diff < epsilon {
			return xNew, iterations, nil
		}

		copy(x, xNew)
	}

	return nil, iterations, fmt.Errorf("достигнуто максимальное количество итераций (%d)", maxIterations)
}

func seidelMethod(A [][]float64, b []float64, epsilon float64) ([]float64, int, error) {
	n := len(b)
	x := make([]float64, n)
	iterations := 0

	for iterations < maxIterations {
		iterations++
		xNew := make([]float64, n)
		copy(xNew, x)

		for i := 0; i < n; i++ {
			sum1 := 0.0
			for j := 0; j < i; j++ {
				sum1 += A[i][j] * xNew[j]
			}

			sum2 := 0.0
			for j := i + 1; j < n; j++ {
				sum2 += A[i][j] * x[j]
			}

			xNew[i] = (b[i] - sum1 - sum2) / A[i][i]
		}

		diff := 0.0
		for i := 0; i < n; i++ {
			diff += math.Pow(xNew[i]-x[i], 2)
		}
		diff = math.Sqrt(diff)

		if diff < epsilon {
			return xNew, iterations, nil
		}

		copy(x, xNew)
	}

	return nil, iterations, fmt.Errorf("достигнуто максимальное количество итераций (%d)", maxIterations)
}

func readInput(filename string) (A [][]float64, b []float64, epsilon float64, err error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	lines := make([]string, 0)

	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}

	if len(lines) < 2 {
		return nil, nil, 0, fmt.Errorf("файл должен содержать минимум 2 строки")
	}

	A = make([][]float64, len(lines)-1)
	for i := 0; i < len(lines)-1; i++ {
		parts := strings.Fields(lines[i])
		A[i] = make([]float64, len(parts))
		for j, p := range parts {
			val, err := strconv.ParseFloat(p, 64)
			if err != nil {
				return nil, nil, 0, fmt.Errorf("ошибка в строке %d: %v", i+1, err)
			}
			A[i][j] = val
		}
	}

	parts := strings.Fields(lines[len(lines)-1])
	b = make([]float64, len(parts)-1)
	for i := 0; i < len(parts)-1; i++ {
		val, err := strconv.ParseFloat(parts[i], 64)
		if err != nil {
			return nil, nil, 0, fmt.Errorf("ошибка в строке правых частей: %v", err)
		}
		b[i] = val
	}

	epsilon, err = strconv.ParseFloat(parts[len(parts)-1], 64)
	if err != nil {
		return nil, nil, 0, fmt.Errorf("ошибка при чтении точности: %v", err)
	}

	return A, b, epsilon, nil
}

func verifySolution(A [][]float64, x, b []float64) float64 {
	n := len(A)
	maxError := 0.0
	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += A[i][j] * x[j]
		}
		error := math.Abs(sum - b[i])
		if error > maxError {
			maxError = error
		}
	}
	return maxError
}

func writeResults(filename string, A [][]float64, b []float64, epsilon float64, xSimple, xSeidel []float64, iterSimple, iterSeidel int) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	defer writer.Flush()

	_, _ = writer.WriteString("Исходная матрица A:\n")
	for _, row := range A {
		for _, val := range row {
			_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
		}
		_, _ = writer.WriteString("\n")
	}

	_, _ = writer.WriteString("\nВектор правых частей b:\n")
	for _, val := range b {
		_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
	}

	_, _ = writer.WriteString(fmt.Sprintf("\n\nТочность вычислений: %.0e\n\n", epsilon))

	_, _ = writer.WriteString("Метод простых итераций:\n")
	_, _ = writer.WriteString(fmt.Sprintf("Количество итераций: %d\n", iterSimple))
	_, _ = writer.WriteString("Решение:\n")
	for _, val := range xSimple {
		_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
	}
	errorSimple := verifySolution(A, xSimple, b)
	_, _ = writer.WriteString(fmt.Sprintf("\nМаксимальная невязка: %.4e\n\n", errorSimple))

	_, _ = writer.WriteString("Метод Зейделя:\n")
	_, _ = writer.WriteString(fmt.Sprintf("Количество итераций: %d\n", iterSeidel))
	_, _ = writer.WriteString("Решение:\n")
	for _, val := range xSeidel {
		_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
	}
	errorSeidel := verifySolution(A, xSeidel, b)
	_, _ = writer.WriteString(fmt.Sprintf("\nМаксимальная невязка: %.4e\n", errorSeidel))

	return nil
}

func main() {
	A, b, epsilon, err := readInput("input.txt")
	if err != nil {
		fmt.Printf("Ошибка чтения: %v\n", err)
		return
	}

	xSimple, iterSimple, err := simpleIteration(A, b, epsilon)
	if err != nil {
		fmt.Printf("Ошибка метода простых итераций: %v\n", err)
		return
	}

	xSeidel, iterSeidel, err := seidelMethod(A, b, epsilon)
	if err != nil {
		fmt.Printf("Ошибка метода Зейделя: %v\n", err)
		return
	}

	if err := writeResults("output.txt", A, b, epsilon, xSimple, xSeidel, iterSimple, iterSeidel); err != nil {
		fmt.Printf("Ошибка записи: %v\n", err)
		return
	}

	fmt.Println("Результаты успешно записаны в файл output.txt")
}
