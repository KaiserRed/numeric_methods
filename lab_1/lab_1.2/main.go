package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

func thomasAlgorithm(a, b, c, d []float64) ([]float64, error) {
	n := len(d)
	if len(a) != n-1 || len(b) != n || len(c) != n-1 {
		return nil, fmt.Errorf("некорректные размеры входных данных")
	}

	cp := make([]float64, n)
	dp := make([]float64, n)

	cp[0] = c[0] / b[0]
	dp[0] = d[0] / b[0]

	for i := 1; i < n; i++ {
		denom := b[i] - a[i-1]*cp[i-1]
		if denom == 0 {
			return nil, fmt.Errorf("деление на ноль на шаге %d", i)
		}
		if i < n-1 {
			cp[i] = c[i] / denom
		}
		dp[i] = (d[i] - a[i-1]*dp[i-1]) / denom
	}

	x := make([]float64, n)
	x[n-1] = dp[n-1]

	for i := n - 2; i >= 0; i-- {
		x[i] = dp[i] - cp[i]*x[i+1]
	}

	return x, nil
}

func extractDiagonals(matrix [][]float64) (a, b, c []float64, err error) {
	n := len(matrix)
	if n == 0 {
		return nil, nil, nil, fmt.Errorf("матрица пуста")
	}

	for _, row := range matrix {
		if len(row) != n {
			return nil, nil, nil, fmt.Errorf("матрица не квадратная")
		}
	}

	a = make([]float64, n-1)
	b = make([]float64, n)
	c = make([]float64, n-1)

	for i := 0; i < n; i++ {
		b[i] = matrix[i][i]
		if i < n-1 {
			a[i] = matrix[i+1][i]
			c[i] = matrix[i][i+1]
		}
	}

	return a, b, c, nil
}

func isStrictlyTridiagonal(matrix [][]float64) bool {
	n := len(matrix)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if (j < i-1 || j > i+1) && matrix[i][j] != 0 {
				return false
			}
		}
	}
	return true
}

func readInput(filename string) (matrix [][]float64, d []float64, err error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	lines := make([]string, 0)

	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}

	if len(lines) < 2 {
		return nil, nil, fmt.Errorf("файл должен содержать минимум 2 строки")
	}

	matrix = make([][]float64, len(lines)-1)
	for i := 0; i < len(lines)-1; i++ {
		parts := strings.Fields(lines[i])
		matrix[i] = make([]float64, len(parts))
		for j, p := range parts {
			val, err := strconv.ParseFloat(p, 64)
			if err != nil {
				return nil, nil, fmt.Errorf("ошибка в строке %d: %v", i+1, err)
			}
			matrix[i][j] = val
		}
	}

	parts := strings.Fields(lines[len(lines)-1])
	d = make([]float64, len(parts))
	for i, p := range parts {
		val, err := strconv.ParseFloat(p, 64)
		if err != nil {
			return nil, nil, fmt.Errorf("ошибка в строке правых частей: %v", err)
		}
		d[i] = val
	}

	return matrix, d, nil
}

func verifySolution(matrix [][]float64, x, d []float64) error {
	n := len(matrix)
	const tolerance = 1e-8

	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += matrix[i][j] * x[j]
		}
		if math.Abs(sum-d[i]) > tolerance {
			return fmt.Errorf("проверка не пройдена: строка %d, ожидалось %.8f, получилось %.8f",
				i+1, d[i], sum)
		}
	}
	return nil
}

func writeResults(filename string, matrix [][]float64, d, x []float64) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	defer writer.Flush()

	_, _ = writer.WriteString("Исходная матрица:\n")
	for _, row := range matrix {
		for _, val := range row {
			_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
		}
		_, _ = writer.WriteString("\n")
	}
	_, _ = writer.WriteString("\n")

	_, _ = writer.WriteString("Вектор правых частей:\n")
	for _, val := range d {
		_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
	}
	_, _ = writer.WriteString("\n\n")

	_, _ = writer.WriteString("Решение СЛАУ:\n")
	for _, val := range x {
		_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
	}
	_, _ = writer.WriteString("\n\n")

	err = verifySolution(matrix, x, d)
	if err != nil {
		_, _ = writer.WriteString(fmt.Sprintf("Проверка решения: %v\n", err))
	} else {
		_, _ = writer.WriteString("Проверка решения: решение корректно\n")
	}

	return nil
}

func main() {
	matrix, d, err := readInput("input.txt")
	if err != nil {
		fmt.Printf("Ошибка чтения: %v\n", err)
		return
	}

	if !isStrictlyTridiagonal(matrix) {
		fmt.Println("Ошибка: матрица не является строго трёхдиагональной!")
		return
	}

	a, b, c, err := extractDiagonals(matrix)
	if err != nil {
		fmt.Printf("Ошибка извлечения диагоналей: %v\n", err)
		return
	}

	x, err := thomasAlgorithm(a, b, c, d)
	if err != nil {
		fmt.Printf("Ошибка решения: %v\n", err)
		return
	}

	if err := writeResults("output.txt", matrix, d, x); err != nil {
		fmt.Printf("Ошибка записи: %v\n", err)
		return
	}

	fmt.Println("Решение успешно записано в output.txt")
}
