package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/KaiserRed/numeric_methods/internal/lu_decompose"
)

type Matrix [][]float64

func getBaseDir() (string, error) {
	_, filename, _, ok := runtime.Caller(0)
	if !ok {
		return "", fmt.Errorf("не удалось определить путь к файлу")
	}
	return filepath.Dir(filename), nil
}

func readMatrixAndVectorFromFile(filename string) ([][]float64, []float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, fmt.Errorf("ошибка при открытии файла: %w", err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	var matrix [][]float64
	var vector []float64

	for scanner.Scan() {
		line := scanner.Text()
		elements := strings.Fields(line)
		if len(elements) < 2 {
			return nil, nil, fmt.Errorf("недостаточно данных в строке")
		}

		b, err := strconv.ParseFloat(elements[len(elements)-1], 64)
		if err != nil {
			return nil, nil, fmt.Errorf("ошибка при преобразовании строки в число: %w", err)
		}
		vector = append(vector, b)

		row := make([]float64, len(elements)-1)
		for i := 0; i < len(elements)-1; i++ {
			num, err := strconv.ParseFloat(elements[i], 64)
			if err != nil {
				return nil, nil, fmt.Errorf("ошибка при преобразовании строки в число: %w", err)
			}
			row[i] = num
		}
		matrix = append(matrix, row)
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, fmt.Errorf("ошибка при чтении файла: %w", err)
	}

	if len(matrix) != len(matrix[0]) {
		return nil, nil, fmt.Errorf("матрица должна быть квадратной")
	}

	return matrix, vector, nil
}

func writeMatrix(writer *bufio.Writer, title string, matrix [][]float64) error {
	_, err := writer.WriteString(title + ":\n")
	if err != nil {
		return err
	}

	for _, row := range matrix {
		for _, val := range row {
			_, err := writer.WriteString(fmt.Sprintf("%8.4f ", val))
			if err != nil {
				return err
			}
		}
		_, err := writer.WriteString("\n")
		if err != nil {
			return err
		}
	}
	_, err = writer.WriteString("\n")
	return err
}

func writeResultsToFile(filename string, x, b []float64, A [][]float64, det float64, invA, L, U, P [][]float64) error {
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("ошибка при создании файла: %w", err)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	if err := writeMatrix(writer, "Матрица A", A); err != nil {
		return fmt.Errorf("ошибка при записи матрицы A: %w", err)
	}

	_, _ = writer.WriteString("Вектор b:\n")
	for _, val := range b {
		_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
	}
	_, _ = writer.WriteString("\n\n")

	_, _ = writer.WriteString("Решение СЛАУ (вектор x):\n")
	for _, val := range x {
		_, _ = writer.WriteString(fmt.Sprintf("%8.4f ", val))
	}
	_, _ = writer.WriteString("\n\n")

	err = lu_decompose.VerifySolution(A, x, b)
	if err != nil {
		_, _ = writer.WriteString(fmt.Sprintf("Проверка решения: %v\n\n", err))
	} else {
		_, _ = writer.WriteString("Проверка решения: решение корректно (Ax = b)\n\n")
	}

	_, _ = writer.WriteString(fmt.Sprintf("Определитель матрицы A: %.4f\n\n", det))

	if err := writeMatrix(writer, "Матрица перестановок P", P); err != nil {
		return fmt.Errorf("ошибка при записи матрицы P: %w", err)
	}

	if err := writeMatrix(writer, "Матрица L", L); err != nil {
		return fmt.Errorf("ошибка при записи матрицы L: %w", err)
	}

	if err := writeMatrix(writer, "Матрица U", U); err != nil {
		return fmt.Errorf("ошибка при записи матрицы U: %w", err)
	}

	if err := writeMatrix(writer, "Обратная матрица A^{-1}", invA); err != nil {
		return fmt.Errorf("ошибка при записи обратной матрицы: %w", err)
	}

	return writer.Flush()
}

func main() {
	baseDir, err := getBaseDir()
	if err != nil {
		fmt.Println("Ошибка при получении базовой директории:", err)
		return
	}

	inputPath := filepath.Join(baseDir, "input.txt")

	A, b, err := readMatrixAndVectorFromFile(inputPath)
	if err != nil {
		fmt.Println("Ошибка при чтении данных из файла:", err)
		return
	}

	x, err := lu_decompose.SolveLinearSystem(A, b)
	if err != nil {

		fmt.Println("Ошибка при решении СЛАУ:", err)
		return
	}

	L, U, P, err := lu_decompose.LUDecomposition(A)
	if err != nil {
		fmt.Println("Ошибка при LU-разложении:", err)
		return
	}
	det := lu_decompose.Determinant(U, P)

	invA, err := lu_decompose.InverseMatrix(A)
	if err != nil {
		fmt.Println("Ошибка при нахождении обратной матрицы:", err)
		return
	}

	outputPath := filepath.Join(baseDir, "output.txt")
	err = writeResultsToFile(outputPath, x, b, A, det, invA, L, U, P)
	if err != nil {
		fmt.Println("Ошибка при записи результатов в файл:", err)
		return
	}

	fmt.Println("Результаты успешно записаны в файл output.txt")
}
