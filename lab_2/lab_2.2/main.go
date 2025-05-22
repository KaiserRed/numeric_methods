package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

type Vector struct {
	X, Y float64
}

const a = 2.0

var (
	epsilon float64
	maxIter int
	output  *os.File
	writer  *bufio.Writer
)

func write(format string, args ...any) {
	fmt.Fprintf(writer, format, args...)
}

func parseInputFile(filename string) (Vector, error) {
	file, err := os.Open(filename)
	if err != nil {
		return Vector{}, err
	}
	defer file.Close()

	start := Vector{}
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "x0=") {
			start.X, _ = strconv.ParseFloat(strings.TrimPrefix(line, "x0="), 64)
		} else if strings.HasPrefix(line, "y0=") {
			start.Y, _ = strconv.ParseFloat(strings.TrimPrefix(line, "y0="), 64)
		} else if strings.HasPrefix(line, "epsilon=") {
			epsilon, _ = strconv.ParseFloat(strings.TrimPrefix(line, "epsilon="), 64)
		} else if strings.HasPrefix(line, "maxIter=") {
			maxIter, _ = strconv.Atoi(strings.TrimPrefix(line, "maxIter="))
		}
	}

	if err := scanner.Err(); err != nil {
		return Vector{}, err
	}

	return start, nil
}

func f1(x, y float64) float64 {
	return a*x*x - x + y*y - 1
}

func f2(x, y float64) float64 {
	return y - math.Tan(x)
}

func df1dx(x, y float64) float64 {
	return 2*a*x - 1
}

func df1dy(x, y float64) float64 {
	return 2 * y
}

func df2dx(x, y float64) float64 {
	return -1 / (math.Cos(x) * math.Cos(x))
}

func df2dy(x, y float64) float64 {
	return 1
}

func phi(x, y float64) (float64, float64) {
	yNext := math.Tan(x)
	f1Val := f1(x, yNext)
	xNext := x - 0.1*f1Val
	return xNext, yNext
}

func printSystem() {
	write("Система уравнений:\n")
	write("f1(x, y) = 2x^2 - x + y^2 - 1 = 0\n")
	write("f2(x, y) = y - tan(x) = 0\n")
}

func newtonMethod(start Vector) Vector {
	write("\n--- Метод Ньютона ---\n")
	x, y := start.X, start.Y
	for i := 0; i < maxIter; i++ {
		fx := f1(x, y)
		fy := f2(x, y)

		J := [2][2]float64{
			{df1dx(x, y), df1dy(x, y)},
			{df2dx(x, y), df2dy(x, y)},
		}
		det := J[0][0]*J[1][1] - J[0][1]*J[1][0]
		if math.Abs(det) < 1e-12 {
			write("Якобиан вырожден\n")
			break
		}

		dx := (fx*J[1][1] - fy*J[0][1]) / det
		dy := (fy*J[0][0] - fx*J[1][0]) / det

		xNext := x - dx
		yNext := y - dy

		errorVal := math.Hypot(xNext-x, yNext-y)
		write("Итерация %d: x = %.6f, y = %.6f, ошибка = %.6e\n", i+1, xNext, yNext, errorVal)

		if errorVal < epsilon {
			return Vector{xNext, yNext}
		}
		x, y = xNext, yNext
	}
	return Vector{x, y}
}

func simpleIteration(start Vector) Vector {
	write("\n--- Метод простой итерации ---\n")
	x, y := start.X, start.Y
	for i := 0; i < maxIter; i++ {
		xNext, yNext := phi(x, y)

		if math.IsNaN(xNext) || math.IsNaN(yNext) || math.Abs(x) > math.Pi/2-1e-3 {
			write("Получено недопустимое значение (NaN или x слишком близко к точке разрыва tan).\n")
			break
		}

		errorVal := math.Hypot(xNext-x, yNext-y)
		write("Итерация %d: x = %.6f, y = %.6f, ошибка = %.6e\n", i+1, xNext, yNext, errorVal)

		if errorVal < epsilon {
			return Vector{xNext, yNext}
		}

		x, y = xNext, yNext
	}
	return Vector{x, y}
}

func main() {
	var err error
	output, err = os.Create("output.txt")
	if err != nil {
		fmt.Println("Ошибка создания файла вывода:", err)
		return
	}
	defer output.Close()

	writer = bufio.NewWriter(output)
	defer writer.Flush()

	start, err := parseInputFile("input.txt")
	if err != nil {
		write("Ошибка чтения файла: %v\n", err)
		return
	}

	printSystem()
	write("\nНачальное приближение: x0 = %.6f, y0 = %.6f\n", start.X, start.Y)
	write("Заданная точность: epsilon = %.6e\n", epsilon)
	write("Максимальное число итераций: %d\n", maxIter)

	resultNewton := newtonMethod(start)
	write("\nРезультат (Метод Ньютона): x = %.6f, y = %.6f\n", resultNewton.X, resultNewton.Y)

	resultSimple := simpleIteration(start)
	write("\nРезультат (Метод простой итерации): x = %.6f, y = %.6f\n", resultSimple.X, resultSimple.Y)
}
