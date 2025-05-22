package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

type Params struct {
	a       float64
	b       float64
	epsilon float64
	maxIter int
}

func readInput(filename string) (Params, error) {
	file, err := os.Open(filename)
	if err != nil {
		return Params{}, err
	}
	defer file.Close()

	params := Params{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.Fields(line)
		if len(parts) < 2 {
			continue
		}

		key := strings.ToLower(parts[0])
		value := parts[1]

		switch key {
		case "a":
			val, err := strconv.ParseFloat(value, 64)
			if err != nil {
				return Params{}, fmt.Errorf("ошибка чтения a: %v", err)
			}
			params.a = val
		case "b":
			val, err := strconv.ParseFloat(value, 64)
			if err != nil {
				return Params{}, fmt.Errorf("ошибка чтения b: %v", err)
			}
			params.b = val
		case "epsilon":
			val, err := strconv.ParseFloat(value, 64)
			if err != nil {
				return Params{}, fmt.Errorf("ошибка чтения epsilon: %v", err)
			}
			params.epsilon = val
		case "max_iter":
			val, err := strconv.Atoi(value)
			if err != nil {
				return Params{}, fmt.Errorf("ошибка чтения max_iter: %v", err)
			}
			params.maxIter = val
		}
	}

	return params, nil
}

func iterationMethod(phi func(float64) float64, params Params) (float64, []float64, int) {
	x0 := (params.a + params.b) / 2
	x := x0
	errors := []float64{}

	for i := 0; i < params.maxIter; i++ {
		xNew := phi(x)
		err := math.Abs(xNew - x)
		errors = append(errors, err)

		if err < params.epsilon {
			return xNew, errors, i + 1
		}
		x = xNew
	}

	return x, errors, params.maxIter
}

func newtonMethod(f, df func(float64) float64, params Params) (float64, []float64, int) {
	x := (params.a + params.b) / 2
	errors := []float64{}

	for i := 0; i < params.maxIter; i++ {
		fx := f(x)
		dfx := df(x)

		if math.Abs(dfx) < 1e-12 {
			return math.NaN(), errors, i + 1
		}

		xNew := x - fx/dfx
		err := math.Abs(xNew - x)
		errors = append(errors, err)

		if err < params.epsilon {
			return xNew, errors, i + 1
		}
		x = xNew
	}

	return x, errors, params.maxIter
}

func writeResults(filename string, params Params, rootIter, rootNewton float64,
	errorsIter, errorsNewton []float64, iterIter, iterNewton int) error {

	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	defer writer.Flush()

	fmt.Fprintln(writer, "Анализ уравнения: ln(x + 2) - x⁴ + 0.5 = 0")
	fmt.Fprintln(writer, "Представлено как пересечение: ln(x + 2) = x⁴ - 0.5")
	fmt.Fprintf(writer, "\nИнтервал поиска: [%.2f, %.2f]\n", params.a, params.b)
	fmt.Fprintf(writer, "Точность: %.0e\n", params.epsilon)
	fmt.Fprintf(writer, "Макс. итераций: %d\n", params.maxIter)

	fmt.Fprintln(writer, "\nМетод простой итерации:")
	fmt.Fprintf(writer, "Корень: %.8f\n", rootIter)
	fmt.Fprintf(writer, "Итераций: %d\n", iterIter)
	fmt.Fprintln(writer, "Погрешности по итерациям:")
	for i, err := range errorsIter {
		fmt.Fprintf(writer, "%3d: %.3e\n", i+1, err)
	}

	fmt.Fprintln(writer, "\nМетод Ньютона:")
	fmt.Fprintf(writer, "Корень: %.8f\n", rootNewton)
	fmt.Fprintf(writer, "Итераций: %d\n", iterNewton)
	fmt.Fprintln(writer, "Погрешности по итерациям:")
	for i, err := range errorsNewton {
		fmt.Fprintf(writer, "%3d: %.3e\n", i+1, err)
	}

	return nil
}

func main() {
	f := func(x float64) float64 { return math.Log(x+2) - math.Pow(x, 4) + 0.5 }
	df := func(x float64) float64 { return 1/(x+2) - 4*math.Pow(x, 3) }
	phi := func(x float64) float64 { return math.Pow(math.Log(x+2)+0.5, 0.25) }

	params, err := readInput("input.txt")
	if err != nil {
		fmt.Printf("Ошибка чтения: %v\n", err)
		return
	}

	if f(params.a)*f(params.b) > 0 {
		fmt.Println("Ошибка: на концах интервала функция имеет одинаковые знаки")
		return
	}

	rootIter, errorsIter, iterIter := iterationMethod(phi, params)
	rootNewton, errorsNewton, iterNewton := newtonMethod(f, df, params)

	if err := writeResults("output.txt", params, rootIter, rootNewton,
		errorsIter, errorsNewton, iterIter, iterNewton); err != nil {
		fmt.Printf("Ошибка записи: %v\n", err)
		return
	}

	fmt.Println("Вычисления завершены. Результаты в output.txt")
}
