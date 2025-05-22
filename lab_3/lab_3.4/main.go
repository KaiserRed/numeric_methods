package main

import (
	"fmt"
)

func main() {
	// Исходные данные
	x := []float64{1.0, 1.5, 2.0, 2.5, 3.0}
	y := []float64{2.0, 2.1667, 2.5, 2.9, 3.3333}
	Xstar := 2.0

	fmt.Println("Исходные данные:")
	fmt.Println("i\tx_i\ty_i")
	for i := 0; i < len(x); i++ {
		fmt.Printf("%d\t%.2f\t%.5f\n", i, x[i], y[i])
	}
	fmt.Println()

	// Найдем индекс точки X*
	idx := -1
	for i, val := range x {
		if val == Xstar {
			idx = i
			break
		}
	}
	if idx == -1 {
		fmt.Println("Точка X* не найдена в массиве x")
		return
	}

	h := x[1] - x[0] // шаг (предполагается равномерный)

	fmt.Println("Результаты:")
	fmt.Printf("Вычисляем первую производную в точке X = %.2f\n\n", Xstar)

	// 1) Левосторонняя производная (первый порядок)
	leftDiff := (y[idx] - y[idx-1]) / (x[idx] - x[idx-1])
	fmt.Println("Левосторонняя производная (первый порядок):")
	fmt.Printf("Формула: y'(X*) ≈ (y_%d - y_%d) / (x_%d - x_%d) = (%.5f - %.5f) / (%.2f - %.2f) = %.6f\n\n",
		idx, idx-1, idx, idx-1, y[idx], y[idx-1], x[idx], x[idx-1], leftDiff)

	// 2) Правосторонняя производная (первый порядок)
	rightDiff := (y[idx+1] - y[idx]) / (x[idx+1] - x[idx])
	fmt.Println("Правосторонняя производная (первый порядок):")
	fmt.Printf("Формула: y'(X*) ≈ (y_%d - y_%d) / (x_%d - x_%d) = (%.5f - %.5f) / (%.2f - %.2f) = %.6f\n\n",
		idx+1, idx, idx+1, idx, y[idx+1], y[idx], x[idx+1], x[idx], rightDiff)

	// 3) Центральная разность (производная второго порядка точности)
	centerDiff := (y[idx+1] - y[idx-1]) / (x[idx+1] - x[idx-1])
	fmt.Println("Производная второго порядка точности (центральная разность):")
	fmt.Printf("Формула: y'(X*) ≈ (y_%d - y_%d) / (x_%d - x_%d) = (%.5f - %.5f) / (%.2f - %.2f) = %.6f\n\n",
		idx+1, idx-1, idx+1, idx-1, y[idx+1], y[idx-1], x[idx+1], x[idx-1], centerDiff)

	fmt.Println("Проверка через сумму:")
	fmt.Printf("Сумма (левосторонняя + правосторонняя) / 2 = (%.6f + %.6f) / 2 = %.6f\n\n",
		leftDiff, rightDiff, (leftDiff+rightDiff)/2)

	// 4) Вторая производная (второй порядок)
	secondDiff := (y[idx+1] - 2*y[idx] + y[idx-1]) / ((x[idx+1] - x[idx]) * (x[idx] - x[idx-1]))
	fmt.Println("Вторая производная (второй порядок):")
	fmt.Printf("Формула: y''(X*) ≈ (y_%d - 2*y_%d + y_%d) / (h^2) = (%.5f - 2*%.5f + %.5f) / (%.2f^2) = %.6f\n\n",
		idx+1, idx, idx-1, y[idx+1], y[idx], y[idx-1], h, secondDiff)

	// 5) Проверка второй производной через разности
	// вычислим f'(x+h) и f'(x-h) по формулам первого порядка
	fpxh := (y[idx+2] - y[idx+1]) / (x[idx+2] - x[idx+1])
	fmxh := (y[idx] - y[idx-1]) / (x[idx] - x[idx-1])
	secondDiffCheck := (fpxh - fmxh) / (x[idx+1] - x[idx])

	fmt.Println("Проверка второй производной через разности первых производных:")
	fmt.Printf("Формула: y''(X*) ≈ (y'(x+h) - y'(x-h)) / (2h) = (%.6f - %.6f) / (2*%.2f) = %.6f\n",
		fpxh, fmxh, h, secondDiffCheck)
}
