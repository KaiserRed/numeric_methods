package main

import (
	"fmt"
	"math"
)

func f(x float64) float64 {
	return 1 / math.Sqrt((2*x+7)*(3*x+4))
}

type StepResult struct {
	x       float64
	y       float64
	rect    float64
	trap    float64
	simpson float64
}

func rectangleMethod(a, b float64, h float64) []StepResult {
	n := int((b-a)/h) + 1
	results := make([]StepResult, n)
	sum := 0.0

	for i := 0; i < n; i++ {
		x := a + float64(i)*h
		y := f(x)
		if i > 0 {
			mid := (results[i-1].x + x) / 2
			sum += f(mid) * h
		}
		results[i] = StepResult{x: x, y: y, rect: sum}
	}
	return results
}

func trapezoidMethod(a, b float64, h float64) []StepResult {
	n := int((b-a)/h) + 1
	results := make([]StepResult, n)
	sum := 0.0
	prevY := f(a)

	for i := 0; i < n; i++ {
		x := a + float64(i)*h
		y := f(x)
		if i > 0 {
			sum += (prevY + y) * h / 2
		}
		results[i] = StepResult{x: x, y: y, trap: sum}
		prevY = y
	}
	return results
}

func simpsonMethod(a, b float64, h float64) []StepResult {
	n := int(math.Round((b - a) / h))
	if n%2 != 0 {
		n++
	}
	h = (b - a) / float64(n)

	results := make([]StepResult, n+1)
	sum := f(a)
	results[0] = StepResult{x: a, y: f(a), simpson: 0}

	for i := 1; i < n; i++ {
		x := a + float64(i)*h
		y := f(x)

		if i%2 == 1 {
			sum += 4 * y
		} else {
			sum += 2 * y
		}

		if i%2 == 1 {
			currentSum := h / 3 * (sum + f(x+h))
			results[i+1] = StepResult{
				x:       x + h,
				y:       f(x + h),
				simpson: currentSum,
			}
		}
	}

	sum += f(b)
	finalResult := h / 3 * sum

	results[n] = StepResult{
		x:       b,
		y:       f(b),
		simpson: finalResult,
	}

	return results
}

func calculateAllMethods(a, b, h float64) []StepResult {
	rect := rectangleMethod(a, b, h)
	trap := trapezoidMethod(a, b, h)
	simp := simpsonMethod(a, b, h)

	maxLen := len(rect)
	if len(trap) > maxLen {
		maxLen = len(trap)
	}
	if len(simp) > maxLen {
		maxLen = len(simp)
	}

	results := make([]StepResult, maxLen)
	for i := range results {
		if i < len(rect) {
			results[i].x = rect[i].x
			results[i].y = rect[i].y
			results[i].rect = rect[i].rect
		}
		if i < len(trap) {
			results[i].trap = trap[i].trap
		}
		if i < len(simp) && simp[i].x == results[i].x {
			results[i].simpson = simp[i].simpson
		}
	}
	return results
}

func printTable(title string, results []StepResult) {
	fmt.Printf("\n%s\n", title)
	fmt.Println("+----+-------+---------+---------------+---------------+---------------+")
	fmt.Println("| i |  x_i  |   y_i   | Прямоугольники|    Трапеции   |    Симпсон    |")
	fmt.Println("+----+-------+---------+---------------+---------------+---------------+")

	for i, res := range results {
		fmt.Printf("| %2d | %5.2f | %7.5f | %13.5f | %13.5f | %13.5f |\n",
			i, res.x, res.y, res.rect, res.trap, res.simpson)
	}
	fmt.Println("+----+-------+---------+---------------+---------------+---------------+")
}

func main() {
	a, b := 0.0, 4.0
	h1, h2 := 1.0, 0.5
	fmt.Println("Интеграл: F = ∫[1.0, 0.5] 1/√((2*x+7)*(3*x+4)) dx")

	resH1 := calculateAllMethods(a, b, h1)
	printTable(fmt.Sprintf("Таблица 1: Шаг h = %.1f", h1), resH1)

	resH2 := calculateAllMethods(a, b, h2)
	printTable(fmt.Sprintf("Таблица 2: Шаг h = %.2f", h2), resH2)

	// Уточнение методом Рунге-Ромберга
	F_rect_h1 := resH1[len(resH1)-1].rect
	F_rect_h2 := resH2[len(resH2)-1].rect
	F_trap_h1 := resH1[len(resH1)-1].trap
	F_trap_h2 := resH2[len(resH2)-1].trap
	F_simp_h1 := resH1[len(resH1)-1].simpson
	F_simp_h2 := resH2[len(resH2)-1].simpson

	rect_refined := rungeRomberg(F_rect_h1, F_rect_h2, h1/h2, 2)
	trap_refined := rungeRomberg(F_trap_h1, F_trap_h2, h1/h2, 2)
	simp_refined := rungeRomberg(F_simp_h1, F_simp_h2, h1/h2, 4)

	// Вывод уточненных результатов
	fmt.Printf("\nУточнение методом Рунге-Ромберга-Ричардсона:\n")
	fmt.Println("+------------------+------------------+------------------+")
	fmt.Println("|   Метод          | Уточненное значение | Погрешность     |")
	fmt.Println("+------------------+------------------+------------------+")
	fmt.Printf("| Прямоугольники   | %16.5f | %14.5f |\n", rect_refined, math.Abs(rect_refined-F_rect_h2))
	fmt.Printf("| Трапеции         | %16.5f | %14.5f |\n", trap_refined, math.Abs(trap_refined-F_trap_h2))
	fmt.Printf("| Симпсон          | %16.5f | %14.5f |\n", simp_refined, math.Abs(simp_refined-F_simp_h2))
	fmt.Println("+------------------+------------------+------------------+")

	// Проверка 1: Ошибки уменьшаются с уменьшением шага
	error_rect_h1 := math.Abs(rect_refined - F_rect_h1)
	error_rect_h2 := math.Abs(rect_refined - F_rect_h2)

	error_trap_h1 := math.Abs(trap_refined - F_trap_h1)
	error_trap_h2 := math.Abs(trap_refined - F_trap_h2)

	error_simp_h1 := math.Abs(simp_refined - F_simp_h1)
	error_simp_h2 := math.Abs(simp_refined - F_simp_h2)

	check1 := error_rect_h2 < error_rect_h1 &&
		error_trap_h2 < error_trap_h1 &&
		error_simp_h2 < error_simp_h1

	// Проверка 2: Метод Симпсона точнее других
	error_rect_final := math.Abs(rect_refined - F_rect_h2)
	error_trap_final := math.Abs(trap_refined - F_trap_h2)
	error_simp_final := math.Abs(simp_refined - F_simp_h2)

	check2 := error_simp_final < error_rect_final &&
		error_simp_final < error_trap_final

	fmt.Println("\nПроверки:")
	fmt.Printf("1. Ошибки уменьшаются с уменьшением шага: %t\n", check1)
	fmt.Printf("2. Метод Симпсона точнее других методов: %t\n", check2)

}

func rungeRomberg(Fh, Fkh, k, p float64) float64 {
	return Fkh + (Fkh-Fh)/(math.Pow(k, p)-1)
}
