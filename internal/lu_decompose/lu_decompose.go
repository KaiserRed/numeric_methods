package lu_decompose

import (
	"fmt"
	"math"
)

func LUDecomposition(A [][]float64) (L, U, P [][]float64, err error) {
	n := len(A)
	if n == 0 || len(A[0]) != n {
		return nil, nil, nil, fmt.Errorf("матрица должна быть квадратной")
	}

	L = make([][]float64, n)
	U = make([][]float64, n)
	P = make([][]float64, n)
	for i := range L {
		L[i] = make([]float64, n)
		U[i] = make([]float64, n)
		P[i] = make([]float64, n)
		P[i][i] = 1
	}

	for i := range A {
		copy(U[i], A[i])
	}

	for k := 0; k < n; k++ {
		maxRow := k
		maxVal := math.Abs(U[k][k])
		for i := k + 1; i < n; i++ {
			if math.Abs(U[i][k]) > maxVal {
				maxVal = math.Abs(U[i][k])
				maxRow = i
			}
		}

		if maxRow != k {
			U[k], U[maxRow] = U[maxRow], U[k]
			P[k], P[maxRow] = P[maxRow], P[k]
			if k > 0 {
				for j := 0; j < k; j++ {
					L[k][j], L[maxRow][j] = L[maxRow][j], L[k][j]
				}
			}
		}

		if math.Abs(U[k][k]) < 1e-12 {
			return nil, nil, nil, fmt.Errorf("матрица вырождена")
		}

		for i := k + 1; i < n; i++ {
			L[i][k] = U[i][k] / U[k][k]
			for j := k; j < n; j++ {
				U[i][j] -= L[i][k] * U[k][j]
			}
		}
	}

	for i := 0; i < n; i++ {
		L[i][i] = 1
	}

	return L, U, P, nil
}

func SolveLinearSystem(A [][]float64, b []float64) ([]float64, error) {
	L, U, P, err := LUDecomposition(A)
	if err != nil {
		return nil, err
	}

	n := len(b)
	pb := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if P[i][j] == 1 {
				pb[i] = b[j]
				break
			}
		}
	}

	// Решаем Lz = Pb
	z := make([]float64, n)
	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < i; j++ {
			sum += L[i][j] * z[j]
		}
		z[i] = pb[i] - sum
	}

	// Решаем Ux = z
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		sum := 0.0
		for j := i + 1; j < n; j++ {
			sum += U[i][j] * x[j]
		}
		x[i] = (z[i] - sum) / U[i][i]
	}

	return x, nil
}

func InverseMatrix(A [][]float64) ([][]float64, error) {
	n := len(A)
	L, U, P, err := LUDecomposition(A)
	if err != nil {
		return nil, err
	}

	invA := make([][]float64, n)
	for i := range invA {
		invA[i] = make([]float64, n)
	}

	for j := 0; j < n; j++ {
		e := make([]float64, n)
		e[j] = 1

		pe := make([]float64, n)
		for i := 0; i < n; i++ {
			for k := 0; k < n; k++ {
				if P[i][k] == 1 {
					pe[i] = e[k]
					break
				}
			}
		}

		// Решаем Lz = Pe
		z := make([]float64, n)
		for i := 0; i < n; i++ {
			sum := 0.0
			for k := 0; k < i; k++ {
				sum += L[i][k] * z[k]
			}
			z[i] = pe[i] - sum
		}

		// Решаем Ux = z
		x := make([]float64, n)
		for i := n - 1; i >= 0; i-- {
			sum := 0.0
			for k := i + 1; k < n; k++ {
				sum += U[i][k] * x[k]
			}
			x[i] = (z[i] - sum) / U[i][i]
		}

		for i := 0; i < n; i++ {
			invA[i][j] = x[i]
		}
	}

	return invA, nil
}

func Determinant(U [][]float64, P [][]float64) float64 {
	n := len(U)
	det := 1.0

	// Определитель U - произведение диагональных элементов
	for i := 0; i < n; i++ {
		det *= U[i][i]
	}

	// Определитель P - (-1)^число перестановок
	sign := 1.0
	for i := 0; i < n; i++ {
		if P[i][i] != 1 {
			sign *= -1
		}
	}

	return det * sign
}

// Ax = b
func VerifySolution(A [][]float64, x []float64, b []float64) error {
	n := len(A)
	if n == 0 || len(A[0]) != n || len(x) != n || len(b) != n {
		return fmt.Errorf("несовместимые размеры матрицы и векторов")
	}

	const tolerance = 1e-8
	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += A[i][j] * x[j]
		}
		if math.Abs(sum-b[i]) > tolerance {
			return fmt.Errorf("решение не удовлетворяет уравнению: строка %d, ожидалось %.8f, получилось %.8f",
				i, b[i], sum)
		}
	}
	return nil
}
