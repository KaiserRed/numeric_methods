package matrix

import "fmt"

func LUDecomposition(A [][]float64) (L, U [][]float64, err error) {
	n := len(A)
	if n == 0 || len(A[0]) != n {
		return nil, nil, fmt.Errorf("Матрица должна быть квадратной")
	}

	L = make([][]float64, n)
	U = make([][]float64, n)
	for i := range L {
		L[i] = make([]float64, n)
		U[i] = make([]float64, n)
	}

	for i := 0; i < n; i++ {
		for k := 0; k < n; k++ {
			sum := 0.0
			for j := 0; j < i; j++ {
				sum += L[i][j] * U[j][k]
			}
			U[i][k] = A[i][k] - sum
		}

		for k := i; k < n; k++ {
			if i == k {
				L[i][i] = 1
			} else {
				sum := 0.0
				for j := 0; j < i; j++ {
					sum += L[k][j] * U[j][i]
				}
				L[k][i] = (A[k][i] - sum) / U[i][i]
			}
		}
	}
	return L, U, nil

}

// L * z = b
func solveLowerTriangularVector(L [][]float64, b []float64) []float64 {
	n := len(L)
	z := make([]float64, n)

	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += L[i][j] * z[j]
		}
		z[i] = (b[i] - sum) / L[i][i]
	}
	return z
}

// U * x = z
func solveUpperTriangularVector(U [][]float64, z []float64) []float64 {
	n := len(U)
	x := make([]float64, n)

	for i := n - 1; i >= 0; i-- {
		sum := 0.0
		for j := i + 1; j < n; j++ {
			sum += U[i][j] * x[j]
		}
		x[i] = (z[i] - sum) / U[i][i]
	}
	return x
}
func inverseMatrix(A [][]float64) ([][]float64, error) {
	n := len(A)
	L, U, err := LUDecomposition(A)
	if err != nil {
		return nil, err
	}

	E := make([][]float64, n)
	for i := range E {
		E[i] = make([]float64, n)
		E[i][i] = 1
	}

	X := make([][]float64, n)
	for i := range X {
		X[i] = make([]float64, n)
	}

	for j := 0; j < n; j++ {
		b := make([]float64, n)
		for i := 0; i < n; i++ {
			b[i] = E[i][j]
		}
		z := solveLowerTriangularVector(L, b)
		x := solveUpperTriangularVector(U, z)

		for i := 0; i < n; i++ {
			X[i][j] = x[i]
		}
	}
	return X, nil
}

func determinant(U [][]float64) float64 {
	det := 1.0
	for i := 0; i < len(U); i++ {
		det *= U[i][i]
	}
	return det
}

func checkSolution(A [][]float64, x []float64, b []float64) error {
	n := len(A)
	result := make([]float64, n)

	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			sum += A[i][j] * x[j]
		}
		result[i] = sum
	}

	for i := 0; i < n; i++ {
		if fmt.Sprintf("%.4f", result[i]) != fmt.Sprintf("%.4f", b[i]) {
			return fmt.Errorf("решение неверное: A * x != b")
		}
	}

	return nil
}

func solveLinearSystem(A [][]float64, b []float64) ([]float64, error) {
	L, U, err := LUDecomposition(A)
	if err != nil {
		return nil, err
	}

	z := solveLowerTriangularVector(L, b)
	x := solveUpperTriangularVector(U, z)

	return x, nil

}
