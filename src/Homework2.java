
public class Homework2 {

	public static void main(String[] args) {
		int n = 3;
		double[] d = new double[n];
		double[][] A = { { 1, -1, 2 }, { -1, 5, -4 }, { 2, -4, 6 } };
		// { 1, -1, 2 }
		// { -1, 5, -4 }
		// { 2, -4, 6 }
		boolean matrixIsValid = isMatrixValid(A, n);

		System.out.println("Det(A) = " + calcDeterminant(A, 3));
		System.out.println("Matrice valida: " + matrixIsValid);

		if (matrixIsValid) {
			// D
			for (int i = 0; i < n; ++i) {
				d[i] = A[i][i];
				for (int j = 0; j <= i - 1; ++j)
					d[i] -= (A[i][i]/A[i][i])*d[j];
				
			}
			System.out.print("D = ");
			for (double it : d)
				System.out.print(it + " ");
			System.out.println();
		}
	}

	public static boolean isMatrixValid(double[][] a, int n) {
		if (a.length != a[0].length)
			throw new RuntimeException("Matricea nu este patratica!");

		for (int i = 1; i < n; ++i)
			for (int j = 1; j < n; ++j) {
				if (a[i][j] != a[j][i])
					throw new RuntimeException("Matricea nu este simetrica!");
			}

		// for (int k = 1; k <= n; ++k)
		// for (int i = 0; i < n - k + 1; ++i)
		// for (int j = 0; j < n - k + 1; ++j) {
		// double aux[][] = getSubMatrix(a, k, i, j);
		// if (calcDeterminant(aux, aux.length) == 0)
		// throw new RuntimeException("Matricea nu este pozitiv definita!");
		// }

		if (calcDeterminant(a, n) == 0)
			return false;

		return true;
	}

	// public static double[][] getSubMatrix(double[][] a, int n, int posI, int
	// posJ) {
	// double[][] r = new double[n][n];
	//
	// for (int i = 0; i < n; ++i) {
	// for (int j = 0; j < n; ++j) {
	// r[i][j] = a[i + posI][j + posJ];
	// // System.out.print(r[i][j] + " ");
	// }
	// // System.out.println();
	// // if (i == n - 1)
	// // System.out.println();
	// }
	//
	// return r;
	// }

	public static double calcDeterminant(double[][] a, int n) {
		double r = 0;

		if (n == 1) {
			r = a[0][0];
		} else if (n == 2) {
			r = a[0][0] * a[1][1] - a[1][0] * a[0][1];
		} else {
			r = 0;
			for (int j1 = 0; j1 < n; j1++) {
				double[][] m = new double[n - 1][];
				for (int k = 0; k < (n - 1); k++) {
					m[k] = new double[n - 1];
				}
				for (int i = 1; i < n; i++) {
					int j2 = 0;
					for (int j = 0; j < n; j++) {
						if (j == j1)
							continue;
						m[i - 1][j2] = a[i][j];
						j2++;
					}
				}
				r += Math.pow(-1.0, 1.0 + j1 + 1.0) * a[0][j1] * calcDeterminant(m, n - 1);
			}
		}
		return r;
	}

	public double[][] transposeMatrix(double[][] a) {
		double[][] r = new double[a.length][a.length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a.length; j++) {
				r[j][i] = a[i][j];
			}
		}
		return r;
	}

	public static double[][] multiplyMatrix(double[][] a, double[][] b) {
		double[][] r = new double[a.length][a.length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a.length; j++) {
				double s = 0;
				for (int k = 0; k < a.length; k++) {
					s += a[i][k] * b[k][j];
				}
				r[i][j] = s;
			}
		}
		return r;
	}

}
