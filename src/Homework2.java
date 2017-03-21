
public class Homework2 {

	public static void main(String[] args) {
		int n = 3;
		double epsilon = 0.000001d;
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
			for (int p = 0; p < n; ++p) {
				d[p] = A[p][p];
				for (int k = 0; k <= p - 1; ++k)
					d[p] -= (A[p][p] / A[p][p]) * d[k];
				// d[p] -= d[k] * Math.pow(A[p][k], 2);
				if (d[p] == 0)
					System.exit('d');
			}
			System.out.print("D = ");
			for (double it : d)
				System.out.print(it + " ");
			System.out.println();

			double[][] L = new double[n][n];

			// L
			// for (int i = 0; i < n; ++i)
			// L[i][i] = 1;
			// for (int i = 0; i < n; ++i) {
			// for (int j = i + 1; j < n; ++j) {
			// double aux = 0d;
			// for (int k = 0; k < j - 1; ++k) {
			// // if (i == k && j == k)
			// // aux += d[k];
			// // else if (i == k)
			// // aux += L[j][k] * d[k];
			// // else if (j == k)
			// aux += L[i][k] * L[j][k] * d[k];
			// }
			// A[i][j] = A[j][i] = (A[j][i] - aux) / d[i];
			// }
			// }

			for (int p = 0; p < n; ++p) {
				for (int i = p + 1; i < n; ++i) {
					for (int k = 0; k <= p - 1; ++k) {
						A[i][p] -= d[k] * A[i][k] * A[p][k];
					}
					A[i][p] /= d[p];
					A[p][i] = A[i][p];
				}
			}
			System.out.println("L = ");
			for (double it : A[0])
				System.out.print(it + " ");
			System.out.println();
			for (int i = 1; i < n; ++i) {
				for (double it : A[i])
					System.out.print(it + " ");
				System.out.println();
			}

			// L[2][1] = L[1][2] = -0.5d;
			// double[][] rez;
			// double[][] D = { { 1, 0, 0 }, { 0, 4, 0 }, { 0, 0, 1 } };
			// rez = multiplyMatrix(L, D);
			// rez = multiplyMatrix(rez, transposeMatrix(L));
			// System.out.println("rez = ");
			// for (double it : rez[0])
			// System.out.print(it + " ");
			// System.out.println();
			// for (int i = 1; i < n; ++i) {
			// for (double it : rez[i])
			// System.out.print(it + " ");
			// System.out.println();
			// }
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

	public static double[][] getSubMatrix(double[][] a, int n, int posI, int posJ) {
		double[][] r = new double[n][n];

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				r[i][j] = a[i + posI][j + posJ];
				// System.out.print(r[i][j] + " ");
			}
			// System.out.println();
			// if (i == n - 1)
			// System.out.println();
		}

		return r;
	}

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

	public static double[][] transposeMatrix(double[][] a) {
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
