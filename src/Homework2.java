
public class Homework2 {

	public static void main(String[] args) {
		int n = 3;
		double[][] x = { { 2, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
		System.out.println(calcDeterminant(x, 3));
		System.out.println(verifyInputMatrix(x, n));
	}

	public static boolean verifyInputMatrix(double[][] a, int n) {
		if (a.length != a[0].length)
			return false;

		for (int i = 1; i < a.length / 2; ++i)
			for (int j = 1; j < a.length / 2; ++j) {
				if (a[i][j] != a[j][i])
					return false;
			}
		for (int k = 1; k <= n; ++k)
			for (int i = 0; i < n - k + 1; ++i)
				for (int j = 0; j < n - k + 1; ++j) {
					double aux[][] = getSubMatrix(a, k, i, j);
					if (calcDeterminant(aux, aux.length) == 0)
						return false;
				}

		return true;
	}

	public static double[][] getSubMatrix(double[][] a, int n, int posI, int posJ) {
		double[][] r = new double[n][n];

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				r[i][j] = a[i + posI][j + posJ];
				System.out.print(r[i][j] + " ");
			}
			System.out.println();
			if (i == n - 1)
				System.out.println();
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
