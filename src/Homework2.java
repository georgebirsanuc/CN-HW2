import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import Jama.CholeskyDecomposition;
import Jama.Matrix;

public class Homework2 {
	static int n = 3;
	static double epsilon = 0.000001d;
	static double[] d = new double[n];
	static double[] z = new double[n];
	static double[] y = new double[n];
	static double[] x = new double[n];
	static double[] b;// = { 1, 2, 3 };
	public static double[][] bb = new double[n][1];
	static double[][] A;// = { { 1, -1, 2 }, { -1, 5, -4 }, { 2, -4, 6 } };
	static double[][] Ainit;// = { { 1, -1, 2 }, { -1, 5, -4 }, { 2, -4, 6 } };
	// double[][] A = { { 1, 2.5, 3 }, { 2.5, 8.25, 15.5 }, { 3, 15.5, 43 }
	// };
	// { 1, -1, 2 }
	// { -1, 5, -4 }
	// { 2, -4, 6 }

	public static void main(String[] args) throws FileNotFoundException {
		getInput();

		System.out.println("A = ");
		printArray(A);
		boolean matrixIsValid = isMatrixValid(A, n);

		System.out.println("Det(A) = " + calcDeterminant(A, 3));
		System.out.println("Matrice valida: " + matrixIsValid);

		if (matrixIsValid) {
			// // D
			// for (int p = 0; p < n; ++p) {
			// double aux = 0;
			// for (int k = 0; k < p; ++k) {
			//// if (verifyDivision(d[k]))
			//// d[p] -= A[p][k] / A[p][k] * d[k];
			// aux += d[k] * Math.pow(A[p][k], 2);
			// }
			// d[p] = A[p][p] - aux;
			// if (d[p] == 0)
			// System.exit('d');
			// }

			// L
			for (int p = 0; p < n; ++p) {
				d[p] = A[p][p];
				for (int kk = 0; kk < p; ++kk)
					d[p] -= d[kk] * Math.pow(A[p][kk], 2);
				for (int i = p + 1; i < n; ++i) {
					for (int k = 0; k < p; ++k) {
						A[i][p] -= d[k] * A[i][k] * A[p][k];
					}
					if (verifyDivision(d[p]))
						A[i][p] /= d[p];
					A[p][i] = A[i][p];
				}
			}
			System.out.print("D = ");
			for (double it : d)
				System.out.print(it + " ");
			System.out.println();
			System.out.println("L = ");
			for (double it : A[0])
				System.out.print(it + " ");
			System.out.println();
			for (int i = 1; i < n; ++i) {
				for (double it : A[i])
					System.out.print(it + " ");
				System.out.println();
			}

			// det(A) = det(L)*det(D)*det(L^T)
			double detA = 1;
			for (int i = 0; i < n; ++i)
				detA *= d[i];
			System.out.println("Det (A) = " + detA);

			// Ax = b
			for (int i = 0; i < n; ++i) {
				y[i] = b[i];
				for (int j = 0; j < i; ++j) {
					if (i == j)
						y[i] -= 1 * y[j];
					else
						y[i] -= A[i][j] * y[j];
				}
			}
			System.out.print("y = ");
			printArray(y);

			for (int i = 0; i < n; ++i) {
				if (verifyDivision(d[i]))
					z[i] = y[i] / d[i];
			}
			System.out.print("z = ");
			printArray(z);

			for (int i = n - 1; i >= 0; --i) {
				x[i] = z[i];
				for (int j = i + 1; j < n; ++j) {
					if (i == j)
						x[i] -= 1 * x[j];
					else
						x[i] -= A[j][i] * x[j];
				}
			}
			System.out.print("x = ");
			printArray(x);

			// ||Ainit * xChol - b||2
			double norm = 0;
			double[] xNorm;

			xNorm = multiplyMatrixWithVector(Ainit, x);
			xNorm = substractVectors(xNorm, b);

			for (int i = 0; i < n; ++i) {
				norm += Math.pow(Math.abs(xNorm[i]), 2);
			}
			norm = Math.sqrt(norm);
			System.out.println("Norma = " + norm);
			System.out.println("xNorm = ");
			printArray(xNorm);

			for (int i = 0; i < n; ++i)
				bb[i][0] = b[i];

			jamaSol();
		}
	}

	public static void jamaSol() {
		Matrix ma = new Matrix(Ainit);
		Matrix mb = new Matrix(bb);
		Matrix mx = ma.solve(mb);
		Matrix residual = ma.times(mx).minus(mb);

		CholeskyDecomposition mChol = ma.chol();

		Matrix mL = mChol.getL();

		mL.print(3, 6);

		mx.print(3, 6);
	}

	public static double[] substractVectors(double[] a, double[] b) {
		double[] r = new double[a.length];
		for (int i = 0; i < a.length; ++i)
			r[i] = a[i] - b[i];
		return r;
	}

	public static void printArray(double[] a) {
		for (double it : a)
			System.out.print(it + " ");
		System.out.println();
	}

	public static void printArray(double[][] a) {
		for (int i = 0; i < a.length; ++i) {
			for (double it : a[i])
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

	public static double[] multiplyMatrixWithVector(double[][] matrix, double[] vector) {
		int rows = matrix.length;
		int columns = matrix[0].length;

		double[] result = new double[rows];

		for (int row = 0; row < rows; row++) {
			double sum = 0;
			for (int column = 0; column < columns; column++) {
				sum += matrix[row][column] * vector[column];
			}
			result[row] = sum;
		}
		return result;
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

	public static void getInput() throws FileNotFoundException {
		Scanner sc = new Scanner(new File("input.txt"));
		epsilon = sc.nextDouble();
		n = sc.nextInt();

		A = new double[n][n];
		Ainit = new double[n][n];
		b = new double[n];

		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				Ainit[i][j] = A[i][j] = sc.nextDouble();

		for (int i = 0; i < n; ++i)
			b[i] = sc.nextDouble();

		sc.close();
	}

	public static boolean verifyDivision(double x) {
		if (Math.abs(x) < epsilon)
			throw new ArithmeticException("Division by zero (custom)");

		return true;
	}
}
