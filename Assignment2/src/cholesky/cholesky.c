#include <stddef.h>
#include <math.h>

#include "utils.h"

/* Initialize upper triangle values of L to zero. */
static void zero_upper_triangle(struct Matrix *L)
{
	size_t i, j;

	for (i = 0; i < L->m; i++) {
		for (j = i + 1; j < L->n; j++)
			L->entries[i][j] = 0.0;
	}
}

/* Cholesky decomposition
 * 
 * Decomposes an n x n real symmetric positive-definite matrix
 * into its Cholesky decomposition L*L^T where L is lower-triangular.
 *
 * The matrix argument A is overwritten with the result of L.
 * Only the lower half of the matrix is touched, the upper-half
 * remains undefined.
 *
 * Parameters:
 * n - dimension of matrix
 * A - matrix to decompose
 *
 * Returns:
 * 0 if operation was successful
 * -1 if the matrix A is not positive-definite.
 */
static int cholesky_decomposition(double **A, size_t n)
{
	double **L = A;	/* The result overwrites lower half of A */
	size_t i, j, k;

	for (j = 0; j < n; j++) {
		/* Check that the matrix is positive-definite */
		if (A[j][j] <= 0.0)
			return -1;

		L[j][j] = sqrt(A[j][j]);

		/* Check again that L[j][j] > 0 in case that
		 * the sqrt introduced round-off errors... */
		if (L[j][j] <= 0.0)
			return -1;

		for (i = j + 1; i < n; i++) {
			L[i][j] = A[i][j] / L[j][j];

			for (k = j + 1; k <= i; k++)
				A[i][k] = A[i][k] - L[i][j] * L[k][j];
		}
	}

	return 0;
}

/* Cholesky decomposition banded
 * 
 * Decomposes an n x n real symmetric positive-definite matrix
 * into its Cholesky decomposition L*L^T where L is lower-triangular.
 *
 * This version exploits the banded nature of A to speed up computation.
 *
 * The matrix argument A is overwritten with the result of L.
 * Only the lower half of the matrix is touched, the upper-half
 * remains undefined.
 *
 * Parameters:
 * n - dimension of matrix
 * A - matrix to decompose
 * b - half-bandwidth of matrix
 *
 * Returns:
 * 0 if operation was successful
 * -1 if the matrix A is not positive-definite.
 */
static int cholesky_decomposition_banded(double **A, size_t n, size_t b)
{
	double **L = A;	/* The result overwrites lower half of A */
	size_t i, j, k;

	for (j = 0; j < n; j++) {
		/* Check that the matrix is positive-definite */
		if (A[j][j] <= 0.0)
			return -1;

		L[j][j] = sqrt(A[j][j]);

		/* Check again that L[j][j] > 0 in case that
		 * the sqrt introduced round-off errors... */
		if (L[j][j] <= 0.0)
			return -1;

		/* Can skip most of the entries after the half bandwidth because they are 0. */
		for (i = j + 1; i < j + b && i < n; i++) {
			L[i][j] = A[i][j] / L[j][j];

			for (k = j + 1; k <= i; k++)
				A[i][k] = A[i][k] - L[i][j] * L[k][j];
		}
	}

	return 0;
}

/* Forward elimination
 *
 * Performs forward elimination to solve the equation Ly = b,
 * where L is a lower-triangular matrix.
 *
 * The vector argument b is overwritten to the solution vector y,
 * and only the lower-half of the matrix L is touched.
 *
 * Parameters:
 * n - dimension of matrix / vector
 * L - lower-triangular matrix
 * b - vector in the equation Ly = b
 */
static void forward_elimination(double *b, double **L, size_t n)
{
	double *y = b;	/* The result overwrites b */
	size_t i, j;

	for (j = 0; j < n; j++) {
		y[j] = b[j] / L[j][j];

		for (i = j + 1; i < n; i++)
			b[i] = b[i] - L[i][j] * y[j];
	}
}

/* Back substitution
 *
 * Performs substitution to solve the equation (L^T)x = y,
 * where L is a lower-triangular matrix.
 *
 * The vector argument y is overwritten to the solution vector x,
 * and only the lower-half of the matrix L is touched.
 *
 * Parameters:
 * n - dimension of matrix / vector
 * L - lower-triangular matrix
 * y - vector in the equation (L^T)x = y
 */
static void back_substitution(double *y, double **L, size_t n)
{
	double *x = y;	/* The result overwrites y */
	size_t i, j;
	size_t t;

	for (t = 0; t < n; t++) {
		i = n - t - 1;
		x[i] = y[i] / L[i][i];

		for (j = 0; j < i; j++)
			y[j] = y[j] - L[i][j] * x[i];
	}
}

/* See cholesky.h header for documentation */
int cholesky_solve_system(struct Vector **xp, const struct Matrix *A, const struct Vector *b, struct Matrix **Lp)
{
	struct Matrix *L;
	struct Vector *x;

	if (A->m != A->n)
		exit_with_error("Matrix A must be a square matrix.");

	if (b->n != A->m)
		exit_with_error("Matrix A and vector b not compatible for the system of equations.");

	if (!Matrix_is_symmetric(A))
		return -1;

	L = Matrix_copy(A);

	if (cholesky_decomposition(L->entries, L->n) != 0) {
		Matrix_delete(L);
		return -1;
	}

	x = Vector_copy(b);

	forward_elimination(x->entries, L->entries, x->n);
	back_substitution(x->entries, L->entries, x->n);

	if (Lp != NULL) {
		zero_upper_triangle(L);
		*Lp = L;
	} else {
		Matrix_delete(L);
	}

	*xp = x;

	return 0;
}

/* See cholesky.h header for documentation */
int cholesky_solve_system_banded(struct Vector **xp, const struct Matrix *A, const struct Vector *b, struct Matrix **Lp, size_t hb)
{
	struct Matrix *L;
	struct Vector *x;

	if (A->m != A->n)
		exit_with_error("Matrix A must be a square matrix.");

	if (b->n != A->m)
		exit_with_error("Matrix A and vector b not compatible for the system of equations.");

	if (!Matrix_is_symmetric(A))
		return -1;

	L = Matrix_copy(A);

	if (cholesky_decomposition_banded(L->entries, L->n, hb) != 0) {
		Matrix_delete(L);
		return -1;
	}

	x = Vector_copy(b);

	forward_elimination(x->entries, L->entries, x->n);
	back_substitution(x->entries, L->entries, x->n);

	if (Lp != NULL) {
		zero_upper_triangle(L);
		*Lp = L;
	} else {
		Matrix_delete(L);
	}

	*xp = x;

	return 0;
}
