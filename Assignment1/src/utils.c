#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#define SIZET_MAX	((size_t)(-1))

#include "utils.h"

/* PROTOTYPES */
static void print_row(const double *row, size_t n);
static double random_double_in_range(double max_absvalue, double precision);
/* END PROTOTYPES */

static void print_row(const double *row, size_t n)
{
	size_t i;

	printf("[");

	for (i = 0; i < n; i++) {
		printf("%f", row[i]);

		if (i != n - 1)
			printf(", ");
	}

	printf("]");
}

static double random_double_in_range(double max_absvalue, double precision)
{
	int nbuckets = (int)(2.0 * max_absvalue / precision + 1.0);
	int middle = nbuckets / 2;
	int random;

	random = (rand() % nbuckets) - middle;
	return random * precision;
}

void exit_with_error(const char *errmsg)
{
	fprintf(stderr, "%s\n", errmsg);
	exit(EXIT_FAILURE);
}

void *malloc_or_fail(size_t count, size_t size)
{
	void *ptr;

	if (count == 0 || size == 0)
		exit_with_error("Invalid arguments for malloc_or_fail.");
	else if (count > SIZET_MAX / size)
		exit_with_error("Count too large for malloc_or_fail.");

	ptr = malloc(count * size);

	if (ptr == NULL) {
		perror("malloc");
		exit(EXIT_FAILURE);
	}

	return ptr;
}

struct Vector *Vector_new(size_t n)
{
	struct Vector *v;

	v = malloc_or_fail(1, sizeof *v);
	v->entries = malloc_or_fail(n, sizeof *(v->entries));
	v->n = n;

	return v;
}

void Vector_delete(struct Vector *v)
{
	free(v->entries);
	free(v);
}

struct Matrix *Matrix_new(size_t m, size_t n)
{
	struct Matrix *M;
	size_t i;

	M = malloc_or_fail(1, sizeof *M);
	M->entries = malloc_or_fail(m, sizeof *(M->entries));

	for (i = 0; i < m; i++)
		M->entries[i] = malloc_or_fail(n, sizeof *(M->entries[i]));

	M->m = m;
	M->n = n;

	return M;
}

void Matrix_delete(struct Matrix *M)
{
	size_t i;

	for (i = 0; i < M->m; i++)
		free(M->entries[i]);

	free(M->entries);
	free(M);
}

void Vector_print(const struct Vector *v)
{
	print_row(v->entries, v->n);
}

void Matrix_print(const struct Matrix *M)
{
	size_t i;

	printf("[");

	for (i = 0; i < M->m; i++) {
		print_row(M->entries[i], M->n);

		if (i != M->m - 1)
			printf(", ");
	}

	printf("]");
}

struct Vector *Vector_copy(const struct Vector *v)
{
	struct Vector *cv;

	cv = Vector_new(v->n);
	memcpy(cv->entries, v->entries, (cv->n) * sizeof *(cv->entries));

	return cv;
}

struct Matrix *Matrix_zero(size_t m, size_t n)
{
	struct Matrix *M;
	size_t i, j;

	M = Matrix_new(m, n);

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			M->entries[i][j] = 0.0;
	}

	return M;
}

struct Matrix *Matrix_copy(const struct Matrix *M)
{
	struct Matrix *cM;
	size_t i;

	cM = Matrix_new(M->m, M->n);

	for (i = 0; i < cM->m; i++)
		memcpy(cM->entries[i], M->entries[i], (cM->n) * sizeof *(cM->entries[i]));

	return cM;
}

struct Matrix *Matrix_random(size_t m, size_t n, double range, double precision, enum MatrixPattern pattern)
{
	struct Matrix *M;
	size_t i, j;

	M = Matrix_new(m, n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			switch (pattern) {
				case MATRIX_PATTERN_LOWER_TRIANGULAR:
					if (j <= i)
						M->entries[i][j] = random_double_in_range(range, precision);
					else
						M->entries[i][j] = 0.0;

					break;
				case MATRIX_PATTERN_UPPER_TRIANGULAR:
					if (j <= i)
						M->entries[i][j] = 0.0;
					else
						M->entries[i][j] = random_double_in_range(range, precision);

					break;
				case MATRIX_PATTERN_SYMMETRIC:
					if (j <= i)
						M->entries[i][j] = random_double_in_range(range, precision);
					else
						M->entries[i][j] = M->entries[j][i];

					break;

				case MATRIX_PATTERN_NONE:
				default:
					M->entries[i][j] = random_double_in_range(range, precision);
					break;
			}
		}
	}

	return M;
}

struct Vector *Vector_random(size_t n, double range, double precision)
{
	struct Vector *v;
	size_t i;

	v = Vector_new(n);

	for (i = 0; i < n; i++)
		v->entries[i] = random_double_in_range(range, precision);

	return v;
}

struct Vector *Vector_matrix_multiply(const struct Matrix *A, const struct Vector *x)
{
	struct Vector *b;
	size_t i, j;

	if (A->n != x->n)
		exit_with_error("Dimensions of matrix and vector incompatible for multiplication.");

	b = Vector_new(A->m);

	for (i = 0; i < A->m; i++) {
		b->entries[i] = 0.0;

		for (j = 0; j < A->n; j++)
			b->entries[i] += (A->entries[i][j]) * (x->entries[j]);
	}

	return b;
}

struct Matrix *Matrix_multiply(const struct Matrix *A, const struct Matrix *B)
{
	struct Matrix *C;
	size_t i, j, k;

	if (A->n != B->m)
		exit_with_error("Dimensions of matrices incompatible for multiplication.");

	C = Matrix_new(A->m, B->n);

	for (i = 0; i < C->m; i++) {
		for (j = 0; j < C->n; j++) {
			C->entries[i][j] = 0.0;

			for (k = 0; k < A->n; k++)
				C->entries[i][j] += A->entries[i][k] * B->entries[k][j];
		}
	}

	return C;
}

struct Matrix *Matrix_transpose(const struct Matrix *A)
{
	struct Matrix *B;
	size_t i, j;

	B = Matrix_new(A->n, A->m);

	for (i = 0; i < B->m; i++) {
		for (j = 0; j < B->n; j++)
			B->entries[i][j] = A->entries[j][i];
	}

	return B;
}

int Matrix_is_symmetric(const struct Matrix *M, double precision)
{
	size_t i, j;
	size_t n;

	if (M->m != M->n)
		return -1;

	n = M->n;

	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			if (abs(M->entries[i][j] - M->entries[j][i]) > precision)
				return 0;
		}
	}

	return 1;
}

struct Vector *Vector_substract(const struct Vector *u, const struct Vector *v)
{
	struct Vector *result;
	size_t i;

	if (u->n != v->n)
		exit_with_error("Dimension of vectors incompatible for subtraction.\n");

	result = Vector_new(u->n);

	for (i = 0; i < u->n; i++)
		result->entries[i] = u->entries[i] - v->entries[i];

	return result;
}

int Vector_equal(const struct Vector *u, const struct Vector *v, double precision)
{
	size_t i;

	if (u->n != v->n)
		return 0;

	for (i = 0; i < u->n; i++) {
		if (abs(u->entries[i] - v->entries[i]) > precision)
			return 0;
	}

	return 1;
}
