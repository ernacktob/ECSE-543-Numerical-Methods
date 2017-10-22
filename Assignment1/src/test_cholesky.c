#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>

#include "cholesky.h"
#include "utils.h"

#define PRECISION	0.000000001
#define RESOLUTION	0.1
#define RANGE_MAX	100.0

#define NTRIALS		10000000

enum TestResult {
	TEST_SUCCESS = 0,	/* Test passed */
	TEST_NOTSPD,		/* Generated matrix was not symmetric positive-definite (round off errors) */
	TEST_WRONGSOL		/* Solution obtained was wrong */
};

static struct Matrix *random_nonsingular_lower_triangular(size_t m, size_t n, double range, double resolution)
{
	struct Matrix *M;
	size_t i, j;

	M = Matrix_new(m, n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (j <= i)
				M->entries[i][j] = random_double_in_range(range, resolution);
			else
				M->entries[i][j] = 0.0;

			/* Avoid generating a singular matrix (with zero entry in diagonal) */
			if (i == j && M->entries[i][j] == 0.0)
				M->entries[i][j] += resolution;
		}
	}

	return M;
}


static enum TestResult test_solver(void)
{
	size_t sizes[] = {2, 3, 4, 5};
	struct Matrix *L, *Ltranspose, *A;
	struct Vector *b, *x, *found_x;
	struct Matrix *found_L;
	size_t n;
	enum TestResult result;

	n = sizes[rand() % (sizeof sizes / sizeof sizes[0])];
	L = random_nonsingular_lower_triangular(n, n, RANGE_MAX, RESOLUTION);
	x = Vector_random(n, RANGE_MAX, RESOLUTION);

	Ltranspose = Matrix_transpose(L);
	A = Matrix_multiply(L, Ltranspose);
	b = Vector_matrix_multiply(A, x);

	if (cholesky_solve_system(&found_x, A, b, &found_L) != 0) {
		/* Debugging information */
		printf("Matrix A was not symmetric positive definite, or round-off error was introduced.\n");
		printf("A = \n");
		Matrix_print(A);
		printf("\n");
		printf("L = \n");
		Matrix_print(L);
		printf("\n\n");
		result = TEST_NOTSPD;
		goto cleanup_;
	}

	if (!Vector_equal(found_x, x, PRECISION)) {
		/* Debugging information */
		printf("Wrong solution for system.\n");
		printf("found_x = \n");
		Vector_print(found_x);
		printf("\n");
		printf("x = \n");
		Vector_print(x);
		printf("\n");
		printf("A = \n");
		Matrix_print(A);
		printf("\n");
		printf("L = \n");
		Matrix_print(L);
		printf("\n");
		printf("found_L = \n");
		Matrix_print(found_L);
		printf("\n\n");
		result = TEST_WRONGSOL;
		goto cleanup_found_L;
	}

	result = TEST_SUCCESS;

cleanup_found_L:
	Matrix_delete(found_L);
	Vector_delete(found_x);
cleanup_:
	Vector_delete(b);
	Matrix_delete(A);
	Matrix_delete(Ltranspose);
	Matrix_delete(L);
	Vector_delete(x);

	return result;
}

int simple_test(void)
{
	const double testA1[][2] = {{1.0, -1.0}, {-1.0, 5.0}};
	const double testb1[] = {-1.0, 9.0};

	const double testA2[][3] = {{1.0, 0.0, -1.0}, {0.0, 4.0, 0.0}, {-1.0, 0.0, 2.0}};
	const double testb2[] = {-2.0, 8.0, 5.0};

	const double testA3[][4] = {{4.0, 4.0, 4.0, 4.0}, {4.0, 5.0, 5.0, 5.0}, {4.0, 5.0, 14.0, 14.0}, {4.0, 5.0, 14.0, 30.0}};
	const double testb3[] = {40.0, 49.0, 112.0, 176.0};

	struct Matrix *A1, *A2, *A3;
	struct Vector *b1, *b2, *b3;
	struct Matrix *L1, *L2, *L3;
	struct Vector *x1, *x2, *x3;

	size_t i, j;
	int result;

	A1 = Matrix_new(2, 2);
	b1 = Vector_new(2);

	A2 = Matrix_new(3, 3);
	b2 = Vector_new(3);

	A3 = Matrix_new(4, 4);
	b3 = Vector_new(4);

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++)
			A1->entries[i][j] = testA1[i][j];

		b1->entries[i] = testb1[i];
	}

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++)
			A2->entries[i][j] = testA2[i][j];

		b2->entries[i] = testb2[i];
	}

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++)
			A3->entries[i][j] = testA3[i][j];

		b3->entries[i] = testb3[i];
	}

	if (cholesky_solve_system(&x1, A1, b1, &L1) != 0) {
		printf("Failed to solve system.\n");
		result = -1;
		goto cleanup;
	}

	if (cholesky_solve_system(&x2, A2, b2, &L2) != 0) {
		printf("Failed to solve system.\n");
		result = -1;
		goto cleanup;
	}

	if (cholesky_solve_system(&x3, A3, b3, &L3) != 0) {
		printf("Failed to solve system.\n");
		result = -1;
		goto cleanup;
	}

	printf("L1 = ");
	Matrix_print(L1);
	printf("\n");
	printf("x1 = ");
	Vector_print(x1);
	printf("\n\n");

	printf("L2 = ");
	Matrix_print(L2);
	printf("\n");
	printf("x2 = ");
	Vector_print(x2);
	printf("\n\n");

	printf("L3 = ");
	Matrix_print(L3);
	printf("\n");
	printf("x3 = ");
	Vector_print(x3);
	printf("\n\n");

	result = 0;

cleanup:
	Matrix_delete(A1);
	Matrix_delete(A2);
	Matrix_delete(A3);
	Vector_delete(b1);
	Vector_delete(b2);
	Vector_delete(b3);
	Matrix_delete(L1);
	Matrix_delete(L2);
	Matrix_delete(L3);
	Vector_delete(x1);
	Vector_delete(x2);
	Vector_delete(x3);

	return result;
}

int main(void)
{
	int success_count = 0;
	int notspd_count = 0;
	int wrongsol_count = 0;
	int i;

	srand(time(NULL));

	if (simple_test() != 0) {
		printf("Failed the simple test... :(\n");
		return 0;
	}

	for (i = 0; i < NTRIALS; i++) {
		switch (test_solver()) {
			case TEST_SUCCESS:
				++success_count;
				break;
			case TEST_NOTSPD:
				++notspd_count;
				break;
			case TEST_WRONGSOL:
				++wrongsol_count;
				break;
		}
	}

	printf("Success rate:\t\t\t\t%d/%d\n", success_count, NTRIALS);
	printf("Not symmetric positive-definite rate:\t%d/%d\n", notspd_count, NTRIALS);
	printf("Wrong solution rate:\t\t\t%d/%d\n", wrongsol_count, NTRIALS);

	return 0;
}
