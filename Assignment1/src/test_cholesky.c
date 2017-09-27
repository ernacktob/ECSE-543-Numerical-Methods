#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>

#include "cholesky.h"
#include "utils.h"

#define PRECISION	0.00000001
#define RANGE_MAX	100.0

#define NTRIALS		100000

static int test_solver(void)
{
	size_t sizes[] = {2, 3, 4, 5};
	struct Matrix *L, *Ltranspose, *A;
	struct Vector *b, *x, *found_x;
	size_t n;
	int result;

	n = sizes[rand() % (sizeof sizes / sizeof sizes[0])];
	L = Matrix_random(n, n, RANGE_MAX, PRECISION, MATRIX_PATTERN_LOWER_TRIANGULAR);
	x = Vector_random(n, RANGE_MAX, PRECISION);

	Ltranspose = Matrix_transpose(L);
	A = Matrix_multiply(L, Ltranspose);
	b = Vector_matrix_multiply(A, x);

	printf("A = ");
	Matrix_print(A);
	printf("\n");

	printf("b = ");
	Vector_print(b);
	printf("\n");

	if (cholesky_solve_system(&found_x, A, b, PRECISION) != 0) {
		printf("Matrix A was not symmetric positive definite, or round-off error was introduced.\n");
		result = -1;
		goto cleanup_;
	}

	printf("L = ");
	Matrix_print(L);
	printf("\n");

	printf("x = ");
	Vector_print(x);
	printf("\n");

	if (!Vector_equal(found_x, x, PRECISION)) {
		printf("Wrong solution for system.\n");
		result = -1;
		goto cleanup_found_x;
	}

	result = 0;

cleanup_found_x:
	Vector_delete(found_x);
cleanup_:
	Vector_delete(b);
	Matrix_delete(A);
	Matrix_delete(Ltranspose);
	Matrix_delete(L);
	Vector_delete(x);

	return result;
}

int main(void)
{
	int success_count = 0;
	int i;

	srand(time(NULL));

	for (i = 0; i < NTRIALS; i++) {
		if (test_solver() == 0)
			++success_count;
	}

	printf("Success rate: %d/%d\n", success_count, NTRIALS);
	return 0;
}
