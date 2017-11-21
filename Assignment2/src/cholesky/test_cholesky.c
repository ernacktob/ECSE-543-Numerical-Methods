#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "cholesky.h"
#include "utils.h"

#include "THE_MATRIX.h"

int check_the_matrix(void)
{
	struct Matrix *A;
	struct Vector *b;
	struct Matrix *L;
	struct Vector *x;

	size_t i, j;
	int result;

	A = Matrix_new(DIMENSION, DIMENSION);
	b = Vector_new(DIMENSION);

	for (i = 0; i < DIMENSION; i++) {
		for (j = 0; j < DIMENSION; j++)
			A->entries[i][j] = testA[i][j];

		b->entries[i] = testb[i];
	}

	if (cholesky_solve_system(&x, A, b, &L) != 0) {
		printf("Failed to solve system.\n");
		result = -1;
		goto cleanup;
	}

	result = 0;

	printf("x = ");
	Vector_print(x);
	printf("\n");

cleanup:
	Matrix_delete(A);
	Vector_delete(b);
	Matrix_delete(L);
	Vector_delete(x);

	return result;
}

int main(void)
{
	if (check_the_matrix() == 0)
		printf("The matrix is positive-definite, yay!\n");
	else
		printf("The matrix was not positive-definite... :(\n");

	return 0;
}
