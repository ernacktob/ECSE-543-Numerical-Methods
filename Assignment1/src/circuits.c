#include <stdio.h>
#include <stdlib.h>

#include "circuits.h"
#include "cholesky.h"
#include "utils.h"

int circuits_parse_file(struct CircuitDescription *circuit, const char *filename)
{
	FILE *filePtr;
	struct Matrix *A, *Y;
	struct Vector *J, *E;
	size_t nnodes, nbranches;
	size_t i, j;
	int value;
	char c;
	double Jk, Rk, Ek;
	int result = -1;

	filePtr = fopen(filename, "r");

	if (filePtr == NULL) {
		perror("fopen");
		result = -1;
		goto cleanup_;
	}

	/* First row is nodes and branch count */
	if (fscanf(filePtr, "%lu %lu\n", &nnodes, &nbranches) != 2) {
		perror("fscanf");
		result = -1;
		goto cleanup_filePtr;
	}

	if (nnodes == 0 || nbranches == 0) {
		fprintf(stderr, "Node and branch counts cannot be zero.\n");
		result = -1;
		goto cleanup_filePtr;
	}

	A = Matrix_new(nnodes, nbranches);

	/* Read the reduced incidence matrix */
	for (i = 0; i < nnodes; i++) {
		for (j = 0; j < nbranches; j++) {
			if (fscanf(filePtr, "%d", &value) != 1) {
				perror("fscanf");
				result = -1;
				goto cleanup_A;
			}

			/* Check that each incidence matrix entry is -1, 0 or 1. */
			if (value != -1 && value != 0 && value != 1) {
				fprintf(stderr, "Incidence matrix can only have entries of -1, 0 or 1.\n");
				result = -1;
				goto cleanup_A;
			}

			if (fscanf(filePtr, "%c", &c) != 1) {
				perror("fscanf");
				result = -1;
				goto cleanup_A;
			}

			/* Format expects space between each branch for the same node.
			 * After the last branch for a given node, we expect a newline. */
			if (j == nbranches - 1) {
				if (c != '\n') {
					fprintf(stderr, "Expected \\n after end of incidence matrix row.\n");
					result = -1;
					goto cleanup_A;
				}
			} else {
				if (c != ' ') {
					fprintf(stderr, "Expected space after entry of incidence matrix.\n");
					result = -1;
					goto cleanup_A;
				}
			}

			A->entries[i][j] = (double)value;
		}
	}

	Y = Matrix_zero(nbranches, nbranches);
	J = Vector_new(nbranches);
	E = Vector_new(nbranches);

	/* Read the branches (current, resistance, voltage) */
	for (j = 0; j < nbranches; j++) {
		/* Read branch current, resistance and voltage */
		if (fscanf(filePtr, "%lf %lf %lf", &Jk, &Rk, &Ek) != 3) {
			perror("fscanf");
			result = -1;
			goto cleanup_EJY;
		}

		if (fscanf(filePtr, "%c", &c) != 1) {
			perror("fscanf");
			result = -1;
			goto cleanup_EJY;
		}

		/* We expect each branch to be on separate line. */
		if (c != '\n') {
			fprintf(stderr, "Expected \\n after the branch entry.\n");
			result = -1;
			goto cleanup_EJY;
		}

		/* We expect each branch to contain a nonzero resistance. */
		if (Rk == 0.0) {
			fprintf(stderr, "Branch with zero resistance is not supported.\n");
			result = -1;
			goto cleanup_EJY;
		}

		J->entries[j] = Jk;
		Y->entries[j][j] = 1.0 / Rk;
		E->entries[j] = Ek;
	}

	circuit->A = A;
	circuit->J = J;
	circuit->Y = Y;
	circuit->E = E;

	result = 0;
	goto cleanup_filePtr;

cleanup_EJY:
	Vector_delete(E);
	Vector_delete(J);
	Matrix_delete(Y);
cleanup_A:
	Matrix_delete(A);
cleanup_filePtr:
	fclose(filePtr);
cleanup_:
	return result;
}

struct Vector *circuits_solve_voltages(const struct CircuitDescription *circuit)
{
	struct Matrix *M, *Atranspose, *YAtranspose;
	struct Vector *b, *YE, *JminusYE;
	struct Vector *V;

	/* Compute M = AYA^T, which is the matrix that is obtained from KCL. */
	Atranspose = Matrix_transpose(circuit->A);
	YAtranspose = Matrix_multiply(circuit->Y, Atranspose);
	M = Matrix_multiply(circuit->A, YAtranspose);

	Matrix_delete(Atranspose);
	Matrix_delete(YAtranspose);

	/* Compute b = A(J - YE), which is the vector of source currents from KCL. */
	YE = Vector_matrix_multiply(circuit->Y, circuit->E);
	JminusYE = Vector_substract(circuit->J, YE);
	b = Vector_matrix_multiply(circuit->A, JminusYE);

	Vector_delete(YE);
	Vector_delete(JminusYE);

	/* Solve the system (AYA^T)V = A(J - YE) for the node voltages V. */
	if (cholesky_solve_system(&V, M, b, NULL) != 0)
		exit_with_error("The matrix AYA^T was not symmetric positive-definite.");

	Vector_delete(b);
	Matrix_delete(M);

	return V;
}

struct Vector *circuits_solve_voltages_banded(const struct CircuitDescription *circuit, size_t hb)
{
	struct Matrix *M, *Atranspose, *YAtranspose;
	struct Vector *b, *YE, *JminusYE;
	struct Vector *V;

	/* Compute M = AYA^T, which is the matrix that is obtained from KCL. */
	Atranspose = Matrix_transpose(circuit->A);
	YAtranspose = Matrix_multiply(circuit->Y, Atranspose);
	M = Matrix_multiply(circuit->A, YAtranspose);

	Matrix_delete(Atranspose);
	Matrix_delete(YAtranspose);

	/* Compute b = A(J - YE), which is the vector of source currents from KCL. */
	YE = Vector_matrix_multiply(circuit->Y, circuit->E);
	JminusYE = Vector_substract(circuit->J, YE);
	b = Vector_matrix_multiply(circuit->A, JminusYE);

	Vector_delete(YE);
	Vector_delete(JminusYE);

	/* Solve the system (AYA^T)V = A(J - YE) for the node voltages V. */
	if (cholesky_solve_system_banded(&V, M, b, NULL, hb) != 0)
		exit_with_error("The matrix AYA^T was not symmetric positive-definite.");

	Vector_delete(b);
	Matrix_delete(M);

	return V;
}

void circuits_destroy(struct CircuitDescription *circuit)
{
	Vector_delete(circuit->E);
	Matrix_delete(circuit->Y);
	Vector_delete(circuit->J);
	Matrix_delete(circuit->A);
}
