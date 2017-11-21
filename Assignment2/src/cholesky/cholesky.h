#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <stddef.h>

#include "utils.h"

/* Cholesky solve system
 *
 * Solve a system of equations represented by the matrix equation Ax = b,
 * where A is a real symmetric positive-definite matrix. The function
 * finds the Cholesky decomposition of A into L*L^T. It then solves the
 * equation Ly = b for y using forward elimination, and finally solves
 * the equation (L^T)x = y by back substitution.
 *
 * Parameters:
 * xp - pointer to solution vector that will be set after success
 * A - n x n real symmetric positive-definite matrix
 * b - n x 1 real vector in the equation Ax = b
 * Lp - if not NULL, stores the L matrix found during Cholesky decomposition.
 *
 * Returns:
 * 0 if operation successful
 * -1 if A is not symmetric positive-definite.
 *
 */
int cholesky_solve_system(struct Vector **xp, const struct Matrix *A, const struct Vector *b, struct Matrix **Lp);

/* This is the same except it also include half bandwidth hb to speed up computations for question 2. */
int cholesky_solve_system_banded(struct Vector **xp, const struct Matrix *A, const struct Vector *b, struct Matrix **Lp, size_t hb);

#endif
