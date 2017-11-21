#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>

/* utils.h
 * Generic utility functions for Vector and Matrix operations.
 */

struct Vector {
	double *entries;
	size_t n;
};

struct Matrix {
	double **entries;
	size_t m;
	size_t n;
};

enum MatrixPattern {
	MATRIX_PATTERN_NONE = 0,
	MATRIX_PATTERN_LOWER_TRIANGULAR,
	MATRIX_PATTERN_UPPER_TRIANGULAR,
	MATRIX_PATTERN_SYMMETRIC
};

/* Generate a random double between -max_absvalue and max_absvalue, as multiple of resolution. */
double random_double_in_range(double max_absvalue, double resolution);

/* Exit program with errormsg printed on stderr. */
void exit_with_error(const char *errmsg);

/* Attempt to malloc memory for 'count' elements of size 'size', or exit with error. */
void *malloc_or_fail(size_t count, size_t size);

/* Vector operations */
struct Vector *Vector_new(size_t n);
void Vector_delete(struct Vector *v);
void Vector_print(const struct Vector *v);
struct Vector *Vector_copy(const struct Vector *v);
struct Vector *Vector_random(size_t n, double range, double resolution);
struct Vector *Vector_matrix_multiply(const struct Matrix *A, const struct Vector *x);
struct Vector *Vector_substract(const struct Vector *u, const struct Vector *v);
int Vector_equal(const struct Vector *u, const struct Vector *v, double precision);

/* Matrix operations */
struct Matrix *Matrix_new(size_t m, size_t n);
void Matrix_delete(struct Matrix *M);
void Matrix_print(const struct Matrix *M);
struct Matrix *Matrix_zero(size_t m, size_t n);
struct Matrix *Matrix_copy(const struct Matrix *M);
struct Matrix *Matrix_random(size_t m, size_t n, double range, double resolution, enum MatrixPattern pattern);
struct Matrix *Matrix_multiply(const struct Matrix *A, const struct Matrix *B);
struct Matrix *Matrix_transpose(const struct Matrix *A);
int Matrix_is_symmetric(const struct Matrix *M);

#endif
