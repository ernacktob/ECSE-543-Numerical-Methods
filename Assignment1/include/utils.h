#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>

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

void exit_with_error(const char *errmsg);
void *malloc_or_fail(size_t count, size_t size);

struct Vector *Vector_new(size_t n);
void Vector_delete(struct Vector *v);
void Vector_print(const struct Vector *v);
struct Vector *Vector_copy(const struct Vector *v);
struct Vector *Vector_random(size_t n, double range, double precision);
struct Vector *Vector_matrix_multiply(const struct Matrix *A, const struct Vector *x);
struct Vector *Vector_substract(const struct Vector *u, const struct Vector *v);
int Vector_equal(const struct Vector *u, const struct Vector *v, double precision);

struct Matrix *Matrix_new(size_t m, size_t n);
void Matrix_delete(struct Matrix *M);
void Matrix_print(const struct Matrix *M);
struct Matrix *Matrix_zero(size_t m, size_t n);
struct Matrix *Matrix_copy(const struct Matrix *M);
struct Matrix *Matrix_random(size_t m, size_t n, double range, double precision, enum MatrixPattern pattern);
struct Matrix *Matrix_multiply(const struct Matrix *A, const struct Matrix *B);
struct Matrix *Matrix_transpose(const struct Matrix *A);
int Matrix_is_symmetric(const struct Matrix *M, double precision);

#endif
