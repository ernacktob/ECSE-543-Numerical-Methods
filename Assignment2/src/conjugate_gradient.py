import math # Used for sqrt function

def infinity_norm(x):
	return max(x)

def two_norm(x):
	return math.sqrt(sum([xi**2 for xi in x]))

def vector_add(x, y):
	return [xi + yi for xi, yi in zip(x, y)]

def vector_subtract(x, y):
	return [xi - yi for xi, yi in zip(x, y)]

def vector_scalar_multiply(a, x):
	return [a * xi for xi in x]

def scalar_product(x, y):
	return sum([xi * yi for xi, yi in zip(x, y)])

def matrix_vector_product(A, x):
	return [scalar_product(row, x) for row in A]

def matrix_inner_product(x, A, y):
	return scalar_product(x, matrix_vector_product(A, y))

# Use the non-preconditioned conjugate gradient method to find
# the solution to the equation Ax = b with initial guess x0.
def conjugate_gradient(A, b, x0):
	residuals = []
	x = x0
	r = vector_subtract(b, matrix_vector_product(A, x))
	p = r

	residuals.append(r)

	# Dimension
	n = len(x0)

	# We need at most n steps to converge to the solution
	for k in range(n):
		alpha = scalar_product(p, r) / matrix_inner_product(p, A, p)
		x = vector_add(x, vector_scalar_multiply(alpha, p))
		r = vector_subtract(b, matrix_vector_product(A, x))
		beta = -matrix_inner_product(p, A, r) / matrix_inner_product(p, A, p)
		p = vector_add(r, vector_scalar_multiply(beta, p))
		residuals.append(r)

	return x, residuals
