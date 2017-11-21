import conjugate_gradient

# Geometry parameters
OUTER_HEIGHT = 0.2
OUTER_WIDTH = 0.2
INNER_HEIGHT = 0.04
INNER_WIDTH = 0.08

# Fixed potentials
OUTER_POTENTIAL = 0.0
INNER_POTENTIAL = 15.0

# Obtain mapping between node number and coordinates (as integer multiples of h)
def get_nodes(X1, X2, Y1, Y2):
	nodes = {}
	reverse_nodes = {}
	n = 0

	for i in range(1, X1 + X2):
		if i < X1 or i > X2:
			jrange = range(1, Y1 + Y2)
		else:
			jrange = list(range(1, Y1)) + list(range(Y2 + 1, Y1 + Y2))

		for j in jrange:
			nodes[n] = (i, j)
			reverse_nodes[(i, j)] = n
			n += 1

	return nodes, reverse_nodes

# Obtain the coefficients for the linear equation corresponding to the
# Laplace equation satisfied at the point (x, y).
def gen_row_equation(x, y, X1, X2, Y1, Y2, nodes, reverse_nodes):
	coeffs = [0.0 for i in range(len(nodes))]
	rhs = 0.0

	coeffs[reverse_nodes[(x, y)]] = 4.0

	# We check that the neighboring nodes are either
	# in the interior, or on the boundary. In the first
	# case the neighbor corresponds to an unknown variable,
	# otherwise it contributes to the right-hand-side.

	if (x + 1, y) in reverse_nodes:
		coeffs[reverse_nodes[(x + 1, y)]] = -1.0
	else:
		if x + 1 == X1:
			rhs += INNER_POTENTIAL
		else:
			rhs += OUTER_POTENTIAL

	if (x, y + 1) in reverse_nodes:
		coeffs[reverse_nodes[(x, y + 1)]] = -1.0
	else:
		if y + 1 == Y1:
			rhs += INNER_POTENTIAL
		else:
			rhs += OUTER_POTENTIAL

	if (x - 1, y) in reverse_nodes:
		coeffs[reverse_nodes[(x - 1, y)]] = -1.0
	else:
		if x - 1 == X2:
			rhs += INNER_POTENTIAL
		else:
			rhs += OUTER_POTENTIAL

	if (x, y - 1) in reverse_nodes:
		coeffs[reverse_nodes[(x, y - 1)]] = -1.0
	else:
		if y - 1 == Y2:
			rhs += INNER_POTENTIAL
		else:
			rhs += OUTER_POTENTIAL

	return coeffs, rhs

# Produce the matrix and fixed vector for the equation Ax = b
# associated with solving the Laplace finite-difference equations.
def gen_matrix_equation(h):
	# Corner coordinates of inner conductor boundary
	X1 = int(((OUTER_WIDTH - INNER_WIDTH) / 2.0) / h)
	Y1 = int(((OUTER_HEIGHT - INNER_HEIGHT) / 2.0) / h)
	X2 = int(((OUTER_WIDTH + INNER_WIDTH) / 2.0) / h)
	Y2 = int(((OUTER_HEIGHT + INNER_HEIGHT) / 2.0) / h)

	nodes, reverse_nodes = get_nodes(X1, X2, Y1, Y2)

	A = []
	b = []

	for i in sorted(nodes.keys()):
		x, y = nodes[i]
		coeffs, rhs = gen_row_equation(x, y, X1, X2, Y1, Y2, nodes, reverse_nodes)
		A.append(coeffs)
		b.append(rhs)

	return A, b, nodes, reverse_nodes

# Create the THE_MATRIX.h file for C usage for Cholesky decomposition
def CREATE_THE_MATRIX(A, b):
	f = open('cholesky/THE_MATRIX.h', 'w')
	f.write("#ifndef THE_MATRIX_H\n")
	f.write("#define THE_MATRIX_H\n")
	f.write("\n")
	f.write("#define DIMENSION %d\n"%len(b))
	f.write("\n")
	f.write("const double testA[DIMENSION][DIMENSION] = {\n")
	f.write(",\n".join(["\t{" + ", ".join(map(str, row)) + "}" for row in A]))
	f.write("};\n")
	f.write("\n")
	f.write("const double testb[DIMENSION] = {" + ", ".join(map(str, b)) + "};\n")
	f.write("\n")
	f.write("#endif\n")
	f.close()

h = 0.02
A, b, nodes, reverse_nodes = gen_matrix_equation(h)
x, residuals = conjugate_gradient.conjugate_gradient(A, b, [0.0 for i in range(len(b))])

#CREATE_THE_MATRIX(A, b)

for i in sorted(nodes.keys()):
	print "%d: (%f, %f) => %f"%(i, nodes[i][0] * h, nodes[i][1] * h, x[i])

#print ",".join(map(str, map(conjugate_gradient.infinity_norm, residuals)))
#print ",".join(map(str, map(conjugate_gradient.two_norm, residuals)))
