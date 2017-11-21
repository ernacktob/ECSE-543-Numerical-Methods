import conjugate_gradient

# Geometry parameters
OUTER_HEIGHT = 0.2
OUTER_WIDTH = 0.2
INNER_HEIGHT = 0.04
INNER_WIDTH = 0.08

# Fixed potentials
OUTER_POTENTIAL = 0.0
INNER_POTENTIAL = 15.0

def get_nodes(NX1, NX2, NY1, NY2):
	nodes = {}
	reverse_nodes = {}
	n = 0

	for i in range(1, NX1 + NX2 + 1):
		if i < NX1:
			jrange = range(1, NY1 + 1)
		else:
			jrange = range(1, NY2)

		for j in jrange:
			nodes[n] = (i, j)
			reverse_nodes[(i, j)] = n
			n += 1

	return nodes, reverse_nodes

def gen_row_equation(x, y, NX1, NX2, NY1, NY2, nodes, reverse_nodes):
	if x == 0.02 and y == 0.08:
		import pdb
		pdb.set_trace()

	coeffs = [0.0 for i in range(len(nodes))]
	rhs = 0.0

	coeffs[reverse_nodes[(x, y)]] = -4.0

	if (x + 1, y) in reverse_nodes:
		coeffs[reverse_nodes[(x + 1, y)]] += 1.0
	else:
		if x + 1 == NX1:
			rhs -= INNER_POTENTIAL
		else:
			coeffs[reverse_nodes[(x - 1, y)]] += 1.0

	if (x, y + 1) in reverse_nodes:
		coeffs[reverse_nodes[(x, y + 1)]] += 1.0
	else:
		if y + h == NY2:
			rhs -= INNER_POTENTIAL
		else:
			coeffs[reverse_nodes[(x, y - 1)]] += 1.0

	if (x - 1, y) in reverse_nodes:
		coeffs[reverse_nodes[(x - 1, y)]] += 1.0
	else:
		rhs -= OUTER_POTENTIAL

	if (x, y - 1) in reverse_nodes:
		coeffs[reverse_nodes[(x, y - 1)]] += 1.0
	else:
		rhs -= OUTER_POTENTIAL

	return coeffs, rhs

def gen_matrix_equation(h):
	# Number of elements depending on presence of inner conductor boundary
	NX1 = int(((OUTER_WIDTH - INNER_WIDTH) / 2.0) / h)
	NY1 = int((OUTER_HEIGHT / 2.0) / h)
	NX2 = int((INNER_WIDTH / 2.0) / h)
	NY2 = int(((OUTER_HEIGHT - INNER_HEIGHT) / 2.0) / h)

	nodes, reverse_nodes = get_nodes(NX1, NX2, NY1, NY2)

	A = []
	b = []

	for i in sorted(nodes.keys()):
		x, y = nodes[i]
		coeffs, rhs = gen_row_equation(x, y, NX1, NX2, NY1, NY2, nodes, reverse_nodes)
		A.append(coeffs)
		b.append(rhs)

	return A, b, nodes, reverse_nodes

h = 0.02
A, b, nodes, reverse_nodes = gen_matrix_equation(h)
print A
print b

x, residuals = conjugate_gradient.conjugate_gradient(A, b, [0.0 for i in range(len(b))])

for i in sorted(nodes.keys()):
	print "%d: (%f, %f) => %f"%(i, nodes[i][0] * h, nodes[i][1] * h, x[i])
