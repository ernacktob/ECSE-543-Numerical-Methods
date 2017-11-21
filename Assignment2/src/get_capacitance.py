import sys

# Dielectric permittivity of vacuum
eps0 = 8.854e-12

# Geometry parameters
OUTER_HEIGHT = 0.2
OUTER_WIDTH = 0.2
INNER_HEIGHT = 0.04
INNER_WIDTH = 0.08
ELEMENT_HEIGHT = 0.02
ELEMENT_WIDTH = 0.02

# Fixed potentials
OUTER_POTENTIAL = 0.0
INNER_POTENTIAL = 15.0

# Number of elements depending on presence of inner conductor boundary
NX1 = int(((OUTER_WIDTH - INNER_WIDTH) / 2.0) / ELEMENT_WIDTH)
NY1 = int((OUTER_HEIGHT / 2.0) / ELEMENT_HEIGHT)
NX2 = int((INNER_WIDTH / 2.0) / ELEMENT_WIDTH)
NY2 = int(((OUTER_HEIGHT - INNER_HEIGHT) / 2.0) / ELEMENT_HEIGHT)

# Compute the scalar product of two vectors
def scalar_product(x, y):
	return sum([xi * yi for xi, yi in zip(x, y)])

# Compute the product of a matrix times a vector
def matrix_vector_multiply(A, x):
	return [scalar_product(row, x) for row in A]

# Compute the energy/eps0 for a given square mesh with potentials U at the corners
# W = 1/2 (U^T)SU
def get_element_W(U):
	# Global S matrix obtained from Question 1 of this assignment
	S = [[1.0, -0.5, 0.0, -0.5], [-0.5, 1.0, -0.5, 0.0], [0.0, -0.5, 1.0, -0.5], [-0.5, 0.0, -0.5, 1.0]]
	SU = matrix_vector_multiply(S, U)
	UTSU = scalar_product(U, SU)
	W = 0.5 * UTSU

	return W

if len(sys.argv) != 2:
	print "Usage: %s <potentials_file>"%sys.argv[0]
	quit()

with open(sys.argv[1], 'r') as f:
	potentials_file = f.read()

nodes = {}
potentials = {}

# Read the nodes and potentials from file
for line in potentials_file.split("\n"):
	if len(line) == 0:
		continue

	entries = line.split(' ')
	node = int(entries[0])
	x = float(entries[1])
	y = float(entries[2])
	phi = float(entries[3])

	nodes[(x, y)] = node
	potentials[(x, y)] = phi

total_W = 0.0

# Calculate the total W by summing the individual W values of each square mesh.
for i in range(NX1 + NX2):
	if i < NX1:
		jrange = range(NY1)
	else:
		jrange = range(NY2)

	for j in jrange:
		vertices = [(i, j + 1), (i, j), (i + 1, j), (i + 1, j + 1)]
		vertices = [(x * ELEMENT_WIDTH, y * ELEMENT_HEIGHT) for x, y in vertices]
		pots = [potentials[(x, y)] for x, y in vertices]
		total_W += get_element_W(pots)

# Take into account the symmetry of the four quarters
total_W *= 4.0

# Voltage is potential difference between inner and outer conductor
V = INNER_POTENTIAL - OUTER_POTENTIAL

# Use equation eps0 * W = 1/2 CV^2
C = (eps0 * 2.0 * total_W) / (V ** 2)

print "Total energy per unit length: %e J/m"%(eps0 * total_W)
print "Total capacitance per unit length: %e F/m"%C
