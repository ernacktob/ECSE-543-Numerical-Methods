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

labels = {}
n = 1

# Generate the node coordinates and labels using lexicographic ordering,
# taking into account the inner conductor's boundary
for i in range(NX1 + NX2 + 1):
	if i <= NX1:
		jrange = range(NY1 + 1)
	else:
		jrange = range(NY2 + 1)

	for j in jrange:
		labels[(i, j)] = n
		n += 1
		print "  %f     %f"%(i * ELEMENT_WIDTH, j * ELEMENT_HEIGHT)

print "/"

# Generate the mesh elements with zero source currents
for i in range(NX1 + NX2):
	if i < NX1:
		jrange = range(NY1)
	else:
		jrange = range(NY2)

	for j in jrange:
		vertices1 = [(i, j), (i + 1, j), (i + 1, j + 1)]
		vertices2 = [(i, j), (i + 1, j + 1), (i, j + 1)]

		print "  %-2d %-2d %-5d 0.0000"%(labels[vertices1[0]], labels[vertices1[1]], labels[vertices1[2]])
		print "  %-2d %-2d %-5d 0.0000"%(labels[vertices2[0]], labels[vertices2[1]], labels[vertices2[2]])

print "/"

# Generate the fixed potentials on outer and inner conductors
for i in range(NX1 + NX2 + 1):
	if i == 0:
		jrange1 = range(NY1 + 1)
	else:
		jrange1 = [0]

	if i < NX1:
		jrange2 = []
	elif i == NX1:
		jrange2 = range(NY2, NY1 + 1)
	else:
		jrange2 = [NY2]

	for j in jrange1:
		print "  %-5d %f"%(labels[(i, j)], OUTER_POTENTIAL)

	for j in jrange2:
		print "  %-5d %f"%(labels[(i, j)], INNER_POTENTIAL)

print "/"
