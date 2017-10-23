#include <stdio.h>
#include <stdlib.h>

#include "circuits.h"
#include "utils.h"

int main(int argc, const char *argv[])
{
	struct CircuitDescription circuit;
	struct Vector *V;
	size_t N;
	double R;

	if (argc != 3) {
		fprintf(stderr, "Usage: %s <filename> <N>\n", argv[0]);
		return 0;
	}

	if (circuits_parse_file(&circuit, argv[1]) != 0) {
		fprintf(stderr, "Failed to parse circuit file.\n");
		return -1;
	}

	N = strtoul(argv[2], NULL, 10);
	V = circuits_solve_voltages_banded(&circuit, N + 1);

	R = (1000.0 * (V->entries[V->n - 1] / 1.0)) / (1.0 - (V->entries[V->n - 1] / 1.0));

	printf("Resistance of mesh: %f ohms.\n", R);

	Vector_delete(V);
	circuits_destroy(&circuit);

	return 0;
}
