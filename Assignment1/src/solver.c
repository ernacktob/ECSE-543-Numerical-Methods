#include <stdio.h>

#include "circuits.h"
#include "utils.h"

int main(int argc, const char *argv[])
{
	struct CircuitDescription circuit;
	struct Vector *V;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
		return 0;
	}

	if (circuits_parse_file(&circuit, argv[1]) != 0) {
		fprintf(stderr, "Failed to parse circuit file.\n");
		return -1;
	}

	V = circuits_solve_voltages(&circuit);

	printf("V = ");
	Vector_print(V);
	printf("\n");

	Vector_delete(V);
	circuits_destroy(&circuit);

	return 0;
}
