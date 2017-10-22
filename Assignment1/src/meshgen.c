#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/errno.h>

void generate_input_file(size_t N)
{
	size_t nnodes, nbranches;
	size_t node, branch;
	size_t i, j;
	size_t bijN, bijS, bijE, bijW; /* North, South, East, West */

	nnodes = 2 * N * N;
	nbranches = 2 * N * (N - 1) + (2 * N - 1) * N;

	/* Don't print node 0 */
	printf("%lu %lu\n", nnodes - 1, nbranches + 1);

	/* Print the reduced incidence matrix (skip row corresponding to node 0) */
	for (node = 1; node < nnodes; node++) {
		i = node / N;
		j = node % N;

		/* -1 means there is no branch in that direction. */
		bijN = bijS = bijE = bijW = -1; 

		if (i > 0)
			bijS = (i - 1) * (2 * N - 1) + (N - 1) + j;

		if (j > 0)
			bijW = i * (2 * N - 1) + j - 1;

		if (j < N - 1)
			bijE = i * (2 * N - 1) + j;

		if (i < 2 * N - 1)
			bijN = i * (2 * N - 1) + (N - 1) + j;

		for (branch = 0; branch < nbranches; branch++) {
			if (branch == bijS || branch == bijW)
				printf("-1");
			else if (branch == bijE || branch == bijN)
				printf("1");
			else
				printf("0");

			printf(" ");
		}

		/* Add the external branch for voltage source */
		if (node == nnodes - 1)
			printf("-1");
		else
			printf("0");

		printf("\n");
	}

	for (branch = 0; branch < nbranches; branch++)
		printf("0 1000 0\n");

	/* Add voltage source with series 1k resistor */
	printf("0 1000 1\n");
}

int main(int argc, const char *argv[])
{
	unsigned long temp;
	size_t N;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s <N>\n", argv[0]);
		return 0;
	}

	temp = strtoul(argv[1], NULL, 10);

	if (temp == ULONG_MAX || (temp == 0 && errno == EINVAL)) {
		perror("strtoul");
		return -1;
	}

	N = (size_t)temp;
	generate_input_file(N);

	return 0;
}
