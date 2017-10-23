#include <stdio.h>
#include <stdlib.h>

#define OUTER_L	0.2
#define INNER_L 0.04
#define INNER_W 0.08

#define INNER_V 15.0
#define OUTER_V 0.0

double **alloc_grid(size_t N)
{
	double **phi;
	size_t i;

	phi = malloc(N * sizeof *phi);

	if (phi == NULL) {
		perror("malloc");
		exit(-1);
	}

	for (i = 0; i < N; i++) {
		phi[i] = malloc(N * sizeof *phi[i]);

		if (phi[i] == NULL) {
			perror("malloc");
			exit(-1);
		}
	}

	return phi;
}

void free_grid(double **phi, size_t N)
{
	size_t i;

	for (i = 0; i < N; i++)
		free(phi[i]);

	free(phi);
}

void print_grid(double **phi, size_t N)
{
	size_t i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			printf("%f", phi[i][j]);

			if (j != N - 1)
				printf("\t");
		}

		printf("\n");
	}

	printf("\n");
}

unsigned int successive_over_relaxation(double ***phip, size_t *Np, double h, double w, double r)
{
	double **phi;
	unsigned int iterations;
	double max_r, current_r;
	size_t N;
	size_t inner_Nx, inner_Ny;
	size_t i, j;

	/* Need false boundary for Neumann condition at planes of symmetry. */
	N = (OUTER_L / 2.0) / h + 2.0;
	phi = alloc_grid(N);

	/* Nodes at which we have the inner boundary */
	inner_Nx = ((OUTER_L - INNER_W) / 2.0) / h;
	inner_Ny = ((OUTER_L - INNER_L) / 2.0) / h;

	/* Initialize everything to 0. */
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			phi[i][j] = 0.0;
	}

	/* Set boundary conditions on inner and outer conductors */
	for (i = 0; i < N; i++)
		phi[i][0] = OUTER_V;

	for (j = 0; j < N; j++)
		phi[0][j] = OUTER_V;

	for (i = inner_Nx; i < N; i++)
		phi[i][inner_Ny] = INNER_V;

	for (j = inner_Ny; j < N; j++)
		phi[inner_Nx][j] = INNER_V;

	iterations = 0;

	do {
		max_r = 0.0;

		for (i = 1; i < inner_Nx; i++) {
			/* Apply Laplace equation */
			for (j = 1; j < N - 1; j++)
				phi[i][j] = (1.0 - w) * phi[i][j] + (w / 4.0) * (phi[i - 1][j] + phi[i][j - 1] + phi[i + 1][j] + phi[i][j + 1]);

			/* Apply Neumann boundary condition at false boundary */
			phi[i][N - 1] = phi[i][N - 3];
		}

		for (i = inner_Nx; i < N - 1; i++) {
			for (j = 1; j < inner_Ny; j++)
				phi[i][j] = (1.0 - w) * phi[i][j] + (w / 4.0) * (phi[i - 1][j] + phi[i][j - 1] + phi[i + 1][j] + phi[i][j + 1]);
		}

		/* Neumann boundary on half plane */
		for (j = 1; j < inner_Ny; j++)
			phi[N - 1][j] = phi[N - 3][j];

		/* Compute residuals at each node */
		for (i = 1; i < inner_Nx; i++) {
			for (j = 1; j < N - 1; j++) {
				current_r = phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1] - 4.0 * phi[i][j];

				if (current_r > max_r)
					max_r = current_r;
			}
		}

		for (i = inner_Nx; i < N - 1; i++) {
			for (j = 1; j < inner_Ny; j++) {
				current_r = phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1] - 4.0 * phi[i][j];

				if (current_r > max_r)
					max_r = current_r;
			}
		}

		++iterations;
	} while (max_r >= r);

	*phip = phi;
	*Np = N;

	return iterations;
}

int main(void)
{
	double **phi;
	size_t N;
	unsigned int iterations;
	double x = 0.06, y = 0.04;
	double h = 0.02;
	double r = 1e-5;
	double w;
	int i;

	printf("w parameter\titerations\tphi at (%f, %f)\n", x, y);

	for (i = 0; i < 10; i++) {
		w = 1.0 + 0.1 * i;
		iterations = successive_over_relaxation(&phi, &N, 0.02, w, r);
		printf("%f\t%u\t\t%f\n", w, iterations, phi[(int)(x / h)][(int)(y / h)]);
		printf("\n");
		free_grid(phi, N);
	}

	return 0;
}
