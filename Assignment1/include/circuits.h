#ifndef CIRCUITS_H
#define CIRCUITS_H

#include "utils.h"

/* Contains the description of a circuit in terms
 * of the reduced incidence matrix A, the conductance
 * matrix Y, the current vector J and the voltage
 * vector E. */
struct CircuitDescription {
	struct Matrix *A;
	struct Matrix *Y;
	struct Vector *J;
	struct Vector *E;
};

/* Fill out a CircuitDescription from an input file. */
int circuits_parse_file(struct CircuitDescription *circuit, const char *filename);

/* Solve for the node voltages in a circuit described by CircuitDescription. */
struct Vector *circuits_solve_voltages(const struct CircuitDescription *circuit);
struct Vector *circuits_solve_voltages_banded(const struct CircuitDescription *circuit, size_t hb);

/* Release memory allocated internally for CircuitDescription. */
void circuits_destroy(struct CircuitDescription *circuit);

#endif
