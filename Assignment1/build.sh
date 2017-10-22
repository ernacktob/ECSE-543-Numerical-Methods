#!/bin/sh

mkdir -p bin
gcc -O2 -Wall -Wextra -pedantic -std=c89 -Iinclude src/test_cholesky.c src/utils.c src/cholesky.c -o bin/test_cholesky
gcc -O2 -Wall -Wextra -pedantic -std=c89 -Iinclude src/solver.c src/circuits.c src/utils.c src/cholesky.c -o bin/circuit_solver
gcc -O2 -Wall -Wextra -pedantic -std=c89 -Iinclude src/meshgen.c -o bin/meshgen
gcc -O2 -Wall -Wextra -pedantic -stc=c89 -Iinclude src/meshsolve.c src/circuits.c src/utils.c src/cholesky.c -o bin/meshsolve
