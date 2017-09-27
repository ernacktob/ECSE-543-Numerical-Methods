#!/bin/sh

gcc -O2 -Wall -Wextra -pedantic -std=c89 -Iinclude src/test_cholesky.c src/utils.c src/cholesky.c -o bin/test_cholesky
gcc -O2 -Wall -Wextra -pedantic -std=c89 -Iinclude src/circuit_solver.c src/utils.c src/cholesky.c -o bin/circuit_solver
