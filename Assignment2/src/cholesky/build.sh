#!/bin/sh

gcc -O2 -Wall -Wextra -pedantic cholesky.c utils.c test_cholesky.c -o test_cholesky
