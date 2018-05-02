#ifndef  DECOMP_HPP_
#define DECOMP_HPP_

#define mdim 8 
#define  EPSILON  2.2e-16
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

int decomp(int n, int ndim,
	double *a, double *cond,
	int pivot[], int *flag);

int solve(int n, int ndim,
	double *a, double b[],
	int pivot[]);

#endif // ! DECOMP_HPP_
