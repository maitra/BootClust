#ifndef STATUTILS_H
#define STATUTILS_H

#include <stdio.h>		/* for fprintf() */
#include <stdlib.h>		/* for malloc() */
#include <math.h>		/* for sqrt() */
#include <assert.h>
#include "array.h"
#define MATHLIB_STANDALONE 
#include <Rmath.h>

int srswor(int n, int k, int *y);
void SampleVarCov(double **X, int N, int p, int K, double **c, int *ic1, int *nc, double **VarCov);
double trW(double *wss, int K);
double range(double *x, size_t n);
double sd(double *x, size_t n);
#endif /* STATUTILS_H */
