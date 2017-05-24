#ifndef HAVE_LEASTSQUAREERROR_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "qrMethods.h"

typedef struct {
	int n;
	double *c;
	double *dc;
	double (*fitFun)(int i, double x);
} leastSquareError;

leastSquareError *leastSquareErrorAlloc(int n);
void leastSquareErrorFree(leastSquareError *ls);
void leastSquareErrorCalcConst(leastSquareError *ls, gsl_vector *x, gsl_vector *y, gsl_vector *dy);
double leastSquareErrorEval(leastSquareError *ls, double x);
double leastSquareErrorEvalMinus(leastSquareError *ls, double x);
double leastSquareErrorEvalPlus(leastSquareError *ls, double x);
void leastSquareErrorSingularValueDecomp(gsl_matrix *A, gsl_matrix *V, gsl_vector *S);
void leastSquareErrorSingularValueCalcConst(leastSquareError *ls, gsl_vector *x, gsl_vector *y, gsl_vector *dy);

#define HAVE_LEASTSQUAREERROR_H
#endif
