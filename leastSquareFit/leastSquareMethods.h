#ifndef HAVE_LEASTSQUARE_H

typedef struct {
	int n;
	double *c;
	double (*fitFun)(int i, double x);
} leastSquare;

leastSquare *leastSquareAlloc(int n);
void leastSquareFree(leastSquare *ls);
void leastSquareCalcConst(leastSquare *ls, gsl_vector *x, gsl_vector *y, gsl_vector *dy);
double leastSquareEval(leastSquare *ls, double x);

#define HAVE_LEASTSQUARE_H
#endif
