#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "leastSquareMethods.h"


double fitFun(int i, double x) {
	switch(i) {
	case 0: return 1.0/x;	break;
	case 1: return 1.0;	break;
	case 2: return x;	break;
	default: {fprintf(stderr,"fitFun: invalide i: %d\n",i); return NAN;}
	}
}

int main(void) {
	int n = 10;
	double xTemp[] = {0.100, 0.145, 0.211, 0.307, 0.447, 0.649, 0.944, 1.372, 1.995, 2.900};
	double yTemp[] = {12.644, 9.235, 7.377, 6.460, 5.555, 5.896, 5.673, 6.964, 8.896, 11.355};
	double dyTemp[] = {0.858, 0.359, 0.505, 0.403, 0.683, 0.605, 0.856, 0.351, 1.083, 1.002};

	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector *y = gsl_vector_alloc(n);
	gsl_vector *dy = gsl_vector_alloc(n);

	printf("# x \t y \t dy\n");
	for(int i = 0; i < n; i++) {
		gsl_vector_set(x,i,xTemp[i]);
		printf("%.3g \t",xTemp[i]);
		gsl_vector_set(y,i,yTemp[i]);
		printf("%.3g \t",yTemp[i]);
		gsl_vector_set(dy,i,dyTemp[i]);
		printf("%.3g \n",dyTemp[i]);
	}
	printf("\n\n");

	int noOfFuns = 3;
	leastSquare *ls = leastSquareAlloc(noOfFuns);
	(*ls).fitFun = &fitFun;
	leastSquareCalcConst(ls,x,y,dy);

	double stepSize = 0.05;
	double start = xTemp[0]-stepSize, end = xTemp[n-1]+stepSize;
	printf("# xFit \t yFit\n");
	for(double i = start; i < end; i += stepSize) {
		printf("%g \t %g\n",i,leastSquareEval(ls,i));
	}

	leastSquareFree(ls);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);

	return 0;
}
