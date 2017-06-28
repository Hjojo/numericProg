#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "ODEMethods.h"

void f(double t, gsl_vector *y, gsl_vector *dydt) {

	double y0 = gsl_vector_get(y,0);

	gsl_vector_set(dydt,0,y0*(1-y0));
}


int main(void) {

	int n = 1, maxSteps = 100;
	double a = 0, b = 2, h, abs = 1e-6, eps = 1e-6;

	gsl_vector *tPath = gsl_vector_alloc(maxSteps);
	gsl_matrix *yPath = gsl_matrix_alloc(n,maxSteps);

	gsl_vector_set(tPath,0,a);
	gsl_matrix_set(yPath,0,0,0.5);
	h = copysign(0.1,b-a);
	int lastStep = driverPathStoring(tPath,b,&h,yPath,abs,eps,&rkStep45,&f);

	if( lastStep == DRIVER_FAIL ) {
		fprintf(stderr,"Driver couldn't find a solution within %d steps\n",maxSteps);
		return 0;
	}

	printf("# Solve y' = y*(1-y) for y(%g) = %g\n",a,0.5);
	printf("# t \t yODE \t yExact\n");
	for(int i = 0; i < lastStep; i++) {
		printf("%.3g \t %.3g",gsl_vector_get(tPath,i),gsl_matrix_get(yPath,0,i));
		printf("\t %.3g\n",1.0/(1+exp(-gsl_vector_get(tPath,i))));
	}

	gsl_vector_free(tPath);
	gsl_matrix_free(yPath);

	return 0;
}
