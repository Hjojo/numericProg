#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "ODEMethods.h"

void f(double t, gsl_vector *y, gsl_vector *dydt) {

	double y0 = gsl_vector_get(y,0);

	gsl_vector_set(dydt,0,y0*(1-y0));
}


int main(void) {

	int n = 1;
	double a = 0, b = 2, h, abs = 1e-6, eps = 1e-6;

	gsl_vector *y = gsl_vector_alloc(n);

	gsl_vector_set(y,0,0.5);
	h = copysign(0.1,b-a);
	driver(&a,b,&h,y,abs,eps,&rkStep45,&f);

	printf("Solve y' = y*(1-y) for y(%g) = %g\n",0.,0.5);
	printf("y(%g) = %g\n",a,gsl_vector_get(y,0));
	printf("Exact solution is %g\n",1.0/(1+exp(-a)));

	gsl_vector_free(y);

	return 0;
}
