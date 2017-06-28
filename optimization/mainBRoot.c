#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "newtonMethods.h"


double rosenbrock(gsl_vector *x) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	return pow(1-x1,2)+100*pow(y1-pow(x1,2),2);
}

void rosenbrockGradient(gsl_vector *x, gsl_vector* df) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	gsl_vector_set(df,0,-2*(1-x1)-400*x1*(y1-pow(x1,2)));
	gsl_vector_set(df,1,200*(y1-pow(x1,2)));
}

int mainRosenbrock(void) {
	int n = 2;
	double dx = 1e-4, epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(n);
	gsl_vector_set(xStart,0,5);
	gsl_vector_set(xStart,1,-7);

	int iterations = newton(&rosenbrockGradient,n,xStart,dx,epsilon);

	gsl_vector *df = gsl_vector_alloc(n);
	rosenbrockGradient(xStart,df);

	for(int i = 0; i < n; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}
	printf("f(xStart) = %g\n",rosenbrock(xStart));
	printf("|df(xStart)| = %g\n",gsl_blas_dnrm2(df));
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(df);

	return 0;
}



double himmelblau(gsl_vector *x) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	return pow( pow(x1,2)+y1-11,2 ) + pow( x1+pow(y1,2)-7,2 );
}

void himmelblauGradient(gsl_vector *x, gsl_vector* df) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	gsl_vector_set(df,0,4*x1*(pow(x1,2)+y1-11)+2*(x1+pow(y1,2)-7));
	gsl_vector_set(df,1,2*(pow(x1,2)+y1-11)+4*y1*(x1+pow(y1,2)-7));
}

int mainHimmelblau(void) {
	int n = 2;
	double dx = 1e-4, epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(n);
	gsl_vector_set(xStart,0,5);
	gsl_vector_set(xStart,1,-7);

	int iterations = newton(&himmelblauGradient,n,xStart,dx,epsilon);

	gsl_vector *df = gsl_vector_alloc(n);
	rosenbrockGradient(xStart,df);

	for(int i = 0; i < n; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}
	printf("f(xStart) = %g\n",rosenbrock(xStart));
	printf("|df(xStart)| = %g\n",gsl_blas_dnrm2(df));
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(df);

	return 0;
}


int main(void) {
	printf("Root finding: Rosenbrock minimization\n");
	mainRosenbrock();
	printf("\n");

	printf("Root finding: Himmelblau minimization\n");
	mainHimmelblau();

	return 0;
}
