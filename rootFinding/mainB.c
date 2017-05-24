#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include "newtonMethods.h"

void fRosenbrock(gsl_vector *xVec, gsl_vector *fx) {
	double x = gsl_vector_get(xVec,0);
	double y = gsl_vector_get(xVec,1);

	// Minimum of (1-x)^2 + 100(y-x^2)^2
	gsl_vector_set(fx,0,-2*(1-x)-400*x*(y-x*x));
	gsl_vector_set(fx,1,200*(y-x*x));
}

void jacobianRosenbrock(gsl_vector *xVec, gsl_matrix *J) {
	double x = gsl_vector_get(xVec,0);
	double y = gsl_vector_get(xVec,1);

	gsl_matrix_set(J,0,0,2-400*y+1200*x*x);
	gsl_matrix_set(J,0,1,-400*x);
	gsl_matrix_set(J,1,0,-400*x);
	gsl_matrix_set(J,1,1,200);

}

int mainRosenbrockAnalytic(void) {
	int noOfVar = 2;
	int noOfFuns = 2;
	double epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(noOfVar);
	gsl_vector_set(xStart,0,3);
	gsl_vector_set(xStart,1,7);

	int iterations = newtonWithJacobian(&fRosenbrock,noOfFuns,&jacobianRosenbrock,xStart,epsilon);

	for(int i = 0; i < noOfVar; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	fRosenbrock(xStart,fx);
	for(int i = 0; i < noOfFuns; i++) {
		printf("Function %d = %g\n",i,gsl_vector_get(fx,i));
	}
	printf("Did it in %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(fx);

	return 0;
}

int mainRosenbrockDiff(void) {
	int noOfVar = 2;
	int noOfFuns = 2;
	double dx = 1e-4, epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(noOfVar);
	gsl_vector_set(xStart,0,3);
	gsl_vector_set(xStart,1,7);

	int iterations = newton(&fRosenbrock,noOfFuns,xStart,dx,epsilon);

	for(int i = 0; i < noOfVar; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	fRosenbrock(xStart,fx);
	for(int i = 0; i < noOfFuns; i++) {
		printf("Function %d = %g\n",i,gsl_vector_get(fx,i));
	}
	printf("Did it in %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(fx);

	return 0;
}


int main(void) {
	printf("Rosenbrock Analytic jacobian\n");
	mainRosenbrockAnalytic();
	printf("\n");

	printf("Rosenbrock diff jacobian\n");
	mainRosenbrockDiff();
	printf("\n");

	return 0;
}
