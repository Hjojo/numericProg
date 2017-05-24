#include <gsl/gsl_vector.h>
#include <math.h>
#include "newtonMethods.h"

void fRandom(gsl_vector *xVec, gsl_vector *fx) {
	double A = 1e4;
	double x = gsl_vector_get(xVec,0);
	double y = gsl_vector_get(xVec,1);

	gsl_vector_set(fx,0, A*x*y - 1 );
	gsl_vector_set(fx,1, exp(-x) + exp(-y) - 1 - 1.0/A);
}

int mainRandomEquations(void) {

	int noOfVar = 2;
	int noOfFuns = 2;
	double dx = 1e-4, epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(noOfVar);
	gsl_vector_set(xStart,0,1e-1);
	gsl_vector_set(xStart,1,0.1e-3);

	int iterations = newton(&fRandom,noOfFuns,xStart,dx,epsilon);

	for(int i = 0; i < noOfVar; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	fRandom(xStart,fx);
	for(int i = 0; i < noOfFuns; i++) {
		printf("Function %d = %g\n",i,gsl_vector_get(fx,i));
	}
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(fx);

	return 0;
}

void fRosenbrock(gsl_vector *xVec, gsl_vector *fx) {
	double x = gsl_vector_get(xVec,0);
	double y = gsl_vector_get(xVec,1);

	gsl_vector_set(fx,0,-2*(1-x)-400*x*(y-x*x));
	gsl_vector_set(fx,1,200*(y-x*x));
}

int mainRosenbrock(void) {
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
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(fx);

	return 0;
}

void fHimmelblau(gsl_vector *xVec, gsl_vector *fx) {
	double x = gsl_vector_get(xVec,0);
	double y = gsl_vector_get(xVec,1);

	gsl_vector_set(fx,0, 4*x*(x*x+y-11) + 2*(x+y*y-7) );
	gsl_vector_set(fx,1, 2*(x*x+y-11) + 4*y*(x+y*y-7) );
}

int mainHimmelblau(void) {
	int noOfVar = 2;
	int noOfFuns = 2;
	double dx = 1e-4, epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(noOfVar);
	gsl_vector_set(xStart,0,2);
	gsl_vector_set(xStart,1,3);

	int iterations = newton(&fHimmelblau,noOfFuns,xStart,dx,epsilon);

	for(int i = 0; i < noOfVar; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	fHimmelblau(xStart,fx);
	for(int i = 0; i < noOfFuns; i++) {
		printf("Function %d = %g\n",i,gsl_vector_get(fx,i));
	}
	printf("Did %d iterations.\n",iterations);


	gsl_vector_free(xStart);
	gsl_vector_free(fx);

	return 0;
}


void fExample1(gsl_vector *xVec, gsl_vector *fx) {
	double x = gsl_vector_get(xVec,0);
	double y = gsl_vector_get(xVec,1);

	// Minumum of 5*x*x + 7*y - 30*x*y
	gsl_vector_set(fx,0, 10*x - 30*y );
	gsl_vector_set(fx,1, 7 - 30*x );
}

int mainExample1(void) {
	int noOfVar = 2;
	int noOfFuns = 2;
	double dx = 1e-4, epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(noOfVar);
	gsl_vector_set(xStart,0,2);
	gsl_vector_set(xStart,1,3);

	int iterations = newton(&fExample1,noOfFuns,xStart,dx,epsilon);

	for(int i = 0; i < noOfVar; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	fExample1(xStart,fx);
	for(int i = 0; i < noOfFuns; i++) {
		printf("Function %d = %g\n",i,gsl_vector_get(fx,i));
	}
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(fx);

	return 0;
}


void fExample2(gsl_vector *xVec, gsl_vector *fx) {
	double x = gsl_vector_get(xVec,0);
	double y = gsl_vector_get(xVec,1);
	double z = gsl_vector_get(xVec,2);

	// Minimum of pow(x-5,2) + pow(y-0.1,2) + pow(z+38,2) - 88
	gsl_vector_set(fx,0, 2*(x-5) );
	gsl_vector_set(fx,1, 2*(y-0.1) );
	gsl_vector_set(fx,2, 2*(z+38) );
}

int mainExample2(void) {
	int noOfVar = 3;
	int noOfFuns = 3;
	double dx = 1e-4, epsilon = 1e-6;
	gsl_vector *xStart = gsl_vector_alloc(noOfVar);
	gsl_vector_set(xStart,0,2);
	gsl_vector_set(xStart,1,3);
	gsl_vector_set(xStart,2,0.1);

	int iterations = newton(&fExample2,noOfFuns,xStart,dx,epsilon);

	for(int i = 0; i < noOfVar; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	fExample2(xStart,fx);
	for(int i = 0; i < noOfFuns; i++) {
		printf("Function %d = %g\n",i,gsl_vector_get(fx,i));
	}
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(fx);

	return 0;
}


int main(void) {
	printf("Random equations\n");
	mainRandomEquations();
	printf("\n");
	printf("Rosenbrock\n");
	mainRosenbrock();
	printf("\n");
	printf("Himmelblau\n");
	mainHimmelblau();
	printf("\n");
	printf("Example 1\n");
	mainExample1();
	printf("\n");
	printf("Example 2\n");
	mainExample2();
	return 0;
}
