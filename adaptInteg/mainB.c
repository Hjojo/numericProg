#include <stdio.h>
#include <math.h>
#include "adaptInteg.h"


int main(void) {
	double a, b, abs = 1e-3, eps = 1e-3, err;
	int noEval;


	double fA(double x) {
		noEval++;
		return exp(-x);
	}
	a = 0; b = INFINITY; noEval = 0;
	double Q1 = adaptIntegWithInf(fA,a,b,abs,eps,&err);
	printf("Calculating int_0^INF exp(-x) dx\n");
	printf("Integral = %g with error = %g\n",Q1,err);
	printf("Should be %g\n",1.);
	printf("Did a total of %d evaluations\n",noEval);
	printf("\n");
	/*
	double fB(double x) {
		noEvals++;
		return
	}
	a = 0; b = INFINITY; noEval = 0;
	double Q2 = adaptInteg(fB,a,b,abs,eps,&err);
	printf("Calculating int_0^1 x^3 dx\n");
	printf("Integral = %g with error = %g\n",Q2,err);
	printf("Should be %g\n",1./4);
	printf("Did a total of %d evaluations\n",noEval);
	printf("\n");

	a = 0; b = INFINITY; noEval = 0;
	double Q3 = adaptInteg(fC,a,b,abs,eps,&err);
	printf("Calculating int_0^1 cos(x) dx\n");
	printf("Integral = %g with error = %g\n",Q3,err);
	printf("Should be %g\n",sin(b));
	printf("Did a total of %d evaluations\n",noEval);
	printf("\n");


	double fPi(double x)
	{
		noEval++;
		return 4*sqrt(1-pow(1-x,2));
	};
	a = 0; b = INFINITY; noEval = 0;
	double QPi = adaptInteg(fPi,a,b,abs,eps,&err);
	printf("Calculating int_0^1 4*sqrt(1-(1-x)^2) dx\n");
	printf("Integral = %.12g with error = %lg\n",QPi,err);
	printf("Should be %.12g\n",M_PI);
	printf("Did a total of %d evaluations\n",noEval);
	*/

	return 0;
}
