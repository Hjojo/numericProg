#include <stdio.h>
#include <math.h>
#include "adaptInteg.h"


int main(void) {
	double a, b, abs = 1e-6, eps = 1e-6, err;
	int noEval;


	double fA(double x) {
		noEval++;
		return exp(-x);
	}
	a = 0; b = INFINITY; noEval = 0;
	double Q1 = adaptIntegOpenWithInf(fA,a,b,abs,eps,&err);
	printf("Calculating int_0^INF exp(-x) dx\n");
	printf("Integral = %g with error = %g\n",Q1,err);
	printf("Should be %g\n",1.);
	printf("Did a total of %d evaluations\n",noEval);
	printf("\n");

	double fB(double x) {
		noEval++;
		return x*exp(-pow(x,2));
	}
	a = -INFINITY; b = 0; noEval = 0;
	double Q2 = adaptIntegOpenWithInf(fB,a,b,abs,eps,&err);
	printf("Calculating int_-INF^0 x*exp(-x^2) dx\n");
	printf("Integral = %g with error = %g\n",Q2,err);
	printf("Should be %g\n",-1./2);
	printf("Did a total of %d evaluations\n",noEval);
	printf("\n");

	double fC(double x) {
		noEval++;
		return 1/(pow(x,2)+25);
	}
	a = -INFINITY; b = INFINITY; noEval = 0;
	double Q3 = adaptIntegOpenWithInf(fC,a,b,abs,eps,&err);
	printf("Calculating int_-INF^INF 1/(x^2+25) dx\n");
	printf("Integral = %g with error = %g\n",Q3,err);
	printf("Should be %g\n",M_PI/5);
	printf("Did a total of %d evaluations\n",noEval);

	return 0;
}
