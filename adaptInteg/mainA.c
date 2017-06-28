#include <stdio.h>
#include <math.h>
#include "adaptInteg.h"

double fA(double x) {
	return sqrt(x);
}

double fB(double x) {
	return 1/sqrt(x);
}

double fC(double x) {
	return log(x)/sqrt(x);
}


int main(void) {
	double a = 0, b = 1, abs = 1e-4, eps = 1e-4, err;

	double Q1 = adaptIntegOpen(fA,a,b,abs,eps,&err);
	printf("Calculating int_0^1 sqrt(x) dx\n");
	printf("Integral = %g with error = %g\n",Q1,err);
	printf("Should be %g\n",2./3);
	printf("\n");

	double Q2 = adaptIntegOpen(fB,a,b,abs,eps,&err);
	printf("Calculating int_0^1 1/sqrt(x) dx\n");
	printf("Integral = %g with error = %g\n",Q2,err);
	printf("Should be %g\n",2.);
	printf("\n");

	double Q3 = adaptIntegOpen(fC,a,b,abs,eps,&err);
	printf("Calculating int_0^1 ln(x)/sqrt(x) dx\n");
	printf("Integral = %g with error = %g\n",Q3,err);
	printf("Should be %g\n",-4.);
	printf("\n");

	abs = 1e-8; eps = 1e-8;
	int noEval = 0;
	double fPi(double x)
	{
		noEval++;
		return 4*sqrt(1-pow(1-x,2));
	};
	double QPi = adaptIntegOpen(fPi,a,b,abs,eps,&err);
	printf("Calculating int_0^1 4*sqrt(1-(1-x)^2) dx\n");
	printf("Integral = %.15g with error = %lg\n",QPi,err);
	printf("Should be  %.15g\n",M_PI);
	printf("Did a total of %d evaluations\n",noEval);


	return 0;
}
