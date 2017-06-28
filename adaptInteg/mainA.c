#include <stdio.h>
#include <math.h>
#include "adaptInteg.h"

double fA(double x) {
	return sqrt(x);
}

double fB(double x) {
	return pow(x,3);
}

double fC(double x) {
	return cos(x);
}

/*double fPi(double x) {
	return 4*sqrt(1-pow(1-x,2));
}*/


int main(void) {
	double a = 0, b = 1, abs = 1e-2, eps = 1e-2, err;

	double Q1 = adaptInteg(fA,a,b,abs,eps,&err);
	printf("Calculating int_0^1 sqrt(x) dx\n");
	printf("Integral = %g with error = %g\n",Q1,err);
	printf("Should be %g\n",2./3);
	printf("\n");

	double Q2 = adaptInteg(fB,a,b,abs,eps,&err);
	printf("Calculating int_0^1 x^3 dx\n");
	printf("Integral = %g with error = %g\n",Q2,err);
	printf("Should be %g\n",1./4);
	printf("\n");

	double Q3 = adaptInteg(fC,a,b,abs,eps,&err);
	printf("Calculating int_0^1 cos(x) dx\n");
	printf("Integral = %g with error = %g\n",Q3,err);
	printf("Should be %g\n",sin(b));
	printf("\n");

	abs = 1e-3; eps = 1e-3;
	int noEval = 0;
	double fPi(double x)
	{
		noEval++;
		return 4*sqrt(1-pow(1-x,2));
	};
	double QPi = adaptInteg(fPi,a,b,abs,eps,&err);
	printf("Calculating int_0^1 4*sqrt(1-(1-x)^2) dx\n");
	printf("Integral = %.12g with error = %lg\n",QPi,err);
	printf("Should be %.12g\n",M_PI);
	printf("Did a total of %d evaluations\n",noEval);


	return 0;
}
