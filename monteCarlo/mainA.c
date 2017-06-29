#include <stdio.h>
#include <math.h>


void plainMonteCarlo(
	int n,
	double *a,
	double *b,
	double f(double *x),
	int N,
	double *result,
	double *error
);


double fSimple1(double *x) {
	return pow(x[0],3)+x[0]*x[2];
}

void simpleIntegral1(void) {
	int n = 3;
	double a[3] = {0, 0, -5};
	double b[3] = {3, 2*M_PI, 5};

	double result, error;
	int N = 1e4;

	plainMonteCarlo(n,a,b,&fSimple1,N,&result,&error);

	printf("Simple integral #1\n");
	printf("Result = %g\n",result);
	printf("Error = %g\n",error);
	printf("Should be %g\n",405*M_PI);
}



double fSimple2(double *x) {
	return pow(x[0],2)+4*x[1];
}

void simpleIntegral2(void) {
	int n = 2;
	double a[2] = {11, 7};
	double b[2] = {14, 10};

	double result, error;
	int N = 1e4;

	plainMonteCarlo(n,a,b,&fSimple2,N,&result,&error);

	printf("Simple integral #2\n");
	printf("Result = %g\n",result);
	printf("Error = %g\n",error);
	printf("Should be %g\n",1719.);
}





double fHard(double *x) {
	return 1/( 1-cos(x[0])*cos(x[1])*cos(x[2]) )/pow(M_PI,3);
}

void hardIntegral(void) {
	int n = 3;
	double a[3] = {0};
	double b[3] = {M_PI, M_PI, M_PI};

	double result, error;
	int N = 1e4;

	plainMonteCarlo(n,a,b,&fHard,N,&result,&error);

	printf("Hard integral\n");
	printf("Result = %g\n",result);
	printf("Error = %g\n",error);
	printf("Should be %g\n",1.3932039296856768591842462603255);
}


int main(void) {
	simpleIntegral1();
	printf("\n");
	simpleIntegral2();
	printf("\n");
	hardIntegral();

	return 0;
}

