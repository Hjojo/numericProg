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



double fHard(double *x) {
	return 1/( 1-cos(x[0])*cos(x[1])*cos(x[2]) )/pow(M_PI,3);
}

int main(void) {
	int n = 3;
	double a[3] = {0};
	double b[3] = {M_PI, M_PI, M_PI};

	double result, error;

	printf("# Hard integral\n");
	printf("# N \t error\n");
	for(int N = 10; N <= 1e7; N *= 10) {
		plainMonteCarlo(n,a,b,&fHard,N,&result,&error);
		printf("%d \t %g\n",N,error);
	}

	return 0;
}
