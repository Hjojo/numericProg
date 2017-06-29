#include <math.h>
#include <stdlib.h>

#define RND ((double)rand()/RAND_MAX)

void randomx(int n, double *a, double *b, double *x) {
	for(int i = 0; i < n; i++) {
		x[i] = a[i] + RND*(b[i]-a[i]);
	}
}

void plainMonteCarlo(
	int n,
	double *a,
	double *b,
	double f(double *x),
	int N,
	double *result,
	double *error
) {
	double volume = 1;
	for(int i = 0; i < n; i++) {
		volume *= b[i]-a[i];
	}

	double sum = 0, sum2 = 0;
	double x[n];
	for(int i = 0; i < N; i++) {
		randomx(n,a,b,x);
		double fx = f(x);
		sum += fx;
		sum2 += pow(fx,2);
	}
	double mean = sum/N;
	double sigma = sqrt(sum2/N - pow(mean,2));

	*result = mean*volume;
	*error = sigma/sqrt(N)*volume;
}
