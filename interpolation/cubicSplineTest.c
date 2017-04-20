#include <stdio.h>
#include "cubicSpline.h"
#include <stdlib.h>
#define RND ((double)rand()/RAND_MAX-0.5)*2

int main() {
	int n = 20;
	double x[n], y[n];
	double xStart = 2;

	printf("# x \t y\n");
	for(int i = 0; i<n; i++) {
		x[i] = xStart+i;
		y[i] = 2*RND;
		printf("%g \t %g\n",x[i],y[i]);
	}
	printf("\n\n");


	cspline *s = cspline_alloc(n,x,y);

	printf("# xSpline \t ySpline \t derivative \t integral\n");
	for(double i = xStart; i<=n-1+xStart; i += 0.05) {
		printf("%g \t %g \t %g \t %g\n",
			i,
			cspline_evaluate(s,i),
			cspline_derivative(s,i),
			cspline_integral(s,i) );
	}

	cspline_free(s);

	return 0;
}
