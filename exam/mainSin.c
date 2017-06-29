#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quadSpline.h"
#include "cubicSpline.h"
#include "cubicSubSpline.h"

int main(void) {

	int n = 20;
	double x[n], y[n];

	printf("# x \t y\n");
	for(int i = 0; i < n; i++) {
		x[i] = i;
		y[i] = sin(i);
		printf("%g \t %g\n",x[i],y[i]);
	}
	printf("\n\n");

	cspline *cSub = cSubSpline_alloc(n,x,y);
	cspline *c = cspline_alloc(n,x,y);
	qspline *qs = qspline_alloc(n,x,y);

	printf("#Cubic sub-spline\n");
	printf("# x \t ySub \t yCubic \t yQuad\n");
	for(double i = 0; i < n; i += 0.05) {
		printf("%g \t %g \t %g \t %g\n",
			i,
			cspline_evaluate(cSub,i),
			cspline_evaluate(c,i),
			qspline_evaluate(qs,i) );
	}

	cspline_free(cSub);
	cspline_free(c);
	qspline_free(qs);

	return 0;
}
