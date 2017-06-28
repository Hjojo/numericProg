#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#define RND ((double)rand()/RAND_MAX-0.5)*2

double linterp(int n, double *x, double *y, double z);

double linterp_integ(int n, double *x, double *y, double z) {
	assert(z >= x[0] && z <= x[n-1]);

	double area = 0, p;
	int i;
	for(i = 0; x[i+1] <= z; i++) {
		double deltax = x[i+1]-x[i];
		p = (y[i+1]-y[i])/deltax;
		area += y[i]*deltax + 1.0/2 * p * pow(deltax,2);
	}
	p = (y[i+1]-y[i])/(x[i+1]-x[i]);
	area += y[i]*(z-x[i]) + 1.0/2 * p * pow(z-x[i],2);
	return area;
}

int main(void) {

	double x0 = 2, xEnd = 20;
	int n = xEnd-x0+1;
	double *x = malloc(n*sizeof(double));
	double *y = malloc(n*sizeof(double));
	double *ysin = malloc(n*sizeof(double));

	printf("# x \t y \t sin(x)\n");
	for(int i = 0; i < n; i++) {
		x[i] = i+x0;
		y[i] = RND;
		ysin[i] = sin(x[i]);
		printf("%g \t %g \t %g\n",x[i],y[i],ysin[i]);
	}
	printf("\n\n");

	double dx = 0.2;
	printf("# xi \t yi \t yint \t ysini \t ysinInt\n");
	for(double i=x0; i<=xEnd; i+=dx) {
		printf("%g \t %g \t %g \t %g \t %g\n",
			i,
			linterp(n,x,y,i),
			linterp_integ(n,x,y,i),
			linterp(n,x,ysin,i),
			linterp_integ(n,x,ysin,i));
	}

	free(x);
	free(y);
	free(ysin);

	return 0;

}
