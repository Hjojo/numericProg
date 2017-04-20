#include <stdlib.h>
#include <math.h>
#include "quadSpline.h"
#include <stdio.h>

qspline * qspline_alloc(int n, double *x, double *y) {
	qspline *s = malloc(sizeof(qspline));
	(*s).n = n;
	(*s).b = malloc((n-1)*sizeof(double));
	(*s).c = malloc((n-1)*sizeof(double));
	(*s).x = malloc(n*sizeof(double));
	(*s).y = malloc(n*sizeof(double));
	for(int i = 0; i<n; i++) {
		(*s).x[i] = x[i];
		(*s).y[i] = y[i];
	}

	double p[n-1], deltax[n-1];
	for(int i = 0; i<n-1; i++) {
		deltax[i] = (*s).x[i+1]-(*s).x[i];
		double deltay = (*s).y[i+1]-(*s).y[i];
		p[i] = deltay/deltax[i];
	}

	(*s).c[0] = 0;
	for(int i = 0; i<n-2; i++) {
		(*s).c[i+1] = (p[i+1]-p[i]-(*s).c[i]*deltax[i])/deltax[i+1];
	}

	(*s).c[n-2] = (*s).c[n-2]/2;
	for(int i = n-3; i>=0; i--) {
		(*s).c[i] = (p[i+1]-p[i]-(*s).c[i+1]*deltax[i+1])/deltax[i];
	}

	for(int i = 0; i<n-1; i++) {
		(*s).b[i] = p[i]-(*s).c[i]*deltax[i];
	}

	return s;
}

double qspline_evaluate(qspline *s, double z) {
	int i = 0, j = (*s).n-1;
	while(j-i>1) {
		int m = (i+j)/2;
		if(z>(*s).x[m]) {
			i = m;
		} else {
			j = m;
		}
	}

	return (*s).y[i] +
		(*s).b[i] * (z-(*s).x[i]) +
		(*s).c[i] * pow( z-(*s).x[i] , 2 );


}

double qspline_derivative(qspline *s, double z) {
	int i = 0, j = (*s).n-1;
	while(j-i>1) {
		int m = (i+j)/2;
		if(z>(*s).x[m]) {
			i = m;
		} else {
			j = m;
		}
	}
	double derivative = (*s).b[i] + 2 * (*s).c[i] * (z-(*s).x[i]);
	return derivative;
}

double qspline_integralSi(double xi, double z, double yi, double b, double c) {
	double deltax = z-xi;
	double integralSi = yi * deltax +
		b/2*pow(deltax,2) +
		c/3*pow(deltax,3);
	return integralSi;
}

double qspline_integral(qspline *s, double z) {
	double integ = 0;
	int i;
	for(i = 0; (*s).x[i+1] <= z; i++) {
		integ += qspline_integralSi(
			(*s).x[i], (*s).x[i+1], (*s).y[i],
			(*s).b[i], (*s).c[i] );
	}
	integ += qspline_integralSi(
			(*s).x[i], z, (*s).y[i],
			(*s).b[i], (*s).c[i] );

	return integ;
}


void qspline_free(qspline *s) {
	free((*s).b);
	free((*s).c);
	free((*s).x);
	free((*s).y);
	free(s);
}
