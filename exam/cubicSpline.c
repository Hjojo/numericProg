#include <stdlib.h>
#include <math.h>
#include "cubicSpline.h"

cspline * cspline_alloc(int n, double *x, double *y) {
	cspline *s = malloc(sizeof(cspline));
	(*s).n = n;
	(*s).b = malloc(n*sizeof(double));
	(*s).c = malloc((n-1)*sizeof(double));
	(*s).d = malloc((n-1)*sizeof(double));
	(*s).x = malloc(n*sizeof(double));
	(*s).y = malloc(n*sizeof(double));
	for(int i = 0; i<n; i++) {
		(*s).x[i] = x[i];
		(*s).y[i] = y[i];
	}

	double h[n-1], p[n-1];
	for(int i = 0; i<n-1; i++) {
		h[i] = x[i+1]-x[i];
		p[i] = (y[i+1]-y[i])/h[i];
	}

	double B[n], D[n], Q[n-1];
	B[0] = 3*p[0];
	B[n-1] = 3*p[n-2];
	D[0] = 2;
	D[n-1] = 2;
	Q[0] = 1;
	for(int i = 1; i<n-1; i++) {
		double hfrac = h[i-1]/h[i];
		B[i] = 3 * ( p[i-1] + p[i]*hfrac );
		D[i] = 2 * hfrac + 2;
		Q[i] = hfrac;
	}

	for(int i = 1; i<n; i++) {
		B[i] -= B[i-1]/D[i-1];
		D[i] -= Q[i-1]/D[i-1];
	}

	(*s).b[n-1] = B[n-1]/D[n-1];
	for(int i = n-2; i>=0; i--) {
		(*s).b[i] = ( B[i] - Q[i]*(*s).b[i+1] ) / D[i];
	}

	for(int i = 0; i<n-1; i++) {
		(*s).c[i] = ( -2*(*s).b[i] - (*s).b[i+1] + 3*p[i] ) / h[i];
		(*s).d[i] = ( (*s).b[i] + (*s).b[i+1] - 2*p[i] ) / pow(h[i],2);
	}

	return s;
}

double cspline_evaluate(cspline *s, double z) {
	int i = 0, j = (*s).n-1;
	while(j-i>1) {
		int m = (i+j)/2;
		if(z>(*s).x[m]) {
			i = m;
		} else {
			j = m;
		}
	}

	double deltax = z-(*s).x[i];
	return (*s).y[i] +
		(*s).b[i] * (z-(*s).x[i]) +
		(*s).c[i] * pow( deltax, 2 ) +
		(*s).d[i] * pow( deltax, 3 );
}

double cspline_derivative(cspline *s, double z) {
	int i = 0, j = (*s).n-1;
	while(j-i>1) {
		int m = (i+j)/2;
		if(z>(*s).x[m]) {
			i = m;
		} else {
			j = m;
		}
	}

	double deltax = z-(*s).x[i];
	return (*s).b[i] +
		2 * (*s).c[i] * deltax +
		3 * (*s).d[i] * pow( deltax, 2 );
}

double cspline_integralSi(cspline *s, int i, double z) {
	double deltax = z-(*s).x[i];
	double integralSi = (*s).y[i] * deltax +
		(*s).b[i]/2*pow(deltax,2) +
		(*s).c[i]/3*pow(deltax,3) +
		(*s).d[i]/4*pow(deltax,4);
	return integralSi;
}

double cspline_integral(cspline *s, double z) {
	double integ = 0;
	int i;
	for(i = 0; (*s).x[i+1] <= z; i++) {
		integ += cspline_integralSi(s,i,(*s).x[i+1]);
	}
	integ += cspline_integralSi(s,i,z);

	return integ;
}


void cspline_free(cspline *s) {
	free((*s).b);
	free((*s).c);
	free((*s).d);
	free((*s).x);
	free((*s).y);
	free(s);
}
