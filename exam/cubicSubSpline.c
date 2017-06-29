#include <stdlib.h>
#include <math.h>
#include "cubicSubSpline.h"
#include "quadSpline.h"

cspline * cSubSpline_alloc(int n, double *x, double *y) {
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

	int nTemp = 3;
	double xTemp[nTemp], yTemp[nTemp];
	for(int i = 1; i < n-1; i++) {
		xTemp[0] = x[i-1]; xTemp[1] = x[i]; xTemp[2] = x[i+1];
		xTemp[0] = y[i-1]; yTemp[1] = y[i]; yTemp[2] = y[i+1];
		qspline *qs = qspline_alloc(nTemp,xTemp,yTemp);

		(*s).b[i] = (*qs).b[1];
		if( i == 1 ) {
			(*s).b[0] = (*qs).b[0];
		}
		if( i == n-2 ) {
			(*s).b[n-1] = (*qs).b[nTemp-1];
		}

		qspline_free(qs);
	}

	double h[n-1], p[n-1];
	for(int i = 0; i<n-1; i++) {
		h[i] = x[i+1]-x[i];
		p[i] = (y[i+1]-y[i])/h[i];
	}

	for(int i = 0; i<n-1; i++) {
		(*s).c[i] = ( -2*(*s).b[i] - (*s).b[i+1] + 3*p[i] ) / h[i];
		(*s).d[i] = ( (*s).b[i] + (*s).b[i+1] - 2*p[i] ) / pow(h[i],2);
	}

	return s;
}
