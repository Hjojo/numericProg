#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>
#include <stdlib.h>
#include "leastSquareMethods.h"
#include "qrMethods.h"

leastSquare *leastSquareAlloc(int n) {
	leastSquare *ls = malloc(sizeof(leastSquare));
	(*ls).n = n;
	(*ls).c = malloc(n);
	return ls;
}

void leastSquareFree(leastSquare *ls) {
	free((*ls).c);
	free(ls);
}

void leastSquareCalcConst(leastSquare *ls, gsl_vector *x, gsl_vector *y, gsl_vector *dy) {
	int n = (*x).size, m = (*ls).n;
	assert( n == (*y).size && n == (*dy).size );

	gsl_matrix *A = gsl_matrix_alloc(n,m);
	gsl_matrix *R = gsl_matrix_alloc(m,m);

	gsl_vector_div(y,dy);
	for(int i = 0; i < n; i++) {
		double dyi = gsl_vector_get(dy,i);
		for(int k = 0; k < m; k++) {
			double Aik = (*ls).fitFun(k,gsl_vector_get(x,i)) / dyi;
			gsl_matrix_set(A,i,k,Aik);
		}
	}

	qr_gs_decomp(A,R);
	qr_gs_solve(A,R,y);

	for(int i = 0; i < m; i++) {
		(*ls).c[i] = gsl_vector_get(y,i);
	}

	gsl_matrix_free(A);
	gsl_matrix_free(R);
}


double leastSquareEval(leastSquare *ls, double x) {
	double eval = 0;
	int m = (*ls).n;
	for(int k = 0; k < m; k++) {
		eval += (*ls).c[k] * (*ls).fitFun(k,x);
	}
	return eval;
}
