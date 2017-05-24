#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "qrMethods.h"
#include "leastSquareErrorMethods.h"

int jacobiEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal);

leastSquareError *leastSquareErrorAlloc(int n) {
	leastSquareError *ls = malloc(sizeof(leastSquareError));
	(*ls).n = n;
	(*ls).c = malloc(n);
	(*ls).dc = malloc(n);
	return ls;
}

void leastSquareErrorFree(leastSquareError *ls) {
	free((*ls).c);
	free((*ls).dc);
	free(ls);
}

void leastSquareErrorCalcConst(leastSquareError *ls, gsl_vector *x, gsl_vector *y, gsl_vector *dy) {
	int n = (*x).size, m = (*ls).n;
	assert( n == (*y).size && n == (*dy).size );

	gsl_matrix *A = gsl_matrix_alloc(n,m);
	gsl_matrix *R = gsl_matrix_alloc(m,m);
	gsl_matrix *RTR = gsl_matrix_alloc(m,m);
	gsl_matrix *sigma = gsl_matrix_alloc(m,m);

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

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,R,0.0,RTR);
	gsl_matrix_set_identity(R);
	qr_gs_inverse(R,RTR,sigma);

	for(int i = 0; i < m; i++) {
		(*ls).c[i] = gsl_vector_get(y,i);
		(*ls).dc[i] = sqrt(gsl_matrix_get(sigma,i,i));
	}

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(RTR);
	gsl_matrix_free(sigma);
}


double leastSquareErrorEval(leastSquareError *ls, double x) {
	double eval = 0;
	int m = (*ls).n;
	for(int k = 0; k < m; k++) {
		eval += (*ls).c[k] * (*ls).fitFun(k,x);
	}
	return eval;
}

double leastSquareErrorEvalPlus(leastSquareError *ls, double x) {
	double eval = 0;
	int m = (*ls).n;
	for(int k = 0; k < m; k++) {
		eval += ( (*ls).c[k] + (*ls).dc[k] ) * (*ls).fitFun(k,x);
	}
	return eval;
}

double leastSquareErrorEvalMinus(leastSquareError *ls, double x) {
	double eval = 0;
	int m = (*ls).n;
	for(int k = 0; k < m; k++) {
		eval += ( (*ls).c[k] - (*ls).dc[k] ) * (*ls).fitFun(k,x);
	}
	return eval;
}

void leastSquareErrorSingularValueDecomp(gsl_matrix *A, gsl_matrix *V, gsl_vector *S) {

	int n = (*A).size1, m = (*A).size2;

	gsl_matrix *ATA = gsl_matrix_alloc(m,m);
	gsl_matrix *AV = gsl_matrix_alloc(n,m);

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,ATA);

	jacobiEigenvalue(ATA,V,S);

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,V,0.0,AV);

	for(int i = 0; i < m; i++) {
		gsl_vector_set(S,i,sqrt( gsl_vector_get(S,i) ));
	}

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			double Djj = gsl_vector_get(S,j);
			double AVij = gsl_matrix_get(AV,i,j);
			gsl_matrix_set(A,i,j, AVij/Djj );
		}
	}

	gsl_matrix_free(ATA);
	gsl_matrix_free(AV);
}

void leastSquareErrorSingularValueSolve(gsl_matrix *V, gsl_matrix *U, gsl_vector *S, gsl_vector *b) {
	int n = (*V).size1;

	gsl_vector *UTb = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasTrans,1.0,U,b,0.0,UTb);

	for(int i = 0; i < n; i++) {
		double UTbi = gsl_vector_get(UTb,i);
		gsl_vector_set(UTb,i,UTbi/gsl_vector_get(S,i));
	}

	for(int i = 0; i < n; i++) {
		double ci = 0;
		for(int j = 0; j < n; j++) {
			ci += gsl_matrix_get(V,i,j)*gsl_vector_get(UTb,j);
		}
		gsl_vector_set(b,i,ci);
	}

	gsl_vector_free(UTb);
}

void leastSquareErrorSingularValueCalcConst(leastSquareError *ls, gsl_vector *x, gsl_vector *y, gsl_vector *dy) {
	int n = (*x).size, m = (*ls).n;
	assert( n == (*y).size && n == (*dy).size );

	gsl_matrix *A = gsl_matrix_alloc(n,m);
	gsl_matrix *V = gsl_matrix_alloc(m,m);
	gsl_vector *S = gsl_vector_alloc(m);
	gsl_matrix *VSMin2 = gsl_matrix_alloc(m,m);
	gsl_matrix *sigma = gsl_matrix_alloc(m,m);

	gsl_vector_div(y,dy);
	for(int i = 0; i < n; i++) {
		double dyi = gsl_vector_get(dy,i);
		for(int k = 0; k < m; k++) {
			double Aik = (*ls).fitFun(k,gsl_vector_get(x,i)) / dyi;
			gsl_matrix_set(A,i,k,Aik);
		}
	}

	leastSquareErrorSingularValueDecomp(A,V,S);
	leastSquareErrorSingularValueSolve(V,A,S,y);

	for(int j = 0; j < m; j++) {
		double Sjj = gsl_vector_get(S,j);
		for(int i = 0; i < m; i++) {
			double Vij = gsl_matrix_get(V,i,j);
			gsl_matrix_set(VSMin2,i,j,Vij/Sjj/Sjj);
		}
	}
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,VSMin2,V,0.0,sigma);

	for(int i = 0; i < m; i++) {
		(*ls).c[i] = gsl_vector_get(y,i);
		(*ls).dc[i] = gsl_matrix_get(sigma,i,i);
	}

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_matrix_free(VSMin2);
	gsl_matrix_free(sigma);
}
