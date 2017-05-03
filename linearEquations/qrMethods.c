#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <assert.h>
#include "qrMethods.h"

void qr_gs_decomp(gsl_matrix *A, gsl_matrix *R) {
	int n = (*A).size1, m = (*A).size2;
	assert( n >= m );
	assert( m == (*R).size1 && m == (*R).size2 );

	gsl_vector *qR = gsl_vector_alloc(n);

	for(int i = 0; i < m; i++) {
		gsl_vector_view columnI = gsl_matrix_column(A,i);
		gsl_matrix_set(R,i,i,gsl_blas_dnrm2(&columnI.vector));
		gsl_vector_scale(&columnI.vector,1/gsl_matrix_get(R,i,i));

		for(int j = i+1; j<m; j++) {
			gsl_vector_view columnJ = gsl_matrix_column(A,j);
			double ddot;
			gsl_blas_ddot(&columnI.vector,&columnJ.vector,&ddot);
			gsl_matrix_set(R,i,j,ddot);

			gsl_vector_memcpy(qR,&columnI.vector);
			gsl_vector_scale(qR,gsl_matrix_get(R,i,j));
			gsl_vector_sub(&columnJ.vector,qR);
		}
	}

	gsl_vector_free(qR);
}

void upperTriangular_solve(const gsl_matrix *U, gsl_vector *c) {
	int n = (*c).size;
	for(int i = n-1; i >= 0; i--) {
		double yi = gsl_vector_get(c,i);
		for(int k = i+1; k < n; k++) {
			yi -= gsl_matrix_get(U,i,k) * gsl_vector_get(c,k);
		}
		gsl_vector_set(c,i,yi/gsl_matrix_get(U,i,i));
	}
}

/*
   Only works for square matrices. This means that the matrix A,
   that is decomposed into QR, must be square.
   The reason for this is that the solution is substituted into b.
*/
void qr_gs_solve(const gsl_matrix *Q, const gsl_matrix *R, gsl_vector *b) {
	int n = (*Q).size1, m = (*Q).size2;
	assert( n == m );

	gsl_vector *QTb = gsl_vector_alloc(m);
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,QTb);

	gsl_vector_memcpy(b,QTb);
	upperTriangular_solve(R,b);

	gsl_vector_free(QTb);
}

void qr_gs_inverse(const gsl_matrix *Q, const gsl_matrix *R, gsl_matrix *B) {
	int n = (*Q).size1;
	assert( (*Q).size1 == (*Q).size2 &&
		(*R).size1 == n && (*R).size2 == n &&
		(*B).size1 == n && (*B).size2 == n );

	for(int j = 0; j < n; j++) {
		for(int i = 0; i < n; i++) {
			gsl_matrix_set(B,i,j,(i == j) ? 1 : 0);
		}
		gsl_vector_view v = gsl_matrix_column(B,j);
		qr_gs_solve(Q,R,&v.vector);
	}
}

void qr_givens_decomp(gsl_matrix *A) {

	int n = (*A).size1, m = (*A).size2;

	for(int p = 0; p < n || p < m; p++) {
		for(int q = p+1; q < n; q++) {
			double theta = atan2( gsl_matrix_get(A,q,p) , gsl_matrix_get(A,p,p) );
			for(int j = p; j < m; j++) {
				double xp = gsl_matrix_get(A,p,j), xq = gsl_matrix_get(A,q,j);
				gsl_matrix_set(A,p,j, xp*cos(theta)+xq*sin(theta));
				gsl_matrix_set(A,q,j,-xp*sin(theta)+xq*cos(theta));
			}
			gsl_matrix_set(A,q,p,theta);
		}
	}
}

void qr_givens_QTb(const gsl_matrix *QR, gsl_vector *b) {
	int n = (*QR).size1, m = (*QR).size2;

	assert( (*b).size == n );

	for(int p = 0; p < n || p < m; p++) {
		for(int q = p+1; q < n; q++) {
			double theta = gsl_matrix_get(QR,q,p);
			double bp = gsl_vector_get(b,p), bq = gsl_vector_get(b,q);
			gsl_vector_set(b,p, bp*cos(theta)+bq*sin(theta));
			gsl_vector_set(b,q,-bp*sin(theta)+bq*cos(theta));
		}
	}
}

void qr_givens_solve(const gsl_matrix *QR, gsl_vector *b) {
	qr_givens_QTb(QR,b);
	upperTriangular_solve(QR,b);
}

void qr_givens_inverse(const gsl_matrix *QR, gsl_matrix *B) {
	int n = (*QR).size1, m = (*QR).size2;

	assert( n == m && (*B).size1 == n && (*B).size2 == m );

	for(int j = 0; j < m; j++) {
		for(int i = 0; i < n; i++) {
			gsl_matrix_set(B,i,j, (i == j) ? 1 : 0 );
		}
		gsl_vector_view v = gsl_matrix_column(B,j);
		qr_givens_solve(QR,&v.vector);
	}
}
