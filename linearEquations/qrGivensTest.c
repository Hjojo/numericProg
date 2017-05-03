#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include "qrMethods.h"

#define RND (double)rand()/RAND_MAX

int equal(double a, double b, double tau, double epsilon);

int main(void) {

	int n = 15, m = 15;

	gsl_matrix *A = gsl_matrix_alloc(n,m);
	gsl_vector *b = gsl_vector_alloc(n);
	gsl_matrix *AOriginal = gsl_matrix_alloc(n,m);
	gsl_vector *bOriginal = gsl_vector_alloc(n);

	fprintf(stderr,"A x = b\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			double matIJ = RND;
			gsl_matrix_set(A,i,j,matIJ);
			fprintf(stderr,"%.2g\t", matIJ);
		}
		double vecI = RND;
		gsl_vector_set(b,i,vecI);
		fprintf(stderr,"x%d\t%.2g\n",i,vecI);
	}
	fprintf(stderr,"\n");


	gsl_matrix_memcpy(AOriginal,A);
	gsl_vector_memcpy(bOriginal,b);

	printf("before decomp\n");
	qr_givens_decomp(A);

	fprintf(stderr,"A decomposed\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			fprintf(stderr,"%.2g\t",gsl_matrix_get(A,i,j));
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");

	qr_givens_solve(A,b);

	printf("between decomp and solve\n");
	gsl_vector *result = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasNoTrans,1.0,AOriginal,b,0.0,result);
	printf("Testing if A x = b...\n");
	fprintf(stderr,"result\tb\n");
	for(int i = 0; i < n; i++) {
		double resultI = gsl_vector_get(result,i), bI = gsl_vector_get(bOriginal,i);
		fprintf(stderr,"%.2g\t%.2g\n",resultI,bI);
		assert( equal(resultI,bI,1e-12,1e-6) );
	}
	printf("test successful\n\n");

	gsl_matrix *B = gsl_matrix_alloc(n,m);
	qr_givens_inverse(A,B);

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,AOriginal,B,0.0,A);

	printf("Testing if A^-1 A = I...\n");
	fprintf(stderr,"A^-1 A\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			fprintf(stderr,"%.1g\t",gsl_matrix_get(A,i,j));
			assert( equal(gsl_matrix_get(A,i,j), (i == j) ? 1 : 0, 1e-12, 1e-6) );
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	printf("Test successful\n\n");


	gsl_matrix_free(A);
	gsl_matrix_free(AOriginal);
	gsl_vector_free(b);
	gsl_vector_free(bOriginal);
	gsl_vector_free(result);
	gsl_matrix_free(B);



	return 0;
}
