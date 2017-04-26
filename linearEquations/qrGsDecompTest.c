#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <assert.h>

#define RND (double)rand()/RAND_MAX


void qr_gs_decomp(gsl_matrix *A, gsl_matrix *R);
void qr_gs_solve(const gsl_matrix *Q, const gsl_matrix *R, gsl_vector *b);
int equal(double a, double b, double tau, double epsilon);

int main(void) {

	int n = 10, m = 10;
	double tau = 1e-12, epsilon = 1e-12;
	gsl_matrix *A = gsl_matrix_alloc(n,m);
	gsl_matrix *R = gsl_matrix_alloc(m,m);

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			gsl_matrix_set(A,i,j,RND);
		}
	}

	gsl_matrix *AOriginal = gsl_matrix_alloc(n,m);
	gsl_matrix_memcpy(AOriginal,A);

	qr_gs_decomp(A,R);

	printf("Checking if R is upper triangular...\n");
	fprintf(stderr,"All values under the diagonal in R should be zero.\n");
	fprintf(stderr,"R matrix\n");
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < m; j++) {
			fprintf(stderr,"%.2g\t",gsl_matrix_get(R,i,j));
			if( i > j ) {
				assert(gsl_matrix_get(R,i,j) == 0);
			}
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	printf("R is triangular.\n\n");

	printf("Checking that Q^T Q = 1...\n");
	fprintf(stderr,"Values for QTQ. Should be equal to the identity matrix.\n");
	gsl_matrix *QTQ = gsl_matrix_alloc(m,m);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0.0,QTQ);
	fprintf(stderr,"QTQ matrix\n");
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < m; j++) {
			fprintf(stderr,"%.1g\t",gsl_matrix_get(QTQ,i,j));
			double expectedValue = (i == j) ? 1 : 0;
			//printf("Exp %g\n",expectedValue);
			assert( equal(gsl_matrix_get(QTQ,i,j),expectedValue,tau,epsilon) );
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	printf("The check for Q^T Q = 1 turned out positive.\n\n");

	printf("Checking if QR = A...\n");
	fprintf(stderr,"Elementens in QR and A should be equal.\n");
	gsl_matrix *QR = gsl_matrix_alloc(n,m);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,R,0.0,QR);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			fprintf(stderr,"QR(%d,%d) = %g, A(%d,%d) = %g\n",
				i,j,gsl_matrix_get(QR,i,j),
				i,j,gsl_matrix_get(AOriginal,i,j));
			assert(
				equal( gsl_matrix_get(QR,i,j),
					gsl_matrix_get(AOriginal,i,j),
					tau,epsilon )
				);
		}
	}
	fprintf(stderr,"\n");
	printf("The check for QR = A turned out positive.\n\n");



	gsl_vector *b = gsl_vector_alloc(n);
	for(int i = 0; i < n; i++) {
		gsl_vector_set(b,i,RND);
	}
	gsl_vector *bOriginal = gsl_vector_alloc(n);
	gsl_vector_memcpy(bOriginal,b);
	qr_gs_solve(A,R,b);

	gsl_vector *bSolved = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasNoTrans,1.0,AOriginal,b,0.0,bSolved);

	printf("Testing that Ax = b...\n");
	fprintf(stderr,"bOriginal \t bSolved\n");
	for(int i = 0; i < n; i++) {
		fprintf(stderr,"%g \t %g\n",
			gsl_vector_get(bOriginal,i),
			gsl_vector_get(bSolved,i));
		assert(
			equal(gsl_vector_get(bOriginal,i),
				gsl_vector_get(bSolved,i),
				tau,epsilon)
			);
	}
	printf("Check turned out successful\n");
	printf("\n");

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(AOriginal);
	gsl_matrix_free(QTQ);
	gsl_matrix_free(QR);
	gsl_vector_free(b);
	gsl_vector_free(bOriginal);
	gsl_vector_free(bSolved);

	return 0;
}
