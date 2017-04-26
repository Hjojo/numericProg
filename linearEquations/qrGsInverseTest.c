#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#define RND (double)rand()/RAND_MAX

void qr_gs_decomp(gsl_matrix *A, gsl_matrix *R);
void qr_gs_inverse(const gsl_matrix *Q, const gsl_matrix *R, gsl_matrix *B);
int equal(double a, double b, double tau, double epsilon);

int main(void) {

	int n = 9;
	gsl_matrix *A = gsl_matrix_alloc(n,n);
	gsl_matrix *AOriginal = gsl_matrix_alloc(n,n);
	gsl_matrix *R = gsl_matrix_alloc(n,n);
	gsl_matrix *B = gsl_matrix_alloc(n,n);

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			gsl_matrix_set(A,i,j,RND);
		}
	}
	gsl_matrix_memcpy(AOriginal,A);

	qr_gs_decomp(A,R);
	qr_gs_inverse(A,R,B);

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,AOriginal,B,0.0,R);
	printf("Testing if AB = I...\n");
	fprintf(stderr,"Matrix AB (should be equal to I)\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			fprintf(stderr,"%.1g\t",gsl_matrix_get(R,i,j));
			double expectedValue = (i == j) ? 1 : 0;
			assert( equal( gsl_matrix_get(R,i,j), expectedValue,
					1e-12, 1e-6 ) );
		}
		fprintf(stderr,"\n");
	}
	printf("Test successful!\n");

	return 0;
}
