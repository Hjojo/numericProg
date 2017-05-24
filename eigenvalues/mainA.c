#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

#define RND (double)rand()/RAND_MAX

int jacobiEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal);
int equal(double a, double b, double tau, double epsilon);

int main(int argc, char *argv[]) {
	assert( argc > 1 );
	int n = atoi(argv[1]);
	printf("Matrix size is (%d,%d)\n",n,n);

	gsl_matrix *A = gsl_matrix_alloc(n,n);
	gsl_matrix *eigVec = gsl_matrix_alloc(n,n);
	gsl_vector *eigVal = gsl_vector_alloc(n);

	for(int i = 0; i < n; i++) {
		gsl_matrix_set(A,i,i,RND);
		for(int j = n-1; j > i; j--) {
			double Aij = RND;
			gsl_matrix_set(A,i,j,Aij);
			gsl_matrix_set(A,j,i,Aij);
		}
	}

	int noOfRotations = jacobiEigenvalue(A,eigVec,eigVal);
	printf("Number of rotations: %d\n\n",noOfRotations);

	for(int i = 0; i < n; i++) {
		for(int j = n-1; j > i; j--) {
			double Aji = gsl_matrix_get(A,j,i);
			gsl_matrix_set(A,i,j,Aji);
		}
	}

	gsl_matrix *VTA = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,eigVec,A,0.0,VTA);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,VTA,eigVec,0.0,A);

	/*printf("VTAT Matrix. Should be a diagonal matrix equal the eigenvalues.\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			printf("%.1g\t",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}

	printf("\nEigenvalue vector.\n");
	for(int i = 0; i < n; i++) {
		printf("%.1g\n",gsl_vector_get(eigVal,i));
	}
	printf("\n");*/

	printf("Testing if VTAV is is equal to a diagonal matrix with the eigenvalues...\n");
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			double diagI = gsl_vector_get(eigVal,i);
			assert(
				equal(gsl_matrix_get(A,i,j), (i == j) ? diagI : 0, 1e-6, 1e-6)
			);
		}
	}
	printf("Test successful.\n");

	gsl_matrix_free(A);
	gsl_matrix_free(eigVec);
	gsl_matrix_free(VTA);
	gsl_vector_free(eigVal);

	return 0;
}
