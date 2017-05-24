#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>

#define RND (double)rand()/RAND_MAX

int jacobiEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal);
int jacobiClassicEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal);

/*
	First argument decides the size of the symmetric square matrix
	and second argument should be 1 for cyclic mehtod and 2 for classic method.
*/
int main(int argc, char *argv[]) {
	assert(argc > 2);
	assert( atoi(argv[2]) != 1 || atoi(argv[2]) != 2 );

	int n = atoi(argv[1]);

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

	int noOfRotations;
	printf("Full diagonalization of the matrix of size (%d,%d).\n",n,n);
	if( atoi(argv[2]) == 1 ) {
		noOfRotations = jacobiEigenvalue(A,eigVec,eigVal);
		printf("Number of rotations for cyclic: ");
	} else if( atoi(argv[2]) == 2 ) {
		noOfRotations = jacobiClassicEigenvalue(A,eigVec,eigVal);
		printf("Number of rotations for classic method: ");
	}
	printf("%d\n",noOfRotations);

	gsl_matrix_free(A);
	gsl_matrix_free(eigVec);
	gsl_vector_free(eigVal);

	return 0;
}
