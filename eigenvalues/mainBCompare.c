#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <assert.h>

#define RND (double)rand()/RAND_MAX

int jacobiEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal);
int jacobiNumberLowestEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal, int numberEigVal);


void restoreA(gsl_matrix *A) {
	int n = (*A).size1;
	for(int i = 0; i < n; i++) {
		for(int j = n-1; j > i; j--) {
			double Aji = gsl_matrix_get(A,j,i);
			gsl_matrix_set(A,i,j,Aji);
		}
	}
}


int main(int argc, char *argv[]) {
	assert(argc > 1);

	int n = atoi(argv[1]);
	printf("Matrix size is (%d,%d).\n",n,n);

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
	printf("Full diagonalization of the matrix of size (%d,%d).\n",n,n);
	int noOfRotationsVal = jacobiNumberLowestEigenvalue(A,eigVec,eigVal,n);
	printf("Number of rotations for ""value-by-value"": %d\n",noOfRotationsVal);

	restoreA(A);
	int noOfRotationsCyclic = jacobiEigenvalue(A,eigVec,eigVal);
	/*int noOfRotationsCyclic = noOfSweeps * ( (n-1)*(n-2)/2 );*/
	printf("Number of rotations for cyclic: %d\n\n",noOfRotationsCyclic);

	restoreA(A);
	int noOfRotationsLowest = jacobiNumberLowestEigenvalue(A,eigVec,eigVal,1);
	printf("Number of rotations for only the lowest eigenvalue with the ""value-by-value"" method: ");
	printf("%d\n\n",noOfRotationsLowest);


	gsl_matrix_free(A);
	gsl_matrix_free(eigVec);
	gsl_vector_free(eigVal);

	return 0;
}
