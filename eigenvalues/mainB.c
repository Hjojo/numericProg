#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define RND (double)rand()/RAND_MAX

int jacobiNumberLowestEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal, int numberEigVal);

int main(int argc, char *argv[]) {
	int n;
	int numberEigVal;

	if(argc < 3) {
		n = 12;
		numberEigVal = n/2;
		printf("Since no argument was passed matrix size if (%dx%d) and %d eigenvalues will be found.\n",n,n,numberEigVal);
	} else {
		n = atoi(argv[1]);
		numberEigVal = atoi(argv[2]);
		printf("Size of matrix is (%d,%d) and number of eigenvalues to be found is %d\n",n,n,numberEigVal);
	}

	printf("The largest eigenvalues can be found by adding a phase of pi/2 to the rotation angle\n\n");

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

	int noOfRotations = jacobiNumberLowestEigenvalue(A,eigVec,eigVal,numberEigVal);
	printf("Number of rotations: %d\n\n",noOfRotations);

	printf("Eigenvalue vector - %d smallest elements.\n",numberEigVal);
	for(int i = 0; i < numberEigVal; i++) {
		printf("%.1g\n",gsl_vector_get(eigVal,i));
	}
	printf("\n");

	gsl_matrix_free(A);
	gsl_matrix_free(eigVec);
	gsl_vector_free(eigVal);

	return 0;
}
