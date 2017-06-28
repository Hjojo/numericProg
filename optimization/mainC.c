#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define RND (double)rand()/RAND_MAX

int downhillSimplex(
	double f(gsl_vector *x),
	gsl_matrix *simplex,
	int n,
	double minDist
);



void printMatrix(gsl_matrix *A) {
	int n = (*A).size1, m = (*A).size2;

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			printf("%.2g \t",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}

}


double rosenbrock(gsl_vector *x) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	return pow(1-x1,2)+100*pow(y1-pow(x1,2),2);
}

int mainRosenbrock(void) {
	int n = 2;

	gsl_matrix *simplexStart = gsl_matrix_alloc(n,n+1);

	for( int i = 0; i < n; i++) {
		for( int j = 0; j < n+1; j++) {
			//gsl_matrix_set(simplexStart,i,j,20*(RND-0.5));
			gsl_matrix_set(simplexStart,i,j,5+(i==j ? 6 : 0));
		}
	}

	double minDist = 1e-6;

	int iterations = downhillSimplex(
		&rosenbrock,
		simplexStart,
		n,
		minDist
	);

	printf("Did %d iterations.\n",iterations);
	printf("Simplex:\n");
	printMatrix(simplexStart);

	gsl_matrix_free(simplexStart);

	return 0;
}


double himmelblau(gsl_vector *x) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	return pow( pow(x1,2)+y1-11,2 ) + pow( x1+pow(y1,2)-7,2 );
}

int mainHimmelblau(void) {
	int n = 2;

	gsl_matrix *simplexStart = gsl_matrix_alloc(n,n+1);

	for( int i = 0; i < n; i++) {
		for( int j = 0; j < n+1; j++) {
			//gsl_matrix_set(simplexStart,i,j,20*(RND-0.5));
			gsl_matrix_set(simplexStart,i,j,-2+(i==j ? 6 : 0));
		}
	}

	double minDist = 1e-6;

	int iterations = downhillSimplex(
		&himmelblau,
		simplexStart,
		n,
		minDist
	);

	printf("Did %d iterations.\n",iterations);
	printf("Simplex:\n");
	printMatrix(simplexStart);

	gsl_matrix_free(simplexStart);

	return 0;
}

int main(void) {
	printf("Simplex: Rosenbrock minimization\n");
	mainRosenbrock();
	printf("\n");
	printf("Simplex: Himmelblau minimization\n");
	mainHimmelblau();

	return 0;
}
