#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int newton(
	double f(gsl_vector *x),
	void gradient(gsl_vector *x, gsl_vector *df),
	void hessian(gsl_vector *x, gsl_matrix *H),
	gsl_vector *xStart,
	double eps
);

int quasiNewton(
	double f(gsl_vector *x),
	void gradient(gsl_vector *x, gsl_vector *df),
	gsl_vector *xStart,
	double eps
);



double rosenbrock(gsl_vector *x) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	return pow(1-x1,2)+100*pow(y1-pow(x1,2),2);
}

void rosenbrockGradient(gsl_vector *x, gsl_vector* df) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	gsl_vector_set(df,0,-2*(1-x1)-400*x1*(y1-pow(x1,2)));
	gsl_vector_set(df,1,200*(y1-pow(x1,2)));
}

void rosenbrockHessian(gsl_vector *x, gsl_matrix *H) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	double H12And21 = -400*x1;
	gsl_matrix_set(H,0,0,2-400*y1+1200*pow(x1,2));
	gsl_matrix_set(H,1,1,200);
	gsl_matrix_set(H,0,1,H12And21);
	gsl_matrix_set(H,1,0,H12And21);
}

int mainRosenbrockNewton(void) {
	int n = 2;

	gsl_vector *xStart = gsl_vector_alloc(n);
	gsl_vector_set(xStart,0,5);
	gsl_vector_set(xStart,1,-7);

	double eps = 1e-6;

	int iterations = newton(
		&rosenbrock,
		&rosenbrockGradient,
		&rosenbrockHessian,
		xStart,
		eps
	);

	gsl_vector *df = gsl_vector_alloc(n);
	rosenbrockGradient(xStart,df);

	for(int i = 0; i < n; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}
	printf("f(xStart) = %g\n",rosenbrock(xStart));
	printf("|df(xStart)| = %g\n",gsl_blas_dnrm2(df));
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(df);

	return 0;
}

int mainRosenbrockQuasiNewton(void) {
	int n = 2;

	gsl_vector *xStart = gsl_vector_alloc(n);
	gsl_vector_set(xStart,0,5);
	gsl_vector_set(xStart,1,-7);

	double eps = 1e-6;

	int iterations = quasiNewton(
		&rosenbrock,
		&rosenbrockGradient,
		xStart,
		eps
	);

	gsl_vector *df = gsl_vector_alloc(n);
	rosenbrockGradient(xStart,df);

	for(int i = 0; i < n; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}
	printf("f(xStart) = %g\n",rosenbrock(xStart));
	printf("|df(xStart)| = %g\n",gsl_blas_dnrm2(df));
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(df);

	return 0;
}


double himmelblau(gsl_vector *x) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	return pow( pow(x1,2)+y1-11,2 ) + pow( x1+pow(y1,2)-7,2 );
}

void himmelblauGradient(gsl_vector *x, gsl_vector* df) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	gsl_vector_set(df,0,4*x1*(pow(x1,2)+y1-11)+2*(x1+pow(y1,2)-7));
	gsl_vector_set(df,1,2*(pow(x1,2)+y1-11)+4*y1*(x1+pow(y1,2)-7));
}

void himmelblauHessian(gsl_vector *x, gsl_matrix *H) {
	double x1 = gsl_vector_get(x,0);
	double y1 = gsl_vector_get(x,1);

	double H12And21 = 4*(x1+y1);
	gsl_matrix_set(H,0,0,4*y1+12*pow(x1,2)-42);
	gsl_matrix_set(H,1,1,4*x1+12*pow(y1,2)-28);
	gsl_matrix_set(H,0,1,H12And21);
	gsl_matrix_set(H,1,0,H12And21);
}


int mainHimmelblauNewton(void) {
	int n = 2;

	gsl_vector *xStart = gsl_vector_alloc(n);
	gsl_vector_set(xStart,0,5);
	gsl_vector_set(xStart,1,-7);

	double eps = 1e-6;

	int iterations = newton(
		&himmelblau,
		&himmelblauGradient,
		&himmelblauHessian,
		xStart,
		eps
	);

	gsl_vector *df = gsl_vector_alloc(n);
	himmelblauGradient(xStart,df);

	for(int i = 0; i < n; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}
	printf("f(xStart) = %g\n",himmelblau(xStart));
	printf("|df(xStart)| = %g\n",gsl_blas_dnrm2(df));
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(df);

	return 0;
}

int mainHimmelblauQuasiNewton(void) {
	int n = 2;

	gsl_vector *xStart = gsl_vector_alloc(n);
	gsl_vector_set(xStart,0,5);
	gsl_vector_set(xStart,1,-7);

	double eps = 1e-6;

	int iterations = quasiNewton(
		&himmelblau,
		&himmelblauGradient,
		xStart,
		eps
	);

	gsl_vector *df = gsl_vector_alloc(n);
	himmelblauGradient(xStart,df);

	for(int i = 0; i < n; i++) {
		printf("xStart(%d) = %g\n",i,gsl_vector_get(xStart,i));
	}
	printf("f(xStart) = %g\n",himmelblau(xStart));
	printf("|df(xStart)| = %g\n",gsl_blas_dnrm2(df));
	printf("Did %d iterations.\n",iterations);

	gsl_vector_free(xStart);
	gsl_vector_free(df);

	return 0;
}


int main(void) {
	printf("Quasi Newton: Rosenbrock minimization\n");
	mainRosenbrockQuasiNewton();
	printf("Newton: Rosenbrock minimization\n");
	mainRosenbrockNewton();
	printf("\n");

	printf("Quasi Newton: Himmelblau minimization\n");
	mainHimmelblauQuasiNewton();
	printf("Newton: Himmelblau minimization\n");
	mainHimmelblauNewton();


	return 0;
}
