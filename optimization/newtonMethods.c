#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "qrMethods.h"


/*
	Makes the jacobian matrix by the finite difference method and
	save it in the matrix J.
	Note that fx must be calculated beforehand.
*/
void jacobianMatrixDiff(
	void f(gsl_vector *x, gsl_vector *fx),
	gsl_vector *x,
	gsl_vector *fx,
	gsl_matrix *J,
	double dx
) {
	int n = (*fx).size, m = (*x).size;

	gsl_vector *dfxk = gsl_vector_alloc(n);

	for(int k = 0; k < m; k++) {
		gsl_vector_set(x,k,gsl_vector_get(x,k)+dx);

		f(x,dfxk);
		gsl_vector_sub(dfxk,fx);
		gsl_vector_scale(dfxk,1.0/dx);

		for(int i = 0; i < n; i++) {
			gsl_matrix_set(J,i,k,gsl_vector_get(dfxk,i));
		}

		gsl_vector_set(x,k,gsl_vector_get(x,k)-dx);
	}

	gsl_vector_free(dfxk);
}

double findOptimalLambda(
		void f(gsl_vector *x, gsl_vector *fx),
		gsl_vector *x,
		gsl_vector *fx,
		gsl_vector *deltax
) {

	double lambda = 1, epsilon = 1e-2;
	gsl_vector *xPlusLambdadx = gsl_vector_alloc((*x).size);
	gsl_vector *lambdadx = gsl_vector_alloc((*x).size);
	gsl_vector_memcpy(xPlusLambdadx,x);
	gsl_vector_memcpy(lambdadx,deltax);
	gsl_vector_scale(lambdadx,lambda);
	gsl_vector_add(xPlusLambdadx,lambdadx);

	while( gsl_blas_dnrm2(xPlusLambdadx) > (1-lambda/2)*gsl_blas_dnrm2(fx)
		&& lambda > epsilon )
	{
		lambda /= 2;
		gsl_vector_memcpy(xPlusLambdadx,x);
		gsl_vector_memcpy(lambdadx,deltax);
		gsl_vector_scale(lambdadx,lambda);
		gsl_vector_add(xPlusLambdadx,lambdadx);
	}

	gsl_vector_free(xPlusLambdadx);
	gsl_vector_free(lambdadx);

	return lambda;
}



int newton(
	void f(gsl_vector *x, gsl_vector *fx),
	int noOfFuns,
	gsl_vector *xStart,
	double dx,
	double epsilon
) {
	int m = (*xStart).size;

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	gsl_matrix *J = gsl_matrix_alloc(noOfFuns,m);
	gsl_matrix *R = gsl_matrix_alloc(m,m);
	gsl_vector *deltax = gsl_vector_alloc(m);

	f(xStart,fx);

	int iterations = 0;
	while( gsl_blas_dnrm2(fx)  > epsilon ) {
		iterations++;

		jacobianMatrixDiff(f,xStart,fx,J,dx);
		gsl_vector_scale(fx,-1.0);
		gsl_vector_memcpy(deltax,fx);

		qr_gs_decomp(J,R);
		qr_gs_solve(J,R,deltax);

		double lambda = findOptimalLambda(f,xStart,fx,deltax);

		gsl_vector_scale(deltax,lambda);
		gsl_vector_add(xStart,deltax);
		f(xStart,fx);
	}

	gsl_vector_free(fx);
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(deltax);

	return iterations;
}

int newtonWithJacobian(
	void f(gsl_vector *x, gsl_vector *fx),
	int noOfFuns,
	void jacobian(gsl_vector *x, gsl_matrix *J),
	gsl_vector *xStart,
	double epsilon
) {
	int m = (*xStart).size;

	gsl_vector *fx = gsl_vector_alloc(noOfFuns);
	gsl_matrix *J = gsl_matrix_alloc(noOfFuns,m);
	gsl_matrix *R = gsl_matrix_alloc(m,m);
	gsl_vector *deltax = gsl_vector_alloc(m);

	f(xStart,fx);

	int iterations = 0;
	while( gsl_blas_dnrm2(fx)  > epsilon ) {
		iterations++;

		jacobian(xStart,J);
		gsl_vector_scale(fx,-1.0);
		gsl_vector_memcpy(deltax,fx);

		qr_gs_decomp(J,R);
		qr_gs_solve(J,R,deltax);

		double lambda = findOptimalLambda(f,xStart,fx,deltax);

		gsl_vector_scale(deltax,lambda);
		gsl_vector_add(xStart,deltax);
		f(xStart,fx);
	}

	gsl_vector_free(fx);
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(deltax);

	return iterations;
}
