#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "qrMethods.h"

double lineSearchBacktracking(
		double f(gsl_vector *x),
		gsl_vector *x,
		gsl_vector *gradf,
		gsl_vector *deltax
) {

	double lambda = 1, alpha = 1e-4, fx = f(x), deltaxdx;
	gsl_vector *xPlusLambdadx = gsl_vector_alloc((*x).size);
	gsl_vector *lambdadx = gsl_vector_alloc((*x).size);

	gsl_vector_memcpy(xPlusLambdadx,x);
	gsl_vector_memcpy(lambdadx,deltax);
	gsl_vector_scale(lambdadx,lambda);
	gsl_vector_add(xPlusLambdadx,lambdadx);

	gsl_blas_ddot(deltax,gradf,&deltaxdx);


	while( f(xPlusLambdadx) >= fx + alpha * lambda * deltaxdx )
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
	double f(gsl_vector *x),
	void gradient(gsl_vector *x, gsl_vector *df),
	void hessian(gsl_vector *x, gsl_matrix *H),
	gsl_vector *xStart,
	double eps
	)
{
	int n = (*xStart).size;

	gsl_vector *df = gsl_vector_alloc(n);
	gsl_matrix *H = gsl_matrix_alloc(n,n);
	gsl_vector *deltax = gsl_vector_alloc(n);
	gsl_matrix *R = gsl_matrix_alloc(n,n);

	gradient(xStart,df);

	int iterations = 0;
	while( gsl_blas_dnrm2(df) > eps ) {
		iterations++;

		hessian(xStart,H);
		gsl_matrix_scale(H,-1.0);
		qr_gs_decomp(H,R);
		gsl_vector_memcpy(deltax,df);
		qr_gs_solve(H,R,deltax);

		double lambda = lineSearchBacktracking(f,xStart,df,deltax);

		gsl_vector_scale(deltax,lambda);
		gsl_vector_add(xStart,deltax);

		gradient(xStart,df);
	}

	gsl_vector_free(df);
	gsl_matrix_free(H);
	gsl_vector_free(deltax);
	gsl_matrix_free(R);

	return iterations;
}


double lineSearchBacktrackingBroyden(
		double f(gsl_vector *x),
		gsl_vector *x,
		gsl_vector *gradf,
		gsl_vector *deltax,
		gsl_matrix *H
) {

	double lambda = 1, alpha = 1e-4, dx = 1e-6, fx = f(x), deltaxdx;
	gsl_vector *xPlusLambdadx = gsl_vector_alloc((*x).size);
	gsl_vector *lambdadx = gsl_vector_alloc((*x).size);

	gsl_vector_memcpy(xPlusLambdadx,x);
	gsl_vector_memcpy(lambdadx,deltax);
	gsl_vector_scale(lambdadx,lambda);
	gsl_vector_add(xPlusLambdadx,lambdadx);

	gsl_blas_ddot(deltax,gradf,&deltaxdx);


	while( f(xPlusLambdadx) >= fx + alpha * lambda * deltaxdx )
	{
		lambda /= 2;
		gsl_vector_memcpy(xPlusLambdadx,x);
		gsl_vector_memcpy(lambdadx,deltax);
		gsl_vector_scale(lambdadx,lambda);
		gsl_vector_add(xPlusLambdadx,lambdadx);

		if( lambda < dx ) {
			gsl_matrix_set_identity(H);
			break;
		}
	}

	gsl_vector_free(xPlusLambdadx);
	gsl_vector_free(lambdadx);

	return lambda;
}


int quasiNewton(
	double f(gsl_vector *x),
	void gradient(gsl_vector *x, gsl_vector *df),
	gsl_vector *xStart,
	double eps
	)
{
	int n = (*xStart).size;

	gsl_vector *df = gsl_vector_alloc(n);
	gsl_vector *y = gsl_vector_alloc(n);
	gsl_vector *s = gsl_vector_alloc(n);
	gsl_matrix *H = gsl_matrix_alloc(n,n);
	gsl_vector *Hy = gsl_vector_alloc(n);
	gsl_matrix_set_identity(H);

	double HyTs;

	gradient(xStart,df);

	int iterations = 0;
	while( gsl_blas_dnrm2(df) > eps ) {
		iterations++;

		gsl_blas_dgemv(CblasNoTrans,-1.0,H,df,0.0,s);

		double lambda = lineSearchBacktrackingBroyden(f,xStart,df,s,H);

		gsl_vector_scale(s,lambda);
		gsl_vector_add(xStart,s);
		gradient(xStart,y);
		gsl_vector_sub(y,df);

		gsl_blas_dgemv(CblasNoTrans,1.0,H,y,0.0,Hy);
		gsl_blas_dgemv(CblasNoTrans,1.0,H,s,0.0,df);
		gsl_blas_ddot(Hy,s,&HyTs);
		gsl_vector_sub(s,Hy);
		gsl_blas_dger(1.0/HyTs,s,df,H);

		gradient(xStart,df);
	}

	gsl_vector_free(df);
	gsl_vector_free(y);
	gsl_vector_free(s);
	gsl_matrix_free(H);
	gsl_vector_free(Hy);

	return iterations;
}

double simplexDiff(gsl_matrix *simplex, int n) {
	double greatestDiff = 0, currDiff;

	for(int j = 1; j < n+1; j++) {
		currDiff = 0;
		for(int i = 0; i < n; i++) {
			currDiff += pow(
				gsl_matrix_get(simplex,i,j)-
				gsl_matrix_get(simplex,i,0),
				2);
		}
		if(currDiff > greatestDiff) {
			greatestDiff = currDiff;
		}
	}
	return sqrt(greatestDiff);
}

void reflection(gsl_vector *highest, gsl_vector *centroid, gsl_vector *reflected) {
	gsl_vector_memcpy(reflected,centroid);
	gsl_vector_scale(reflected,2);
	gsl_vector_sub(reflected,highest);
}

void expansion(gsl_vector *highest, gsl_vector *centroid, gsl_vector *expanded) {
	gsl_vector_memcpy(expanded,centroid);
	gsl_vector_sub(expanded,highest);
	gsl_vector_scale(expanded,2);
	gsl_vector_add(expanded,centroid);
}

void contraction(gsl_vector *highest, gsl_vector *centroid, gsl_vector *contracted) {
	gsl_vector_memcpy(contracted,centroid);
	gsl_vector_add(contracted,highest);
	gsl_vector_scale(contracted,0.5);
}

void reduction(gsl_vector *lowest, gsl_matrix *simplex, int lowIndex, int n) {
	for(int j = 0; j < n+1; j++) {
		if( j == lowIndex ) {
			continue;
		}
		gsl_vector_view colJ = gsl_matrix_column(simplex,j);
		gsl_vector_add(&colJ.vector,lowest);
		gsl_vector_scale(&colJ.vector,0.5);
	}
}

void findHighestLowest(double f(gsl_vector *p), gsl_matrix *simplex, gsl_vector *highest, gsl_vector *lowest, int *highIndex, int *lowIndex, int n) {
	double currF, highF = -INFINITY, lowF = INFINITY;

	for(int j = 0; j < n+1; j++) {
		gsl_vector_view colJ = gsl_matrix_column(simplex,j);
		currF = f(&colJ.vector);
		if( currF > highF ) {
			highF = currF;
			*highIndex = j;
		}
		if( currF < lowF ) {
			lowF = currF;
			*lowIndex = j;
		}
	}
	gsl_vector_view highCol = gsl_matrix_column(simplex,*highIndex);
	gsl_vector_memcpy(highest,&highCol.vector);

	gsl_vector_view lowCol = gsl_matrix_column(simplex,*lowIndex);
	gsl_vector_memcpy(lowest,&lowCol.vector);
}


void findCentroid(gsl_matrix *simplex, gsl_vector *centroid, int highIndex, int n) {

	for(int j = 0; j < n+1; j++) {
		if( j == highIndex ) {
			continue;
		}
		gsl_vector_view colJ = gsl_matrix_column(simplex,j);
		gsl_vector_add(centroid,&colJ.vector);
	}
	gsl_vector_scale(centroid,1.0/n);

}

int downhillSimplex(
	double f(gsl_vector *p),
	gsl_matrix *simplex,
	int n,
	double minDist
	)
{

	gsl_vector *reflected = gsl_vector_alloc(n);
	gsl_vector *other = gsl_vector_alloc(n);
	gsl_vector *centroid = gsl_vector_calloc(n);
	gsl_vector *highest = gsl_vector_alloc(n);
	gsl_vector *lowest = gsl_vector_alloc(n);
	int lowIndex, highIndex, iterations = 0;

	while( simplexDiff(simplex,n) > minDist ) {
		iterations++;

		findHighestLowest(f,simplex,highest,lowest,&highIndex,&lowIndex,n);
		findCentroid(simplex,centroid,highIndex,n);

		reflection(highest,centroid,reflected);
		if( f(reflected) < f(lowest) ) {
			expansion(highest,centroid,other);
			gsl_vector_view highVec = gsl_matrix_column(simplex,highIndex);
			if( f(other) < f(reflected) ) {
				gsl_vector_memcpy(&highVec.vector,other);
				printf("Updated expansion\n");
			} else {
				gsl_vector_memcpy(&highVec.vector,reflected);
				printf("Updated reflection 1\n");
			}
		} else {
			if( f(reflected) < f(highest) ) {
				gsl_vector_view highVec = gsl_matrix_column(simplex,highIndex);
				gsl_vector_memcpy(&highVec.vector,reflected);
				printf("Updated reflection 2\n");
			} else {
				contraction(highest,centroid,other);
				if( f(other) < f(highest) ) {
					gsl_vector_view highVec = gsl_matrix_column(simplex,highIndex);
					gsl_vector_memcpy(&highVec.vector,other);
					printf("Updated contraction\n");
				} else {
					reduction(lowest,simplex,lowIndex,n);
					printf("Updated reduction\n");
				}
			}
		}
	}

	return iterations;
}
