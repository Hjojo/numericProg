#include <gsl/gsl_blas.h>
#include <math.h>
#include "ODEMethods.h"

void rkStep45(
	double t,
	double h,
	gsl_vector *y,
	void f(double t, gsl_vector *y, gsl_vector *dydt),
	gsl_vector *yh,
	gsl_vector *err
) {
	int n = (*y).size;

	gsl_vector *K1 = gsl_vector_alloc(n);
	gsl_vector *K2 = gsl_vector_alloc(n);
	gsl_vector *K3 = gsl_vector_alloc(n);
	gsl_vector *K4 = gsl_vector_alloc(n);
	gsl_vector *K5 = gsl_vector_alloc(n);
	gsl_vector *K6 = gsl_vector_alloc(n);
	gsl_vector *yt = gsl_vector_alloc(n);

	f(t,y,K1);
	gsl_vector_memcpy(yt,y);
	gsl_blas_daxpy(h/4,K1,yt);

	f(t+1./4*h,yt,K2);
	gsl_vector_memcpy(yt,y);
	gsl_blas_daxpy(h*3./32,K1,yt);
	gsl_blas_daxpy(h*9./32,K2,yt);

	f(t+3./8*h,yt,K3);
	gsl_vector_memcpy(yt,y);
	gsl_blas_daxpy(h*1932./2197,K1,yt);
	gsl_blas_daxpy(-h*7200./2197,K2,yt);
	gsl_blas_daxpy(h*7296./2197,K3,yt);

	f(t+12./13*h,yt,K4);
	gsl_vector_memcpy(yt,y);
	gsl_blas_daxpy(h*439./216,K1,yt);
	gsl_blas_daxpy(-h*8,K2,yt);
	gsl_blas_daxpy(h*3680./513,K3,yt);
	gsl_blas_daxpy(-h*845./4104,K4,yt);

	f(t+h,yt,K5);
	gsl_vector_memcpy(yt,y);
	gsl_blas_daxpy(-h*8./27,K1,yt);
	gsl_blas_daxpy(h*2,K2,yt);
	gsl_blas_daxpy(-h*3544./2565,K3,yt);
	gsl_blas_daxpy(h*1859./4104,K4,yt);
	gsl_blas_daxpy(-h*11./40,K5,yt);

	f(t+1./2*h,yt,K6);

	gsl_vector_memcpy(yh,y);
	gsl_blas_daxpy(h*16./135,K1,yh);
	gsl_blas_daxpy(h*6656./12825,K3,yh);
	gsl_blas_daxpy(h*28561./56430,K4,yh);
	gsl_blas_daxpy(-h*9./50,K5,yh);
	gsl_blas_daxpy(h*2./55,K6,yh);

	gsl_vector_memcpy(err,y);
	gsl_blas_daxpy(h*25./216,K1,err);
	gsl_blas_daxpy(h*1408./2565,K3,err);
	gsl_blas_daxpy(h*2197./4104,K4,err);
	gsl_blas_daxpy(-h/5,K5,err);

	gsl_vector_sub(err,yh);
	gsl_vector_scale(err,-1);


	gsl_vector_free(K1);
	gsl_vector_free(K2);
	gsl_vector_free(K3);
	gsl_vector_free(K4);
	gsl_vector_free(K5);
	gsl_vector_free(K6);
	gsl_vector_free(yt);
}



void driver(
	double *t,
	double b,
	double *h,
	gsl_vector *y,
	double abs,
	double eps,
	void stepper(
		double t, double h, gsl_vector *y,
		void f(double t, gsl_vector *y, gsl_vector *dydt),
		gsl_vector *yh, gsl_vector *err
		),
	void f(double t, gsl_vector *y, gsl_vector *dydt)
) {

	int n = (*y).size;
	double a = *t, tol, normErr;
	gsl_vector *yh = gsl_vector_alloc(n);
	gsl_vector *err = gsl_vector_alloc(n);

	while( *t < b ) {
		if(*t+*h>b) {
			*h = b-*t;
		}
		stepper(*t,*h,y,f,yh,err);

		tol = ( eps*gsl_blas_dnrm2(yh) + abs ) * sqrt(*h/(b-a));
		normErr = gsl_blas_dnrm2(err);

		if( normErr < tol ) {
			*t += *h;
			gsl_vector_memcpy(y,yh);
		}

		*h *= pow(tol/normErr,0.25)*0.95;
	}

	gsl_vector_free(yh);
	gsl_vector_free(err);
}

int driverPathStoring(
	gsl_vector *tPath,
	double b,
	double *h,
	gsl_matrix *yPath,
	double abs,
	double eps,
	void stepper(
		double t, double h, gsl_vector *y,
		void f(double t, gsl_vector *y, gsl_vector *dydt),
		gsl_vector *yh, gsl_vector *err
		),
	void f(double t, gsl_vector *y, gsl_vector *dydt)
) {
	int n = (*yPath).size1;
	int maxSteps = (*tPath).size, step = 0;
	double a = gsl_vector_get(tPath,0), tol, normErr, t;
	gsl_vector *yh = gsl_vector_alloc(n);
	gsl_vector *err = gsl_vector_alloc(n);
	gsl_vector_view y, yNext;

	while( gsl_vector_get(tPath,step) < b ) {

		t = gsl_vector_get(tPath,step);
		y = gsl_matrix_column(yPath,step);

		if(t+*h>b) {
			*h = b-t;
		}
		stepper(t,*h,&y.vector,f,yh,err);

		tol = ( eps*gsl_blas_dnrm2(yh) + abs ) * sqrt(*h/(b-a));
		normErr = gsl_blas_dnrm2(err);

		if( normErr < tol ) {
			step++;
			if( step+1 > maxSteps ) {
				return DRIVER_FAIL;
			}
			gsl_vector_set(tPath,step,t+*h);
			yNext = gsl_matrix_column(yPath,step);
			gsl_vector_memcpy(&yNext.vector,yh);
		}

		*h *= pow(tol/normErr,0.25)*0.95;
	}

	gsl_vector_free(yh);
	gsl_vector_free(err);

	return step+1;
}
