#include <math.h>
#include <assert.h>
#include "adaptInteg.h"

double adaptClosedRectTrapz(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err,
	double f1,
	double f3,
	int noRecursions
) {
	assert(noRecursions < 1e6);

	double f2 = f((a+b)/2);

	double Q = (f1+2*f2+f3)*(b-a)/4;
	double q = (f1+f2+f3)*(b-a)/2;
	*err = fabs(Q-q);
	double tol = abs + eps*fabs(Q);

	if(*err < tol) {
		return Q;
	} else {
		double err1, err2;
		double Q1 = adaptClosedRectTrapz(f,a,(a+b)/2,abs/sqrt(2),eps,&err1,f1,f2,noRecursions+1);
		double Q2 = adaptClosedRectTrapz(f,(a+b)/2,b,abs/sqrt(2),eps,&err2,f2,f3,noRecursions+1);
		*err = sqrt( pow(err1,2) + pow(err2,2) );
		return Q1+Q2;
	}
}


double adaptOpenRectTrapz(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err,
	double f2,
	double f3,
	int noRecursions
) {
	assert(noRecursions < 1e6);

	double f1 = f(a+(b-a)/6);
	double f4 = f(a+(b-a)*5/6);

	double Q = (2*f1+f2+f3+2*f4)*(b-a)/6;
	double q = (f1+f2+f3+f4)*(b-a)/4;
	*err = fabs(Q-q);
	double tol = abs + eps*fabs(Q);

	if(*err < tol) {
		return Q;
	} else {
		double err1, err2;
		double Q1 = adaptOpenRectTrapz(f,a,(a+b)/2,abs/sqrt(2),eps,&err1,f1,f2,noRecursions+1);
		double Q2 = adaptOpenRectTrapz(f,(a+b)/2,b,abs/sqrt(2),eps,&err2,f3,f4,noRecursions+1);
		*err = sqrt( pow(err1,2) + pow(err2,2) );
		return Q1+Q2;
	}
}


double adaptIntegClosed(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
) {
	double f1 = f(a);
	double f3 = f(b);

	int noRecursions = 0;
	return adaptClosedRectTrapz(f,a,b,abs,eps,err,f1,f3,noRecursions);
}

double adaptIntegOpen(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
) {
	double f2 = f(a+(b-a)*2/6);
	double f3 = f(a+(b-a)*4/6);

	int noRecursions = 0;
	return adaptOpenRectTrapz(f,a,b,abs,eps,err,f2,f3,noRecursions);
}

double adaptIntegClosedWithInf(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
) {

	if( (isinf(a) == -1) & (isinf(b) == 1) ) {
		double f(double t) {
			return f( t/(1-pow(t,2)) ) * (1+pow(t,2))/pow(1-pow(t,2),2);
		}
		a = -1; b = 1;
	} else if( isinf(a) == -1 ) {
		double f(double t) {
			return f( b-(1-t)/t )/pow(t,2);
		}
		a = 0; b = 1;
	} else if( isinf(b) == 1 ) {
		double f(double t) {
			return f( a+(1-t)/t )/pow(t,2);
		}
		a = 0; b = 1;
	}

	return adaptIntegClosed(f,a,b,abs,eps,err);

}

double adaptIntegOpenWithInf(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
) {

	if( (isinf(a) == -1) & (isinf(b) == 1) ) {
		double pmInf(double t) {
			return f( t/(1-pow(t,2)) ) * (1+pow(t,2))/pow(1-pow(t,2),2);
		}
		double aNew = -1, bNew = 1;
		return adaptIntegOpen(pmInf,aNew,bNew,abs,eps,err);
	} else if( isinf(a) == -1 ) {
		double plusInf(double t) {
			return f( b-(1-t)/t )/pow(t,2);
		}
		double aNew = 0, bNew = 1;
		return adaptIntegOpen(plusInf,aNew,bNew,abs,eps,err);
	} else if( isinf(b) == 1 ) {
		double minInf(double t) {
			return f( a+(1-t)/t )/pow(t,2);
		}
		double aNew = 0, bNew = 1;
		return adaptIntegOpen(minInf,aNew,bNew,abs,eps,err);
	}

	return adaptIntegOpen(f,a,b,abs,eps,err);
}
