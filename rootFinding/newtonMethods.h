#ifndef HAVE_NEWTONMETHODS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int newton(
	void f(gsl_vector *x, gsl_vector *fx),
	int noOfFuns,
	gsl_vector *xStart,
	double dx,
	double epsilon
);

int newtonWithJacobian(
	void f(gsl_vector *x, gsl_vector *fx),
	int noOfFuns,
	void jacobian(gsl_vector *x, gsl_matrix *J),
	gsl_vector *xStart,
	double epsilon
);

#define HAVE_NEWTONMETHODS_H
#endif
