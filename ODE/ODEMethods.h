#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define DRIVER_FAIL (int) -1

void rkStep45(
	double t,
	double h,
	gsl_vector *y,
	void f(double t, gsl_vector *y, gsl_vector *dydt),
	gsl_vector *yh,
	gsl_vector *err
);

void driver(
	double *t,
	double b,
	double *h,
	gsl_vector *y,
	double acc,
	double eps,
	void stepper(
		double t, double h, gsl_vector *y,
		void f(double t, gsl_vector *y, gsl_vector *dydt),
		gsl_vector *yh, gsl_vector *err
		),
	void f(double t, gsl_vector *y, gsl_vector *dydt)
);

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
);
