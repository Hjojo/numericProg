#ifndef HAVE_CSPLINE_H
typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;

cspline * cspline_alloc(int n, double *x, double *y);
double cspline_evaluate(cspline *s, double z);
double cspline_derivative(cspline *s, double z);
double cspline_integral(cspline *s, double z);
void cspline_free(cspline *s);

#define HAVE_CSPLINE_H
#endif
