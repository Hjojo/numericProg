typedef struct {int n; double *x,*y,*b,*c;} qspline;

qspline * qspline_alloc(int n, double *x, double *y);
double qspline_evaluate(qspline *s, double z);
double qspline_derivative(qspline *s, double z);
double qspline_integral(qspline *s, double z);
void qspline_free(qspline *s);
