

double adaptIntegClosed(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
);

double adaptIntegOpen(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
);

double adaptIntegClosedWithInf(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
);

double adaptIntegOpenWithInf(
	double f(double),
	double a,
	double b,
	double abs,
	double eps,
	double *err
);
