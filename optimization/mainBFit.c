#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int quasiNewton(
	double f(gsl_vector *x),
	void gradient(gsl_vector *x, gsl_vector *df),
	gsl_vector *xStart,
	double eps
);


double function(double t, gsl_vector *params) {
	double A = gsl_vector_get(params,0);
	double T = gsl_vector_get(params,1);
	double B = gsl_vector_get(params,2);

	return A*exp(-t/T)+B;
}


double objective(gsl_vector *x) {
	double A = gsl_vector_get(x,0);
	double T = gsl_vector_get(x,1);
	double B = gsl_vector_get(x,2);

	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int n = sizeof(t)/sizeof(t[0]);

	double result = 0;
	for( int i = 0; i < n; i++) {
		result += pow( (A*exp(-t[i]/T)+B-y[i])/e[i] , 2 );
	}

	return result;
}

void objectiveGradient(gsl_vector *x, gsl_vector* df) {
	double A = gsl_vector_get(x,0);
	double T = gsl_vector_get(x,1);
	double B = gsl_vector_get(x,2);

	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int n = sizeof(t)/sizeof(t[0]);


	double result1 = 0, result2 = 0, result3 = 0, temp;
	for( int i = 0; i < n; i++) {
		temp = 2*(A*exp(-t[i]/T)+B-y[i])/pow(e[i],2);
		result1 += exp(-t[i]/T)*temp;
		result2 += A*t[i]/pow(T,2)*exp(-t[i]/T)*temp;
		result3 += temp;
	}

	gsl_vector_set(df,0,result1);
	gsl_vector_set(df,1,result2);
	gsl_vector_set(df,2,result3);
}


int main(void) {
	int noParams = 3;

	gsl_vector *xStart = gsl_vector_alloc(noParams);
	gsl_vector_set(xStart,0,5);
	gsl_vector_set(xStart,1,10);
	gsl_vector_set(xStart,2,0.5);

	double eps = 1e-6;

	quasiNewton(
		&objective,
		&objectiveGradient,
		xStart,
		eps
	);

	double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77};
	double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46};
	double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24};
	int n = sizeof(t)/sizeof(t[0]);
	printf("# t \t y \t e\n");
	for(int i = 0; i < n; i++) {
		printf("%g \t %g \t %g \n",t[i],y[i],e[i]);
	}

	printf("\n\n# t \t y \t with A = %.2g, T = %.2g, B = %.2g\n",
			gsl_vector_get(xStart,0),
			gsl_vector_get(xStart,1),
			gsl_vector_get(xStart,2));
	for(double t = 0; t < 10; t += 0.01) {
		printf("%g \t %g\n",t,function(t,xStart));
	}

	gsl_vector_free(xStart);
	return 0;
}
