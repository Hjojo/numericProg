#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>


int jacobiEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal) {
	int n = (*A).size1;
	assert( n == (*A).size2 );

	for(int i = 0; i < n; i++) {
		gsl_vector_set(eigVal,i,gsl_matrix_get(A,i,i));
		for(int j = 0; j < n; j++) {
			gsl_matrix_set(eigVec,i,j,(i == j) ? 1 : 0);
		}
	}

	double c, s, theta, App, Aqq, Apq;

	int iterate = 1, noOfRotations = 0;
	while( iterate ) {
		iterate = 0;
		for(int p = 0; p < n-1; p++) {
			for(int q = p+1; q < n; q++) {
				App = gsl_vector_get(eigVal,p);
				Aqq = gsl_vector_get(eigVal,q);
				Apq = gsl_matrix_get(A,p,q);
				theta = atan2(2*Apq, Aqq - App) / 2;
				c = cos(theta);
				s = sin(theta);
				double AppCalc = c*c*App-2*s*c*Apq+s*s*Aqq;
				double AqqCalc = s*s*App+2*s*c*Apq+c*c*Aqq;


				if( App != AppCalc || Aqq != AqqCalc ) {
					noOfRotations++;
					iterate = 1;
					gsl_vector_set(eigVal,p,AppCalc);
					gsl_vector_set(eigVal,q,AqqCalc);
					gsl_matrix_set(A,p,q,0.0);

					for(int i = 0; i < p; i++) {
						double Aip = gsl_matrix_get(A,i,p);
						double Aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,i,p,c*Aip-s*Aiq);
						gsl_matrix_set(A,i,q,s*Aip+c*Aiq);
					}

					for(int i = p + 1; i < q; i++) {
						double Api = gsl_matrix_get(A,p,i);
						double Aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,p,i,c*Api-s*Aiq);
						gsl_matrix_set(A,i,q,s*Api+c*Aiq);
					}

					for(int i = q + 1; i < n; i++) {
						double Api = gsl_matrix_get(A,p,i);
						double Aqi = gsl_matrix_get(A,q,i);
						gsl_matrix_set(A,p,i,c*Api-s*Aqi);
						gsl_matrix_set(A,q,i,s*Api+c*Aqi);
					}

					for(int i = 0; i < n; i++) {
						double Vip = gsl_matrix_get(eigVec,i,p);
						double Viq = gsl_matrix_get(eigVec,i,q);
						gsl_matrix_set(eigVec,i,p,c*Vip-s*Viq);
						gsl_matrix_set(eigVec,i,q,s*Vip+c*Viq);
					}
				}
			}
		}
	}

	return noOfRotations;
}
