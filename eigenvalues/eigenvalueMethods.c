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

int jacobiEliminateRowP(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal, int p, double phase) {
	int n = (*A).size1;
	assert( p >= 0 && p < n-1 );

	double c, s, theta, App, Aqq, Apq;

	int iterate = 1, noOfRotations = 0;
	while( iterate ) {
		iterate = 0;
		for(int q = p+1; q < n; q++) {
			App = gsl_vector_get(eigVal,p);
			Aqq = gsl_vector_get(eigVal,q);
			Apq = gsl_matrix_get(A,p,q);
			theta = phase+atan2(2*Apq, Aqq - App) / 2;
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

	return noOfRotations;
}

int jacobiNumberEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal, int noEigVal, double phase) {
	int n = (*A).size1;
	assert( n == (*A).size2 );
	assert( noEigVal > 0 && noEigVal <= n );

	for(int i = 0; i < n; i++) {
		gsl_vector_set(eigVal,i,gsl_matrix_get(A,i,i));
		for(int j = 0; j < n; j++) {
			gsl_matrix_set(eigVec,i,j,(i == j) ? 1 : 0);
		}
	}

	int noOfRotations = 0;
	for(int i = 0; i < noEigVal && i < n-1; i++) {
		noOfRotations += jacobiEliminateRowP(A,eigVec,eigVal,i,phase);
	}
	return noOfRotations;
}

int jacobiNumberLowestEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal, int noEigVal) {
	double phase = 0;
	return jacobiNumberEigenvalue(A,eigVec,eigVal,noEigVal,phase);
}

int jacobiNumberHighestEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal, int noEigVal) {
	double phase = M_PI/2;
	return jacobiNumberEigenvalue(A,eigVec,eigVal,noEigVal,phase);
}

void largestElementInEachUpperTriangularRow(const gsl_matrix *A, gsl_vector *L) {
	int n = (*L).size;
	assert( n == (*A).size1-1 );
	int m = (*A).size2;

	double max = -INFINITY, curr;
	for(int i = 0; i < n; i++) {
		//max = gsl_matrix_get(A,i,i+1);
		//gsl_vector_set(L,i,i+1);
		for(int j = i+1; j < m; j++) {
			curr = gsl_matrix_get(A,i,j);
			if(curr > max) {
				max = curr;
				gsl_vector_set(L,i,j);
			}
		}
	}
}


void findLargestElementIndex(int *p, int *q, const gsl_matrix *A, const gsl_vector *L) {
	int n = (*L).size;

	double max = -INFINITY, curr;
	for(int i = 0; i < n; i++) {
		curr = gsl_matrix_get(A,i,gsl_vector_get(L,i));
		if( curr > max ) {
			max = curr;
			*p = i;
			*q = gsl_vector_get(L,i);
		}
	}
}

void updateLargestElementInRow(int i, int j, const gsl_matrix *A, gsl_vector *L) {
	double currMax = gsl_matrix_get(A,i,gsl_vector_get(L,i));
	double newElement = gsl_matrix_get(A,i,j);
	if( newElement > currMax ) {
		gsl_vector_set(L,i,j);
	}
}

int jacobiClassicEigenvalue(gsl_matrix *A, gsl_matrix *eigVec, gsl_vector *eigVal) {
	int n = (*A).size1;
	assert( n == (*A).size2 );

	gsl_vector *L = gsl_vector_alloc(n-1);
	largestElementInEachUpperTriangularRow(A,L);

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

		int p, q;
		findLargestElementIndex(&p,&q,A,L);

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
				updateLargestElementInRow(i,p,A,L);
				updateLargestElementInRow(i,q,A,L);
			}

			for(int i = p + 1; i < q; i++) {
				double Api = gsl_matrix_get(A,p,i);
				double Aiq = gsl_matrix_get(A,i,q);
				gsl_matrix_set(A,p,i,c*Api-s*Aiq);
				gsl_matrix_set(A,i,q,s*Api+c*Aiq);
				updateLargestElementInRow(p,i,A,L);
				updateLargestElementInRow(i,q,A,L);
			}

			for(int i = q + 1; i < n; i++) {
				double Api = gsl_matrix_get(A,p,i);
				double Aqi = gsl_matrix_get(A,q,i);
				gsl_matrix_set(A,p,i,c*Api-s*Aqi);
				gsl_matrix_set(A,q,i,s*Api+c*Aqi);
				updateLargestElementInRow(p,i,A,L);
				updateLargestElementInRow(q,i,A,L);
			}

			for(int i = 0; i < n; i++) {
				double Vip = gsl_matrix_get(eigVec,i,p);
				double Viq = gsl_matrix_get(eigVec,i,q);
				gsl_matrix_set(eigVec,i,p,c*Vip-s*Viq);
				gsl_matrix_set(eigVec,i,q,s*Vip+c*Viq);
			}
		}
	}

	gsl_vector_free(L);

	return noOfRotations;
}
