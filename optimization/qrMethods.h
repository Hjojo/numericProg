#ifndef HAVE_QRMETHODS_H

void qr_gs_decomp(gsl_matrix *A, gsl_matrix *R);
void upperTriangular_solve(const gsl_matrix *U, gsl_vector *c);
void qr_gs_solve(const gsl_matrix *Q, const gsl_matrix *R, gsl_vector *b);

#define HAVE_QRMETHODS_H
#endif
