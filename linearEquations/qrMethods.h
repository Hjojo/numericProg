#ifndef HAVE_QRMETHODS_H

void qr_gs_decomp(gsl_matrix *A, gsl_matrix *R);
void upperTriangular_solve(const gsl_matrix *U, gsl_vector *C);
void qr_gs_solve(const gsl_matrix *Q, const gsl_matrix *R, gsl_vector *b);
void qr_gs_inverse(const gsl_matrix *Q, const gsl_matrix *R, gsl_matrix *B);
void qr_givens_decomp(gsl_matrix *A);
void qr_givens_QTb(const gsl_matrix *QR, gsl_vector *b);
void qr_givens_solve(const gsl_matrix *QR, gsl_vector *b);
void qr_givens_inverse(const gsl_matrix *QR, gsl_matrix *B);

#define HAVE_QRMETHODS_H
#endif
