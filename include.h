#ifndef INCLUDE_H
#define INCLUDE_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <immintrin.h>

#define NDIM 4
#define TOLERANCE 1e-5
#define DELTA 1e-7
#define NX 100
#define NY 100
#define NZ 100

/* Matrix functions */

double determinant3x3(double mat[3][3]);
double determinant4x4(double mat[NDIM][NDIM]);
void cofactor(double mat[NDIM][NDIM], double cofactorMat[NDIM][NDIM]);
void transpose(double mat[NDIM][NDIM], double transposed[NDIM][NDIM]);
int inverse_matrix(double mat[NDIM][NDIM], double inverse[NDIM][NDIM]);
void check_inverse(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);

/* Metric functions */

void verify_metric(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
void calculate_metric(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]); 
void calculate_christoffel(double X[NDIM], double h, double gamma[NDIM][NDIM][NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM], char *metric);
void check_symmetry_christoffel(double gamma[NDIM][NDIM][NDIM]);
void calculate_sphere_metric(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]);

/* Tensor functions */

void initialize_riemann_tensor(double R[NDIM][NDIM][NDIM][NDIM]); 
void calculate_riemann(double Gamma[NDIM][NDIM][NDIM], double Gamma_plus_h[NDIM][NDIM][NDIM], double Gamma_minus_h[NDIM][NDIM][NDIM], double Riemann[NDIM][NDIM][NDIM][NDIM], double h);
void calculate_ricci(double Riemann[NDIM][NDIM][NDIM][NDIM], double g_inv[NDIM][NDIM], double Ricci[NDIM][NDIM]); 
void calculate_einstein_tensor(double Ricci[NDIM][NDIM], double g[NDIM][NDIM], double Ricci_scalar, double G[NDIM][NDIM]);
void check_riemann_symmetries(double R[NDIM][NDIM][NDIM][NDIM], double tol); 
double centered_difference(double Gamma_plus_h[NDIM][NDIM][NDIM], double Gamma_minus_h[NDIM][NDIM][NDIM], int rho, int mu, int nu, double h);
void contract_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM], double Ricci[NDIM][NDIM]); 
/*utils*/


double calculate_kretschmann(double Riemann[NDIM][NDIM][NDIM][NDIM], double g_inv[NDIM][NDIM]);

/* Print functions */

void print_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM]);
void print_christoffel_matrix(double gamma[NDIM][NDIM][NDIM]);
void print_matrix(const char* name, double mat[NDIM][NDIM]);
void Kerr_surfaces(double x[NDIM], double a);

/* Grid functions */
void init_grid();
void compute_on_grid();
#endif

