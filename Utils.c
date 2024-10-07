#include "include.h"


void check_symmetry_christoffel(double gamma[NDIM][NDIM][NDIM]) {
    printf("\nChristoffel symbols verif with delta = %e:\n", TOLERANCE);
    for (int lambda = 0; lambda < NDIM; lambda++) {
        for (int mu = 0; mu < NDIM; mu++) {
            for (int nu = mu; nu < NDIM; nu++) { 
                double diff = fabs(gamma[lambda][mu][nu] - gamma[lambda][nu][mu]);
                if (diff > TOLERANCE) {
                    printf("!symmetry at Gamma^%d_%d%d: |%f - %f| = %f\n", 
                           lambda, mu, nu, gamma[lambda][mu][nu], gamma[lambda][nu][mu], diff);
                }
            }
        }
    }
    printf("Symmetric !\n");
}


void initialize_riemann_tensor(double R[NDIM][NDIM][NDIM][NDIM]) {
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    R[mu][nu][rho][sigma] = 0.0; 
                }
            }
        }
    }
}


void print_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM]) {
    printf("\nRiemann Tensor:\n"); 
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            printf("\nRiemann[%d][%d]:\n", rho, sigma);
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    printf("%12.6f\t", Riemann[rho][sigma][mu][nu]);
                }
                printf("\n");
            }
        }
    }
}


void print_christoffel_matrix(double gamma[NDIM][NDIM][NDIM]) {
    printf("\nChristoffel Symbols:\n");
    for (int lambda = 0; lambda < NDIM; lambda++) {
        printf("\nGamma^%d:\n", lambda);
        for (int mu = 0; mu < NDIM; mu++) {
            for (int nu = 0; nu < NDIM; nu++) {
                printf("%12.6f\t", gamma[lambda][mu][nu]);
            }
            printf("\n");
        }
    }
}

