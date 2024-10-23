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
    const double threshold = 1e-10; 
    printf("\nRiemann Tensor (Non-zero components):\n"); 
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    double value = Riemann[rho][sigma][mu][nu];
                    if (fabs(value) > threshold) {
                        printf("Riemann[%d][%d][%d][%d] = %12.6f\n", rho, sigma, mu, nu, value);
                    }
                }
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

void generate_3sphere_points(double a, int num_steps) {
    double dChi = M_PI / num_steps;
    double dTheta = M_PI / num_steps;
    double dPhi = 2 * M_PI / num_steps;

    FILE *fp = fopen("3sphere_points.txt", "w");

    for (double chi = 0; chi <= M_PI; chi += dChi) {
        for (double theta = 0; theta <= M_PI; theta += dTheta) {
            for (double phi = 0; phi <= 2 * M_PI; phi += dPhi) {
                double x = a * sin(chi) * sin(theta) * cos(phi);
                double y = a * sin(chi) * sin(theta) * sin(phi);
                double z = a * sin(chi) * cos(theta);
                double w = a * cos(chi);

                double denom = 1 - (w / a);
                if (fabs(denom) > 1e-6) { 
                    double xp = x / denom;
                    double yp = y / denom;
                    double zp = z / denom;

                    fprintf(fp, "%f %f %f\n", xp, yp, zp);
                }
            }
        }
    }

    fclose(fp);
}
