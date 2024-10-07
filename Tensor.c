#include "include.h"

#define NDIM 4
#define TOLERANCE 1e-5



double centered_difference(double Gamma_plus_h[NDIM][NDIM][NDIM], 
                           double Gamma_minus_h[NDIM][NDIM][NDIM], 
                           int rho, int mu, int nu, double h) {
    return (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / ( h);
}

void check_riemann_symmetries(double R[NDIM][NDIM][NDIM][NDIM], double tol) {
    int symmetric = 1;

    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    if (fabs(R[mu][nu][rho][sigma] + R[nu][mu][rho][sigma]) > tol) {
                        printf("Non-antisymmetry R[%d][%d][%d][%d] et R[%d][%d][%d][%d]\n",
                               mu, nu, rho, sigma, nu, mu, rho, sigma);
						printf("R[%d][%d][%d][%d] = %f\n", mu, nu, rho, sigma, R[mu][nu][rho][sigma]);
						printf("R[%d][%d][%d][%d] = %f\n", nu, mu, rho, sigma, R[nu][mu][rho][sigma]);
                        symmetric = 0;
                    }
                }
            }
        }
    }

    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    if (fabs(R[mu][nu][rho][sigma] + R[mu][nu][sigma][rho]) > tol) {
                        printf("Non-antisymmetry R[%d][%d][%d][%d] et R[%d][%d][%d][%d]\n",
                               mu, nu, rho, sigma, mu, nu, sigma, rho);
                        symmetric = 0;
                    }
                }
            }
        }
    }

    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    if (fabs(R[mu][nu][rho][sigma] - R[rho][sigma][mu][nu]) > tol) {
                        printf("Non-symmetry R[%d][%d][%d][%d] et R[%d][%d][%d][%d]\n",
                               mu, nu, rho, sigma, rho, sigma, mu, nu);
						printf("R[%d][%d][%d][%d] = %f\n", mu, nu, rho, sigma, R[mu][nu][rho][sigma]);
						printf("R[%d][%d][%d][%d] = %f\n", rho, sigma, mu, nu, R[rho][sigma][mu][nu]);
                        symmetric = 0;
                    }
                }
            }
        }
    }

    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    double bianchi_identity = R[mu][nu][rho][sigma] +
                                              R[mu][rho][sigma][nu] +
                                              R[mu][sigma][nu][rho];
                    if (fabs(bianchi_identity) > tol) {
                        printf("Violation of Bianchi identity : (%d, %d, %d, %d)\n", 
                               mu, nu, rho, sigma);
                        symmetric = 0;
                    }
                }
            }
        }
    }

    if (symmetric) {
        printf("all Riemann Tensor symmetries are good with tolerance : %e.\n", tol);
    }
}




void calculate_riemann(double Gamma[NDIM][NDIM][NDIM], 
                       double Gamma_plus_h[NDIM][NDIM][NDIM], 
                       double Gamma_minus_h[NDIM][NDIM][NDIM], 
                       double Riemann[NDIM][NDIM][NDIM][NDIM], double h) {
    memset(Riemann, 0, sizeof(double) * NDIM * NDIM * NDIM * NDIM);

    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    double dGamma_mu_nu = centered_difference(Gamma_plus_h, Gamma_minus_h, rho, mu, nu, h * 2.0);
                    double dGamma_nu_mu = centered_difference(Gamma_plus_h, Gamma_minus_h, rho, nu, mu, h * 2.0);

                    double Gamma_terms = 0.0;
                    for (int lambda = 0; lambda < NDIM; lambda++) {
                        Gamma_terms += Gamma[mu][rho][lambda] * Gamma[nu][lambda][sigma]
                                     - Gamma[nu][rho][lambda] * Gamma[mu][lambda][sigma];
                    }

                    Riemann[rho][sigma][mu][nu] = dGamma_mu_nu - dGamma_nu_mu + Gamma_terms;
                }
            }
        }
    }
}


double calculate_kretschmann(double Riemann[NDIM][NDIM][NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    double K = 0.0;
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    double R_contravariant = 0.0;
                    for (int alpha = 0; alpha < NDIM; alpha++) {
                        for (int beta = 0; beta < NDIM; beta++) {
                            R_contravariant += g_inv[mu][alpha] * g_inv[nu][beta] * Riemann[alpha][beta][rho][sigma];
                        }
                    }
                    K += Riemann[mu][nu][rho][sigma] * R_contravariant;
                }
            }
        }
    }
    return K;
}

void calculate_ricci(double Riemann[NDIM][NDIM][NDIM][NDIM], double g_inv[NDIM][NDIM], double Ricci[NDIM][NDIM]) {
    memset(Ricci, 0, sizeof(double) * NDIM * NDIM);
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            double sum = 0.0;
            for (int lambda = 0; lambda < NDIM; lambda++) {
                sum += g_inv[lambda][lambda] * Riemann[lambda][mu][nu][lambda];
            }
            Ricci[mu][nu] = sum;
        }
    }
}

void calculate_einstein_tensor(double Ricci[NDIM][NDIM], double g[NDIM][NDIM], double Ricci_scalar, double G[NDIM][NDIM]) {
    memset(G, 0, sizeof(double) * NDIM * NDIM);
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            G[mu][nu] = Ricci[mu][nu] - 0.5 * g[mu][nu] * Ricci_scalar;
        }
    }
}


