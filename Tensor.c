#include "include.h"

#define NDIM 4
#define TOLERANCE 1e-5




double centered_difference(double Gamma_plus_h[NDIM][NDIM][NDIM], 
                           double Gamma_minus_h[NDIM][NDIM][NDIM], 
                           int rho, int mu, int nu, double h) {
    double diff = (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / h;
    printf("Centered Difference for Gamma[%d][%d][%d]: Plus = %e, Minus = %e, Diff = %e\n",
           rho, mu, nu, Gamma_plus_h[rho][mu][nu], Gamma_minus_h[rho][mu][nu], diff);
    return diff;
}

double richardson_difference(double Gamma_plus_h[NDIM][NDIM][NDIM], 
                             double Gamma_minus_h[NDIM][NDIM][NDIM], 
                             int rho, int mu, int nu, double h) {
    double diff1 = (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / h;
    double half_h = h / 2.0;

    double Gamma_plus_half_h = (Gamma_plus_h[rho][mu][nu] + Gamma_minus_h[rho][mu][nu]) / 2.0;
    double Gamma_minus_half_h = (Gamma_minus_h[rho][mu][nu] + Gamma_plus_h[rho][mu][nu]) / 2.0;
    double diff2 = (Gamma_plus_half_h - Gamma_minus_half_h) / half_h;
	/* printf("Richardson Difference for Gamma[%d][%d][%d]: Plus = %e, Minus = %e, Diff1 = %e, Diff2 = %e\n", */
		   /* rho, mu, nu, Gamma_plus_h[rho][mu][nu], Gamma_minus_h[rho][mu][nu], diff1, diff2); */
    return (4.0 * diff2 - diff1) / 3.0;
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
                       double Riemann[NDIM][NDIM][NDIM][NDIM], 
                       double h, double x[NDIM]) {
    memset(Riemann, 0, sizeof(double) * NDIM * NDIM * NDIM * NDIM);

    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    // Différences centrées pour les dérivées
                    double dGamma_mu_nu = richardson_difference(Gamma_plus_h, Gamma_minus_h, rho, mu, nu, h); 
                    double dGamma_nu_mu = richardson_difference(Gamma_plus_h, Gamma_minus_h, sigma, nu, mu, h);

                    // Terme croisé
                    double Gamma_terms = 0.0;
                    for (int lambda = 0; lambda < NDIM; lambda++) {
                        double term1 = Gamma[mu][rho][lambda] * Gamma[nu][lambda][sigma];
                        double term2 = Gamma[nu][rho][lambda] * Gamma[mu][lambda][sigma];
                        Gamma_terms += term1 - term2;
                    }

                    // Calcul final
                    Riemann[rho][sigma][mu][nu] = dGamma_mu_nu - dGamma_nu_mu + Gamma_terms;

                    // Validation analytique (pour Schwarzschild)
                    double analytical = 0.0;
                    if (rho == 1 && sigma == 0 && mu == 0 && nu == 1) {
                        analytical = 2.0 / pow(x[1], 3.0);
                    } else if (rho == 2 && sigma == 1 && mu == 1 && nu == 2) {
                        analytical = -1.0 / pow(x[1], 2.0);
                    } else if (rho == 3 && sigma == 2 && mu == 2 && nu == 3) {
                        analytical = -sin(x[2]) * cos(x[2]);
                    }

                    double error = fabs(Riemann[rho][sigma][mu][nu] - analytical);
                    if (analytical != 0 && error > TOLERANCE) {
                        printf("Discrepancy in Riemann component R[%d][%d][%d][%d]:\n", 
                               rho, sigma, mu, nu);
                        printf("Computed: %e, Analytical: %e, Error: %e\n", 
                               Riemann[rho][sigma][mu][nu], analytical, error);
                    }
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


