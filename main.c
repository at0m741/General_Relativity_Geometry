#include "include.h"

#define NDIM 4
#define N_r 40    
#define N_theta 50 
#define N_phi 50 

void contract_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM], double Ricci[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    memset(Ricci, 0, sizeof(double) * NDIM * NDIM);

    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    Ricci[mu][nu] += g_inv[rho][sigma] * Riemann[rho][sigma][mu][nu];
                }
            }
        }
    }
}


void calculate_energy_momentum_tensor(double g_inv[NDIM][NDIM], double rho, double p, double u[NDIM], double T[NDIM][NDIM]) {
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            T[mu][nu] = (rho + p) * u[mu] * u[nu] + p * g_inv[mu][nu];
        }
    }
}


int verify_normalization(double g[NDIM][NDIM], double u[NDIM]) {
    double result = 0.0;
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            result += g[mu][nu] * u[mu] * u[nu];
        }
    }
    if (fabs(result + 1.0) < TOLERANCE) {
		printf("La condition de normalisation est satisfaite : g_uv u^u u^v = %f\n", result);	
		return 1;
    } else {
        printf("La condition de normalisation n'est pas satisfaite : g_uv u^u u^v = %f\n", result);
        return 0;
    }
}

void calculate_quadrivector_orbit(double r, double u[NDIM], double g[NDIM][NDIM]) {
    double M = 1.0;
    double Omega = sqrt(M / (r * r * r)); 
    double g_tt = g[0][0];
    double g_phiphi = g[3][3];
    
    u[0] = 1.0 / sqrt(- (g_tt + g_phiphi * Omega * Omega));
    u[1] = 0.0; 
    u[2] = 0.0; 
    u[3] = Omega * u[0];
	for (int i = 0; i < NDIM; i++) {
		printf("u[%d] = %f\n", i, u[i]);
	}
}

double calculate_dT_partial(double T[NDIM][NDIM], int mu, int nu, double delta) {
    double T_plus = (T[mu + 1][nu]); 
    double T_minus = (T[mu][nu]); 
    
    return (T_plus - T_minus) / (2 * DELTA);
}



void calculate_covariant_divergence(double T[NDIM][NDIM], double Gamma[NDIM][NDIM][NDIM], double div_T[NDIM], double delta) {
    memset(div_T, 0, sizeof(double) * NDIM);
    
    for (int nu = 0; nu < NDIM; nu++) {
        for (int mu = 0; mu < NDIM; mu++) {
            double dT_partial = calculate_dT_partial(T, mu, nu, DELTA);

            div_T[nu] += dT_partial;
            for (int lambda = 0; lambda < NDIM; lambda++) {
                div_T[nu] += Gamma[mu][mu][lambda] * T[lambda][nu]; 
                div_T[nu] += Gamma[nu][mu][lambda] * T[mu][lambda];
				printf("Gamma[%d][%d][%d] = %f\n", mu, mu, lambda, Gamma[mu][mu][lambda]);
				printf("div_T[%d] = %f\n", nu, div_T[nu]);
            }
        }
    }
}




int main() {
    double x[NDIM] = {0.0, 10.0, M_PI / 2.0, 0.0};
    double g[NDIM][NDIM], g_inv[NDIM][NDIM];
    double gamma[NDIM][NDIM][NDIM], Riemann[NDIM][NDIM][NDIM][NDIM];
    double Gamma_plus_h[NDIM][NDIM][NDIM], Gamma_minus_h[NDIM][NDIM][NDIM];
    double Ricci[NDIM][NDIM], G[NDIM][NDIM];
    double h = 1e-9, Lambda = 1.0;

    calculate_metric(x, g, g_inv);
    verify_metric(g, g_inv);

    printf("Metric g:\n");
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            printf("g[%d][%d] = %f\t", i, j, g[i][j]);
        }
        printf("\n");
    }
    printf("\nInverse Metric g_inv:\n");
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            printf("g_inv[%d][%d] = %f\t", i, j, g_inv[i][j]);
        }
        printf("\n");
    }

    calculate_christoffel(x, h, gamma, g, g_inv, "Kerr");
    print_christoffel_matrix(gamma);

    initialize_riemann_tensor(Riemann);
    calculate_christoffel(x, h, Gamma_plus_h, g, g_inv, "Kerr");
    calculate_christoffel(x, -h, Gamma_minus_h, g, g_inv, "Kerr");
    calculate_riemann(gamma, Gamma_plus_h, Gamma_minus_h, Riemann, h);
    print_riemann(Riemann);
    check_riemann_symmetries(Riemann, TOLERANCE);
    check_symmetry_christoffel(gamma);
	
	double Ricci2[NDIM][NDIM];
	contract_riemann(Riemann, Ricci2, g_inv);
	printf("\nRicci Tensor:\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			printf("%12.6f\t", Ricci2[i][j]);
		}
		printf("\n");
	}
	double Ricci_scalar = 0.0;
	for (int i = 0; i < NDIM; i++) {
		Ricci_scalar += Ricci2[i][i] * g_inv[i][i];
	}
	printf("\nRicci Scalar: %f\n", Ricci_scalar);
	calculate_einstein_tensor(Ricci2, g, Ricci_scalar, G);
	printf("\nEinstein Tensor:\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			printf("%12.6f\t", G[i][j]);
		}
		printf("\n");
	}

    double K = calculate_kretschmann(Riemann, g_inv);
	double u[NDIM];

	double rho = 1.0;
	double p = 0.1;
	calculate_quadrivector_orbit(x[1], u, g);
	calculate_energy_momentum_tensor(g_inv, rho, p, u, G);
	printf("\nEnergy-Momentum Tensor:\n");
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			printf("%12.6f\t", G[i][j]);
		}
		printf("\n");
	}
	double T[NDIM][NDIM];
	double div_T[NDIM];
	verify_normalization(g, u);
	calculate_covariant_divergence(T, gamma, div_T, DELTA);
	printf("\nDivergence of the Energy-Momentum Tensor:\n");
	for (int i = 0; i < NDIM; i++) {
		printf("%12.6f\t", div_T[i]);
	}
    printf("\nKretschmann Scalar K = %f\n", K);

    return 0;
}
