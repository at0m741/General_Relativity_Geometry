#include "include.h"

void contract_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM], double Ricci[NDIM][NDIM]) {
    memset(Ricci, 0, sizeof(double) * NDIM * NDIM);

    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int lambda = 0; lambda < NDIM; lambda++) {
                Ricci[mu][nu] += Riemann[lambda][mu][lambda][nu] ;
            }
        }
    }

    printf("Vérification des symétries du tenseur de Ricci avec tolérance = %e :\n", TOLERANCE);
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            if (fabs(Ricci[mu][nu] - Ricci[nu][mu]) > TOLERANCE) {
                printf("Symétrie violée à Ricci[%d][%d] et Ricci[%d][%d]: %e vs %e\n", 
                       mu, nu, nu, mu, Ricci[mu][nu], Ricci[nu][mu]);
            }

            if (fabs(Ricci[mu][nu]) < TOLERANCE) {
                printf("Ricci[%d][%d] est proche de zéro : %e\n", mu, nu, Ricci[mu][nu]);
            }
        }
    }
    printf("Vérification du tenseur de Ricci terminée.\n");
}

int main() {
    double x[NDIM] = {0.0, 1.3, M_PI / 2.0, 0.0};
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

    calculate_christoffel(x, h, gamma, g, g_inv);
    print_christoffel_matrix(gamma);

    initialize_riemann_tensor(Riemann);
    calculate_christoffel(x, h, Gamma_plus_h, g, g_inv);
    calculate_christoffel(x, -h, Gamma_minus_h, g, g_inv);
    calculate_riemann(gamma, Gamma_plus_h, Gamma_minus_h, Riemann, h);
    print_riemann(Riemann);
    check_riemann_symmetries(Riemann, TOLERANCE);
    check_symmetry_christoffel(gamma);

	double Ricci2[NDIM][NDIM];
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
    printf("\nKretschmann Scalar K = %f\n", K);

    return 0;
}
