#include "include.h"


int main() {
    double x[NDIM] = {0.0, 1.3, M_PI / 2.0, 0.0};
    double g[NDIM][NDIM];
    double g_inv[NDIM][NDIM];
    double gamma[NDIM][NDIM][NDIM];

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
    double h = 1e-9;
    calculate_christoffel(x, h, gamma, g, g_inv);
    print_christoffel_matrix(gamma);
    double Riemann[NDIM][NDIM][NDIM][NDIM];
	double Gamma_plus_h[NDIM][NDIM][NDIM];
	double Gamma_minus_h[NDIM][NDIM][NDIM];
	initialize_riemann_tensor(Riemann);
	calculate_christoffel(x, h, Gamma_plus_h, g, g_inv);
	calculate_christoffel(x, -h, Gamma_minus_h, g, g_inv);
	calculate_riemann(gamma, Gamma_plus_h, Gamma_minus_h, Riemann, h);
    print_riemann(Riemann);
	check_riemann_symmetries(Riemann, TOLERANCE);
    double Ricci[NDIM][NDIM];
    calculate_ricci(Riemann, g_inv, Ricci);
    
    printf("\nRicci Tensor:\n");
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            printf("%12.6f\t", Ricci[i][j]);
        }
        printf("\n");
    }

    double Ricci_scalar = 0.0;
    for (int mu = 0; mu < NDIM; mu++) {
        Ricci_scalar += g_inv[mu][mu] * Ricci[mu][mu];
    }
    printf("\nRicci Scalar R = %f\n", Ricci_scalar);

    double G[NDIM][NDIM];
    calculate_einstein_tensor(Ricci, g, Ricci_scalar, G);
    printf("\nEinstein Tensor:\n");
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            printf("%12.6f\t", G[i][j]);
        }
        printf("\n");
    }
	check_symmetry_christoffel(gamma);
	double K = calculate_kretschmann(Riemann, g_inv);
	printf("\nKretschmann Scalar K = %f\n", K);


    return 0;
}
