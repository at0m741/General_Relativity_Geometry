#include "include.h"


void verify_metric(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]) {
    int i, j, k;
    double identity[NDIM][NDIM] = {0};
    double delta;

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            identity[i][j] = 0.0;
            for (k = 0; k < NDIM; k++) {
                identity[i][j] += gcon[i][k] * gcov[k][j];
            }
        }
    }

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            if (i == j) {
                delta = 1.0;
            } else {
                delta = 0.0;
            }

            if (fabs(identity[i][j] - delta) > TOLERANCE) {
                printf("Erreur: identity[%d][%d] = %e au lieu de %e\n", i, j, identity[i][j], delta);
            }
        }
    }
    printf("\n");
    check_inverse(gcov, gcon);
    printf("\n");
    print_matrix("gcov", gcov);
    print_matrix("gcon", gcon);
    printf("\n");
    printf("Identity matrix:\n");
    print_matrix("identity", identity);

}


void calculate_metric(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    double r = x[1];
    double theta = x[2];
    double M = 1.0;   
    double sin_theta2 = sin(theta) * sin(theta);

    memset(g, 0, sizeof(double) * NDIM * NDIM);
    memset(g_inv, 0, sizeof(double) * NDIM * NDIM);

    double f = 1.0 - (2.0 * M) / r;

    g[0][0] = -f;  
    g[1][1] = 1.0 / f;
    g[2][2] = r * r; 
    g[3][3] = r * r * sin_theta2; 
	inverse_matrix(g, g_inv);
}

void calculate_minkowski_metric(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM])
{
	memset(g, 0, sizeof(double) * NDIM * NDIM);
	memset(g_inv, 0, sizeof(double) * NDIM * NDIM);

	g[0][0] = -1.0;
	g[1][1] = 1.0;
	g[2][2] = 1.0;
	g[3][3] = 1.0;

	inverse_matrix(g, g_inv);
}

void calculate_sphere_metric(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    double r = x[0];
    double chi = x[1]; 
    double theta = x[2];
    double phi = x[3]; 

    double sin_chi = sin(chi);
    double sin_chi2 = sin_chi * sin_chi;
    double sin_theta = sin(theta);
    double sin_theta2 = sin_theta * sin_theta;

    memset(g, 0, sizeof(double) * NDIM * NDIM);
    memset(g_inv, 0, sizeof(double) * NDIM * NDIM);

    g[0][0] = 1.0;       
    g[1][1] = r * r;     
    g[2][2] = r * r * sin_chi2;            
    g[3][3] = r * r * sin_chi2 * sin_theta2; 
	
	inverse_matrix(g, g_inv);
}


void calculate_christoffel(double X[NDIM], double h, double gamma[NDIM][NDIM][NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM], char *metric) {
    double tmp[NDIM][NDIM][NDIM];
    double Xh[NDIM], Xl[NDIM]; 
    double gh[NDIM][NDIM], gl[NDIM][NDIM]; 

    memset(gamma, 0, sizeof(double) * NDIM * NDIM * NDIM);

    for (int mu = 0; mu < NDIM; mu++) {
        for (int kap = 0; kap < NDIM; kap++) {
            Xh[kap] = X[kap];
            Xl[kap] = X[kap];
        }

        Xh[mu] += DELTA;
        Xl[mu] -= DELTA;
		if (strcmp(metric, "Kerr") == 0) {
			calculate_metric(Xh, gh, g_inv);
			calculate_metric(Xl, gl, g_inv);
		}
		else if (strcmp(metric, "3-Sphere") == 0) {
			calculate_sphere_metric(Xh, gh, g_inv);
			calculate_sphere_metric(Xl, gl, g_inv);
		}
		else if (strcmp(metric, "Minkowski") == 0) {
			calculate_minkowski_metric(Xh, gh, g_inv);
			calculate_minkowski_metric(Xl, gl, g_inv);
		}
		else {
			printf("Unknown metric\n");
			return;
		}

        for (int lam = 0; lam < NDIM; lam++) {
            for (int nu = 0; nu < NDIM; nu++) {
                gamma[lam][nu][mu] = (gh[lam][nu] - gl[lam][nu]) / (2 * DELTA);
            }
        }
    }

    for (int lam = 0; lam < NDIM; lam++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int mu = 0; mu < NDIM; mu++) {
                tmp[lam][nu][mu] = 0.5 * (gamma[nu][lam][mu] + 
                                          gamma[mu][lam][nu] - 
                                          gamma[mu][nu][lam]);
            }
        }
    }

    for (int lam = 0; lam < NDIM; lam++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int mu = 0; mu < NDIM; mu++) {
                gamma[lam][nu][mu] = 0.0;
                for (int kap = 0; kap < NDIM; kap++) {
                    gamma[lam][nu][mu] += g_inv[lam][kap] * tmp[kap][nu][mu];
                }
            }
        }
    }

    
	double gamma_rtt_expected = (1.0 * (X[1] - 2.0 * 1.0)) / (X[1] * X[1] * X[1]);
    double gamma_rrr_expected = -(1.0 / (X[1] * (X[1] - 2.0)));
    double gamma_rthth_expected = -(X[1] - 2.0);
    double gamma_rphph_expected = -(X[1] - 2.0) * sin(X[2]) * sin(X[2]);
    double gamma_trt_expected = (1.0 / (X[1] * (X[1] - 2.0)));
    double gamma_thrth_expected = (1.0 / X[1]);
    double gamma_thphph_expected = -sin(X[2]) * cos(X[2]);
    double gamma_phrph_expected = (1.0 / X[1]);
    double gamma_phthph_expected = (cos(X[2]) / sin(X[2]));

    // Verifications
    if (fabs(gamma[1][0][0] - gamma_rtt_expected) < TOLERANCE) {
        printf("Gamma^r_{tt} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][0][0], gamma_rtt_expected);
    } else {
        printf("Gamma^r_{tt} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][0][0], gamma_rtt_expected);
    }

    if (fabs(gamma[1][1][1] - gamma_rrr_expected) < TOLERANCE) {
        printf("Gamma^r_{rr} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][1][1], gamma_rrr_expected);
    } else {
        printf("Gamma^r_{rr} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][1][1], gamma_rrr_expected);
    }

    if (fabs(gamma[1][2][2] - gamma_rthth_expected) < TOLERANCE) {
        printf("Gamma^r_{\u03b8\u03b8} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][2][2], gamma_rthth_expected);
    } else {
        printf("Gamma^r_{\u03b8\u03b8} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][2][2], gamma_rthth_expected);
    }

    if (fabs(gamma[1][3][3] - gamma_rphph_expected) < TOLERANCE) {
        printf("Gamma^r_{\u03c6\u03c6} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][3][3], gamma_rphph_expected);
    } else {
        printf("Gamma^r_{\u03c6\u03c6} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[1][3][3], gamma_rphph_expected);
    }

    if (fabs(gamma[0][1][0] - gamma_trt_expected) < TOLERANCE) {
        printf("Gamma^t_{rt} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[0][1][0], gamma_trt_expected);
    } else {
        printf("Gamma^t_{rt} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[0][1][0], gamma_trt_expected);
    }

    if (fabs(gamma[2][1][2] - gamma_thrth_expected) < TOLERANCE) {
        printf("Gamma^\u03b8_{r\u03b8} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[2][1][2], gamma_thrth_expected);
    } else {
        printf("Gamma^\u03b8_{r\u03b8} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[2][1][2], gamma_thrth_expected);
    }

    if (fabs(gamma[2][3][3] - gamma_thphph_expected) < TOLERANCE) {
        printf("Gamma^\u03b8_{\u03c6\u03c6} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[2][3][3], gamma_thphph_expected);
    } else {
        printf("Gamma^\u03b8_{\u03c6\u03c6} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[2][3][3], gamma_thphph_expected);
    }

    if (fabs(gamma[3][1][3] - gamma_phrph_expected) < TOLERANCE) {
        printf("Gamma^\u03c6_{r\u03c6} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[3][1][3], gamma_phrph_expected);
    } else {
        printf("Gamma^\u03c6_{r\u03c6} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[3][1][3], gamma_phrph_expected);
    }

    if (fabs(gamma[3][2][3] - gamma_phthph_expected) < TOLERANCE) {
        printf("Gamma^\u03c6_{\u03b8\u03c6} is correct\n");
        printf("Computed: %e, Expected: %e\n", gamma[3][2][3], gamma_phthph_expected);
    } else {
        printf("Gamma^\u03c6_{\u03b8\u03c6} is incorrect\n");
        printf("Computed: %e, Expected: %e\n", gamma[3][2][3], gamma_phthph_expected);
    }}

