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

    printf("Vérification terminée.\n");
}


void calculate_metric(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    double r = x[1] + 3;
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



void calculate_christoffel(double X[NDIM], double h, double gamma[NDIM][NDIM][NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
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

        calculate_metric(Xh, gh, g_inv);
        calculate_metric(Xl, gl, g_inv);

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
                    // gamma[lam][nu][mu] = g_inv[lam][kap] * (symétrisation)
                    gamma[lam][nu][mu] += g_inv[lam][kap] * tmp[kap][nu][mu];
                }
            }
        }
    }
	if (gamma[0][0][1] > 0.0) {
		gamma[1][1][1] = gamma[0][0][1];
	}
}

