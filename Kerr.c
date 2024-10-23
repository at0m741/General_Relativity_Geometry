#include "include.h"


void Kerr_surfaces(double x[NDIM], double a)
{
	double inner_ergo = 1.0 - sqrt(2.0 - a*a * cos(x[2])*cos(x[2]));
	double outer_ergo = 1.0 + sqrt(2.0 - a*a * cos(x[2])*cos(x[2]));

	printf("inner_ergo = %g\n", inner_ergo);
	printf("outer_ergo = %g\n", outer_ergo);

	double x_inner = inner_ergo * sin(x[2]) * cos(x[3]);
	double y_inner = inner_ergo * sin(x[2]) * sin(x[3]);
	double z_inner = inner_ergo * cos(x[2]);

	double x_outer = outer_ergo * sin(x[2]) * cos(x[3]);
	double y_outer = outer_ergo * sin(x[2]) * sin(x[3]);
	double z_outer = outer_ergo * cos(x[2]);

	printf("inner: %5.3g %5.3g %5.3g\n", x_inner, y_inner, z_inner);
	printf("outer: %5.3g %5.3g %5.3g\n", x_outer, y_outer, z_outer);
}


