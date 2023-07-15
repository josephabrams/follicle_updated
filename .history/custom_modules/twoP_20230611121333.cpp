#include "./twoP.h"

//two parameter and three parameter equations for cryobiology
//includes a number of conversions for solutions

std::vector <double> cell_exterior_densities;
std::vector <double> cell_interior_densities;
std::vector <std::vector <double> > virial_coefficients;
double cell_volume=0.0;
void two_parameter_water_osmolality(const double *Lp, const double *R, const double *T, double *surface_area, double *exterior_osmolality, double *internal_osmolality )
{
    //dVw/dt=-Lp*R*T
    return;
}