#include "./twoP.h"

//two parameter and three parameter equations for cryobiology
//includes a number of conversions for solutions

std::vector <double> cell_exterior_densities;
std::vector <double> cell_interior_densities;
std::vector <std::vector <double> > virial_coefficients;
double cell_volume=0.0;
std::vector <double> dVw{0.0,0.0};//current and previous dVw_dt 
void two_parameter_water_osmolality(const double *Lp, const double *R, const double *T, double *surface_area, double *exterior_osmolality, double *internal_osmolality, double *cell_volume, double *prev_volume )
{
    //function=dVw/dt=-Lp*A*R*T*(osmo^ext-osmo^int)
    //dVw[0]=previous dVw_dt
    dVw[1]=-1*(*Lp)*(*Lp)*(*R)*(*T);
    return;
}