#include "./twoP.h"

//two parameter and three parameter equations for cryobiology
//includes a number of conversions for solutions

std::vector <double> cell_exterior_densities;
std::vector <double> cell_interior_densities;
std::vector <std::vector <double> > virial_coefficients;
double cell_volume=0.0;
std::vector<double> water_volume{0.0,0.0};
std::vector <double> dVw{0.0,0.0};//current and previous dVw_dt 
std::vector <double> dN{0.0,0.0};
void dVw_osmolality(double const *Lp, double const *R, double const *T, double *surface_area, double *exterior_osmolality, double *internal_osmolality, std::vector<double> *dVw)
{
    //function=dVw/dt=-Lp*A*R*T*(osmo^ext-osmo^int)
    //dVw[0]=previous dVw_dt
    (*dVw)[1]=-1*(*Lp)*(*surface_area)*(*R)*(*T)*((*exterior_osmolality)-(*internal_osmolality));
    return;
}
void dN_molarity(double *Ps, double *surface_area, double *exterior_molarity, double *internal_molarity, std::vector<double> *dN)
{
    //function=dN/dt=Ps*A(mol^ext-mol^int)
    //dN[0]=previous dN_dt
    (*dN)[1]=(*Ps)*(*surface_area)*((*exterior_molarity)-(*internal_molarity));
    return;
}
void two_p(double const *Lp, double const *R, double const *T, double *surface_area, double *exterior_osmolality, double *internal_osmolality, std::vector<double> *dVw, double *exterior_molarity, double *Ps, double *internal_molarity, std::vector<double> *dN, double dt, std::vector <double> *water_volume )
{
    (*water_volume)[0]=(*water_volume)[1];//current
    (*dN)[0]=(*dN)[1];//current
    if(PhysiCell_globals.current_time<dt)//forward euler
    {
        dVw_osmolality((*Lp), (*R, *T, *surface_area, *exterior_osmolality, *internal_osmolality, *dVw);
        (*water_volume)[1]=(*water_volume)[0]+((*dVw)[1]*dt);
        dN_molarity(*Ps, *surface_area, *exterior_molarity, *internal_molarity, *dN);


    }
    else
    {

    }
    return;
}