#include "./twoP.h"

//two parameter and three parameter equations for cryobiology
//includes a number of conversions for solutions

std::vector <double> cell_exterior_densities;
std::vector <double> cell_interior_densities;
std::vector <std::vector <double> > virial_coefficients;
double cell_volume=0.0;
std::vector<double> water_volume{0.0,0.0};
std::vector<double> moles{0.0,0.0};
std::vector <double> dVw{0.0,0.0};//current and previous dVw_dt 
std::vector <double> dN{0.0,0.0};
void dVw_osmolality(double *Lp, double *R, double *T, double *surface_area, double *exterior_osmolality, double *internal_osmolality, std::vector<double> *dVw)
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
void Adams_Bashforth_ODE_2nd_Order(std::vector<double> *Y_new, std::vector<double> std::vector<double> *prev_df_dts, std::vector<double> *df_dts, double *step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
  for (unsigned int i = 0; i < Y_values.size(); i++) {
    Ynext[i] = Y_values[i] + step_size / 2 * (3 * (df_dts[i] - prev_df_dts[i]));
  }
  return;
}
void two_p(double const *Lp, double const *R, double const *T, double *surface_area, double *exterior_osmolality, double *internal_osmolality, std::vector<double> *dVw, double *exterior_molarity, double *Ps, double *internal_molarity, std::vector<double> *dN, double dt, std::vector <double> *water_volume, std::vector <double> *moles )
{
    (*water_volume)[0]=(*water_volume)[1];//current to previous
    (*moles)[0]=(*moles)[1];//current to previous
    dVw_osmolality((*Lp), (*R), (*T), (*surface_area), (*exterior_osmolality), (*internal_osmolality), (*dVw));//curent
    dN_molarity(*Ps, *surface_area, *exterior_molarity, *internal_molarity, *dN);//current
    if(PhysiCell_globals.current_time<dt)//forward euler
    {
        (*water_volume)[1]=(*water_volume)[0]+((*dVw)[1]*dt);
        (*moles)[1]=(*moles)[0]+((*dN)[1]*dt);
    }
    else //adams bashhforth
    {

    }
    return;
}