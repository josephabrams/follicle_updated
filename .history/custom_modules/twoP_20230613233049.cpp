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
void dN_molarity(std::vector<double> *Ps, double *surface_area, std::vector<double> *exterior_molarity, std::vector<double> *internal_molarity, std::vector<double> *dN)
{
    //multisolute
    //function=dN/dt=Ps*A(mol^ext-mol^int) for the ith solute, Ps=0 for non-permeating
    //internal molarity depends on Vw
    for (size_t i = 0; i < (*dN).size(); i++)
    {
        (*dN)[i]=(*Ps)[i]*(*surface_area)*((*exterior_molarity)[i]-(*internal_molarity)[i]);
    }
    

    return;
}
void Adams_Bashforth_ODE_2nd_Order(std::vector<double> *Y_next, std::vector<double> *Y_current, std::vector<double> *df_dts, std::vector<double> *previous_df_dts, double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
  for (unsigned int i = 0; i < (*Y_current).size(); i++) {
    (*Y_next)[i] = (*Y_current)[i] + (step_size / 2) * (3 * ((*df_dts)[i] - (*previous_df_dts)[i]));
  }
  return;
}
void Adams_Bashforth_ODE_2nd_Order(std::vector<double> *Y, std::vector<double> *df_dt, double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
    (*Y)[1] = (*Y)[0] + (step_size / 2) * (3 * ((*df_dt)[1] - (*df_dt)[0]));
  return;
}
void Adams_Bashforth_ODE_2nd_Order(std::vector< std::vector<double>> *Y, std::vector< std::vector<double>> *df_dt, double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
  //Y_{n+1} is Y[1], F(x_{n},t_{n}) is df_dt[1]
  for (unsigned int i = 0; i < (*Y)[0].size(); i++) {
    (*Y)[1][i] = (*Y)[0][i] + (step_size / 2) * (3 * ((*df_dt)[1][i] - (*df_dt)[0][i]));
  }
  return;
}
// two_p model updates vectors attached to cell object and passes out the concentration
// 
void two_p(double* Lp, double* R, double* T, double* surface_area, double* exterior_osmolality, double* internal_osmolality, double* dVw, double* previous_dVw, std::vector<double> *exterior_molarity, std::vector<double> *Ps, std::vector <double> *internal_molarity, std::vector<double> *dN, std::vector<double> *previous_dN, double dt, double *water_volume, double *previous_water_volume, std::vector <double> *moles, std::vector <double> *previous_moles, std::vector <double> *uptake_concentrations )
{

    dVw_osmolality((Lp), (R), (T), (surface_area), (exterior_osmolality), (internal_osmolality), (dVw));//curent set dVw[1]
    dN_molarity(Ps, surface_area, exterior_molarity, internal_molarity, dN);//current set dVw[1]
    if(PhysiCell_globals.current_time<dt)//forward euler
    {
        (*water_volume)=(*previous_water_volume)+((*dVw)[1]*dt);
        for (size_t i = 0; i < (*moles).size(); i++)
        {
            (*moles)[i]=(*moles)[i]+((*dN)[i]*dt);//
            (*previous_dN)[i]=(*dN)[i];//(*dN)[0]=(*dN)[1];
            (*uptake_concentrations)[i]=(*moles)[i]/(*water_volume);
        }
    }
    else //adams bashhforth second order
    {
        double current_water_volume=(*water_volume);
        Adams_Bashforth_ODE_2nd_Order(water_volume,dVw, 0.01);
        (*previous_dVw)=(*dVw);// set previous for next step (*dVw)[0]=(*dVw)[1];
        (*previous_water_volume)=(*water_volume);//current to previous
        for (size_t i = 0; i < (*moles).size(); i++)
        {
            (internal_molarity)[i]=moles[i]/water_volume
            Adams_Bashforth_ODE_2nd_Order(moles,dN, 0.01);
            (*previous_dN)[i]=(*dN)[i];// set previous for next step (*dN)[0]=(*dN)[1];
        }

    }



    return;
}