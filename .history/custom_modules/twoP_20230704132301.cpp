#include "./twoP.h"

//two parameter and three parameter equations for cryobiology
//includes a number of conversions for solutions

std::vector <double> cell_exterior_densities;
std::vector <double> cell_interior_densities;
//std::vector <std::vector <double> > virial_coefficients{{0,0,0},{1,},{1},{1.678,0.044,0},{}};//same index as solutes {kdiss,B,C}: empty, EG, GLY, PBS, HM, SUC
double cell_volume=0.0;
std::vector<double> water_volume{0.0,0.0};
std::vector<double> moles{0.0,0.0};
std::vector <double> dVw{0.0,0.0};//current and previous dVw_dt 
std::vector <double> dN{0.0,0.0};
void dVw_osmolality(double *Lp, double *R, double *T, double *surface_area, double *exterior_osmolality, double *internal_osmolality, double *dVw )
{
    //function=dVw/dt=-Lp*A*R*T*(osmo^ext-osmo^int)
    (*dVw)=-1*(*Lp)*(*surface_area)*(*R)*(*T)*((*exterior_osmolality)-(*internal_osmolality));
    return;
}
void get_molarities(Cell* pCell,std::vector<double> *exterior_molarity, std::vector<double> *internal_molarity, std::vector<double> *moles, double *water_volume)
{
    for (size_t i = 0; i < pCell->get_microenvironment()->nearest_gradient_vector(pCell->position).size(); i++) //for all solutes so only include the ones needed in a simulation
    {
        (*exterior_molarity)[i]=concentration_at_boundary(pCell, i);//external function is located in follicle_utilities, BioFVM diffusion has external values in microenvironment
        (*internal_molarity)[i]=(*moles)[i]/(*water_volume);//moles and water volume updated from 2p model and initial conditions
    }
    return;
}
/* The following equation computes the osmolality of a ternary solution with m ininitial molalilities using B and C coefficients,with kdiss for ions*/
double get_viral_osmolality(std::vector<double> m, std::vector <double> B, std::vector <double> C, std::vector <double> kdiss)
{
    std::vector temp_m=m[i];//use copy for parallel

    for (size_t i = 0; i < kdiss.size(); i++)
    {
        m[i]=temp_m[i]*kdiss[i]; //kdiss is 1 for non-osmolyte
    }

    //compute virial equation
    double osmolality=0.0;
    #pragma omp reduction(+:osmolality)
    for (size_t i = 0; i < m.size(); i++)
    {
        osmolality+=m[i];
        for (size_t j = 0; j < B.size(); j++)
        {
            osmolality+=((B[i]+B[j])/2)*m[i]*m[j];
            for (size_t k = 0; k < C.size(); k++)
            {
                osmolality+=((C[i]*C[j]*C[k])^(1/3))*m[i]*m[j]*m[k];
            }
            
        }
    }
    
    return osmolality;
}
void osmolality_from_component_moles(double *osmolality, std::vector<double> *moles,std::vector<double> *molar_mass)//component moles->component molality->solution molality(virial equation or dilute addition)
{
  double um_to_m=1.0e-18;
  double water_density=.997; //g/m^3
  double water_mass=water_volume*um_to_m*water_density; // um^3*(m^3/um^3*(g/m^3)1 cubic micrometer = 1e-18 cubic meter
  double solution_mass=water_mass;//grams solution
  #pragma omp reduction(+:solution_mass)
  
  for(size_t i=0; i<(*moles).size(); i++)
  {
    solution_mass+=(*moles)[i]*(*molar_mass)[i];
  }
  
  return;
}
void two_p_update_volume(Cell* pCell, double *water_volume, double *solid_volume, std::vector<double> *solute_specific_volume)
{
    double total_solute_volume=0.0;
    #pragma omp reduction(+:total_solute_volume)
    for (size_t i = 0; i < pCell->get_microenvironment()->nearest_gradient_vector(pCell->position).size(); i++)
    {
        total_solute_volume+=(*solute_specific_volume)[i]*(*moles)[i];
    }
    

    pCell->phenotype.volume.total=(*water_volume)+(*solid_volume)+total_solute_volume;   
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
void Adams_Bashforth_ODE_2nd_Order(double *Y_next, double *Y_current, double *df_dts, double *previous_df_dts, double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))

    (*Y_next) = (*Y_current) + (step_size / 2) * (3 * ((*df_dts) - (*previous_df_dts)));
  return;
}
// void Adams_Bashforth_ODE_2nd_Order(std::vector<double> *Y, std::vector<double> *df_dt, double step_size) {
//   // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
//     (*Y)[1] = (*Y)[0] + (step_size / 2) * (3 * ((*df_dt)[1] - (*df_dt)[0]));
//   return;
// }
// void Adams_Bashforth_ODE_2nd_Order(std::vector< std::vector<double>> *Y, std::vector< std::vector<double>> *df_dt, double step_size) {
//   // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
//   //Y_{n+1} is Y[1], F(x_{n},t_{n}) is df_dt[1]
//   for (unsigned int i = 0; i < (*Y)[0].size(); i++) {
//     (*Y)[1][i] = (*Y)[0][i] + (step_size / 2) * (3 * ((*df_dt)[1][i] - (*df_dt)[0][i]));
//   }
//   return;
// }
// two_p model updates vectors attached to cell object and passes out the concentration
// 
void two_p(double* Lp, double* R, double* T, double* surface_area, double* exterior_osmolality, double* internal_osmolality, double* dVw, double* previous_dVw, std::vector<double> *exterior_molarity, std::vector<double> *Ps, std::vector <double> *internal_molarity, std::vector<double> *dN, std::vector<double> *previous_dN, double dt, double *water_volume, double *previous_water_volume, std::vector <double> *moles, std::vector <double> *previous_moles, std::vector <double> *uptake_concentrations )
{
    dVw_osmolality((Lp), (R), (T), (surface_area), (exterior_osmolality), (internal_osmolality), (dVw));//curent set dVw[1]
    dN_molarity(Ps, surface_area, exterior_molarity, internal_molarity, dN);//already loops through all solutes
    if(PhysiCell_globals.current_time<dt)//forward euler for the first time step
    {
        (*water_volume)=(*previous_water_volume)+((*dVw)*dt);
        for (size_t i = 0; i < (*moles).size(); i++)
        {
            (*moles)[i]=(*moles)[i]+((*dN)[i]*dt);//update moles
        }
    }
    else //adams bashhforth second order
    {    
        Adams_Bashforth_ODE_2nd_Order(water_volume,previous_water_volume,dVw,previous_dVw, dt);
        Adams_Bashforth_ODE_2nd_Order(moles,previous_moles,dN,previous_dN, dt);//update moles uses vector version with built in loop
    }
    for (size_t i = 0; i < (*moles).size(); i++)
    {
        (*uptake_concentrations)[i]=((*previous_moles)[i]/(*previous_water_volume))-((*moles)[i]/(*water_volume));
        (*previous_dN)[i]=(*dN)[i];// set previous for next step (*dN)[0]=(*dN)[1];
    }
    (*previous_dVw)=(*dVw);// set previous for next step (*dVw)[0]=(*dVw)[1];
    (*previous_water_volume)=(*water_volume);//current to previous


    return;
}
