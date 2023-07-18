
//This file and the associated header contain a custom class to enfold the standard cell class in physicell
//Author Joseph Abrams
// Note this was built for PhysiCell v1.10 but it should work with future versions
// TODO Someday: this class could be improved with the addition of some vector <int> for passing voxel indicies around


#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <memory>
#include <stdexcept>

#include "follicle_utilities.h"
#include "springs.h"
using namespace BioFVM;
using namespace PhysiCell;
//double V_wbar=18.015/1000; //ml/mol->l/mol partial molar water volume
//double V_nbar=27.1/1000; //ml/mol partial molar salt volume
//double V_egbar=55.766/1000; //ml/mol partial molar eg volume not currently used
//double V_glybar=0.0730903; //ml/mol	

namespace Springs{

std::vector <Spring_Cell*> all_spring_cells;
std::vector <Spring_Cell*> spring_cell_by_pCell_index;

/* eventually make its own class objects*/
// Cryo_State_Variables::Cryo_State_Variables(int solute_num)
// {
//     temperature=25;
//     gas_constant=8.314;
//     surface_area=0;
//     next_water_volume=0;//for ABM 2nd order
//     water_volume=0;
//     dVw=0;//cell water flux
//     previous_dVw=0;
//     solute_moles.resize(solute_num,0.0);
//     dN.resize(solute_num,0.0);//cell mole flux of solutes

//     exterior_osmolality=0.0;//total exterior osmolality salt+CPA (mole/kg)
//     interior_osmolality=0.0;// total internal osmolality salt+CPA
//     internal_molarity.resize(solute_num,0.0);
//     external_molarity.resize(solute_num,0.0);
//     Ps.resize(solute_num,0.0);
//     Lp=0.0;
//   return;
// }

void Spring_Cell::initialize()
{
  //initializing stuff after constructor is called
  return;
}
/* IMPORTANT NOTE: solute vectors are all of the form {0,HM,EG,GLY,PBS}
 * temp is currently fixed
 * cryo state variables are part of spring cell, a cell type encapsulating class
 * virial equation is used for ternary solutions and binary solutions, quaternary
 * solution osmolality is from component addition and is thus likely and underestimate
 * since it is missing some cross terms
 */

const double Spring_Cell::temperature=293.15;//parameters.doubles("temperature");//293.15;
const double Spring_Cell::gas_constant=0.08205;//parameters.doubles("gas_constant");//0.08205;
const std::vector <double> Spring_Cell::solute_specific_volume={0.01661,0.0557414,0.0730903,0.01661};//{HM,EG,GLY,PBS} in um^3/femptomole (l/mole)
//eventually pull values from xml for solution thermodynamic constants
// const std::vector <double> Spring_Cell::solute_specific_volume={0.0,0.01661,0.0557414,0.0730903,0.01661};//{0,HM,EG,GLY,PBS} in um^3/femptomole (l/mole)
const std::vector <double> Spring_Cell::kdiss={1.678,1.0,1.0,1.678};//{HM,EG,GLY,PBS} virial ionic constant
const std::vector <double> Spring_Cell::virial_B={0.044,0.037,0.023,0.044};
const std::vector <double> Spring_Cell::virial_C={0.0,0.0,-0.001,0.0};
const std::vector <double> Spring_Cell::molar_mass={58.44,62.07,92.09,58.44};// note HM and PBS can have different molar masses depending on fomulation (especially with disodium phosphate) but we treat them as NaCl
Spring_Cell::Spring_Cell(Cell* pCell)//(Cell* pCell,Phenotype& phenotype)
{
  simulation_selected=parameters.ints("selected_simulation");
  is_basement_connected=false;
  basement_length=0.0;
  m_my_pCell=pCell;
// state parameters for Cryo
  surface_area=pCell->custom_data["surface_area"];
  next_water_volume=0.0;//for ABM 2nd order
  solid_volume=pCell->custom_data["initial_volume"]*pCell->custom_data["Vb_fraction"];
  solute_volume=1449.82;
  water_volume=pCell->custom_data["initial_volume"]-solid_volume-solute_volume;
  use_virial=true;//ternary use, quatenary dont
  dVw=0.0;//cell water flux
  previous_dVw=0.0;
  solute_moles.resize(microenvironment.number_of_densities(),0.0);
  next_solute_moles.resize(microenvironment.number_of_densities(),0.0);
  dN.resize(microenvironment.number_of_densities(),0.0); //cell mole flux of solutes
  previous_dN.resize(microenvironment.number_of_densities(),0.0);//previous mole flux of solutes

  exterior_osmolality=0.255;//total exterior osmolality salt+CPA (mole/kg)
  interior_osmolality=0.255;// total internal osmolality salt+CPA
  interior_molarity.resize(microenvironment.number_of_densities(),0.0);//1.68 is the kdiss for NaCl from virial
  exterior_molarity=interior_molarity;
  interior_component_molality=interior_molarity;
  exterior_component_molality=exterior_molarity;
  Ps.resize(microenvironment.number_of_densities(),0.0);
  Lp= parameters.doubles("oocyte_Lp_EG");
  solute_uptake=solute_moles;
  solute_volume=0.0;
  index=0;
  water_uptake=0.0;
 return;
}
//     Cryo_State_Variables(int solute_num);
//     void get_external_molarity(int component, Cell* pCell);
//     double compute_EG_molality(double EG_molarity);
//     double compute_NaCl_molality(double EG_molarity);
//     double compute_Gly_molality(double EG_molarity);
//  virial_solution_osmolality
//  compute_component_internal_molarity //moles/water_volume
//  component_virial_osmolality
//  dN 
//  dVw
//  Two_P_forward_step

Spring_Cell::~Spring_Cell()
{
  delete_spring_cell(this->index);
  return;
}
void Spring_Cell::add_spring(Spring_Cell* other_pSCell, double spring_length)
{
  #pragma omp critical//could speed up with private merge
  {
    std::cout<<"new spring made"<<std::endl;
    Spring* new_spring;
    new_spring=new Spring(other_pSCell->m_my_pCell,spring_length);
    this->m_springs.push_back(new_spring);
  }
  return;
}

void Spring_Cell::remove_spring(Spring_Cell* other_pSCell)
{
  //copied from removing default spring attachments
  #pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < this->m_springs.size() ) 
		{
			
			if( this->m_springs[i]->m_pNeighbor == other_pSCell->m_my_pCell )
			{
        int n = this->m_springs.size();
        Spring* spring_ptr= this->m_springs[n-i];
				// copy last entry to current position 
				this->m_springs[i] = this->m_springs[n-1]; 
        
				// shrink by one 
				this->m_springs.pop_back(); 
        delete spring_ptr;
				found = true; 
			}
			i++; 
		}
	}
	return; 
}
void Spring_Cell::set_2p_initial_conditions()
{
  output_cell_and_voxels(get_exterior_voxels(this->m_my_pCell), this->m_my_pCell);
  // 1 mol/L = 1 fmol/um^3 
  // double HM_density=1.0045; //NaCl density g/cm^3 extrapolated from CRC 
  // double solvent_volume= this->m_my_pCell->phenotype.volume.total-this->solid_volume;
  //set initial stat vectors

  this->interior_molarity[0]=default_microenvironment_options.initial_condition_vector[0];
  // set up permeability vector 
  this->interior_molarity[1]=default_microenvironment_options.initial_condition_vector[1];
  this->Ps[0]=0.0;//hm
  this->solute_moles[0]=(this->interior_molarity[0])*(this->m_my_pCell->phenotype.volume.total-this->solid_volume);//HM osmolality/kdiss*kg/L solvent*density= femptomoles HM
  std::cout<<"Initial HM moles"<<this->solute_moles[0]<<"\n";
  if(this->simulation_selected==1)
  {//EG
    // std::cout<<"FLAG!!: "<<"\n"; 
    if(this->m_my_pCell->type==1)
    { 
      this->Ps[1]=parameters.doubles("oocyte_Ps_EG");
      this->Lp=parameters.doubles("oocyte_Lp_EG");
    }
    else
    {
      this->Ps[1]=parameters.doubles("gran_Ps_EG");
      this->Lp=parameters.doubles("gran_Lp_EG");
    }
    this->solute_moles[1]=0;
      this->interior_component_molality[0]=molarity_to_molality(this->interior_molarity[0], "NaCl");
      this->interior_component_molality[1]=molarity_to_molality(this->interior_molarity[1], "EG");
      double vir=ternary_virial(this->interior_component_molality[0],this->interior_component_molality[1],"NaCl","EG"); 
      std::cout<<"Initial OSMOLALITY: "<< vir<<"\n";   
      this->interior_osmolality=ternary_virial(this->interior_component_molality[0],this->interior_component_molality[1],"NaCl","EG");
  }
  else if(this->simulation_selected==2)
  {//GLY
    if(this->m_my_pCell->type==1)
    { 
      this->Ps[1]=parameters.doubles("oocyte_Ps_GLY");
      this->Lp=parameters.doubles("oocyte_Lp_GLY");
    }
    else
    {
      this->Ps[1]=parameters.doubles("gran_Ps_GLY");
      this->Lp=parameters.doubles("gran_Lp_GLY");
    }
    this->solute_moles[1]=0;
      this->interior_component_molality[0]=molarity_to_molality(this->interior_molarity[0], "NaCl");
      this->interior_component_molality[1]=molarity_to_molality(this->interior_molarity[1], "GLY");
      double vir=ternary_virial(this->interior_component_molality[0],this->interior_component_molality[1],"NaCl","GLY"); 
      std::cout<<"Initial OSMOLALITY: "<< vir<<"\n";   
      this->interior_osmolality=ternary_virial(this->interior_component_molality[0],this->interior_component_molality[1],"NaCl","GLY");
  }
  else if(this->simulation_selected==3)
  {//PBS
    if(this->m_my_pCell->type==1)
    { 
      this->Ps[1]=parameters.doubles("oocyte_Ps_NaCl");
      this->Lp=parameters.doubles("oocyte_Lp_NaCl");
    }
    else
    {
      this->Ps[1]=parameters.doubles("gran_Ps_NaCl");
      this->Lp=parameters.doubles("gran_Lp_NaCl");
    }
    this->solute_moles[1]=0;
      this->interior_component_molality[0]=molarity_to_molality(this->interior_molarity[0], "NaCl");
      this->interior_component_molality[1]=molarity_to_molality(this->interior_molarity[1], "NaCl");
      double vir=ternary_virial(this->interior_component_molality[0],this->interior_component_molality[1],"NaCl","NaCl"); 
      std::cout<<"Initial OSMOLALITY: "<< vir<<"\n";   
      this->interior_osmolality=ternary_virial(this->interior_component_molality[0],this->interior_component_molality[1],"NaCl","NaCl");
  }
  else if(this->simulation_selected==4)
  {//EG & GLY
    if(this->m_my_pCell->type==1)
    { 
      this->Ps[1]=parameters.doubles("oocyte_Ps_EG");
      this->Lp=parameters.doubles("oocyte_Lp_EG_and_GLY");
      this->Ps[2]=parameters.doubles("oocyte_Ps_GLY");
    }
    else
    {
      this->Ps[1]=parameters.doubles("gran_Ps_EG");
      this->Lp=parameters.doubles("gran_Lp_EG_and_GLY");
      this->Ps[2]=parameters.doubles("gran_Ps_GLY");
    }
    this->solute_moles[1]=0;
    this->solute_moles[2]=0;
      this->interior_component_molality[0]=molarity_to_molality(this->interior_molarity[0], "NaCl");
      this->interior_component_molality[1]=molarity_to_molality(this->interior_molarity[1], "EG");
      this->interior_component_molality[2]=molarity_to_molality(this->interior_molarity[2], "GLY");
      double vir=ternary_virial(this->interior_component_molality[1],this->interior_component_molality[2],"EG","GLY")+binary_virial(this->interior_component_molality[0], "NaCl"); 
      std::cout<<"Initial OSMOLALITY: "<< vir<<"\n";   
      this->interior_osmolality=ternary_virial(this->interior_component_molality[1],this->interior_component_molality[2],"EG","GLY")+binary_virial(this->interior_component_molality[0], "NaCl");
  }
  return;
}
Spring::Spring(Cell* pNeighbor, double length)
{
  
  m_pNeighbor=pNeighbor;
  m_length=length;
  std::cout<< "spring constructed"<<std::endl;
  return;
}

Spring_Cell* create_spring_cell(Cell* pCell)//( Cell* pCell,Phenotype& phenotype )
{
	//enfold already existing Cell into super class spring cell
	//auto pNewSpring_Cell=std::make_shared<Spring_Cell>;
  Spring_Cell* pNew; 
  #pragma omp critical //could be made multithread probably if needed but shouldn't be an issue
  { 
    
	  pNew = new Spring_Cell(pCell);//(pCell,phenotype);		
	  all_spring_cells.push_back( pNew ); 
	  pNew->index=all_spring_cells.size()-1;
    spring_cell_by_pCell_index[pNew->m_my_pCell->index]=pNew;
  	
  }
  return pNew; 
}
void delete_spring_cell( int index)
{
  #pragma omp critical
  {
	  Spring_Cell* pDeleteMe = all_spring_cells[index]; 
	  // move last item to index location  
	  all_spring_cells[ all_spring_cells.size()-1 ]->index=index;
	  all_spring_cells[index] = all_spring_cells[ all_spring_cells.size()-1 ];
	  // shrink the vector
	  all_spring_cells.pop_back();	
	  delete pDeleteMe; 
	  
  }
  return; 
}


void Spring_Cell::two_p_update_volume()
{
  // std::cout<<"exterior_osmolality: "<<this->exterior_osmolality<<"\n";
  // std::cout<<"interior_osmolality: "<<this->interior_osmolality<<"\n";
    double total_solute_volume=0.0;
    #pragma omp reduction(+:total_solute_volume)
    for (size_t i = 0; i < this->solute_moles.size(); i++)
    {
        total_solute_volume+=(this->solute_specific_volume)[i]*(this->solute_moles)[i];
    }
    
  #pragma omp critical
  {
    // std::cout<<"Solute volume: "<<total_solute_volume<<"\n";
    // std::cout<<"water_volume: "<<this->water_volume<<"\n";
    this->solute_volume=total_solute_volume;
    // std::cout<<"solid volume: "<<this->solid_volume<<"\n";
    // std::cout<<"Total volume: "<<this->m_my_pCell->phenotype.volume.total<<"\n";
    this->m_my_pCell->phenotype.volume.total=this->water_volume+this->solid_volume+total_solute_volume;  
    // std::cout<<"Total volume: "<<this->m_my_pCell->phenotype.volume.total<<"\n";
  }
  return;
}

void Spring_Cell::dVw_Osmolality()
{ //-LpART(m_e-m_i)
  this->previous_dVw=this->dVw;
  this->dVw=-1*this->Lp*(this->surface_area)*this->temperature*this->gas_constant*(this->exterior_osmolality-this->interior_osmolality);//total exterior osmolality salt+CPA (mole/kg)
  // std::cout<<"exterior_osmolality: "<<this->exterior_osmolality<<"\n";
  // std::cout<<"interior_osmolality: "<<this->interior_osmolality<<"\n";
  // std::cout<<"difference: "<< this->exterior_osmolality-this->interior_osmolality<<"\n";
  // std::cout<<"dVw: "<<this->dVw<<"\n";
  // std::cout<<"Lp: "<< this->Lp<<"\n";
  // std::cout<<"temperature: "<< this->temperature<<"\n";
  // std::cout<<"gas_constant: "<< this->gas_constant<<"\n";
  return;
}
void Spring_Cell::dN_molarity()
{
    //multisolute
    //function=dN/dt=Ps*A(mol^ext-mol^int) for the ith solute, Ps=0 for non-permeating
    //internal molarity depends on Vw
    for (size_t i = 0; i < (dN).size(); i++)
    {
      this->previous_dN=this->dN;
      this->dN[i]=this->Ps[i]*(this->surface_area)*((this->exterior_molarity[i])-(this->interior_molarity[i]));
    }
  // std::cout<<"Ps: "<< this->Ps<<"\n";
  // std::cout<<"dN: "<<this->dN<<"\n";
  // std::cout<<"exterior_molarity: "<<this->exterior_molarity[0]<<"\n";
  // std::cout<<"interior_molarity: "<<this->interior_molarity[0]<<"\n";
  // std::cout<<"exterior_molarity 1: "<<this->exterior_molarity[1]<<"\n";
  // std::cout<<"interior_molarity 1: "<<this->interior_molarity[1]<<"\n";
    return;
}
void Spring_Cell::two_p_forward_step()
{
  // std::cout<<"TEST RUN FLAG!"<<"\n";
  if (PhysiCell_globals.current_time<=0.01) {
  // std::cout<<"Second FLAG!"<<"\n";
    this->next_water_volume=this->water_volume+(this->dVw*0.01);//forward_euler
    for (size_t i = 0; i < (dN).size(); i++)
    {
      this->next_solute_moles[i]=this->solute_moles[i]+(0.01*this->dN[i]);//forward_euler
      // std::cout<<"current moles: "<< this->solute_moles[1]<<"\n";
      // std::cout<<"osmotic volume: "<< this->m_my_pCell->phenotype.volume.total-this->solid_volume<<"\n";
      // std::cout<<"interior concentration: "<< (this->solute_moles[0]+this->solute_moles[1])/(this->m_my_pCell->phenotype.volume.total-this->solid_volume)<<"\n";
      // std::cout<<"sum of molarity:  "<<this->interior_molarity[0]+this->interior_molarity[1]<<"\n"; 
      // std::cout<<"interior osmolality: "<< this->interior_osmolality<<"\n"; 
      this->solute_uptake[i]=this->next_solute_moles[i]-this->solute_moles[i];//forward_euler
      this->solute_moles[i]=this->next_solute_moles[i];//update for next step
      // std::cout<<"uptake: "<< this->solute_uptake[i]<<"\n";
    }
  }
  else {
    
  
    this->next_water_volume=this->water_volume+(this->dVw*0.01);//forward_euler
  // std::cout<<"Third FLAG!"<<"\n";
  // this->next_water_volume=this->water_volume+(0.01/2)*((3*this->dVw)-(this->previous_dVw));
    for (size_t i = 0; i < (dN).size(); i++)
    {
    
      this->next_solute_moles[i]=this->solute_moles[i]+(0.01*this->dN[i]);//forward_euler
      // this->next_solute_moles[i]=this->solute_moles[i]+(0.01/2)*((3*this->dN[i])-(this->previous_dN[i]));
      // std::cout<<"current moles: "<< this->solute_moles[1]<<"\n";
      // std::cout<<"osmotic volume: "<< this->m_my_pCell->phenotype.volume.total-this->solid_volume<<"\n";
      // std::cout<<"interior concentration: "<< (this->solute_moles[0]+this->solute_moles[1])/(this->m_my_pCell->phenotype.volume.total-this->solid_volume)<<"\n";
      // std::cout<<"sum of molarity:  "<<this->interior_molarity[0]+this->interior_molarity[1]<<"\n"; 
      // std::cout<<"interior osmolality: "<< this->interior_osmolality<<"\n"; 
      // std::cout<<"exterior osmolality"<< this->exterior_osmolality<<"\n"; 
      //std::cout<<"current moles: "<< this->solute_moles[1]<<"\n";
      // std::cout<<"next moles: "<< this->next_solute_moles[1]<<"\n";
      this->solute_uptake[i]=this->next_solute_moles[i]-this->solute_moles[i];//forward_euler
      this->solute_moles[i]=this->next_solute_moles[i];//update for next step
      // std::cout<<"MOLES: "<< this->solute_moles[i]<<"\n";
    }
  }
  this->water_uptake=this->next_water_volume-this->water_volume;
  this->water_volume=this->next_water_volume;//update for next step
  return;
}
/*
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
*/

};
//std::vector <Spring_Cell*> all_spring_cells{};
