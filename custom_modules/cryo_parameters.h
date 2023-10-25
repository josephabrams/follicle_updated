/* Author: Joseph E.S. Abrams
 * Please cite if you use this code
 * This is a library of functions and paramaters for simualting cells and tissue undergoing cryopreservation*/

#ifndef __CRYO_PARAMETERS_H__
#define __CRYO_PARAMETERS_H__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <vector>
using namespace BioFVM;
using namespace PhysiCell;

namespace cryo_parameters{

class Cryo_Parameters
{
 private:
 public:
  int simulation_selected;// if using XML value to specify a condition
  //constants available to every spring cell, future version could have changing temp
  double temperature;
  const double gas_constant=0.08205;
  /////paramaters for unit conversions using the virial equation
    ///
  const std::vector <double> solute_specific_volume={0.01661,0.0557414,0.0730903,0.01661};//{HM,EG,GLY,PBS} in um^3/femptomole (l/mole)
  const std::vector <double> kdiss={1.678,1.0,1.0,1.678};//{HM,EG,GLY,PBS} virial ionic constant
  const std::vector <double> virial_B={0.044,0.037,0.023,0.044};
  const std::vector <double> virial_C={0.0,0.0,-0.001,0.0};
  const std::vector <double> molar_mass={58.44,62.07,92.09,58.44};// note HM and PBS can have different molar masses depending on fomulation (especially with disodium phosphate) but we treat them as NaCl
  /////////////////////////////////////////////////////////////  
  double solute_volume; //total volume taken up by the free moles of solute in the cell
  double surface_area; 
  double next_water_volume;//for ABM 2nd order stores the next value of water volume
  double water_volume;//current water volume 
  double solid_volume;//non-osmotically active volume of cell
  double dVw;//cell water flux
  double previous_dVw;
  double previous_radius;
  std::vector<double> solute_moles;
  std::vector<double> next_solute_moles;
  std::vector<double> uptake;//molar uptake/secretion of solutes
  std::vector<double> dN;//cell mole flux of solutes
  std::vector<double> previous_dN;//previous mole flux of solutes
  double exterior_osmolality;//total exterior osmolality salt+CPA (mole/kg)
  double interior_osmolality;// total internal osmolality salt+CPA
  std::vector<double> interior_molarity;
  std::vector<double> exterior_molarity;
  std::vector<double> interior_component_molality;
  std::vector<double> exterior_component_molality;
  std::vector<double> Ps;
  double Lp;
  double toxicity;
  std::vector <double> solute_uptake;
  double water_uptake;
  std::vector <int> uptake_voxels;
  bool use_virial;
	Cryo_Parameters(Cell* pCell); 
}; 





};//namespace cryo_parameters











#endif //__CRYO_PARAMETERS_H__

