#include "./cryo_parameters.h"

using namespace cryo_parameters;
Cryo_Parameters::Cryo_Parameters(Cell* pCell)
{
  // if using XML value to specify a condition
  simulation_selected=parameters.ints("selected_simulation");
  // state parameters for Cryo, should be defined in custom_data for cell
  surface_area=pCell->custom_data["surface_area"];
  next_water_volume=0.0;//for ABM 2nd order
  solid_volume=pCell->custom_data["initial_volume"]*pCell->custom_data["Vb_fraction"];
  solute_volume=pCell->custom_data["initial_solute_volume"];//1449.82;
  water_volume=pCell->custom_data["initial_volume"]-solid_volume-solute_volume;
  use_virial=true;//ternary use, quatenary dont
  dVw=0.0;//cell water flux
  previous_dVw=0.0;
  solute_moles.resize(microenvironment.number_of_densities(),0.0);
  next_solute_moles.resize(microenvironment.number_of_densities(),0.0);
  dN.resize(microenvironment.number_of_densities(),0.0); //cell mole flux of solutes
  previous_dN.resize(microenvironment.number_of_densities(),0.0);//previous mole flux of solutes

  exterior_osmolality=binary_virial(parameters.doubles("initial_EG_only_molarity"),"EG");//total exterior osmolality salt+CPA (mole/kg)
  interior_osmolality=exterior_osmolality;// total internal osmolality salt+CPA
  interior_molarity.resize(microenvironment.number_of_densities(),0.0);//1.68 is the kdiss for NaCl from virial
  exterior_molarity=interior_molarity;
  interior_component_molality=interior_molarity;
  exterior_component_molality=exterior_molarity;
  Ps.resize(microenvironment.number_of_densities(),0.0);
  Lp=0.0;
  solute_uptake.resize(microenvironment.number_of_densities(),0.0);
  solute_volume=0.0;
  water_uptake=0.0;
  toxicity=0.0; 
  previous_radius=0; 
  //constants available to every spring cell, future version could have changing temp
  uptake_voxels={};
	return; 
} 
