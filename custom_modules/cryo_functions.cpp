#include "cryo_functions.h"
#include "cryo_utilities.h"
#include <cstddef>
#include <cstdlib>
namespace cryo_parameters{
//solute names and Lp Ps names in XML must all match or function will fail; HM=NaCl
//XML parameters should be of the form {type_name}_{Lp or Ps}_{solute_name} and should be defined for all cell types in the simulation
//solute name should be the same in the XML,microenvironment and cryo_utilities
void set_2p_initial_conditions(Cell* pCell, Cryo_Parameters &cryo_p)
{
  // output_cell_and_voxels(get_exterior_voxels(this->m_my_pCell), this->m_my_pCell);
  // 1 mol/L = 1 fmol/um^3 
  // double HM_density=1.0045; //NaCl density g/cm^3 extrapolated from CRC 
  // double solvent_volume= this->m_my_pCell->phenotype.volume.total-this->solid_volume;
  //set initial stat vectors
  if(cryo_p.Ps.size()==default_microenvironment_options.initial_condition_vector.size())
  {
    
    std::string Lp_name="";
    cryo_p.interior_molarity=default_microenvironment_options.initial_condition_vector;
    for(size_t i=0;i<default_microenvironment_options.initial_condition_vector.size();i++)
    {

      std::string solute_name=pCell->get_microenvironment()->density_names[i];
      std::string Ps_name=pCell->type_name+"_Ps_"+solute_name;
      Lp_name=pCell->type_name+"_Lp_"+solute_name;

      cryo_p.Ps[i]=parameters.doubles(Ps_name);
      
      if(cryo_p.Ps[i]==0)
      {
        cryo_p.solute_moles[i]=(cryo_p.interior_molarity[i])*(pCell->phenotype.volume.total-cryo_p.solid_volume);//HM osmolality/kdiss*kg/L solvent*density= femptomoles HM
      }
      else {
        cryo_p.solute_moles[i]=0.0;    
      }
      cryo_p.interior_component_molality[i]=molarity_to_molality(cryo_p.interior_molarity[i], solute_name);
      
    }
     
    cryo_p.Lp=parameters.doubles(Lp_name);
    if(pCell->get_microenvironment()->number_of_densities()==1)
    {
      cryo_p.interior_osmolality=binary_virial(cryo_p.interior_component_molality[0],pCell->get_microenvironment()->density_names[0]);
    }
    else if(pCell->get_microenvironment()->number_of_densities()==2)
    {
        cryo_p.interior_osmolality=ternary_virial(cryo_p.interior_component_molality[0],cryo_p.interior_component_molality[1],pCell->get_microenvironment()->density_names[0],pCell->get_microenvironment()->density_names[1]);
    }
    else if(pCell->get_microenvironment()->number_of_densities()==3)
    {
      cryo_p.interior_osmolality=binary_virial(cryo_p.interior_component_molality[0],pCell->get_microenvironment()->density_names[0])+ternary_virial(cryo_p.interior_component_molality[1],cryo_p.interior_component_molality[2],pCell->get_microenvironment()->density_names[1],pCell->get_microenvironment()->density_names[2]);
    }
    else {
      std::cout<<"WARNING!! UNKNOWN SOLUTE OR INVALID NUMBER OF DENSITIES!!"<<"\n";
      exit(1);
    }
  }
  else {
      std::cout<<"WARNING!! UNKNOWN SOLUTE OR INVALID NUMBER OF DENSITIES!!"<<"\n";
      exit(1);
  
  }
  return;
}
// using namespace Springs;
//  #include <stdef.h>
//   Global Variables - will probably add to XML eventually to avoid macros
};
