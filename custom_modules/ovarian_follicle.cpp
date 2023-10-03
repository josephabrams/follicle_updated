#include "./ovarian_follicle.h"
#include "./springs.h"
#include "./follicle_utilities.h"
#include <iostream>
#include <omp.h>
#include <algorithm>
// Global Variables
#define PI 3.14159265
#define R 0.08205//granulosa (10^-3 J/mole*k)
using namespace Springs;
Cell_Definition oocyte_cell;
Cell_Definition granulosa_cell;
Cell_Definition hepato_cell;

std::string output_filename;

double k_granulosa=0.3;
double k_oocyte=0.3;
double k_basement=0.3;
std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }
std::vector<std::string> follicle_coloring_function( Cell* pCell )//function for coloring cells in SVG output
{
	std::vector< std::string > output( 4, "black" );
	//std::vector< int > output( 4, 0.0 );
		output[0] = "white"; // orangered // "purple"; // 128,0,128
		output[1] = "black"; // "magenta";
		output[2] = "white"; // "magenta"; //255,0,255
		output[3] = "white" ;


	if( pCell->position[1] > 53 && pCell->position[1] < 103 && pCell->position[0]>-2 && pCell->position[0]<5 &&  pCell->position[2]<2 && pCell->position[2]>-2 )
	{
		output[0] = "green"; // orangered // "purple"; // 128,0,128
		output[1] = "black"; // "magenta";
		output[2] = "green"; // "magenta"; //255,0,255
		output[3] = "green" ;
		std::ofstream ofs;
		ofs.open ("output/color_point.csv", std::ofstream::out | std::ofstream::app);
		ofs <<"pCell names: "<<pCell<< " Position "<< pCell->position[0]<<","<<pCell->position[1]<<","<<pCell->position[2] <<"\n";
		ofs.close();
		
		return output;
	}
/*
	if(pCell->type==2)	
	{
		output[0] = "rgb(145,191,219)"; // orangered // "purple"; // 128,0,128
		output[1] = "black"; // "magenta";
		output[2] = "rgb(145,191,219)";
		output[3] = "rgb(145,191,219)";// "magenta"; //255,0,255
	}
	if(pCell->type==1)	
	{
		output[0] = "rgb(252,141,89)"; // orangered // "purple"; // 128,0,128
		output[1] = "black"; // "magenta";
		output[2] = "rgb(252,141,89)"; // "magenta"; //255,0,255
		output[3] = "rgb(252,141,89)";
	}
*/
	if( pCell->phenotype.death.dead == false )
	{


	}

	return output;
}
void spring_cell_functions(Cell* pCell, Phenotype& phenotype, double dt)
{
  
  Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index];
  //update_external_concentration and spring cell values
  //2P reaction terms and volume change
  //update internal concentations and spring cell values
  //pass it all back to pCell
  // two_parameter_single_step(pCell,phenotype,dt);
  //sum velocities from non_connected_neighbor_pressure
  //sum velocities from intercellular connections
  //sum velocities from basement membrane 
  update_all_forces(pCell,dt,pCell->custom_data["spring_constant"]);
  
  // if oocyte break connections that are past the max breakage distance
  // if(pCell->type==1){
    //function checks distance
    // break_TZPs(SPcell, pCell->custom_data["max_breakage_distance"]);
  // }
  return;
}
void output_tzp_score(Cell* pCell, Phenotype& phenotype, double dt)
{
  // check that cell is the oocyte
  if(pCell->type==1){
    Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index];
    double tzp_score=(double)(SPcell->m_springs.size())/(double)(SPcell->initial_number_of_connections);
    if (SPcell->is_outside==true)
    {tzp_score=-1.0;}
    if (PhysiCell_globals.current_time>299.9) {
      std::ofstream ofs;
      ofs.open ("./output/TZP_score.csv", std::ofstream::out | std::ofstream::app);
      ofs << parameters.ints("selected_simulation")<<", "<< default_microenvironment_options.Dirichlet_condition_vector<<PhysiCell_globals.current_time<<", "<<pCell->type<<", "<<tzp_score<<", "<< k_oocyte<<", "<<k_granulosa<<", "<< k_basement<<"\n";
      ofs.close();
    }     
  }

  return;
}
void create_cell_types( void )
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs
	SeedRandom( parameters.ints("random_seed") );
	// housekeeping
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	// turn the default cycle model to live,
	// so it's easier to turn off proliferation
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live );
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;
	// cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = false; // true;
	// set to no motility for cancer cells
	cell_defaults.phenotype.motility.is_motile = false;
	// use default proliferation and death
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
	// set default uptake and secretion
	// oxygen
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 0;
	// immunostimulatory
	cell_defaults.phenotype.secretion.saturation_densities[1] = 0;
	// set the default cell type to o2-based proliferation with the effect of the
	// on oncoprotein, and secretion of the immunostimulatory factor
	cell_defaults.functions.update_phenotype = NULL;
	// add the extra bit of "attachment" mechanics
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.name = "the default cell";
	cell_defaults.type = 0;
	// add custom data
	//Parameter<double> paramD;
  // "1/min" , 0.01 );  /* param */
  cell_defaults.custom_data.add_variable("spring_constant","unitless",0.0);
  cell_defaults.custom_data.add_variable("neighborhood_radius","um",0);
  cell_defaults.custom_data.add_variable("initial_volume","unitless",0);
  cell_defaults.custom_data.add_variable("allowed_overlap","um",2);//radius beyond cell surface where connections occur
  // cell_defaults.custom_data.add_variable("Lp","unitless",0);
  // cell_defaults.custom_data.add_variable("Ps","unitless",0);
  // cell_defaults.custom_data.add_variable("solute_1_permeability",0.0);
  cell_defaults.custom_data.add_variable("surface_area","unitless",0);
  cell_defaults.custom_data.add_variable("R","unitless",0);
  cell_defaults.custom_data.add_variable("Temperature","unitless",0);
  // cell_defaults.custom_data.add_variable("initial_water_volume","unitless",0);
  // cell_defaults.custom_data.add_variable("Partial_molar_volume","unitless",0);
  cell_defaults.custom_data.add_variable("Vb_fraction","unitless",0);
  // cell_defaults.custom_data.add_variable( "kill rate" , "1/min" , 0 ); // how often it tries to kill
  // cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached
  // cell_defaults.custom_data.add_variable( "attachment rate" , "1/min" ,0 ); // how long it wants to wander before attaching
  // cell_defaults.custom_data.add_variable("osmotically_active_water_volume","unitless",0.0);
  // cell_defaults.custom_data.add_variable("prev_osmotically_active_water_volume","unitless",0.0);
  cell_defaults.custom_data.add_variable("max_interaction_distance","unitless",0.5);
  //cell_defaults.custom_data.add_variable("total_internal_concentration","unitless",0.0);
  //cell_defaults.custom_data.add_variable("total_external_concentration", "unitless", 0.0);
  // cell_defaults.custom_data.add_variable("dVw_dt", "unitless", 0.0);
  //cell_defaults.custom_data.add_variable("dS_dt_1", "unitless", 0.0);
  //cell_defaults.custom_data.add_variable("external_concentration_of_solute_2","unitless",0);
  //cell_defaults.custom_data.add_variable("Prev_total_volume","unitless",0);
  //ell_defaults.custom_data.add_variable("test_uptake","unitless",0);
  //cell_defaults.custom_data.add_variable("test_current_voxel","unitless",0);
  // create the cell types
	create_oocyte_cell_type();
	create_granulosa_cell_type();
  create_hepato_cell_type();
	return;
}
//oocyte////////////////////
void create_hepato_cell_type(void)
{
	hepato_cell = cell_defaults;
  hepato_cell.name = "hepato cell";
  hepato_cell.type = 3;
	// turn off proliferation but turn on live cell;
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
	hepato_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;
  hepato_cell.phenotype.death.rates={0,0,0};
	hepato_cell.phenotype.secretion.uptake_rates[0] *=0;
	hepato_cell.phenotype.secretion.secretion_rates[0] *=0;
	hepato_cell.phenotype.secretion.uptake_rates[1] *=0;
	hepato_cell.phenotype.secretion.secretion_rates[1] *=0;
	hepato_cell.phenotype.secretion.uptake_rates[2] *=0;
	hepato_cell.phenotype.secretion.secretion_rates[2] *=0;
	// set functions to use custom functions note that position, and radius are still updated within the core code
  hepato_cell.functions.update_phenotype = hepato_phenotype_rule;
  hepato_cell.functions.custom_cell_rule = hepato_cell_rule;
  hepato_cell.functions.update_velocity = hepato_velocity_rule;
	// set custom data values
	hepato_cell.custom_data.add_variable("cell_k","unitless",k_oocyte);//oocyte spring constant
  hepato_cell.custom_data.add_variable("neighborhood_radius","um",15);//radius beyond cell surface where connections occur
  hepato_cell.custom_data.add_variable("allowed_overlap","um",2);//radius beyond cell surface where connections occur
	hepato_cell.custom_data[ "initial_volume" ] = 14137.17;//parameters.doubles("oocyte_isotonic_volume");//817283 um^3
	/*
		DMSO_Lp 1.68 Immature 1.01 Mature
		DMSO_Ps .15 immature .24 mature
		PROH_Lp .72          .86
		PROH_Ps .31          .56
	*/
	hepato_cell.custom_data["surface_area"]=2727.43;//43411.1;
	hepato_cell.custom_data["R"]= R;
	hepato_cell.custom_data["Temperature"]= 293.15;// uh I should get all the temperature stuff lined up
	hepato_cell.custom_data["Vb_fraction"]=.288;
  hepato_cell.custom_data.add_variable("basement_k","unitless",k_basement);
	hepato_cell.custom_data.add_variable("spring_constant","N/m",k_oocyte);
	hepato_cell.custom_data.add_variable("max_breakage_distance","um",10);
  hepato_cell.custom_data.add_variable("max_interaction_distance","unitless",15);
	return;
}
int counting_time=0;
void hepato_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    // std::cout<<"VOLUME:"<<pCell->phenotype.volume.total<<"\n"; 
  //std::vector <double> basement_membrane_center={0.0,0.0,0.0};
	return;
}
void hepato_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
  Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index];
  two_parameter_single_step(pCell,phenotype, dt);
  update_all_forces(pCell,dt,pCell->custom_data["spring_constant"]);
  // if(PhysiCell_globals.current_time<dt){
    // SPcell->initial_number_of_connections=SPcell->m_springs.size();
  // }
  // output_tzp_score(pCell, phenotype, dt);
  int thread_id=omp_get_thread_num();
  // spring_cell_functions(pCell,phenotype,dt);	
  std::ofstream ofs;
	ofs.open ("./output/hepato.csv", std::ofstream::out | std::ofstream::app);
  ofs <<parameters.ints("selected_simulation")<<", "<< default_microenvironment_options.Dirichlet_condition_vector<<", "<<PhysiCell_globals.current_time<<", "<<pCell->type_name<<", "<< pCell->phenotype.volume.total<<", "<<pCell->position[0]<<", "<<pCell->position[1]<<", "<<pCell->position[2]<<", "<<SPcell->m_springs.size() <<", "<<k_oocyte<<", "<<k_granulosa<<", "<<k_basement<<"\n";
	ofs.close();
  // std::cout<<"VOLUME:"<<pCell->phenotype.volume.total<<"\n"; 
  
	return;

}
void hepato_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	return;
}
void create_oocyte_cell_type(void)
{
	oocyte_cell = cell_defaults;
	oocyte_cell.name = "oocyte cell";
	oocyte_cell.type = 1;
	// turn off proliferation but turn on live cell;
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
	oocyte_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;
  int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	//
	oocyte_cell.phenotype.secretion.uptake_rates[0] *=0;
	oocyte_cell.phenotype.secretion.secretion_rates[0] *=0;
	oocyte_cell.phenotype.secretion.uptake_rates[1] *=0;
	oocyte_cell.phenotype.secretion.secretion_rates[1] *=0;
	oocyte_cell.phenotype.secretion.uptake_rates[2] *=0;
	oocyte_cell.phenotype.secretion.secretion_rates[2] *=0;
	// set functions to use custom functions note that position, and radius are still updated within the core code
	oocyte_cell.functions.update_phenotype = oocyte_phenotype_rule;
	oocyte_cell.functions.custom_cell_rule = oocyte_cell_rule;
	oocyte_cell.functions.update_velocity = oocyte_velocity_rule;
	// set custom data values
	oocyte_cell.custom_data.add_variable("cell_k","unitless",k_oocyte);//oocyte spring constant
  oocyte_cell.custom_data.add_variable("neighborhood_radius","um",15);//radius beyond cell surface where connections occur
  oocyte_cell.custom_data.add_variable("allowed_overlap","um",2);//radius beyond cell surface where connections occur
	oocyte_cell.custom_data[ "initial_volume" ] = parameters.doubles("oocyte_isotonic_volume");//817283 um^3
	/*
		DMSO_Lp 1.68 Immature 1.01 Mature
		DMSO_Ps .15 immature .24 mature
		PROH_Lp .72          .86
		PROH_Ps .31          .56
	*/
	oocyte_cell.custom_data["surface_area"]=42273.3;//43411.1;
	oocyte_cell.custom_data["R"]= R;
	oocyte_cell.custom_data["Temperature"]= 293.15;// uh I should get all the temperature stuff lined up
	oocyte_cell.custom_data["Vb_fraction"]=.288;
  oocyte_cell.custom_data.add_variable("basement_k","unitless",k_basement);
	oocyte_cell.custom_data.add_variable("spring_constant","N/m",k_oocyte);
	oocyte_cell.custom_data.add_variable("max_breakage_distance","um",10);
  oocyte_cell.custom_data.add_variable("max_interaction_distance","unitless",15);
	return;
}

void oocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    // std::cout<<"VOLUME:"<<pCell->phenotype.volume.total<<"\n"; 
  //std::vector <double> basement_membrane_center={0.0,0.0,0.0};
	return;
}
void oocyte_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
  // Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index];
  // two_parameter_single_step(pCell,phenotype, dt);
  // if(PhysiCell_globals.current_time<dt){
    // SPcell->initial_number_of_connections=SPcell->m_springs.size();
  // }
  // output_tzp_score(pCell, phenotype, dt);
  // int thread_id=omp_get_thread_num();
  // spring_cell_functions(pCell,phenotype,dt);	
 //  std::ofstream ofs;
	// ofs.open ("./output/oocyte_output.csv", std::ofstream::out | std::ofstream::app);
 //  ofs <<parameters.ints("selected_simulation")<<", "<< default_microenvironment_options.Dirichlet_condition_vector<<", "<<PhysiCell_globals.current_time<<", "<<pCell->type_name<<", "<< pCell->phenotype.volume.total<<", "<<pCell->position[0]<<", "<<pCell->position[1]<<", "<<pCell->position[2]<<", "<<SPcell->m_springs.size() <<", "<<k_oocyte<<", "<<k_granulosa<<", "<<k_basement<<"\n";
	// ofs.close();
  // std::cout<<"VOLUME:"<<pCell->phenotype.volume.total<<"\n"; 
  
	return;

}
void oocyte_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	return;
}
//!granulosas
void create_granulosa_cell_type( void )
{
	granulosa_cell = cell_defaults;
	granulosa_cell.name = "granulosa cell";
	granulosa_cell.type = 2;
	// turn off proliferation and turn on live cell;
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
	granulosa_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );
	granulosa_cell.phenotype.secretion.uptake_rates[0] *=0;
	granulosa_cell.phenotype.secretion.secretion_rates[0] *=0;
	granulosa_cell.phenotype.secretion.uptake_rates[1] *=0;
	granulosa_cell.phenotype.secretion.secretion_rates[1] *=0;
	granulosa_cell.phenotype.secretion.uptake_rates[2] *=0;
	granulosa_cell.phenotype.secretion.secretion_rates[2] *=0;
	//turn off apoptosis
	granulosa_cell.phenotype.death.rates[apoptosis_index] = 0;
	// turn off motility;
	granulosa_cell.phenotype.motility.is_motile = false;
  //turn off built in mechanics
	granulosa_cell.phenotype.mechanics.cell_cell_adhesion_strength *=0;
  granulosa_cell.phenotype.mechanics.cell_cell_repulsion_strength *=0;
	//
	granulosa_cell.phenotype.secretion.uptake_rates[0] *=0;
	granulosa_cell.phenotype.secretion.secretion_rates[0] *=0;
	granulosa_cell.phenotype.secretion.uptake_rates[1] *=0;
	granulosa_cell.phenotype.secretion.secretion_rates[1] *=0;
	granulosa_cell.phenotype.secretion.uptake_rates[2] *=0;
	granulosa_cell.phenotype.secretion.secretion_rates[2] *=0;
	// vector mass transport parameters inserted into Phenotype class by Joseph S Abrams, vectors are the same size as the number of densities in BioFVM
	//external concentrations of densities
    /*
	granulosa_cell.phenotype.secretion.Mse[0]=0.0;
	granulosa_cell.phenotype.secretion.Mse[1]=0.0;
	//permeating moles
	granulosa_cell.phenotype.secretion.s[0]=0.0;
	granulosa_cell.phenotype.secretion.s[1]=0.0;
	//permeating internal concentration
	granulosa_cell.phenotype.secretion.Msi[0]=0.0;
	granulosa_cell.phenotype.secretion.Msi[1]=0.0;
	granulosa_cell.phenotype.secretion.Msi[2]=0.0;
	//nonpermeating external concentration
	granulosa_cell.phenotype.secretion.Mne[0]=0.2534;//holding media
	//nonpermeating internal concentration
	granulosa_cell.phenotype.secretion.Mni[0]=0.2534;//isosmotic
	granulosa_cell.phenotype.secretion.ni[0]=0.2534*418.88;//femptomoles of np solute initial concentration/initial osmoticly active water volume
	//permeabilities
	granulosa_cell.phenotype.secretion.Ps[1]=Ps_granulosa;//EG
	granulosa_cell.phenotype.secretion.Ps[2]=Ps_granulosa_gly;//GLY
	//fluxes of moles for passing back to bioFVM in updating concentration functions
	granulosa_cell.phenotype.secretion.flux[1]=0.0;
	granulosa_cell.phenotype.secretion.prev_s[0] =0.0;
	granulosa_cell.phenotype.secretion.prev_s[1] =0.0;
	granulosa_cell.phenotype.secretion.prev_s[2] =0.0;
	granulosa_cell.phenotype.secretion.prev_dS_dt[0] =0.0;
	granulosa_cell.phenotype.secretion.prev_dS_dt[1] =0.0;
	granulosa_cell.phenotype.secretion.prev_dS_dt[2] =0.0;
	granulosa_cell.phenotype.secretion.delta[0]=1;
	granulosa_cell.phenotype.secretion.delta[1]=1;
	granulosa_cell.phenotype.secretion.delta[2]=1;
	granulosa_cell.phenotype.secretion.partial_molar_volumes[0]=0;
	granulosa_cell.phenotype.secretion.partial_molar_volumes[1]=0.0557414;//EG
	granulosa_cell.phenotype.secretion.partial_molar_volumes[2]=0.0730903;//GLY
	// set functions
  */
	granulosa_cell.functions.update_phenotype = granulosa_phenotype_rule;
	granulosa_cell.functions.custom_cell_rule = granulosa_cell_rule;
	granulosa_cell.functions.update_velocity = granulosa_velocity_rule;
	// set custom data values
	// Parameter<double> paramD;
	granulosa_cell.custom_data.add_variable("basement_k","unitless",k_basement);
  granulosa_cell.custom_data.add_variable("neighborhood_radius","unitless",7.5);//radius beyond cell surface where connections occur
	granulosa_cell.custom_data[ "initial_volume" ] = parameters.doubles("gran_isotonic_volume");//523.6;
	granulosa_cell.custom_data.add_variable("connection_length","micometers",0.5);
	granulosa_cell.custom_data["surface_area"]= 314.159;//parameters.doubles("granulosa_Area");
	granulosa_cell.custom_data["R"]= R;//R;//parameters.doubles("granulosa_R");
	granulosa_cell.custom_data["Temperature"]=293.15;//parameters.doubles("granulosa_Temperature");
	granulosa_cell.custom_data["Vb_fraction"]=.2;
	granulosa_cell.custom_data.add_variable("cell_k","unitless",k_granulosa);
	granulosa_cell.custom_data.add_variable("spring_constant","N/m",k_granulosa);
  granulosa_cell.custom_data.add_variable("allowed_overlap","um",2);//radius beyond cell surface where connections occur
  granulosa_cell.custom_data.add_variable("max_breakage_distance","um",10);
  granulosa_cell.custom_data.add_variable("max_interaction_distance","unitless",10);
	return;
}
void granulosa_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{

}
int count=0;

void granulosa_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
  Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index];
  int thread_id=omp_get_thread_num();
  two_parameter_single_step(pCell,phenotype, dt);
  update_all_forces(pCell,dt,pCell->custom_data["spring_constant"]);
 //  // std::cout<<" FORCE PARAM: "<< k_oocyte<<", "<<k_granulosa<<", "<<k_basement<<"\n";
  // spring_cell_functions(pCell,phenotype,dt);
  // if(count%100==0){
    std::ofstream ofs;
    ofs.open ("./output/2p_Test_granulosa.csv", std::ofstream::out | std::ofstream::app);
    ofs << PhysiCell_globals.current_time<<", "<<pCell->index<<", "<<SPcell->m_my_pCell->index<<", "<< pCell<<", "<< pCell->phenotype.volume.total<<", "<<SPcell->exterior_osmolality<<", "<<SPcell->interior_osmolality <<", "<<SPcell->exterior_molarity[1]<<", "<<SPcell->interior_molarity[1]<<", "<<PhysiCell_settings.omp_num_threads<<", "<< thread_id<<"\n";
    ofs.close();
  //   count++;
  // }
  // else{
  //   count++;
  // }
 //  std::ofstream ofs;
	// ofs.open ("./output/granulosa_output.csv", std::ofstream::out | std::ofstream::app);
	// ofs << parameters.ints("selected_simulation")<<", "<< default_microenvironment_options.Dirichlet_condition_vector<<", "<<PhysiCell_globals.current_time<<", "<<pCell->type_name<<", "<<pCell<<", "<< pCell->phenotype.volume.total<<", "<<pCell->position[0]<<", "<<pCell->position[1]<<", "<<pCell->position[2]<<SPcell->m_springs.size() <<", "<<k_oocyte<<", "<<k_granulosa<<", "<<k_basement<<", "<<SPcell->is_basement_connected<<"\n";
	// ofs.close();
  return;
}
void granulosa_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	//std::vector <double> basement_membrane_center={0.0,0.0,0.0};
	return;
}
//! /////////////////////////////PhysiCell///////////////////////////////////////////////////////////////////
void setup_microenvironment( void )
{
  default_microenvironment_options.use_oxygen_as_first_field = false;
	// set domain parameters
  /*
	default_microenvironment_options.X_range = {-1000, 1000};
	default_microenvironment_options.Y_range = {-1000, 1000};
	default_microenvironment_options.Z_range = {-1000, 1000};
  */
  /*
	// now in XML
	default_microenvironment_options.X_range = {-750, 750};
	default_microenvironment_options.Y_range = {-750, 750};
	default_microenvironment_options.Z_range = {-750, 750};
  */
  // resizing is key it fixes the initial oxygen problem correctly and makes sure everything lines up correctly!!! 
  microenvironment.resize_densities(0);
	// if( default_microenvironment_options.simulate_2D == true )
	// {
  // std::cout << "Warning: overriding 2D setting to return to 3D" << "\n";
  // default_microenvironment_options.simulate_2D = false;
	// }
	// gradients are needed for this example
	//default_microenvironment_options.calculate_gradients = true;
  // microenvironment.add_density( "blank", "dimensionless" );
  // microenvironment.diffusion_coefficients[0] = 5555;
  // microenvironment.decay_rates[0] = 0;
	// add the immunostimulatory factor
	// let's do these in XML later
  static int simulation_selected =parameters.ints( "selected_simulation" ); 
  static int selected_PBS_simulation =parameters.ints("selected_PBS_simulation");
  microenvironment.add_density( "HM", "dimensionless" );
	microenvironment.diffusion_coefficients[0] = 510;
	microenvironment.decay_rates[0] = 0;
  if(simulation_selected==1)
  {
    std::cout<<"running EG loading"<<"\n";
    microenvironment.add_density( "EG", "dimensionless" );
    microenvironment.diffusion_coefficients[1] = 530;
    microenvironment.decay_rates[1] = 0;
    // default_microenvironment_options.initial_condition_vector[2]=2.2;
    // default_microenvironment_options.Dirichlet_condition_vector[2] = 2.2;
    // default_microenvironment_options.Dirichlet_activation_vector[2]=false;
    // default_microenvironment_options.Dirichlet_xmax_values[2]=3.4;
	  //solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);
    default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    default_microenvironment_options.initial_condition_vector[1]=0.00;
    default_microenvironment_options.Dirichlet_condition_vector[1] = parameters.doubles("initial_EG_only_molarity");//parameters.doubles("initial_EG_only_molarity");//1.00;//6.66;
    default_microenvironment_options.Dirichlet_activation_vector[1]=true;	
  }
  else if(simulation_selected==2)
  {
    std::cout<<"running GLY loading"<<"\n";
    microenvironment.add_density( "GLY", "dimensionless" );
	  microenvironment.diffusion_coefficients[1] = 530;
	  microenvironment.decay_rates[1] = 0;
    default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    default_microenvironment_options.initial_condition_vector[1]=0.00;
    default_microenvironment_options.Dirichlet_condition_vector[1] = parameters.doubles("initial_GLY_only_molarity");//1.00;//6.66;
    default_microenvironment_options.Dirichlet_activation_vector[1]=true;	
  	//solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);
  }
  else if(simulation_selected==3)
  {
    //select which PBS simulation, default is 1x
    std::cout<<"running PBS loading"<<"\n";
    microenvironment.add_density( "PBS", "dimensionless" );
    microenvironment.diffusion_coefficients[0] = 510;
    microenvironment.decay_rates[0] = 0;
    if (selected_PBS_simulation==0) {
      //0.5x PBS simulation
      default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
      default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_05xPBS_molarity");
      default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    }
    else if (selected_PBS_simulation==2) {
      default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
      default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_2xPBS_molarity");
      default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
      //2x PBS simulation
    }
    else if (selected_PBS_simulation==5) {
      //5x PBS simulation 
      default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
      default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_5xPBS_molarity");
      default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    }
    else{
      //1x PBS simulation 
      default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
      default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_1xPBS_molarity");
      default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    }
    //solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);
  }
  else if(simulation_selected==4)
  { 
    std::cout<<"running EG & GLY loading"<<"\n";
    microenvironment.add_density( "EG", "dimensionless" );
	  microenvironment.diffusion_coefficients[1] = 530;
	  microenvironment.decay_rates[1] = 0;
    microenvironment.add_density( "GLY", "dimensionless" );
	  microenvironment.diffusion_coefficients[2] = 530;
	  microenvironment.decay_rates[2] = 0;
    default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    default_microenvironment_options.initial_condition_vector[1]=0.00;
    default_microenvironment_options.Dirichlet_condition_vector[1] = parameters.doubles("initial_EG_mixed_molarity");//1.00;//6.66;
    default_microenvironment_options.Dirichlet_activation_vector[1]=true;	
    default_microenvironment_options.initial_condition_vector[2]=0.00;
    default_microenvironment_options.Dirichlet_condition_vector[2] = parameters.doubles("initial_GLY_mixed_molarity");//1.00;//6.66;
    default_microenvironment_options.Dirichlet_activation_vector[2]=true;	
  }
  else if(simulation_selected==5)
  {
    std::cout<<"running Karlsson EG loading"<<"\n";
    microenvironment.add_density( "EG", "dimensionless" );
    microenvironment.diffusion_coefficients[1] = 530;
    microenvironment.decay_rates[1] = 0;
    // default_microenvironment_options.initial_condition_vector[2]=2.2;
    // default_microenvironment_options.Dirichlet_condition_vector[2] = 2.2;
    // default_microenvironment_options.Dirichlet_activation_vector[2]=false;
    // default_microenvironment_options.Dirichlet_xmax_values[2]=3.4;
	  //solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);
    default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    default_microenvironment_options.initial_condition_vector[1]=0.00;
    default_microenvironment_options.Dirichlet_condition_vector[1] = parameters.doubles("initial_EG_Karlsson_molarity");//1.00;//6.66;
    default_microenvironment_options.Dirichlet_activation_vector[1]=true;	
  }
  else if(simulation_selected==6)
  {
    std::cout<<"running Karlsson DMSO loading"<<"\n";
    microenvironment.add_density( "DMSO", "dimensionless" );
    microenvironment.diffusion_coefficients[1] = 530;
    microenvironment.decay_rates[1] = 0;
    // default_microenvironment_options.initial_condition_vector[2]=2.2;
    // default_microenvironment_options.Dirichlet_condition_vector[2] = 2.2;
    // default_microenvironment_options.Dirichlet_activation_vector[2]=false;
    // default_microenvironment_options.Dirichlet_xmax_values[2]=3.4;
	  //solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);
    default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    default_microenvironment_options.initial_condition_vector[1]=0.00;
    default_microenvironment_options.Dirichlet_condition_vector[1] = parameters.doubles("initial_DMSO_Karlsson_molarity");//1.00;//6.66;
    default_microenvironment_options.Dirichlet_activation_vector[1]=true;	
  }
  else if(simulation_selected==7)
  {
    std::cout<<"running Karlsson PROH loading"<<"\n";
    microenvironment.add_density( "PROH", "dimensionless" );
    microenvironment.diffusion_coefficients[1] = 530;
    microenvironment.decay_rates[1] = 0;
    // default_microenvironment_options.initial_condition_vector[2]=2.2;
    // default_microenvironment_options.Dirichlet_condition_vector[2] = 2.2;
    // default_microenvironment_options.Dirichlet_activation_vector[2]=false;
    // default_microenvironment_options.Dirichlet_xmax_values[2]=3.4;
	  //solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);
    default_microenvironment_options.initial_condition_vector[0]=parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_condition_vector[0] = parameters.doubles("initial_HM_molarity");
    default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
    default_microenvironment_options.initial_condition_vector[1]=0.00;
    default_microenvironment_options.Dirichlet_condition_vector[1] = parameters.doubles("initial_PROH_Karlsson_molarity");//1.00;//6.66;
    default_microenvironment_options.Dirichlet_activation_vector[1]=true;	
  }
  else
  {
    std::cout<<"WARNING!!! running all solutes test case. THIS IS NOT A SET SIMULATION!!!"<<"\n";
    microenvironment.add_density( "HM (PBS)", "dimensionless" );
    microenvironment.diffusion_coefficients[0] = 510;
    microenvironment.decay_rates[0] = 0;
  
    microenvironment.add_density( "EG", "dimensionless" );
	  microenvironment.diffusion_coefficients[1] = 530;
	  microenvironment.decay_rates[1] = 0;
	
    microenvironment.add_density( "GLY", "dimensionless" );
    microenvironment.diffusion_coefficients[2] = 530;//530; //TODO: check about the diffusivity of this value
    microenvironment.decay_rates[2] = 0;
    
    microenvironment.add_density( "PBS", "dimensionless" );
    microenvironment.diffusion_coefficients[3] = 510;
    microenvironment.decay_rates[3] = 0;

    microenvironment.add_density( "SUC", "dimensionless" );
    microenvironment.diffusion_coefficients[4] = 1075;
    microenvironment.decay_rates[4] = 0;
  }

  // default_microenvironment_options.initial_condition_vector[0]=0.155;
	// default_microenvironment_options.Dirichlet_condition_vector[0] = 0.155;
	// default_microenvironment_options.Dirichlet_activation_vector[0]=true;	
  // default_microenvironment_options.initial_condition_vector[1]=0.00;
	// default_microenvironment_options.Dirichlet_condition_vector[1] = 6.66;
	// default_microenvironment_options.Dirichlet_activation_vector[1]=true;	
	// set Dirichlet conditions

	//default_microenvironment_options.Dirichlet_condition_vector[2] = 0;

	//solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);
    // std::cout<<"Size initial "<< default_microenvironment_options.initial_condition_vector.size()<<"\n"; 
    // std::cout<<"desnity size initial "<< microenvironment.number_of_densities()<<"\n"; 
	default_microenvironment_options.outer_Dirichlet_conditions = false;
	initialize_microenvironment();
  //dirichlet_nodes_radius is the CPA diffusion to tissue/cell distance, set close to follicle or oocyte for rapid change in exterior conditions
  double dirichlet_nodes_radius=400.0;

  std::vector<std::vector<double>> hole_centers={{0,0,0},{250,300,0},{-250,300,0},
    {300,0,0},{-300,0,0},{250,-300,0},{-250,-300,0}};
  std::vector<double>hole_radii={50,25,25,25,25,25,25};
  // std::vector<std::vector<double>> hole_dirchlet={default_microenvironment_options.Dirichlet_condition_vector,{0.155,0},{0.155,0},{0.155,0},{0.155,0},{0.155,0},{0.155,0}};
  std::vector<std::vector<double>> hole_dirchlet={{0.155,0},{0.155,0},{0.155,0},{0.155,0},{0.155,0},{0.155,0},{0.155,0}};
  // activate_nodes(dirichlet_nodes_radius);//place the Dirichlet nodes at a radius around the follicle 
 for (int j=0;j<hole_dirchlet.size();j++)
  {
    std::cout<<"HOLE DIRCHLET: "<< hole_dirchlet[j]<<"\n";
    set_dirchlet_nodes_cylinder(hole_radii[j],50,hole_centers[j],hole_dirchlet[j]); 

  }
  activate_edges(400);	
  //total_uptake.resize(microenvironment.mesh.voxels.size());
	//all_total_uptakes.resize(number_of_permeating_solutes,total_uptake);
	return;
}
void setup_secondary_stage_follicle(void)
{
  double initial_granulosa_radius=parameters.doubles("gran_radius"); 
  double initial_oocyte_radius=parameters.doubles("oocyte_radius");
	double initial_overlap=.1;//representation of the initial packing density -- taken from physicell simulations of cancer spheroids
	double follicle_radius =initial_oocyte_radius+3*(initial_granulosa_radius*2)+initial_granulosa_radius; // 3 to 4 layers
	
  double cell_spacing = initial_granulosa_radius-initial_overlap;//slight overlap to represent cells up against each other better variable name
	//creating the oocyte
	std::cout << "\tPlacing oocyte cells ... "<< "\n";
	// place the oocyte cell at origin
	Cell* pCell_oocyte;
	pCell_oocyte = create_cell(oocyte_cell);//sets cell type
	pCell_oocyte->assign_position(0,0,0); //assign oocyte position
	//pC->set_total_volume(817283); //sets initial total volume manually
	pCell_oocyte->set_total_volume(pCell_oocyte->custom_data["initial_volume"]);//will set cell to custom value
	//create granulosa
	std::cout << "\tPlacing granulosa cells ... " << "\n";
	std::vector<std::vector<double>> granulosa_positions = create_cell_sphere_positions(cell_spacing,follicle_radius,initial_oocyte_radius);//hollow spheroid around centered oocyte
	std::cout << "creating " << granulosa_positions.size() << " granulosa cells " << "\n";
	for( int i=0; i < granulosa_positions.size(); i++ )
	{
		Cell* pCell_granulosa;//placed this inside of the loop since it makes more sense to me
		pCell_granulosa = create_cell(granulosa_cell); //
		pCell_granulosa->assign_position( granulosa_positions[i] );
		pCell_granulosa->set_total_volume(pCell_granulosa->custom_data["initial_volume"]);
	}
  std::cout<<"Follicle is... "<< follicle_radius<<" um in radius."<<"\n";
	return;
}

void setup_n_layer_follicle(int n_layers)
{
    /*
		follicle_radius =initial_oocyte_radius+n_layers*(initial_granulosa_radius*2)+initial_granulosa_radius; // n to n+1 layers
		double initial_overlap=.1;//representation of the initial packing density from config eventually
		double cell_spacing = initial_granulosa_radius-initial_overlap;//slight overlap to represent cells up against each other better variable name
		//creating the oocyte
		std::cout << "\tPlacing oocyte cells ... "<< "\n";
		// place the oocyte cell at origin
		Cell* pCell_oocyte;
		pCell_oocyte = create_cell(oocyte_cell);//sets cell type
		pCell_oocyte->assign_position(0,0,0); //assign oocyte position
		//pC->set_total_volume(817283); //sets initial total volume manually
		pCell_oocyte->set_total_volume(pCell_oocyte->custom_data["initial_volume"]);//will set cell to custom value
		//create granulosa
		std::cout << "\tPlacing granulosa cells ... " << "\n";
		std::vector<std::vector<double>> granulosa_positions = create_cell_sphere_positions(cell_spacing,follicle_radius,initial_oocyte_radius);//hollow spheroid around centered oocyte
		std::cout << "creating " << granulosa_positions.size() << " granulosa cells " << "\n";
		//ofs.open ("initial_granulosa_size.csv", std::ofstream::out | std::ofstream::trunc);
		for( int i=0; i < granulosa_positions.size(); i++ )
		{
			Cell* pCell_granulosa;//placed this inside of the loop since it makes more sense to me
			pCell_granulosa = create_cell(granulosa_cell); //
			pCell_granulosa->assign_position( granulosa_positions[i] );
			pCell_granulosa->set_total_volume(pCell_granulosa->custom_data["initial_volume"]);
		}
   */
	return;
}
void setup_toy_granulosa_model(void)
{
    /*
	std::cout << "\tPlacing granulosa cell 1 ... " << "\n";
	Cell* pCell_granulosa1;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa1 = create_cell(granulosa_cell); //
	pCell_granulosa1->assign_position(-25,5,0 );
	pCell_granulosa1->set_total_volume(pCell_granulosa1->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 2 ... " << "\n";
	Cell* pCell_granulosa2;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa2 = create_cell(granulosa_cell); //
	pCell_granulosa2->assign_position(-15,5,0 );
	pCell_granulosa2->set_total_volume(pCell_granulosa2->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 3 ... " << "\n";
	Cell* pCell_granulosa3;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa3 = create_cell(granulosa_cell); //
	pCell_granulosa3->assign_position(-5,5,0 );
	pCell_granulosa3->set_total_volume(pCell_granulosa3->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 4 ... " << "\n";
	Cell* pCell_granulosa4;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa4 = create_cell(granulosa_cell); //
	pCell_granulosa4->assign_position(5,5,0 );
	pCell_granulosa4->set_total_volume(pCell_granulosa4->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 5 ... " << "\n";
	Cell* pCell_granulosa5;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa5 = create_cell(granulosa_cell); //
	pCell_granulosa5->assign_position(15,5,0 );
	pCell_granulosa5->set_total_volume(pCell_granulosa5->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 6 ... " << "\n";
	Cell* pCell_granulosa6;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa6 = create_cell(granulosa_cell); //
	pCell_granulosa6->assign_position(25,5,0 );
	pCell_granulosa6->set_total_volume(pCell_granulosa6->custom_data["initial_volume"]);
*/
return;
}

void setup_just_oocyte(void)
{
		//double initial_overlap=.1;//representation of the initial packing density from config eventually
		//double cell_spacing = initial_granulosa_radius-initial_overlap;//slight overlap to represent cells up against each other better variable name
		//creating the oocyte
		std::cout << "\tPlacing oocyte cells ... "<< "\n";
		// place the oocyte cell at origin
		Cell* pCell_oocyte;
		pCell_oocyte = create_cell(oocyte_cell);//sets cell type
		pCell_oocyte->assign_position(0,0,0); //assign oocyte position
		//pC->set_total_volume(817283); //sets initial total volume manually
		pCell_oocyte->set_total_volume(pCell_oocyte->custom_data["initial_volume"]);//will set cell to custom value
    std::cout<<"Initial VOLUME: "<< pCell_oocyte->phenotype.volume.total;
  return;
}

void setup_4_granulosa_test_case(void)
{
	//set up test follicle with one follicle and 4 granulosa\
	//creating the oocyte
	std::cout << "\tPlacing oocyte cells ... "<< "\n";
	// place the oocyte cell at origin
	Cell* pCell_oocyte;
	pCell_oocyte = create_cell(oocyte_cell);//sets cell type
	pCell_oocyte->assign_position(0,0,0); //assign oocyte position
	pCell_oocyte->set_total_volume(pCell_oocyte->custom_data["initial_volume"]);//will set cell to custom value
	//8 granulosa

  std::cout << "\tPlacing granulosa cell 1a ... " << "\n";
	Cell* pCell_granulosa1a;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa1a = create_cell(granulosa_cell); //
	pCell_granulosa1a->assign_position(5,20,0 );
	pCell_granulosa1a->set_total_volume(pCell_granulosa1a->custom_data["initial_volume"]);
	
	std::cout << "\tPlacing granulosa cell 1b ... " << "\n";
	Cell* pCell_granulosa1b;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa1b = create_cell(granulosa_cell); //
	pCell_granulosa1b->assign_position(-7,28,0 );
	pCell_granulosa1b->set_total_volume(pCell_granulosa1b->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 2a ... " << "\n";
	Cell* pCell_granulosa2a;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa2a = create_cell(granulosa_cell); //
	pCell_granulosa2a->assign_position(40,3.5,0 );
	pCell_granulosa2a->set_total_volume(pCell_granulosa2a->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 2b ... " << "\n";
	Cell* pCell_granulosa2b;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa2b = create_cell(granulosa_cell); //
	pCell_granulosa2b->assign_position(100,-7.8,0 );
	pCell_granulosa2b->set_total_volume(pCell_granulosa2b->custom_data["initial_volume"]);
	
	std::cout << "\tPlacing granulosa cell 3a ... " << "\n";
	Cell* pCell_granulosa3a;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa3a = create_cell(granulosa_cell); //
	pCell_granulosa3a->assign_position(17,-37,0 );
	pCell_granulosa3a->set_total_volume(pCell_granulosa3a->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 3b ... " << "\n";
	Cell* pCell_granulosa3b;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa3b = create_cell(granulosa_cell); //
	pCell_granulosa3b->assign_position(-5,-62,0 );
	pCell_granulosa3b->set_total_volume(pCell_granulosa3b->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 4a ... " << "\n";
	Cell* pCell_granulosa4a;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa4a = create_cell(granulosa_cell); //
	pCell_granulosa4a->assign_position(-20,2,0 );
	pCell_granulosa4a->set_total_volume(pCell_granulosa4a->custom_data["initial_volume"]);

	std::cout << "\tPlacing granulosa cell 4b ... " << "\n";
	Cell* pCell_granulosa4b;//placed this inside of the loop since it makes more sense to me
	pCell_granulosa4b = create_cell(granulosa_cell); //
	pCell_granulosa4b->assign_position(-100,-6,0 );
	pCell_granulosa4b->set_total_volume(pCell_granulosa4b->custom_data["initial_volume"]);
	
	//Cell* pCell_to_spring_cell=(*all_cells)[0];
	//Spring_Cell oocyte(pCell_to_spring_cell);
	//print_voxels_for_quick_plotting(pCell_granulosa4b,basement_membrane_voxels, basement_membrane_voxels);
	return;
}
std::vector<std::vector<double>> create_polygon_positions(double cell_radius,std::vector<std::vector<double>> polygon_points, std::vector<std::vector<double>> hole_centers, std::vector<double> hole_radius)//wedge between y=2x and y=-2x
{
 //  std::vector<double> x_points(polygon_points.size(),0.0);
 //  std::vector<double> y_points(polygon_points.size(),0.0);
 //  std::vector<double> z_points(polygon_points.size(),0.0);
 //  for (int i=0; i<polygon_points.size(); i++)
 //  {
 //     x_points[i]=polygon_points[i][0];
 //     y_points[i]=polygon_points[i][1];
 //     z_points[i]=polygon_points[i][2];
 //     }
 //  #pragma omp critical
 //  {
 //    std::sort(x_points.begin(),x_points.end());
 //    std::sort(y_points.begin(),y_points.end());
 //    std::sort(z_points.begin(),z_points.end());
 //  }
 //  int x_middle=x_points.size()/2;
 //  int y_middle=y_points.size()/2;
 //  int z_middle=z_points.size()/2;
	//
 //  double sphere_radius=100;
 //  double inner_radius=0;
	std::vector<std::vector<double>> cells;
  std::vector<double> v1={0,0,0};
  cells.push_back(v1);	
  // int xc=0,yc=0,zc=0;
	// double x_spacing= cell_radius*sqrt(3);
	// double y_spacing= cell_radius*2;
	// double z_spacing= cell_radius*sqrt(3);
	// double z_bottom=-20;
	// double z_top=20;
	// std::vector<double> tempPoint(3,0.0);
	// // std::vector<double> cylinder_center(3,0.0);
 //  //non-cannoical for loops - I take it setup is not parrallel 
	// for(double z=z_points.front(); z<z_points.back(); z+=z_spacing, zc++)
	// {
	// 	for(double x=x_points.front();x<x_points.back();x+=x_spacing, xc++)
	// 	{
	// 		for(double y=y_points.front();y<y_points.back();y+=y_spacing, yc++)
	// 		{
	// 			tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
	// 			tempPoint[1]=y + (xc%2) * cell_radius;
	// 			tempPoint[2]=z;
	//
	// 			if(sqrt(norm_squared(tempPoint))>)
	// 			{
	// 				if(sqrt(norm_squared(tempPoint))> 0)//inner radius set to 0 for now
	// 				{
	// 					if(tempPoint[0]<(tempPoint[1]/2) && tempPoint[0]>(-tempPoint[1]/2) && tempPoint[1]>0)
	// 					{
	// 						// if(norm(hole_center-tempPoint)>hole_radius)
	// 						// {
	// 							cells.push_back(tempPoint);	
	// 						// }
	// 						
	// 					}
	// 					/*
 //                        			std::ofstream ofs;
	// 					ofs.open ("temp_points_location.csv", std::ofstream::out | std::ofstream::app);
	// 					ofs <<(sqrt(norm_squared(tempPoint))-30)<<"," <<"\n";
	// 					ofs.close();
	// 					*/
	// 						
	// 				}
	// 			}
	// 		}
	//
	// 	}
	// }
	return cells;

}
std::vector<std::vector<double>> create_hepato_positions(double cell_radius,
                                                         double sphere_radius, 
                                                         std::vector<std::vector<double>> hole_centers, 
                                                         std::vector<double> hole_radii,
                                                         std::vector<Cell*> sinusoids
                                                         )
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	double z_bottom=-20;
	double z_top=20;
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);

	for(double z=z_bottom; z<z_top; z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				cells.push_back(tempPoint);		
							
	
			}

		}
	}
  for(int i=0; i<cells.size();i++ )
  {
      std::cout<<cells.size()<<"\n";
    for(int j=0; j<hole_centers.size(); j++)
    {
      while(norm(cells[i]-hole_centers[j])<=(hole_radii[j]))
      {
        cells.erase(cells.begin()+i);
      }
    }
  
  }
  for(int k=0; k<cells.size();k++ )
  {
    for(int m=0; m<sinusoids.size(); m++)
    {

      while(norm(cells[k]-sinusoids[m]->position)<=
        (sinusoids[m]->phenotype.geometry.radius+cell_radius-2))
      {

        cells.erase(cells.begin()+k);
      }
    }
  }
	return cells;

}



std::vector<std::vector<double>> create_cell_disk_wedge_positions(double cell_radius, double sphere_radius, double inner_radius, std::vector<double> hole_center, double hole_radius)//wedge between y=2x and y=-2x
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	double z_bottom=-20;
	double z_top=20;
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);

	for(double z=z_bottom; z<z_top; z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;

				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{
					if(sqrt(norm_squared(tempPoint))> 0)//inner radius set to 0 for now
					{
						if(tempPoint[0]<(tempPoint[1]/2) && tempPoint[0]>(-tempPoint[1]/2) && tempPoint[1]>0)
						{
							if(norm(hole_center-tempPoint)>hole_radius)
							{
								cells.push_back(tempPoint);	
							}
							
						}
						/*
                        			std::ofstream ofs;
						ofs.open ("temp_points_location.csv", std::ofstream::out | std::ofstream::app);
						ofs <<(sqrt(norm_squared(tempPoint))-30)<<"," <<"\n";
						ofs.close();
						*/
							
					}
				}
			}

		}
	}
	return cells;

}
void setup_wedge(void)
{
    double initial_granulosa_radius= 5;
    double initial_oocyte_radius=54;
		double initial_overlap=.1;//representation of the initial packing density from config eventually
		double cell_spacing = initial_granulosa_radius-initial_overlap;//slight overlap to represent cells up against each other better variable name
		
		double wedge_radius=250;
		double hole_radius=30;
    std::vector<double> hole={0.0,5.0,3.0};
		//create granulosa wedge
		std::cout << "\tPlacing granulosa cells ... " << "\n";
		std::vector<std::vector<double>> granulosa_positions = create_cell_disk_wedge_positions(cell_spacing,wedge_radius,initial_oocyte_radius,hole,hole_radius);//wedge hopefully
		std::cout << "creating " << granulosa_positions.size() << " wedge " << "\n";
		//ofs.open ("initial_granulosa_size.csv", std::ofstream::out | std::ofstream::trunc);
		for( int i=0; i < granulosa_positions.size(); i++ )
		{
			Cell* pCell_granulosa;//placed this inside of the loop since it makes more sense to me
			pCell_granulosa = create_cell(granulosa_cell); //
			pCell_granulosa->assign_position( granulosa_positions[i] );
			pCell_granulosa->set_total_volume(pCell_granulosa->custom_data["initial_volume"]);
		}
	return;
}

void setup_tissue( void )
{
	//colored_cell_list.resize(0);
  setup_liver();
  // setup_curve();
  //select tissue to set up:
	// setup_secondary_stage_follicle();
	//setup_wedge();
	//setup_toy_granulosa_model();
  // setup_just_oocyte();	
  // setup_4_granulosa_test_case();
	initialize_spring_cells();
	//setup_n_layer_follicle(7);

	return;
}


void setup_curve(){
  bez_curve_3 test;
  test.initalize(0);
  std::vector<double> p1={0,0};
  std::vector<double> p2={150,50};
  std::vector<double> p3={300,90};
  std::vector<double> p4={210,150};
  std::vector<double> p5={0,0,0.0};
  test.set_end_points(p1,p4);
  test.set_control_points(p2,p3);
  std::vector<std::vector<double>> test_points={};
  double t_ratio=0.05; 
  test.points_on_me(t_ratio,&test_points);
  std::vector<std::vector<double>> test_locations={};
  for (int j=0; j<test_points.size(); j++) {
    // std::cout<<test_points[j]<<"\n";
  }
  std::vector<double>t_points={0.01};
  int number_of_ts=50;
  for (int k=0; k<number_of_ts; k++) {
    t_points.push_back(t_points[k]+(1/double(number_of_ts)));
    std::cout<<t_points[k]<<"\n";
  }
  Cell* a;
	a = create_cell(granulosa_cell); //
  a->assign_position(p5);
  a->set_total_volume(a->custom_data["initial_volume"]);
  #pragma omp critical
  {
    cells_on_bez(t_points,&test,&test_locations,5);
    // cells_on_path(test_points,&test_locations,a);
  }
  for(int i=0;i<test_locations.size();i++)
  {
    std::cout<<"locations: "<< test_locations[0]<<"\n";
    Cell* gran;
    gran= create_cell(granulosa_cell);
    gran->assign_position(test_locations[i]);
    gran->set_total_volume(gran->custom_data["initial_volume"]);
  }
  return;
}





poly_corner_2D::poly_corner_2D(std::vector<double> position)
{
  std::vector<double> my_position;
  std::vector<double> connected_point_1;
  std::vector<double> connected_point_2;
  std::vector<double> segment_1; 
  std::vector<double> segment_2; 
  void calculate_connections(std::vector<std::vector<double>> &points);
  void calculate_segments();
  bool is_within(std::vector<double> polygon_center_point);

}

polygon::polygon(int number_of_sides){
return;
}
void polygon::populate(std::vector<poly_corner_2D*> corner_points)
{
  return;
}
void bez_curve_3::initalize(int number){
  index=number;
  end_points.resize(2,std::vector<double>(2,0.0));
  control_points.resize(2,std::vector<double>(2,0.0));
  ratios={1,1,1,1};
  return;
};
void bez_curve_3::set_end_points(std::vector<double> p1,std::vector<double> p2){
  this->end_points={p1,p2};  
  return;
}
void bez_curve_3::set_control_points(std::vector<double> p1,std::vector<double> p2){
  this->control_points={p1,p2}; 
  return;
}
void bez_curve_3::set_ratios(double r1,double r2, double r3, double r4){
  this->ratios={r1,r2,r3,r4}; 
  return;
}
void bez_curve_3::get_point(double t, std::vector<double> *return_point){
    double X=this->end_points[0][0]*(pow((1-t),3))*this->ratios[0]+
              3*this->control_points[0][0]*(pow((1-t),2))*t*this->ratios[1]+
              3*this->control_points[1][0]*(1-t)*(pow(t,2))*this->ratios[2]+
              this->end_points[1][0]*(pow(t,3))*this->ratios[3];
    double Y=this->end_points[0][1]*(pow((1-t),3))*this->ratios[0]+
              3*this->control_points[0][1]*(pow((1-t),2))*t*this->ratios[1]+
              3*this->control_points[1][1]*(1-t)*(pow(t,2))*this->ratios[2]+
              this->end_points[1][1]*(pow(t,3))*this->ratios[3];
    std::vector<double> point={X,Y};
    *return_point=point;
  return;
}
void bez_curve_3::points_on_curve(double &t_values_as_ratio,std::vector<std::vector<double>> &end_points, 
                     std::vector<std::vector<double>> &control_points, 
                     std::vector<double> &ratios, 
                     std::vector<std::vector<double>> *return_points){
 for(size_t i=0;i<100;i+=int(t_values_as_ratio*100))
  {
    double t=double(i)/100;
    double X=end_points[0][0]*(pow((1-t),3))*ratios[0]+3*control_points[0][0]*(pow((1-t),2))*t*ratios[1]
              +3*control_points[1][0]*(1-t)*(pow(t,2))*ratios[2]+end_points[1][0]*(pow(t,3))*ratios[3];
    double Y=end_points[0][1]*(pow((1-t),3))*ratios[0]+3*control_points[0][1]*(pow((1-t),2))*t*ratios[1]
              +3*control_points[1][1]*(1-t)*(pow(t,2))*ratios[2]+end_points[1][1]*(pow(t,3))*ratios[3];
    std::vector<double> point={X,Y};
    return_points->push_back(point);
  }
  return;
}
void bez_curve_3::points_on_me(double& t_values_as_ratio,std::vector<std::vector<double>> *return_points){
  this->points_on_curve(t_values_as_ratio,this->end_points, this->control_points, 
                     this->ratios, return_points);
  
  return;
}
void cells_on_bez(std::vector<double> t_values,bez_curve_3 *bez, std::vector<std::vector<double>>* cell_positions, double radius)
{
  double spacing=1;
  std::vector<double> xy_point={}; 
  std::vector<double> next_point={}; 
  double t=0.01;
  double increment=0.1;
  int count=0;
  while(t<1 && count <1000)
  {
    bez->get_point( t,&xy_point);
    bez->get_point((t+increment), &next_point);
    int count=0;
    if(norm(next_point-xy_point)>((radius)+spacing))
    {
      increment=increment-(increment/2); 
    }
    else if(norm(next_point-xy_point)<((radius-spacing)))
    {
      increment=increment+(increment/2);
    }
    else {
      cell_positions->push_back(xy_point);
      t=t+increment;
    }

    count++; 
  }
//approximate curvy path with line segments, if segment is smaller than a cell, place cells at end points
  std::cout<<"FINISHED WITH "<<bez->index<<"\n";
  return;
}

void cells_on_path(std::vector<std::vector<double>> path_points, std::vector<std::vector<double>>* cell_positions, Cell* pCell)
{
  std::vector<double> line_start(2,0.0); 
  std::vector<double> line_end(2,0.0); 
  for(size_t i=0; i<(path_points.size()-1);i+=2)
  {
    std::cout<<"TEST: "<<i<<"\n";
    line_start=path_points[i];
    line_end=path_points[i+1];
    std::cout<<i<<": "<<line_start<<line_end<<"\n";
    cell_positions->push_back(line_start);
    cell_positions->push_back(line_end);
    // cells_on_line(line_start,line_end,cell_positions,pCell);  
  }
//approximate curvy path with line segments, if segment is smaller than a cell, place cells at end points
  return;
}
void cells_on_line(std::vector<double>&line_start,std::vector<double>&line_end, std::vector<std::vector<double>>* cell_positions, Cell* pCell)//place cells along line segements
{
  double length=norm(line_end-line_start);
  if (length<pCell->phenotype.geometry.radius) {
    cell_positions->push_back(line_start);
    cell_positions->push_back(line_end);
    return;
  }
  else{
    double cutoff=(pCell->phenotype.geometry.radius-pCell->custom_data["allowed_overlap"]);
    double slope=(line_end[1]-line_start[1])/(line_end[0]-line_start[0]);
    double b= line_start[1]-(line_start[0]*slope);
    double spacing=pCell->phenotype.geometry.radius-pCell->custom_data["allowed_overlap"];
    double number_of_centers=((line_end[0]-line_start[0])/spacing)+1;
    double length_remainder= number_of_centers-std::floor(number_of_centers);
    int number_of_cells= int(number_of_centers);
    if(length_remainder>cutoff)
    {
      number_of_cells=int(number_of_centers)+1;
    }
    double new_x=line_start[0];
    double new_y=line_start[1];
    std::vector<double> point(2,0.0);
    for(int i=0; i<number_of_cells;i++){
      point={new_x,new_y};
      cell_positions->push_back(point);
      double dx=new_x;
      new_x+=(length/number_of_cells)*i;  
      new_y= new_x*slope+b;
    }
  
  }

  ///cell positions += with push_back
  
  return;
}
void subtract_overlapping_cells(){
  //can almost certainly use voxels to speed up// if center-to_center with other cell < (radius-overlap) remove it from the cell_list
  return;
}
void add_venule(std::vector<double> position, double radius){return;}
void add_portals(std::vector<std::vector<double>> positions, double radius){return;}
void add_sinusoids(double density){return;}
void setup_liver(){
  std::vector<Cell*> sinusoids={};
  std::vector<std::vector<double>> hole_centers={{0,0,0},{250,300,0},{-250,300,0},
    {300,0,0},{-300,0,0},{250,-300,0},{-250,-300,0}};
  std::vector<double>hole_radii={50,25,25,25,25,25,25};
  bez_curve_3 one;
  one.initalize(1);
  one.set_end_points({50,50}, {250,300});
  one.set_control_points({0,100},{25,0});
  bez_curve_3 two;
  two.set_end_points({-50,-50}, {-250,-300});
  two.set_control_points({0,-80},{-200,0});
  two.initalize(2);
  bez_curve_3 three;
  three.set_end_points({-50,50}, {-250,300});
  three.set_control_points({0,100},{-25,0});
  three.initalize(3);
  bez_curve_3 four;
  four.set_end_points({50,-50}, {250,-300});
  four.set_control_points({0,-90},{10,0});
  four.initalize(4);
  bez_curve_3 five; 
  five.set_end_points({50,0}, {300,0});
  five.set_control_points({45,45},{90,-20});
  five.initalize(5);
  bez_curve_3 six;
  six.set_end_points({-50,0}, {-300,0});
  six.set_control_points({-45,-45},{-90,20});
  six.initalize(6);
  std::vector<double>t_points={0.01};
  int number_of_ts=50;
  for (int k=0; k<number_of_ts-1; k++) {
    t_points.push_back(t_points[k]+(1/double(number_of_ts)));
    // std::cout<<t_points[k]<<"\n";
  }
  std::vector<std::vector<double>> one_locations={};
  std::vector<std::vector<double>> two_locations={};
  std::vector<std::vector<double>> three_locations={};
  std::vector<std::vector<double>> four_locations={};
  std::vector<std::vector<double>> five_locations={};
  std::vector<std::vector<double>> six_locations={};
  std::vector<std::vector<double>> seven_locations={};
  #pragma omp critical
  {
    cells_on_bez(t_points,&one,&one_locations,5);
    // cells_on_path(test_points,&test_locations,a);
  }
  for(int i=0;i<one_locations.size();i++)
  {
    Cell* gran;
    gran= create_cell(granulosa_cell);
    gran->assign_position(one_locations[i]);
    gran->set_total_volume(gran->custom_data["initial_volume"]);
  }
  #pragma omp critical
  {
    cells_on_bez(t_points,&two,&two_locations,5);
    // cells_on_path(test_points,&test_locations,a);
  }
  for(int i=0;i<two_locations.size();i++)
  {
    Cell* gran;
    gran= create_cell(granulosa_cell);
    gran->assign_position(two_locations[i]);
    gran->set_total_volume(gran->custom_data["initial_volume"]);
  }
  #pragma omp critical
  {
    cells_on_bez(t_points,&three,&three_locations,5);
    // cells_on_path(test_points,&test_locations,a);
  }
  for(int i=0;i<three_locations.size();i++)
  {
    Cell* gran;
    gran= create_cell(granulosa_cell);
    gran->assign_position(three_locations[i]);
    gran->set_total_volume(gran->custom_data["initial_volume"]);
  }
  #pragma omp critical
  {
    cells_on_bez(t_points,&four,&four_locations,5);
    // cells_on_path(test_points,&test_locations,a);
  }
  for(int i=0;i<four_locations.size();i++)
  {
    Cell* gran;
    gran= create_cell(granulosa_cell);
    gran->assign_position(four_locations[i]);
    gran->set_total_volume(gran->custom_data["initial_volume"]);
  }
  #pragma omp critical
  {
    cells_on_bez(t_points,&five,&five_locations,5);
    // cells_on_path(test_points,&test_locations,a);
  }
  for(int i=0;i<five_locations.size();i++)
  {
    Cell* gran;
    gran= create_cell(granulosa_cell);
    gran->assign_position(five_locations[i]);
    gran->set_total_volume(gran->custom_data["initial_volume"]);
  }
  #pragma omp critical
  {
    cells_on_bez(t_points,&six,&six_locations,5);
    // cells_on_path(test_points,&test_locations,a);
  }
  for(int i=0;i<six_locations.size();i++)
  {
    Cell* gran;
    gran= create_cell(granulosa_cell);
    gran->assign_position(six_locations[i]);
    gran->set_total_volume(gran->custom_data["initial_volume"]);
  }
  for(int i=0; i<(*all_cells).size();i++)
  {
    Cell* temp_cell=(*all_cells)[i];
    sinusoids.push_back(temp_cell);
  }
  std::vector<std::vector<double>> hepato=create_hepato_positions(14.5, 400, hole_centers, hole_radii,sinusoids);
  for( int j=0; j < hepato.size(); j++ )
		{
			Cell* pCell_hepato;//placed this inside of the loop since it makes more sense to me
			pCell_hepato = create_cell(hepato_cell); //
			pCell_hepato->assign_position( hepato[j] );
			pCell_hepato->set_total_volume(pCell_hepato->custom_data["initial_volume"]);
		}

  return;
}
