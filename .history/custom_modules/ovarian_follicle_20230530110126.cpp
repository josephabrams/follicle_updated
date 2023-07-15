
#include "./ovarian_follicle.h"
#include "./springs.h"
#include "./follicle_utilities.h"
//#include "./EG_simulation_parameters.cpp"
//#include "./PBS_simulation_parameters.cpp"
//#include "./GLY_simulation_parameters.cpp"
//#include "./EG_and_GLY_simulation_parameters.cpp"
#include "./test_parameters.h"

// Global Variables
#define PI 3.14159265
#define R 0.08205//granulosa (10^-3 J/mole*k)
using namespace Springs;
Cell_Definition oocyte_cell;
Cell_Definition granulosa_cell;
std::string output_filename;

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

	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );



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
		//cell_defaults.custom_data.add_variable("spring_constant","unitless",0.0);
        cell_defaults.custom_data.add_variable("neighborhood_radius","um",0);
		cell_defaults.custom_data.add_variable("initial_volume","unitless",0);
		cell_defaults.custom_data.add_variable("Lp","unitless",0);
		cell_defaults.custom_data.add_variable("Ps","unitless",0);



		cell_defaults.custom_data.add_variable("surfaceArea","unitless",0);
		cell_defaults.custom_data.add_variable("R","unitless",0);

		cell_defaults.custom_data.add_variable("Temperature","unitless",0);

		cell_defaults.custom_data.add_variable("isosmotic_volume","unitless",0);

		cell_defaults.custom_data.add_variable("Partial_molar_volume","unitless",0);

		cell_defaults.custom_data.add_variable("Vb_fraction","unitless",0);
		cell_defaults.custom_data.add_variable( "kill rate" , "1/min" , 0 ); // how often it tries to kill
		cell_defaults.custom_data.add_variable( "attachment lifetime" , "min" , 0 ); // how long it can stay attached
		cell_defaults.custom_data.add_variable( "attachment rate" , "1/min" ,0 ); // how long it wants to wander before attaching
		cell_defaults.custom_data.add_variable("osmotically_active_water_volume","unitless",0.0);
		cell_defaults.custom_data.add_variable("prev_osmotically_active_water_volume","unitless",0.0);
		cell_defaults.custom_data.add_variable("max interaction distance","unitless",0.5);
		//cell_defaults.custom_data.add_variable("total_internal_concentration","unitless",0.0);
		//cell_defaults.custom_data.add_variable("total_external_concentration", "unitless", 0.0);
		cell_defaults.custom_data.add_variable("dVw_dt", "unitless", 0.0);
		//cell_defaults.custom_data.add_variable("dS_dt_1", "unitless", 0.0);
		//cell_defaults.custom_data.add_variable("external_concentration_of_solute_2","unitless",0);
		//cell_defaults.custom_data.add_variable("Prev_total_volume","unitless",0);
		//ell_defaults.custom_data.add_variable("test_uptake","unitless",0);
		//cell_defaults.custom_data.add_variable("test_current_voxel","unitless",0);
	// create the cell types

	create_oocyte_cell_type();
	create_granulosa_cell_type();
	return;
}
//oocyte////////////////////
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

	// set functions
	oocyte_cell.functions.update_phenotype = oocyte_phenotype_rule;
	oocyte_cell.functions.custom_cell_rule = oocyte_cell_rule;

	oocyte_cell.functions.update_velocity = oocyte_velocity_rule;

	// set custom data values
	oocyte_cell.custom_data.add_variable("cell_k","unitless",-1.0);//oocyte spring constant
   	oocyte_cell.custom_data.add_variable("neighborhood_radius","um",-1.0);//radius beyond cell surface where connections occur
	oocyte_cell.custom_data[ "initial_volume" ] = 817283;//um^3
	oocyte_cell.custom_data["Lp"]=-1.0;//Lp_oocyte;  //0.0166666666;//0.0135;////um/atm/sec //EG 1/60=.016666667;https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4141568/;
	oocyte_cell.custom_data["Ps"]=-1.0;//Ps_oocyte; //0.0; //.12;//0.0;um/sec

	/*
		DMSO_Lp 1.68 Immature 1.01 Mature
		DMSO_Ps .15 immature .24 mature
		PROH_Lp .72          .86
		PROH_Ps .31          .56
	*/
	oocyte_cell.custom_data["surfaceArea"]=42273.3;//43411.1;
	oocyte_cell.custom_data["R"]= R;
	oocyte_cell.custom_data["Temperature"]= -1.0;
	oocyte_cell.custom_data["isosmotic_volume"]=817283;
	//oocyte_cell.custom_data["Partial_molar_volume"]=0.0557414;
	oocyte_cell.custom_data["Vb_fraction"]=.29;



	oocyte_cell.custom_data.add_variable("dVw_dt", "unitless", 0.0);
	oocyte_cell.custom_data.add_variable("prev_dVw_dt", "unitless", 0.0);


	oocyte_cell.custom_data.add_variable("Toxicity","unitless",0);

	oocyte_cell.custom_data.add_variable("Prev_total_volume","unitless",817283);
	oocyte_cell.custom_data.add_variable("test_uptake","unitless",0);
	oocyte_cell.custom_data.add_variable("test_current_voxel","unitless",0);

	oocyte_cell.custom_data.add_variable("osmotically_active_water_volume","unitless",580271);//.29*initial_volume
	oocyte_cell.custom_data.add_variable("prev_osmotically_active_water_volume","unitless",580271);
	oocyte_cell.custom_data.add_variable("external_GLY_concentration","mole/liter",0.0);
	oocyte_cell.custom_data.add_variable("external_EG_concentration","mole/liter",0.0);
	oocyte_cell.custom_data.add_variable("internal_GLY_concentration","mole/liter",0.0);
	oocyte_cell.custom_data.add_variable("internal_EG_concentration","mole/liter",0.0);
	oocyte_cell.custom_data.add_variable("spring_constant","unitless",50.0);
	oocyte_cell.custom_data.add_variable("max interaction distance","unitless",0.5);

	return;
}

void oocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout<<"number of springs: "<<spring_cell_by_pCell_index.size()<<"\n";
	for (size_t i=0; i<all_spring_cells.size();i++)
	{
		std::cout<<"all springs by index: "<< spring_cell_by_pCell_index[i]->m_my_pCell<<"\n";
	 	//std::cout<<"has neighbors "<< spring_cell_by_pCell_index[i]->m_springs.size()<<"\n";
		std::cout<<" IS basement? "<< spring_cell_by_pCell_index[i]->is_basement_connected<<"\n";
		std::cout<<" Position "<< spring_cell_by_pCell_index[i]->m_my_pCell->position<<"\n";
	 	// for (size_t j = 0; j <spring_cell_by_pCell_index[i]->m_springs.size() ; j++)
	 	// {
	 	// 	std::cout<<"at : "<< i<< ", "<< j<<std::endl;
	 	// 	std::cout<<"connected to : "<< i<< ", "<< j<< ", "<< spring_cell_by_pCell_index[i]->m_springs[j]->m_length<<"\n";
	 	// }
		
		
	 }

	// //non_connected_neighbor_pressure(pCell,dt);
	// double test= concentration_at_boundary(pCell,1);
	//std::cout<<"Concentration average: "<< test<< "\n";
	//std::cout<<" Voxel Volume?: "<< pCell->get_container()->underlying_mesh.dV<< "\n";
	rasterize_my_uptake(pCell,1);

	//std::vector<int> test=get_interior_voxels(pCell);

	//std::cout<< "Granulosa Velocity: "<< pCell->velocity<<std::endl;
	//std::cout<<"Oocyte velocity is "<< pCell->velocity<<std::endl;
	//two_parameter_single_step(pCell,phenotype,dt);
	return;

}
void oocyte_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{

	return;

}
void oocyte_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::vector <double> basement_membrane_center={0.0,0.0,0.0};
	basement_membrane_mechanics(spring_cell_by_pCell_index[pCell->index], 100, basement_membrane_center);
	//custom_add_potentials_for_pCells(pCell);
	// std::ofstream ofs;
	// ofs.open ("./output/oocyte_test_file2.csv", std::ofstream::out | std::ofstream::app);
	// ofs << PhysiCell_globals.current_time<<", "<< pCell<<", "<< pCell->position[0]<<", "<<pCell->position[1]<< ", "<< pCell->position[2]<<", "<<pCell->velocity <<", "<<omp_get_thread_num()<<"\n";
	// ofs.close();
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

	//Parameter<double> paramD;
	granulosa_cell.custom_data.add_variable("basement_k","unitless",-1.0);

    granulosa_cell.custom_data.add_variable("neighborhood_radius","unitless",-1.0);//radius beyond cell surface where connections occur
	granulosa_cell.custom_data[ "initial_volume" ] = 523.6;
	granulosa_cell.custom_data["Lp"]=0;
	granulosa_cell.custom_data["Ps"]=0;
	granulosa_cell.custom_data.add_variable("connection_length","micometers",0.5);

	granulosa_cell.custom_data["surfaceArea"]= 314.159;//parameters.doubles("granulosa_Area");
	granulosa_cell.custom_data["R"]= R;//parameters.doubles("granulosa_R");
	granulosa_cell.custom_data["Temperature"]=-1.0; //parameters.doubles("granulosa_Temperature");
	granulosa_cell.custom_data["isosmotic_volume"]=523.6;//parameters.doubles("isosmotic_granulosa_volume");
	//granulosa_cell.custom_data["Partial_molar_volume"]=0.0557414;//parameters.doubles("EG_partial_molar_volume");
	granulosa_cell.custom_data["Vb_fraction"]=.2;
	granulosa_cell.custom_data.add_variable("cell_k","unitless",-1.0);





	//	granulosa_cell.custom_data.add_variable("Water_volume", "unitless", 418.88);

	granulosa_cell.custom_data.add_variable("dVw_dt", "unitless", 0.0);
	granulosa_cell.custom_data.add_variable("prev_dVw_dt", "unitless", 0.0);


	granulosa_cell.custom_data.add_variable("Toxicity","unitless",0);
	granulosa_cell.custom_data.add_variable("external_concentration_of_solute_2","unitless",0);
	granulosa_cell.custom_data.add_variable("Prev_total_volume","unitless",523.6);
	granulosa_cell.custom_data.add_variable("test_uptake","unitless",0);
	granulosa_cell.custom_data.add_variable("test_current_voxel","unitless",0);

	granulosa_cell.custom_data.add_variable("osmotically_active_water_volume","unitless",418.88);//.29*initial_volume
	granulosa_cell.custom_data.add_variable("prev_osmotically_active_water_volume","unitless",418.88);
	granulosa_cell.custom_data.add_variable("external_GLY_concentration","mole/liter",0.0);
	granulosa_cell.custom_data.add_variable("external_EG_concentration","mole/liter",0.0);
	granulosa_cell.custom_data.add_variable("internal_GLY_concentration","mole/liter",0.0);
	granulosa_cell.custom_data.add_variable("internal_EG_concentration","mole/liter",0.0);
	granulosa_cell.custom_data.add_variable("spring_constant","unitless",0.5);
	granulosa_cell.custom_data.add_variable("max interaction distance","unitless",0.5);
	return;
}
void granulosa_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{

}
void granulosa_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{

	return;
}
int count=0;
void granulosa_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::vector <double> basement_membrane_center={0.0,0.0,0.0};
	//custom_add_potentials_for_pCells(pCell);
	basement_membrane_mechanics(spring_cell_by_pCell_index[pCell->index], 100, basement_membrane_center);
	//std::cout<<"velocity"<<pCell->velocity<<"\n";
	// std::ofstream ofs;
	// if(PhysiCell_globals.current_time<90.01 &&PhysiCell_globals.current_time>89.99 )
	// {
	// 	ofs.open ("./output/test_file3.csv", std::ofstream::out | std::ofstream::app);
	// 	ofs << PhysiCell_globals.current_time<<", "<< pCell<<", "<< pCell->position[0]<<", "<<pCell->position[1]<< ", "<< pCell->position[2]<<", "<<pCell->velocity <<", "<<omp_get_thread_num()<<"\n";
	// 	ofs.close();
	// }
	//std::cout<< "Granulosa Velocity: "<< pCell->velocity<<std::endl;
	return;
}
//! /////////////////////////////PhysiCell///////////////////////////////////////////////////////////////////
void setup_microenvironment( void )
{

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

	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding 2D setting to return to 3D" << "\n";
		default_microenvironment_options.simulate_2D = false;
	}

	// gradients are needed for this example

	//default_microenvironment_options.calculate_gradients = true;

	// add the immunostimulatory factor

	// let's do these in XML later
	microenvironment.add_density( "EG", "dimensionless" );
	microenvironment.diffusion_coefficients[1] = 530;
	microenvironment.decay_rates[1] = 0;


	microenvironment.add_density( "GLY", "dimensionless" );
	microenvironment.diffusion_coefficients[2] = 530;//530; //TODO: check about the diffusivity of this value
	microenvironment.decay_rates[2] = 0;

	microenvironment.add_density( "HM", "dimensionless" );
	microenvironment.diffusion_coefficients[3] = 510;
	microenvironment.decay_rates[3] = 0;
	
	microenvironment.add_density( "PBS", "dimensionless" );
	microenvironment.diffusion_coefficients[4] = 510;
	microenvironment.decay_rates[4] = 0;

	microenvironment.add_density( "SUC", "dimensionless" );
	microenvironment.diffusion_coefficients[5] = 1075;
	microenvironment.decay_rates[5] = 0;
	
	// do not let BioFVM use oxygen as the default
	default_microenvironment_options.use_oxygen_as_first_field = false;
	

	// set Dirichlet conditions
	//default_microenvironment_options.outer_Dirichlet_conditions = true;
	//default_microenvironment_options.Dirichlet_condition_vector[0] = 0;
	//default_microenvironment_options.Dirichlet_condition_vector[1] = 444.44;
	//default_microenvironment_options.Dirichlet_condition_vector[2] = 0;
	//default_microenvironment_options.Dirichlet_activation_vector[1]=false;
	//solute_loading(1.11,2.22,3.33,4.44,5.55,6.66);

	
	initialize_microenvironment();

	//total_uptake.resize(microenvironment.mesh.voxels.size());
	//all_total_uptakes.resize(number_of_permeating_solutes,total_uptake);
	return;
}
void setup_secondary_stage_follicle(void)
{
    /*

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
	return;
    */
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
		pCell_oocyte->assign_position(30,30,0); //assign oocyte position
		//pC->set_total_volume(817283); //sets initial total volume manually
		pCell_oocyte->set_total_volume(100000);//will set cell to custom value


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
	//pC->set_total_volume(817283); //sets initial total volume manually
	pCell_oocyte->set_total_volume(3*pCell_oocyte->custom_data["initial_volume"]);//will set cell to custom value
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

	initialize_spring_cells();
	//initialize_spring_connections();
	//Cell* pCell_to_spring_cell=(*all_cells)[0];
	//Spring_Cell oocyte(pCell_to_spring_cell);
	return;
}
void setup_tissue( void )
{

		

	
    
	//granulosa_k=gran_k;
	//std::cout<<"granulosa_k:"<<granulosa_k<<"\n";
	//set_concentrations( 0.2534, 3 );//set external concetratioin of HM to 0.255 
	//colored_cell_list.resize(0);
	//select tissue to set up:
	//setup_secondary_stage_follicle();
	//setup_wedge();
	//setup_toy_granulosa_model();
	//setup_just_oocyte();
	setup_4_granulosa_test_case();

	//setup_n_layer_follicle(7);
	//find and assign neighbors
    /*
	for (int i = 0; i < (*all_cells).size(); i++)// clear all the connections
	{

		Cell* pCell=(*all_cells)[i];
		populate_my_neighborhood(pCell);
		//std::cout<<"pCell: "<<pCell<<"\n";
		//for(int i=0; i<pCell->state.neighbors.size();i++)//go through all of pCells neighbors
		//{

			//std::cout<<"neighbor "<< i<<": "<<	pCell->state.neighbors[i]<< "\n";
			//std::cout<<"length "<< i<<": "<<	pCell->state.rest_length[i]<< "\n";
		//}
		populate_colored_cell_list(pCell);

	}
	//attach outter cells to basal lamina 
	//TODO: TEST OUTTER CELLS WITH THE COLORING FUNCTION

	for (int i = 0; i < (*all_cells).size(); i++)// clear all the connections
	{
		Cell* pCell=(*all_cells)[i];
		double distance_from_center=sqrt(pCell->position[0]*pCell->position[0]+pCell->position[1]*pCell->position[1]+pCell->position[2]*pCell->position[2]);
		double distance_to_basal_lamina=follicle_radius-initial_granulosa_radius*2;//go outside follicle then inward about 1 granulosa layer, attach these

		if(norm(pCell->position)>= distance_to_basal_lamina)
		{


			basal_lamina_list.push_back(pCell);//only for calculating basal distance
			//store attachment point to basement membrane and rest length
			std::vector <double> attachement_point_location=normalize(pCell->position);
			for( int j=0; j<attachement_point_location.size(); j++)
			{
				attachement_point_location[j]=attachement_point_location[j]*pCell->phenotype.geometry.radius;
			}
			//attachement_point*=follicle_radius;
			pCell->state.basement_attachment_point=pCell->position+attachement_point_location;
			pCell->state.basement_distance=0;//redundant remove
			pCell->state.is_basement=true;
			outtermost_granulosa=pCell;
		}
		
	}
	   started building hole mechanics
		for (int i = 0; i < (*all_cells).size(); i++)// clear all the connections
		{
			Cell* pCell=(*all_cells)[i];
			double distance_from_center=sqrt(pCell->position[0]*pCell->position[0]+pCell->position[1]*pCell->position[1]+pCell->position[2]*pCell->position[2]);
			double distance_to_basal_lamina=follicle_radius-initial_granulosa_radius*2;//go outside follicle then inward about 1 granulosa layer, attach these
			std::vector <double> vector_to_hole;
			for(int i=0; i<3;i++)
			{
				vector_to_hole[i]=hole[i]-pCell->position[i];
			}
			if(norm(vector_to_hole)<20)
			{


				basal_lamina_list.push_back(pCell);//only for calculating basal distance
				//store attachment point to basement membrane and rest length
				std::vector <double> attachement_point_location=normalize(vector_to_hole);
				for( int j=0; j<attachement_point_location.size(); j++)
				{
					attachement_point_location[j]=attachement_point_location[j]*pCell->phenotype.geometry.radius;
				}
				//attachement_point*=follicle_radius;
				pCell->state.basement_attachment_point=pCell->position+attachement_point_location;
				pCell->state.basement_distance=norm(vector_to_hole);
				pCell->state.is_basement=true;
				outtermost_granulosa=pCell;
			}
			
		}
	*/

	return;
}
