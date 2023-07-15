#include "./follicle_utilities.h"
#include "follicle_utilities.h"
#include <memory>
// #include "./springs.h"
// using namespace Springs;
//  #include <stdef.h>
//   Global Variables
#define PI 3.14159265
#define R 0.08205 // granulosa (10^-3 J/mole*k)
std::vector<double> displacement_between_membranes(Cell *pCell_1,
                                                    Cell *pCell_2) {
  // from pCell_1 to pCell_2 vector operator defined in bioFVM
  std::vector<double> displacement = pCell_2->position - pCell_1->position;
  displacement = (norm(pCell_2->position - pCell_1->position) -
                  pCell_1->phenotype.geometry.radius -
                  pCell_2->phenotype.geometry.radius) *
                 normalize(displacement);
  return displacement;
}
double distance_between_membranes(Cell *pCell_1, Cell *pCell_2) {
  double distance = std::abs(norm(pCell_2->position - pCell_1->position) -
                        pCell_1->phenotype.geometry.radius -
                        pCell_2->phenotype.geometry.radius);
  return distance;
}
// Initialize custom neighborhood based on distance between membranes
// loop through all neighboors
void initialize_neighboring_spring_connections() {
  //old version probably deprecate
  for (int i = 0; i < (*all_cells).size(); i++) {

    Cell *pCell = (*all_cells)[i];
    attach_neighboring_springs(pCell, 1.0);
  }
  return;
}
void attach_neighboring_springs(Cell *pCell, double max_spring_length) 
{
  // check nearby cells to identify neighboors using built in mechanics voxel
  // if less than max spring length store (pCell, rest length)
  std::vector<Cell *> neighbors = find_nearby_cells(pCell);
#pragma omp critical // make safe for push_back
  {
    for (int i = 0; i < neighbors.size(); i++) {
      Cell *pNeighbor = neighbors[i];
      double distance = distance_between_membranes(pCell, pNeighbor);
      if (distance < max_spring_length) {
        //testing spring class here
        //pCell->state.connection_and_spring_length.push_back(std::make_pair(pNeighbor, distance));
        pCell->state.neighbors.push_back(pNeighbor); // uses built in neighboor storage, default velocity
                        // function must be off or this will get overwritten
      }
    }
  }
  return;
}
std::vector <int> diffusion_bounding_box(Cell* pCell)
{
  //gets the box of voxels cell is contained in, includes the voxels inside
  std::vector<int> bounding_box_by_index={};
  //bounding box of voxels in the microenvironment
 //for speed get center voxel figure out x,y,z offset and only check local voxels
    Voxel center_voxel=microenvironment.nearest_voxel(pCell->position);
    std::vector<double> center_voxels_center= center_voxel.center; // get center voxel
    std::vector<double> offset=center_voxels_center-pCell->position;
    //figure out the corners of my bounding box of voxels to search
    std::vector<double> lower_point(3,0.0);
    std::vector<double> upper_point(3,0.0);
    for (size_t i = 0; i < 3; i++)
    {
      lower_point[i]= center_voxels_center[i]-offset[i]-ceil(pCell->phenotype.geometry.radius);
      upper_point[i]= center_voxels_center[i]+offset[i]+ceil(pCell->phenotype.geometry.radius);
    }
    Voxel starting_voxel=microenvironment.nearest_voxel(lower_point);
    Voxel ending_voxel=microenvironment.nearest_voxel(upper_point);
    bounding_box_by_index=general_voxel_bounding_box(starting_voxel.center,ending_voxel.center,default_microenvironment_options.dx,microenvironment.mesh);
  //std::cout<<"Bounding box"<<bounding_box_by_index<<"\n";
  return bounding_box_by_index;
}
std::vector<std::vector<double>> get_voxel_corners(std::vector<double> voxel_center)
{
  //find all 8 corners of the voxel cube
  double zz=default_microenvironment_options.dz/2.0;
  double yy=default_microenvironment_options.dy/2.0;
  double xx=default_microenvironment_options.dx/2.0;
  std::vector <std::vector <double>> corners(8,std::vector<double>(3,0.0));
  int count=0;
  for (int i = -1; i < 2; i+=2)
  {
    for (int j = -1; j < 2; j+=2)
    {
      for (int k = -1; k < 2; k+=2)
      {
        std::vector<double> temp_point={xx*i,yy*j,zz*k};
        corners[count]=voxel_center+temp_point;
        count++;
      }
    }   
  }
  return corners;
}
std::vector<int> get_intersecting_voxels(Cell* pCell)
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector <int> bounding_voxels=diffusion_bounding_box(pCell);
    std::vector<int> intersecting_voxels={};
    std::vector<int> return_voxel_index={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
    
     
      int sum=0;
      #pragma omp private(sum,exterior_voxels)
      for(size_t j =0; j<test_corners.size();j++)
      {
        //std::cout<<j<<" "<<test_corners[j]<<" distance: "<< norm(test_corners[j]-pCell->position)<<std::endl;
        if(norm(test_corners[j]-pCell->position)<pCell->phenotype.geometry.radius)
        {
          sum+=1;
        }
      }
      
      if(sum!=8 && sum !=0)
      {
        intersecting_voxels.push_back(bounding_voxels[i]);
      }

    }
    #pragma omp critical
    {
      return_voxel_index.insert(return_voxel_index.end(), intersecting_voxels.begin(), intersecting_voxels.end());
    }

    return return_voxel_index;
}
std::vector<int> get_exterior_voxels(Cell* pCell)
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector <int> bounding_voxels=diffusion_bounding_box(pCell);
    std::vector<int> exterior_voxels={};
    std::vector<int> return_voxel_index={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
    
     
      int sum=0;
      #pragma omp private(sum,exterior_voxels)
      for(size_t j =0; j<test_corners.size();j++)
      {
        //std::cout<<j<<" "<<test_corners[j]<<" distance: "<< norm(test_corners[j]-pCell->position)<<std::endl;
        if(norm(test_corners[j]-pCell->position)<pCell->phenotype.geometry.radius)
        {
          sum+=1;
        }
      }
      
      if(sum!=8 && sum !=0 && norm(test_voxel_center-pCell->position)>pCell->phenotype.geometry.radius)
      {
        exterior_voxels.push_back(bounding_voxels[i]);
      }

    }
    #pragma omp critical
    {
      return_voxel_index.insert(return_voxel_index.end(), exterior_voxels.begin(), exterior_voxels.end());
    }

    return return_voxel_index;
}
std::vector<int> get_interior_voxels(Cell* pCell)
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector <int> bounding_voxels=diffusion_bounding_box(pCell);
    std::vector<int> interior_voxels={};
    std::vector<int> return_voxel_index={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
    
     
      int sum=0;
      #pragma omp private(interior_voxels)

      if(norm(test_voxel_center-pCell->position)<=pCell->phenotype.geometry.radius)
      {
        interior_voxels.push_back(bounding_voxels[i]);
      }

    }
    #pragma omp critical
    {
      return_voxel_index.insert(return_voxel_index.end(), interior_voxels.begin(), interior_voxels.end());
    }

    return return_voxel_index;
}

void solute_loading( double oxygen, double final_solute_concentration_1, double final_solute_concentration_2,double final_solute_concentration_3,double final_solute_concentration_4, double final_solute_concentration_5)
{
	//Dirichlet node for all the voxels located outside of the ring
	std::vector<double> dirichlet_solute_end_state( 6 , 0.0 );
	dirichlet_solute_end_state[0]=oxygen;//oxygen
	dirichlet_solute_end_state[1]=final_solute_concentration_1;//eg
	dirichlet_solute_end_state[2]=final_solute_concentration_2;//gly
	dirichlet_solute_end_state[3]=final_solute_concentration_3;//hm
	dirichlet_solute_end_state[4]=final_solute_concentration_4;//pbs
	dirichlet_solute_end_state[5]=final_solute_concentration_5;//suc
	//dirichlet_solute[4]=0.0;//sucrose
	//Dirichlet nodes for all the voxels located inside of the ring
/*	std::vector<double> dirichlet_solute_start_state( 3 , 0 );
	dirichlet_solute_start_state[0]=0.0;//annoying oxygen no idea how to get rid of it but its turned off
	dirichlet_solute_start_state[1]=initial_solute_concentration_1;//final_external_osmolarity;//eg
	dirichlet_solute_start_state[2]=initial_solute_concentration_2;//gly
*/
/*
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		if(dist(microenvironment.voxels(i).center, {0.0,0.0,0.0})>=130)
		{
			microenvironment.update_dirichlet_node( i,dirichlet_solute_end_state );
		}
		else
		{
			microenvironment.remove_dirichlet_node(i);
		}
	}
*/
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
	{
		if(dist(microenvironment.voxels(i).center, {0.0,0.0,0.0})>=130)
		{
			microenvironment.update_dirichlet_node( i,dirichlet_solute_end_state );
		}
		else
		{
			microenvironment.remove_dirichlet_node(i);
		}
	}
	return;
}
double concentration_at_boundary(Cell* pCell, int solute_index)
{
  //previously used microenvironment.nearest_density_vector(boundary_point)[solute]
  //internalized substrates is the vector of substrates inside the cell (inside basic agent)
  double average=0;
  double sum=0;
  std::vector <int> exterior_voxel_index={};
  exterior_voxel_index=get_exterior_voxels(pCell);
  #pragma omp reduction(+:sum)
  for (size_t i = 0; i < exterior_voxel_index.size(); i++)
  {
     //std::cout<<" concentration "<<microenvironment.nearest_density_vector(exterior_voxel_index[i] )[solute_index]<<"\n";
      //std::cout<<" exterior voxel index"<< exterior_voxel_index[i]<<"\n";
      //std::cout<<" concentration for sum: "<<microenvironment.nearest_density_vector(exterior_voxel_index[i] )[solute_index]<<"\n";
      sum+=microenvironment.nearest_density_vector( exterior_voxel_index[i])[solute_index];
  }
  #pragma omp critical
  {
    average=sum/exterior_voxel_index.size();
  }
 
  return average; 
}



std::vector<Cell *> cells_in_me(Cell *pCell) // uses mechanics vectors to search for cells that are within or equal to pCell radius can also be used for bounding boxes
{
  //NOTE: current physicell uses the same cartesian mesh for mechanics and diffusion?? could use functions for diffusion bounding box
  std::vector<Cell *> agents_in_me;
  std::vector<int> search_voxels=diffusion_bounding_box(pCell);
  std::vector<Cell *> agents_in_voxel;
  #pragma omp private(agents_in_voxel)
  for (int i = 0; i < search_voxels.size();i++) 
  { 
      for (int j = 0; j < pCell->get_container()->agent_grid[search_voxels[i]].size(); j++) 
      {
        Cell* temp_agent=pCell->get_container()->agent_grid[search_voxels[i]][j];
        if ( temp_agent!= pCell && norm(temp_agent->position-pCell->position)<=pCell->phenotype.geometry.radius) 
        {
          agents_in_voxel.push_back(pCell->get_container()->agent_grid[search_voxels[i]][j]);
          // std::cout<< i << " Cell: "<<
          // agents_in_voxel[i]->position<<std::endl;
        }

      }

  }
  #pragma omp critical
  {
    agents_in_me.insert(agents_in_me.end(), agents_in_voxel.begin(), agents_in_voxel.end());
  }
      // std::cout<<"voxel " << temp_voxel.mesh_index<< " is inside my radius
      // and located at: "<< temp_voxel.center[0]<<",
      // "<<temp_voxel.center[1]<<", "<<temp_voxel.center[2]<<std::endl;
    
  // std::cout<<agents_in_me.size()<< " agents in me:
  // "<<agents_in_me[0]<<std::endl;
  return agents_in_me;
}
void find_basement_membrane_voxels(std::vector<double> center_of_sphere, std::vector<double> radius_of_sphere) {
  return;
}
std::vector<double> Hookes_law_force(std::vector<double> direction, double rest_length, double current_length, double spring_constant) 
{
  //TODO double check direction
  std::vector<double> force(3,0.0);
  force = (spring_constant * std::abs(rest_length - current_length)) *normalize(direction);//f=-kX
  //std::cout<<"Hookes law: "<< force<<"\n";
  return force;
}

void non_connected_neighbor_pressure(Cell* pCell, double dt) {
  // similar to the new function in PhysiCell 1.10 automatically pushes on cells
  // that press up against the cell; however this version doesnt require the
  // cell to be the same size as mechanics voxel uses same spring value as
  // connections, avoids double pushing on connected cells
  std::vector<Cell *> possible_neighbors = cells_in_me(pCell); // vector containing cells that could be interacting
                                                               // with pCell that aren't in neighboors
  std::vector<Cell *> non_connected_neighbors;
#pragma omp private(possible_neighbors)
  for (int j = 0; j < possible_neighbors.size(); j++) {
    for (int i = 0; i < spring_cell_by_pCell_index[pCell->index]->m_springs.size(); i++) {
      if (possible_neighbors[j] == spring_cell_by_pCell_index[pCell->index]->m_springs[i]->m_pNeighbor) {
        possible_neighbors.erase(possible_neighbors.begin() + j);
      }
    }
  }
#pragma omp critical
  {
    non_connected_neighbors.insert(non_connected_neighbors.end(),
                                   possible_neighbors.begin(),
                                   possible_neighbors.end());
  }
  // std::cout<<"non connected neighbors list"<<
  // non_connected_neighbors.size()<< std::endl;
  double sum_x_velocity = pCell->velocity[0];
  double sum_y_velocity = pCell->velocity[1];
  double sum_z_velocity = pCell->velocity[2];
  Cell *neighbor={};
#pragma omp private(neighbor) for reduction(+:sum_x_velocity,sum_y_velocity,sum_z_velocity)//firstprivate(sum_x_velocity, sum_y_velocity, sum_z_velocity) for//not sure if omp reduction will work with operator overloading so I did it this way
  for (int i = 0; i < non_connected_neighbors.size(); i++) {
    // std::cout<<" I is "<< i<<std::endl;
    neighbor = non_connected_neighbors[i];
    // std::cout<<neighbor<<std::endl;
    std::vector<double> force_on_neighbor_direction =
        neighbor->position - pCell->position;
    // if statement assumes strict interpentration of neighboors membrane into
    // pCells membrane to exert spring pressure
    if ((norm(force_on_neighbor_direction) -
         neighbor->phenotype.geometry.radius) <
        pCell->phenotype.geometry.radius) {
      double spring_stretch = distance_between_membranes(pCell, neighbor);
      std::vector<double> force_on_neighbor =
          Hookes_law_force(force_on_neighbor_direction, 0, spring_stretch,
                           pCell->custom_data["spring_constant"]);
      double neighbor_mass = neighbor->phenotype.volume.total; // mass_of_cell(neighboor);
              //TODO make mass of cell function
      double pCell_mass = pCell->phenotype.volume.total;
      //#pragma omp private(neighbor) for
      for (int coord=0; coord<3; coord++)
      {
              // forward euler velocity for now v=v+dt*dv/vt assume constant -
              // -acceleration and mass over dt
        neighbor->velocity[coord] += (abs(force_on_neighbor[coord] / neighbor_mass * dt) < 1e-12) ? 0.0 : (force_on_neighbor[coord] / neighbor_mass *dt);
      }


      // neighbor->velocity[0] +=
      //     (abs(force_on_neighbor[0] / neighbor_mass * dt) < 1e-16)
      //         ? 0.0
      //         : (force_on_neighbor[0] / neighbor_mass * dt);
      // neighbor->velocity[1] +=
      //     (abs(force_on_neighbor[1] / neighbor_mass * dt) < 1e-16)
      //         ? 0.0
      //         : (force_on_neighbor[1] / neighbor_mass * dt);
      // neighbor->velocity[2] +=
      //     (abs(force_on_neighbor[2] / neighbor_mass * dt) < 1e-16)
      //         ? 0.0
      //         : (force_on_neighbor[2] / neighbor_mass * dt);
      sum_x_velocity += (std::abs(force_on_neighbor[0] / pCell_mass * dt) < 1e-12)
                            ? 0.0
                            : (force_on_neighbor[0] / pCell_mass * dt);
      sum_y_velocity += (std::abs(force_on_neighbor[1] / pCell_mass * dt) < 1e-12)
                            ? 0.0
                            : (force_on_neighbor[1] / pCell_mass * dt);
      sum_z_velocity += (std::abs(force_on_neighbor[2] / pCell_mass * dt) < 1e-12)
                            ? 0.0
                            : (force_on_neighbor[2] / pCell_mass * dt);
    }
    // TODO: possibly, eventually use Adams_Bashforth_ODE_2nd_Order for velocity
  }
#pragma omp critical
  {
    pCell->velocity[0] = (sum_x_velocity < 1e-12) ? 0.0 : sum_x_velocity;
    pCell->velocity[1] = (sum_y_velocity < 1e-12) ? 0.0 : sum_y_velocity;
    pCell->velocity[2] = (sum_z_velocity < 1e-12) ? 0.0 : sum_z_velocity;
  }
}
// connection break- sever spring connection but maintain neighboor list to
// check for membrane overlap two lists, connected_cells[] for cells connected
// and neighboorhood, use default for granulosa and custom for oocyte
// oocyte neighboorhood - get voxel neighboorhood around oocyte
// non-spring mechanics (membrane overlap)

// using custom parameter or vector define initial spring length- function to
// calculate change in velocity Adams-bashforth second order ODE solver 2P Molar
// Concentration linked ODEs - requires storing previous molar concentrations
// and current for solutes and water voxelized secretion voxelized uptake
double Adams_Bashforth_ODE_2nd_Order(double y_value, double prev_df_dt, double df_dt, double step_size) 
{
  double ynew = y_value + step_size / 2 * (3 * (df_dt - prev_df_dt));
  return ynew;
}
double Forward_Euler(double f_value, double step_size, double df_dt) 
{
  double fnew = f_value + step_size * df_dt;
  return fnew;
}
std::vector<double> Adams_Bashforth_ODE_2nd_Order(std::vector<double> Y_values, std::vector<double> prev_df_dts, std::vector<double> df_dts, double step_size) {
  // vector form of Y_{n+1}=Y_{n}+h/2(3*(F(x_{n},t_{n})-F(x_{n-1}, t_{n-})))
std::vector<double> Ynext;
#pragma omp critical // lock thread for safety if running in parrallel requires
                     // include omp.h
  { Ynext.resize(Y_values.size()); }
  for (unsigned int i = 0; i < Y_values.size(); i++) {
    Ynext[i] = Y_values[i] + step_size / 2 * (3 * (df_dts[i] - prev_df_dts[i]));
  }
  return Ynext;
}
std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius,double inner_radius) 
{
  std::vector<std::vector<double>> cells;
  int xc = 0, yc = 0, zc = 0;
  double x_spacing = cell_radius * sqrt(3);
  double y_spacing = cell_radius * 2;
  double z_spacing = cell_radius * sqrt(3);
  std::vector<double> tempPoint(3, 0.0);

  for (double z = -sphere_radius; z < sphere_radius; z += z_spacing, zc++) {
    for (double x = -sphere_radius; x < sphere_radius; x += x_spacing, xc++) {
      for (double y = -sphere_radius; y < sphere_radius; y += y_spacing, yc++) {
        tempPoint[0] = x + (zc % 2) * 0.5 * cell_radius;
        tempPoint[1] = y + (xc % 2) * cell_radius;
        tempPoint[2] = z;

        if (sqrt(norm_squared(tempPoint)) < sphere_radius) {
          if (sqrt(norm_squared(tempPoint)) > inner_radius) {
            /*
std::ofstream ofs;
            ofs.open ("temp_points_location.csv", std::ofstream::out |
std::ofstream::app); ofs <<(sqrt(norm_squared(tempPoint))-30)<<"," <<"\n";
            ofs.close();
            */
            cells.push_back(tempPoint);
          }
        }
      }
    }
  }
  return cells;
}

void two_parameter_single_step(Cell *pCell, Phenotype &phenotype, double dt) // calculates current volume based on
                                          // previous volume and parameters
{
  // TODO update this function to only use custom variables 2023-05-02
  // could probably be improved in the future with some error checking
  // numberofsolutes is the number of solutes being modeled
  // Ms stores osmolality of all external permeating solutes
  // Mn stores osmolality of all internal non-permeating solutes
  // delta stores the deltas for converting femptomoles/liter to osmolality
  // s stores the femptomoles of internal solutes
  // Vw is the osmotically active volume of the cell
  // osmotically_active_water_volume[0] previous Vw,
  // osmotically_active_water_volume[1] current Vw
  //
  // double
  // vb_volume=pCell->custom_data["initial_volume"]*pCell->custom_data["Vb_fraction"];
  double s_volume = 0;  // will be the total volume taken by the moles of solute
  double total_Mse = 0; // total external permeating solute mole/L
  double total_Msi = 0; // total internal permeating solute mole/L
  double total_Mne = 0; // total external non-permeating solute mole/L
  double total_Mni = 0; // total inernatl non-permeating solute mole/L
  double xwe = 0;       // external water mole fraction
  double xwi = 0;       // intenral water mole fraction
  double cwe = 0;
  double cwi = 0;

  double V_wbar = 18.015 / 1000; // ml/mol->l/mol partial molar water volume
  double V_nbar = 27.1 / 1000;   // ml/mol partial molar salt volume
  double V_egbar =
      55.766 / 1000; // ml/mol partial molar eg volume not currently used
  double V_glybar = 0.0730903; // ml/mol
  double cpa_volume_ext =
      pCell->custom_data["external_EG_concentration"] * V_egbar +
      pCell->custom_data["external_GLY_concentration"] *
          V_glybar; // moles/liter*liter/mole use
  double cpa_volume_int =
      pCell->custom_data["internal_EG_concentration"] * V_egbar +
      pCell->custom_data["internal_GLY_concentration"] *
          V_glybar; // moles/liter*liter/mole
  double external_salt_volume =
      pCell->custom_data["external_holding_media_concentration"] * V_nbar +
      pCell->custom_data["external_PBS_concentration"] *
          V_nbar; // HM +PBS assumes both are NACL
  double internal_salt_volume =
      pCell->custom_data["internal_holding_media_concentration"] * V_nbar +
      pCell->custom_data["internal_PBS_concentration"] * V_nbar;
  // total_Mse=sampled_molarity(pCell);//sum the mole/kg of permeating external
  // solute
  // total_Msi=pCell->custom_data["prev_total_permeating_solute"]/pCell->custom_data["prev_osmotically_active_water_volume"]);

  //	cwe=1/V_wbar*(1-external_salt_volume-cpa_volume_ext);//mole/ml
  //	xwe=cwe/(cwe+total_Mse+total_Mne); //xwe=cwe/(cwe+cce+ccn)
  // std::cout<< "xwe: "<< xwe << "\n";
  //	cwi=1/V_wbar*(1-internal_salt_volume-cpa_volume_int);//mole/ml
  //	xwi=cwi/(cwi+total_Msi+total_Mni);//xwi=cwi/(cwi+cci+ccn)
  // external solute concentration
  // internal solute concentration
  //	pCell->custom_data["dVw_dt"]=pCell->custom_data["Lp"]*pCell->custom_data["surfaceArea"]*pCell->custom_data["R"]*pCell->custom_data["Temperature"]/V_wbar*log(xwe/xwi);
  //  pCell->custom_data["EG_dS_dt"]=pCell->custom_data["EG_Ps"]*pCell->custom_data["surfaceArea"]*(pCell->custom_data["Mse"]-pCell->custom_data["Msi"]);
  /*


for(int j=0;  j<number_of_NP_solutes; j++)// nonpermeating solutes
{
  total_Mne+=pCell->phenotype.secretion.Mne[j];//sum the mole/kg of
non-permeating external solute
          pCell->phenotype.secretion.Mni[j]=(pCell->phenotype.secretion.ni[j]/pCell->custom_data["prev_osmotically_active_water_volume"]);//moles/volume
  total_Mni+=pCell->phenotype.secretion.Mni[j];

}


//dVw is change in osmotically_active_water_volume
//dVw_dt=-Lp*surfaceArea*gasConst*temperature*(Mse+Mne-Msi-Mni);
  //std::cout<<(-1)<<", "<<pCell->custom_data["Lp"]<<",
"<<pCell->custom_data["surfaceArea"]<<", "<<pCell->custom_data["R"]<<",
"<<pCell->custom_data["Temperature"]<<", "<<(total_Mse+total_Mne)<<",
"<<(total_Msi+total_Mni)<<"\n";

  //change above to Molarity
  //std::cout<<pCell->phenotype.secretion.Ps[1]<<",
"<<pCell->custom_data["surfaceArea"]<<", "<<
pCell->phenotype.secretion.Mse[1]<<",
"<<pCell->phenotype.secretion.Msi[1]<<"\n";
  //pCell->phenotype.secretion.dS_dt[1]=pCell->phenotype.secretion.Ps[1]*pCell->custom_data["surfaceArea"]*(pCell->phenotype.secretion.Mse[1]-pCell->phenotype.secretion.Msi[1]);
  //std::cout<<"dVw_dt: " <<pCell->custom_data["dVw_dt"]<<"\n";
  for(int k=0; k<number_of_permeating_solutes;k++)//handles the mole flux for
each permeating solute, for non-permeating ps=0
{
  pCell->phenotype.secretion.dS_dt[k]=pCell->phenotype.secretion.Ps[k]*pCell->custom_data["surfaceArea"]*(pCell->phenotype.secretion.Mse[k]-pCell->phenotype.secretion.Msi[k]);

  }

  if(PhysiCell_globals.current_time<dt)//forward euler
      {
          pCell->custom_data["osmotically_active_water_volume"]=pCell->custom_data["prev_osmotically_active_water_volume"]+(pCell->custom_data["dVw_dt"]*dt);
          for(int l=0; l<number_of_permeating_solutes; l++)
          {
              pCell->phenotype.secretion.s[l]=pCell->phenotype.secretion.prev_s[l]+(pCell->phenotype.secretion.dS_dt[l]*dt);//femptomoles
of solute

          }
      }


  if(PhysiCell_globals.current_time>=dt)//adams-bashforth
      {
          //calculate Vw with second order adams bashforth
          pCell->custom_data["osmotically_active_water_volume"]=pCell->custom_data["prev_osmotically_active_water_volume"]+(dt/2)*(3*pCell->custom_data["dVw_dt"]-pCell->custom_data["prev_dVw_dt"]);
          //calculate s[i] with second order adams bashforth
          for(int m=0;  m<number_of_permeating_solutes;m++)
          {
              pCell->phenotype.secretion.s[m]=pCell->phenotype.secretion.prev_s[m]+(dt/2)*((3*pCell->phenotype.secretion.dS_dt[m])-pCell->phenotype.secretion.prev_dS_dt[m]);

          }
      }

//add all the internal permeating solute volumes together
for(int n=0; n<number_of_permeating_solutes;n++)
{
  s_volume+=pCell->phenotype.secretion.s[n]*pCell->phenotype.secretion.partial_molar_volumes[n];
//femptomoles*l/femptomole
          //std::cout<<"solute volume"<<s_volume<<"\n";
}
//set Volume= vb_volume
(liters)+osmotically_active_water_volume(liters)+s[i]*Vbar[i](femptomoles*liters/femptomole=liters)
double
total_volume=vb_volume+pCell->custom_data["osmotically_active_water_volume"]+s_volume;

  pCell->set_total_volume(total_volume);
  //pCell->set_total_volume(pCell->custom_data["Water_volume"]+(Partial_molar_volume*pCell->custom_data["Solute"])+Vb_volume);

//everything is calculated for this step set previous values for next step
pCell->custom_data["prev_dVw_dt"]=pCell->custom_data["dVw_dt"];
pCell->custom_data["prev_osmotically_active_water_volume"]=pCell->custom_data["osmotically_active_water_volume"];
for(int o=0;o<number_of_permeating_solutes;  o++)
{
          pCell->phenotype.secretion.flux[o]=pCell->phenotype.secretion.s[o]-pCell->phenotype.secretion.prev_s[o];//uptake
moles estimate for passing change in density back to bioFVM
  pCell->phenotype.secretion.prev_dS_dt[o]=pCell->phenotype.secretion.dS_dt[o];
  pCell->phenotype.secretion.prev_s[o]=pCell->phenotype.secretion.s[o];
}
*/
  return;
}
void rasterize_my_uptake(Cell* pCell, double solute_index)
{
  double total_voxel_volume=0.0;
  double uptake_per_voxel=0.0;
  //test double
  double custom_cell_uptake=2.22;//kg/um^3 dt already factored in at calculate 2p//TODO add this to custom_data this should be the value adjusted by Compute_2P
  std::vector<int> interior_voxels=get_interior_voxels(pCell);
  //get my interior voxels
  if(interior_voxels.size()==1)
  {
    //smaller than voxel
    total_voxel_volume=pCell->get_container()->underlying_mesh.dV;
    uptake_per_voxel=custom_cell_uptake*total_voxel_volume/(pCell->phenotype.volume.total);
  }
  else
  {
    //if multivoxel divide 
    total_voxel_volume=interior_voxels.size()*pCell->get_container()->underlying_mesh.dV;
    uptake_per_voxel=custom_cell_uptake*total_voxel_volume/(pCell->phenotype.volume.total);
  }
  //for my voxels uptake (secretion=-uptake)
  #pragma omp critical
  {
    for (size_t i = 0; i < interior_voxels.size(); i++)
    {
      microenvironment.density_vector(interior_voxels[i])[solute_index]+=uptake_per_voxel;
      
    }
  }
  return;
  
}
//as of May 23 2023 mechanics voxels seem to equal microenvironment voxels
std::vector <int> basement_membrane_voxels{};//external list for checking if a cell is intersecting the basement membrane
std::vector<int> get_basement_membrane_intersection(std::vector <double> center_point, double radius)
{
  double voxel_size=default_microenvironment_options.dx;
  //double inner_radius=radius-voxel_size;
  double outer_radius=radius+(voxel_size);//~2 voxels thick
     //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector <int> bounding_voxels=spherical_bounding_box(center_point,outer_radius);
    std::vector<int> intersecting_voxels={};
    std::vector<int> return_voxel_index={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=microenvironment.voxels(bounding_voxels[i]).center;
      std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
    
     
      int sum=0;
      #pragma omp private(sum,exterior_voxels)
      for(size_t j =0; j<test_corners.size();j++)
      {
        //std::cout<<j<<" "<<test_corners[j]<<" distance: "<< norm(test_corners[j]-pCell->position)<<std::endl;
        if(norm(test_corners[j]-center_point)<outer_radius && norm(test_corners[j]-center_point)>radius)//include catch zone of thickness
        {
          sum+=1;
        }
      }
      
      //if(sum!=8 && sum !=0)
      if(sum!=0)
      {
        intersecting_voxels.push_back(bounding_voxels[i]);
      }

    }
    #pragma omp critical
    {
      return_voxel_index.insert(return_voxel_index.end(), intersecting_voxels.begin(), intersecting_voxels.end());
    }

    return return_voxel_index;
}

void initialize_basement_membrane(std::vector <double> center_point, double radius)
{
  //calculate voxels that represent the BM
  basement_membrane_voxels=get_basement_membrane_intersection(center_point,radius);
  return;
}
std::vector <int> variable_moore_neighborhood(Cell* pCell, double radius)
{
  //works just like bounding diffusion box to get mechnical voxels around pCell with any radius
  //returns voxel by index
  //NOTE!!! currently PhysiCells microenvironment voxels are used for both mechanics and diffusion
  if(radius<=pCell->phenotype.geometry.radius)
  {
    std::cout<< "WARNING!! Are you sure you want to use a moore neighborhood smaller than the cell?"<<"\n";
  }

  std::vector<int> bounding_box_by_index={};
  //bounding box of voxels in the microenvironment
 //for speed get center voxel figure out x,y,z offset and only check local voxels
    Voxel center_voxel=microenvironment.nearest_voxel(pCell->position);
    std::vector<double> center_voxels_center= center_voxel.center; // get center voxel
    std::vector<double> offset=center_voxels_center-pCell->position;
    //figure out the corners of my bounding box of voxels to search
    std::vector<double> lower_point(3,0.0);
    std::vector<double> upper_point(3,0.0);
    for (size_t i = 0; i < 3; i++)
    {
      lower_point[i]= center_voxels_center[i]-offset[i]-ceil(radius);
      upper_point[i]= center_voxels_center[i]+offset[i]+ceil(radius);
    }
    Voxel starting_voxel=microenvironment.nearest_voxel(lower_point);
    Voxel ending_voxel=microenvironment.nearest_voxel(upper_point);
    if(starting_voxel.mesh_index==ending_voxel.mesh_index)
    {
      //contained entirely in 1 voxel my exterior is my voxel
      bounding_box_by_index.push_back(starting_voxel.mesh_index);
      return bounding_box_by_index;
    }
    //std::vector<int> box_dimensions(3,0.0); probably not needed recasting stuff as int, I don't like it but hopefully safe
   int num_voxels_across_box= int((ceil(ending_voxel.center[0])-floor(starting_voxel.center[0]))/default_microenvironment_options.dx); 
   int resizer=num_voxels_across_box*num_voxels_across_box*num_voxels_across_box;

    //loop through voxel centers, get index and store it
    if (floor(starting_voxel.center[2])==0 && ceil(ending_voxel.center[2])==0)
    {
      #pragma omp critical
      { //to avoid race conditions on resize if multiple cells are resized, unsure if this is an issue
      bounding_box_by_index.resize(num_voxels_across_box*num_voxels_across_box);
      }
      int count=0;
      for (int yi = int(floor(starting_voxel.center[1])); yi <int(ceil(ending_voxel.center[1])) ; yi=yi+int(default_microenvironment_options.dy))
      {

        for (int xi = int(floor(starting_voxel.center[0])); xi <int(ceil(ending_voxel.center[0])) ; xi+=int(default_microenvironment_options.dx))
        {
 
          std::vector<double> position={double(xi),double(yi),0.0};
          bounding_box_by_index[count]=microenvironment.nearest_voxel_index(position);
          count++;
        }
      }
    
    }
    else
    {
    #pragma omp critical
    { //to avoid race conditions on resize if multiple cells are resized, unsure if this is an issue
      bounding_box_by_index.resize(num_voxels_across_box*num_voxels_across_box*num_voxels_across_box);
    }
      int count=0;
      for (int zi = int(floor(starting_voxel.center[2])); zi <int(ceil(ending_voxel.center[2]))  ; zi+=default_microenvironment_options.dz)
      {
        for (int yi = int(floor(starting_voxel.center[1])); yi <int(ceil(ending_voxel.center[1])) ; yi+=default_microenvironment_options.dy)
        {
          for (int xi = int(floor(starting_voxel.center[0])); xi <int(ceil(ending_voxel.center[0])) ; xi+=default_microenvironment_options.dx)
          {
            std::vector<double> position={double(xi),double(yi),double(zi)};
            //Voxel bounding_voxel=microenvironment.nearest_voxel(position);
            bounding_box_by_index[count]=microenvironment.nearest_voxel_index(position);
            count++;
          }
        }
      }
    }
  //std::cout<<"Bounding box"<<bounding_box_by_index<<"\n";
  return bounding_box_by_index;
}
std::vector <int> max_interaction_variable_moore_neighborhood(Cell* pCell, double maximum_interaction_distance)
{
  double radius=pCell->phenotype.geometry.radius+maximum_interaction_distance;
  //works just like bounding diffusion box to get mechnical voxels around pCell with any radius
  //returns voxel by index
  //NOTE!!! currently PhysiCells microenvironment voxels are used for both mechanics and diffusion
  if(radius<=pCell->phenotype.geometry.radius)
  {
    std::cout<< "WARNING!! Are you sure you want to use a moore neighborhood smaller than the cell?"<<"\n";
  }

  
  //bounding box of voxels in the microenvironment
 //for speed get center voxel figure out x,y,z offset and only check local voxels
    Voxel center_voxel=microenvironment.nearest_voxel(pCell->position);
    std::vector<double> center_voxels_center= center_voxel.center; // get center voxel
    std::vector<double> offset=center_voxels_center-pCell->position;
    //figure out the corners of my bounding box of voxels to search
    std::vector<double> lower_point(3,0.0);
    std::vector<double> upper_point(3,0.0);
    for (size_t i = 0; i < 3; i++)
    {
      lower_point[i]= center_voxels_center[i]-offset[i]-ceil(radius);
      upper_point[i]= center_voxels_center[i]+offset[i]+ceil(radius);
 
    }

  std::vector<int> bounding_box_by_index=general_voxel_bounding_box(lower_point,upper_point,default_microenvironment_options.dx,(pCell->get_container()->underlying_mesh));
  //std::cout<<"Bounding box"<<bounding_box_by_index<<"\n";
  return bounding_box_by_index;
}
std::vector <int> general_voxel_bounding_box(std::vector <double> starting_position, std::vector <double> ending_position, double voxel_length, BioFVM::Cartesian_Mesh a_mesh)
{
  //function to find a box of voxels in a cartesian mesh given the lowest corner to the highest, should work in 2D and 3D
  std::vector<int> bounding_box_by_index={-1};
  Voxel starting_voxel=a_mesh.nearest_voxel(starting_position);// should work with pCell-get_container or microenvironment
  Voxel ending_voxel=a_mesh.nearest_voxel(ending_position);
  if(starting_voxel.mesh_index==ending_voxel.mesh_index)
  {
    //contained entirely in 1 voxel my exterior is my voxel
    bounding_box_by_index[0]=starting_voxel.mesh_index;
    return bounding_box_by_index;
  }

   std::vector <int> num_voxels_across(3,0.0);
   //loop through voxel centers, get index and store it
   std::vector <int> dimensions={0,0,0};
   for (size_t i = 0; i < 3; i++)
   {
    num_voxels_across[i]=int((ceil(ending_voxel.center[i])-floor(starting_voxel.center[i]))/voxel_length); 
    if(floor(starting_voxel.center[i])==ceil(ending_voxel.center[i]))
    {
      dimensions[i]=1;//if a dimension is flat mark it with 1
    }
   }
   
    int upper_z=0;
    int upper_y=0;
    int upper_x=0;
    #pragma omp critical
    { //to avoid race conditions on resize if multiple cells are resized, unsure if this is an issue
     
      upper_z=num_voxels_across[2]+dimensions[2];
      upper_y=num_voxels_across[1]+dimensions[1];
      upper_x=num_voxels_across[0]+dimensions[0];
      bounding_box_by_index.resize(upper_x*upper_y*upper_z);
    }
    std::vector<double> position(3,0.0);
    int count=0;
    for(int zi=0; zi<upper_z; zi++)
    { 
      position[2]=double(starting_voxel.center[2]+(double(zi-dimensions[2])*voxel_length));
      for (int yi = 0; yi <upper_y ; yi++)
      {
        position[1]=double(starting_voxel.center[1]+(double(yi-dimensions[1])*voxel_length));
        for (int xi = 0; xi <upper_x ; xi++)
        {
            position[0]=double(starting_voxel.center[0]+(double(xi-dimensions[0])*voxel_length));
            bounding_box_by_index[count]=microenvironment.nearest_voxel_index(position);
            count++;
        }
      }
    }
 
  return bounding_box_by_index;
}
void connect_spring_cells(Spring_Cell* SpCell_1, Spring_Cell* SpCell_2)
{
  double spring_length=distance_between_membranes(SpCell_1->m_my_pCell,SpCell_2->m_my_pCell);
  SpCell_1->add_spring(SpCell_2,spring_length);
  SpCell_2->add_spring(SpCell_1,spring_length);
  return;
}
std::vector <Cell*> cells_in_neighborhood(Cell* pCell, double maximum_interaction_distance)
{
  std::vector <int> voxels_to_search_indicies=max_interaction_variable_moore_neighborhood(pCell,maximum_interaction_distance);
  std::vector <Cell*> cells_found_private{};
  std::vector <Cell*> cells_found{};
  #pragma omp private(cells_found_private)
  for (size_t i = 0; i < voxels_to_search_indicies.size(); i++)
  {
    for (size_t j = 0; j < pCell->get_container()->agent_grid[voxels_to_search_indicies[j]].size(); j++) 
    {
      if(pCell->get_container()->agent_grid[voxels_to_search_indicies[i]][j]!=pCell)
      {
        cells_found_private.push_back(pCell->get_container()->agent_grid[voxels_to_search_indicies[i]][j]);
      }
    }
  }
  #pragma omp critical
  {
       cells_found.insert(cells_found.end(), cells_found_private.begin(), cells_found_private.end());
  }
  return cells_found;
}
void custom_add_potentials(Spring_Cell* SpCell)
{
  double dt=0.01;
  std::vector <double> sum_of_forces_private(3,0.0);
  std::vector <double> sum_of_forces;
  #pragma omp private(sum_of_forces_private)
  for (size_t i = 0; i < SpCell->m_springs.size(); i++)//loop through my springs
  {
    double test_length=2.0;
    double spring_stretch=distance_between_membranes(SpCell->m_my_pCell,SpCell->m_springs[i]->m_pNeighbor);
    std::vector<double> force_direction= displacement_between_membranes(SpCell->m_my_pCell,SpCell->m_springs[i]->m_pNeighbor);//from me to neighbor
    sum_of_forces_private+=Hookes_law_force(force_direction,SpCell->m_springs[i]->m_length,spring_stretch-test_length,SpCell->m_my_pCell->custom_data["cell_k"]);
    //std::cout<<"spring stretch "<<spring_stretch<<"\n";
    //std::cout<<"force_direction "<<force_direction<<"\n";
    //std::cout<<"force_direction "<<force_direction<<"\n";
    // f/m=a=dv/dt
  }
  
  #pragma omp critical
  {
    sum_of_forces.insert(sum_of_forces.end(),sum_of_forces_private.begin(),sum_of_forces_private.end());
    //std::cout<<"sum of forces"<<sum_of_forces<<std::endl;
    sum_of_forces[0]=(dt*sum_of_forces[0]/SpCell->m_my_pCell->phenotype.volume.total); //f/m=a=(dv/dt)
    sum_of_forces[1]=(dt*sum_of_forces[1]/SpCell->m_my_pCell->phenotype.volume.total); //f/m=a=(dv/dt)
    sum_of_forces[2]=(dt*sum_of_forces[2]/SpCell->m_my_pCell->phenotype.volume.total); //f/m=a=(dv/dt)
    SpCell->m_my_pCell->velocity+=sum_of_forces;
  }
  return;
}
void custom_add_potentials_for_pCells(Cell* pCell)
{
  Spring_Cell* current_spring_cell=spring_cell_by_pCell_index[pCell->index];
  custom_add_potentials(current_spring_cell);
  return;
}
std::vector <Spring_Cell*> spring_cells_in_neighborhood(Spring_Cell* SpCell, double maximum_interaction_distance)
{
  //std::cout<<"Third!!"<<"\n";
  //TODO fix redundant code from above function
  std::vector <int> voxels_to_search_indicies=max_interaction_variable_moore_neighborhood(SpCell->m_my_pCell,maximum_interaction_distance);
  //std::cout<< "searching "<< voxels_to_search_indicies.size()<<std::endl;
  std::vector <Cell*> cells_found_private{};
  std::vector <Cell*> cells_found{};
  std::vector <Spring_Cell*> spring_cells_found{};
  #pragma omp private(cells_found_private)
  for (size_t i = 0; i < voxels_to_search_indicies.size(); i++)
  {
   // std::cout<<"checking "<<SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]].size()<<std::endl;
    for (size_t j = 0; j < SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]].size(); j++) 
    {
   //  std::cout<< "found "<< SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]][j]<<std::endl;
      if(SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]][j]!=SpCell->m_my_pCell)
      {
        cells_found_private.push_back(SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]][j]);
      }
    }
  }
  #pragma omp critical
  {
    cells_found.insert(cells_found.end(), cells_found_private.begin(), cells_found_private.end());
    spring_cells_found.resize(cells_found.size());
  }
  for(size_t i =0; i< cells_found.size(); i++)
  {
    spring_cells_found[i]=spring_cell_by_pCell_index[cells_found[i]->index];
   //std::cout<< "spring cells found: " <<spring_cell_by_pCell_index[cells_found[i]->index]<<"\n";
  }
  return spring_cells_found;
}
void initialize_spring_connections()
{
  
  //using modified moore neighborhood search for "spring cell" neighbors and connect them with spring objects
  for(size_t i=0;i<all_spring_cells.size();i++) //loop through all spring cells
  {
    //std::cout<<"Second!!"<<std::endl;
    Spring_Cell* SpCell=all_spring_cells[i];
    //std::cout<<" spcell "<<SpCell<<std::endl;
    std::vector<Spring_Cell*> SpNeighbors=spring_cells_in_neighborhood(SpCell,1);
   
    for (size_t j = 0; j < SpNeighbors.size(); j++)
    {
      //std::cout <<"Connecting "<<SpNeighbors[j]<<"\n";
      connect_spring_cells(SpCell,SpNeighbors[j]);
    }
  }
  return;
}
void create_spring_cell_from_cell(Cell* pCell, Phenotype& phenotype)
{
  //cell types must be created individually to pass phenotype
  
  return;
}
void initialize_spring_cells()//make all cells spring ce
{

  //encapsulate all cells in the super class spring cell
  //find initial spring lengths and connect neighbors
  //set up basement membrane and connected exterior cells
  std::vector <double> basement_membrane_center={0.0,0.0,0.0};
  double basement_membrane_radius=100;
  spring_cell_by_pCell_index.resize((*all_cells).size());
  for (size_t i = 0; i < (*all_cells).size(); i++)
  {
    Spring_Cell* springy=create_spring_cell((*all_cells)[i]);
  }
  

  initialize_spring_connections();
  initialize_basement_membrane(basement_membrane_center,basement_membrane_radius);
  initialize_basement_membrane_connections( basement_membrane_center, basement_membrane_radius);
  
  return;
}
void initialize_basement_membrane_connections( std::vector <double> basement_membrane_center, double basement_membrane_radius)
{
  //TODO: should be thread safe but verify
  for (int i=0; i< basement_membrane_voxels.size();i++)
  {
    for( int j=0; j<(all_spring_cells).size();j++)
    {
      Spring_Cell* SpCell=(all_spring_cells)[j];
      if(is_in_voxel(SpCell->m_my_pCell,basement_membrane_voxels[i]))
      {
        double length=basement_membrane_radius-norm(SpCell->m_my_pCell->position-basement_membrane_center)-SpCell->m_my_pCell->phenotype.geometry.radius;
        SpCell->is_basement_connected=true;
        SpCell->basement_length=length;
      }
    }
  }

  return;
}
bool is_in_voxel(Cell* pCell, Voxel* pVoxel)
{
  int voxel_index=pVoxel->mesh_index;
  bool value=false;
  for (int i = 0; i < pCell->get_container()->agent_grid[voxel_index].size(); i++) 
  {
  // std::cout<< i <<std::endl;
    if (pCell->get_container()->agent_grid[voxel_index][i] == pCell) 
    {
      value=true;
    }
  }
  return value;
}
bool is_in_voxel(Cell* pCell, int voxel_index)
{
  bool value=false;
  for (int i = 0; i < pCell->get_container()->agent_grid[voxel_index].size(); i++) 
  {
  // std::cout<< i <<std::endl;
    if (pCell->get_container()->agent_grid[voxel_index][i] == pCell) 
    {
      value=true;
    }
  }
  return value;
}
void basement_membrane_mechanics(Spring_Cell* SpCell, double basement_membrane_radius, std::vector <double> basement_membrane_center)
{
  //basement membrane radius is the center of the membrane which is 2-3 voxels thick
  std::vector<double>force_direction=-1*normalize(SpCell->m_my_pCell->position);//force direction always inward
  for (int i=0; i< basement_membrane_voxels.size();i++)
  {
    if(is_in_voxel(SpCell->m_my_pCell,basement_membrane_voxels[i]))
    {
      double rest_length=SpCell->m_my_pCell->phenotype.geometry.radius; //not connected
     /// microenvironment.mesh.voxels(basement_membrane_voxels[i])
     if(SpCell->is_basement_connected==true)
     {
        rest_length=SpCell->basement_length;
     }
      double current_length=basement_membrane_radius-norm(SpCell->m_my_pCell->position-basement_membrane_center)-SpCell->m_my_pCell->phenotype.geometry.radius;
     if(current_length+SpCell->m_my_pCell->phenotype.geometry.radius+default_microenvironment_options.dx<0)
     {
      //updated to use outer_radius from basement membrane function
      std::cout<< "WARNING!! CELL PASSED THROUGH BASEMENT MEMBRANE!"<<std::endl;
     } 
      current_length=current_length-2.0;
      std::vector<double> force= Hookes_law_force(force_direction,rest_length,current_length,SpCell->m_my_pCell->custom_data["basement_k"]);
      double dt=0.01;
      #pragma omp critical
      {
        std::cout<<"MASS: "<<SpCell->m_my_pCell->phenotype.volume.total<<"\n";
        force[0]=force[0]*dt/SpCell->m_my_pCell->phenotype.volume.total;//a=F/m
        force[1]=force[1]*dt/SpCell->m_my_pCell->phenotype.volume.total;//a=F/m
        force[2]=force[2]*dt/SpCell->m_my_pCell->phenotype.volume.total;//a=F/m
        std::cout<<"FORCE: "<<force<<"\n";
        SpCell->m_my_pCell->velocity=SpCell->m_my_pCell->velocity+force;
        std::cout<<"VELOCITY: "<<force<<"\n";
      }
      return;
    } 
  }
} 
void basement_membrane_mechanics(Cell* pCell, double basement_membrane_radius) //for future versions in follicle distance test will be faster currently not working
{
  std::vector <int>  test_voxel_indicies{};
  //mechanics and diffusion are using the same mesh as of may 23 2023
  std::vector <int> interior_indicies=get_interior_voxels(pCell);
  std::vector <int> exterior_indicies=get_exterior_voxels(pCell);
  if(exterior_indicies.size()!=1)//allows for testing if multivoxel cell intersects the basement membrane
  {
    #pragma omp critical
    {
      test_voxel_indicies.resize(interior_indicies.size()+exterior_indicies.size());
    }
    int exterior_count=0;
    for(int i=0;i<test_voxel_indicies.size();i++)
    {
      if(i<interior_indicies.size())
      {
        test_voxel_indicies[i]=interior_indicies[i];
      }
      else
      {
        test_voxel_indicies[i]=exterior_indicies[exterior_count];
        exterior_count++;
      }
    }
  }
  for(int j=0; j<test_voxel_indicies.size();j++)
  {
    for(int k=0; k<basement_membrane_voxels.size();k++)
    {
      if(test_voxel_indicies[j]==basement_membrane_voxels[k])
      {
        
      }
    }
  }
  return;
}

std::vector <int> spherical_bounding_box(std::vector <double> center_point, double radius)
{
  //gets voxels of a bounding box around an arbitrary sphere
  //NOTE!!! currently PhysiCells microenvironment voxels are both mechanics and diffusion
  std::vector<int> bounding_box_by_index={};
  //bounding box of voxels in the microenvironment
 //for speed get center voxel figure out x,y,z offset and only check local voxels
    Voxel center_voxel=microenvironment.nearest_voxel(center_point);
    std::vector<double> center_voxels_center= center_voxel.center; // get center voxel
    std::vector<double> offset=center_voxels_center-center_point;
    //figure out the corners of my bounding box of voxels to search
    std::vector<double> lower_point(3,0.0);
    std::vector<double> upper_point(3,0.0);
    for (size_t i = 0; i < 3; i++)
    {
      lower_point[i]= center_voxels_center[i]-offset[i]-(ceil(radius));
      upper_point[i]= center_voxels_center[i]+offset[i]+ceil(radius);
    }
    Voxel starting_voxel=microenvironment.nearest_voxel(lower_point);
    Voxel ending_voxel=microenvironment.nearest_voxel(upper_point);
    bounding_box_by_index=general_voxel_bounding_box(starting_voxel.center,ending_voxel.center,default_microenvironment_options.dx,microenvironment.mesh);
  return bounding_box_by_index;
}

void print_voxels_for_quick_plotting(Cell* pCell,std::vector <int> bounding_voxels, std::vector <int> sub_section)
{
      std::cout<<"bounding voxels centers "<<std::endl;
    std::cout<<"[ ";
    for(int n=0; n<bounding_voxels.size(); n++)
    {  
      if(n<bounding_voxels.size()-1)
      { 
        std::cout<<pCell->get_container()->underlying_mesh.voxels[bounding_voxels[n]].center[0]<<", ";
      }
      else
      {
        std::cout<<pCell->get_container()->underlying_mesh.voxels[bounding_voxels[n]].center[0];
      }
    }
    std::cout<<"],";
    std::cout<<"[ ";
    for(int n=0; n<bounding_voxels.size(); n++)
    { 
      if(n<bounding_voxels.size()-1)
        { 
          std::cout<<pCell->get_container()->underlying_mesh.voxels[bounding_voxels[n]].center[1]<<", ";
        }
      else
      {
        std::cout<<pCell->get_container()->underlying_mesh.voxels[bounding_voxels[n]].center[1];
      }
    }
    std::cout<<"]";
    std::cout<<" ]"<<std::endl;
    std::cout<<"interior voxel centers "<< std::endl;
    std::cout<<"[ ";
    for(int m=0; m<sub_section.size(); m++)
    {
        if(m<sub_section.size()-1)
        {
          std::cout<<pCell->get_container()->underlying_mesh.voxels[sub_section[m]].center[0]<<", ";
        }
        else
        {
          std::cout<<pCell->get_container()->underlying_mesh.voxels[sub_section[m]].center[0];
        }
    }
    
    std::cout<<"],";
    std::cout<<"[ ";
     for(int m=0; m<sub_section.size(); m++)
    {
        if(m<sub_section.size()-1)
        {
          std::cout<<pCell->get_container()->underlying_mesh.voxels[sub_section[m]].center[1]<<", ";
        }
        else
        {
          std::cout<<pCell->get_container()->underlying_mesh.voxels[sub_section[m]].center[1];
        }
    }
    std::cout<<" ]";
    std::cout<<" ]"<<std::endl;
    return;
}
double moles_in_voxel(int voxel_id, std::string solute_name)
{
  int vox=microenvironment.find_density_index(solute_name);
  double moles=0.0;
  //get water volume of voxel
  std::vector <double> solute_specific_volumes={0.0,1.1,2.2,3.3,4.4};//grams/mole//TODO put in real values and move to config

  double total_volume=microenvironment.mesh.voxels[voxel_id].volume;
  double sum=0.0;//get volume taken up by all solutes in a voxel
  #pragma omp private(sum)
  for (size_t i = 0; i < microenvironment.density_vector(vox).size(); i++)
  {
    if(solute_specific_volumes.size()!=microenvironment.density_vector(vox).size())
    {
      std::cout<<"MISSING/EXTRA SPECIFIC VOLUME!"<<"\n";
      return moles;
    }
    //sum+=microenvironment.density_vector(vox)[i]*solute_molar_masses[i];//mole/kg/m^3/kg*10^18um^3/m^3
  }
  

  return moles;
}