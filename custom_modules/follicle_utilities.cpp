#include "./follicle_utilities.h"
// #include "./springs.h"
using namespace BioFVM;
using namespace PhysiCell;
using namespace Springs;
// using namespace Springs;
//  #include <stdef.h>
//   Global Variables - will probably add to XML eventually to avoid macros
#define PI 3.14159265
#define R 0.08205 // granulosa (10^-3 J/mole*k)
//const static std::vector<std::vector<double>> molar_conversion_coeff={{1.1,1.2},{2.1,2.2},{3.1,3.2}};
/* returns the disance vector pointing from cell_1 to cell_2 */
std::vector<double> displacement_between_membranes(Cell *pCell_1, Cell *pCell_2) {
  // from pCell_1 to pCell_2 vector operator defined in bioFVM
  std::vector<double> displacement = pCell_2->position - pCell_1->position;
  displacement = distance_between_membranes(pCell_1, pCell_2) * normalize(displacement);
  return displacement;
}
/* return the distance between spherical cell membranes always positive*/
double distance_between_membranes(Cell *pCell_1, Cell *pCell_2) {
  double distance = std::abs(norm(pCell_2->position - pCell_1->position) - pCell_1->phenotype.geometry.radius - pCell_2->phenotype.geometry.radius);
  //lets avoid floating-point hell
  
  if(std::abs(distance)<1e-16){distance=0;}
  return distance;
}
double signed_distance_between_membranes(Cell *pCell_1, Cell *pCell_2) {
  double distance = norm(pCell_2->position - pCell_1->position) - pCell_1->phenotype.geometry.radius - pCell_2->phenotype.geometry.radius; 
  if(std::abs(distance)<1e-16){distance=0;}
  return distance;
}

// Initialize custom neighborhood based on distance between membranes
// loop through all neighboors
void initialize_neighboring_spring_connections() {
  //old version probably deprecate
//  for (int i = 0; i < (*all_cells).size(); i++) {
//
//    Cell *pCell = (*all_cells)[i];
//    attach_neighboring_springs(pCell, 1.0);
//  }
  return;
}
void attach_neighboring_springs(Cell *pCell, double max_spring_length) 
{
  // this is old code now using Spring_Cell class need to deprecate and fix any dependent functions
  // check nearby cells to identify neighboors using built in mechanics voxel
  // if less than max spring length store (pCell, rest length)
//  std::vector<Cell *> neighbors = find_nearby_cells(pCell);
//#pragma omp critical // make safe for push_back
//  {
//    for (int i = 0; i < neighbors.size(); i++) {
//      Cell *pNeighbor = neighbors[i];
//      double distance = distance_between_membranes(pCell, pNeighbor);
//      if (distance <= max_spring_length) {
//        //testing spring class here
//        //pCell->state.connection_and_spring_length.push_back(std::make_pair(pNeighbor, distance));
//        pCell->state.neighbors.push_back(pNeighbor); // uses built in neighboor storage, default velocity
//                        // function must be off or this will get overwritten
//      }
//    }
//  }
  return;
}
//   // check nearby cells to identify neighboors using built in mechanics voxel
// void dettach_neighboring_springs(Cell *pCell, double max_spring_length) 
// {  // if less than max spring length store (pCell, rest length)
//   std::vector<Cell *> neighbors = find_nearby_cells(pCell);
//   #pragma omp critical // make safe for push_back
//   {
//     for (int i = 0; i < neighbors.size(); i++) {
//       Cell *pNeighbor = neighbors[i];
//       double distance = distance_between_membranes(pCell, pNeighbor);
//       if (distance > max_spring_length) {
//         //testing spring class here
//         //pCell->state.connection_and_spring_length.push_back(std::make_pair(pNeighbor, distance));
//         pCell->state.neighbors.push_back(pNeighbor); // uses built in neighboor storage, default velocity
//                         // function must be off or this will get overwritten
//       }
//     }
//   }
//   return;
//
// }
void diffusion_bounding_box(Cell* pCell, std::vector<int>* bounding_box_by_index)
{
  //gets the box of voxels cell is contained in, includes the voxels inside
  //bounding box of voxels in the microenvironment
  //for speed get center voxel figure out x,y,z offset and only check local voxels
  int center_voxel_index=microenvironment.nearest_voxel_index(pCell->position);
  std::vector<double> center_voxel_position=pCell->get_container()->underlying_mesh.voxels[center_voxel_index].center;
  std::vector<double> offset=center_voxel_position-pCell->position;
  //figure out the corners of my bounding box of voxels to search
  std::vector<double> lower_point(3,0.0);
  std::vector<double> upper_point(3,0.0);
  for (size_t i = 0; i < 3; i++)
  {
    lower_point[i]= center_voxel_position[i]-offset[i]-ceil(pCell->phenotype.geometry.radius);
    upper_point[i]= center_voxel_position[i]+offset[i]+ceil(pCell->phenotype.geometry.radius);
  }
  std::vector<double> starting_voxel_center=microenvironment.nearest_voxel(lower_point).center;
  std::vector<double> ending_voxel_center=microenvironment.nearest_voxel(upper_point).center;
  general_voxel_bounding_box(bounding_box_by_index,starting_voxel_center,ending_voxel_center,default_microenvironment_options.dx,microenvironment.mesh);
  //std::cout<<"Bounding box"<<bounding_box_by_index<<"\n";
  return;
}
void get_voxel_corners(std::vector<double> &voxel_center, std::vector<std::vector<double>> &return_corners )
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
  #pragma omp critical
  {
    return_corners[0].assign(corners[0].begin(),corners[0].end());
    return_corners[1].assign(corners[1].begin(),corners[1].end());
    return_corners[2].assign(corners[2].begin(),corners[2].end());
    return_corners[3].assign(corners[3].begin(),corners[3].end());
    return_corners[4].assign(corners[4].begin(),corners[4].end());
    return_corners[5].assign(corners[5].begin(),corners[5].end());
    return_corners[6].assign(corners[6].begin(),corners[6].end());
    return_corners[7].assign(corners[7].begin(),corners[7].end());
  }
  return;
}
void get_intersecting_voxels(Cell* pCell,std::vector<int>* return_intersecting_voxel_indicies)
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector<int> bounding_voxels={};
    diffusion_bounding_box(pCell,&bounding_voxels);
    std::vector<int> intersecting_voxels={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      get_voxel_corners(test_voxel_center,test_corners);
    
     
      int sum=0;
      // #pragma omp private(sum,intersecting_voxels)
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
      return_intersecting_voxel_indicies->assign(intersecting_voxels.begin(),intersecting_voxels.end());
      //return_intersecting_voxel_indicies->insert(return_intersecting_voxel_indicies->end(), intersecting_voxels.begin(), intersecting_voxels.end());
    }

    return;
}
void get_exterior_voxels(Cell* pCell, std::vector<double>* return_exterior_voxel_indicies)//some but not all the coners of the voxel are within the cell and the voxel center is outside the cell
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector <int> bounding_voxels={};
    diffusion_bounding_box(pCell,&bounding_voxels);
    std::vector<int> exterior_voxels={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      // std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      get_voxel_corners(test_voxel_center,test_corners);
    
     
      int sum=0;
      // #pragma omp private(sum,exterior_voxels)
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
      return_exterior_voxel_indicies->assign(exterior_voxels.begin(), exterior_voxels.end());
    }

    return;
}
void get_interior_voxels(Cell* pCell, std::vector<int>* return_interior_voxel_indicies) //voxel center is inside the cell
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector<int> bounding_voxels={};
    diffusion_bounding_box(pCell,&bounding_voxels);
    
    std::vector<int> interior_voxels={};
    // #pragma omp private(interior_voxels)
    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      // std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      get_voxel_corners(test_voxel_center,test_corners);
     
      int sum=0;
      // #pragma omp private(interior_voxels)

      if(norm(test_voxel_center-pCell->position)<=pCell->phenotype.geometry.radius)
      {
        interior_voxels.push_back(bounding_voxels[i]);
      }

    }
    #pragma omp critical
    {
      return_interior_voxel_indicies->assign(interior_voxels.begin(), interior_voxels.end());
    }

    return;
}

void solute_loading( double oxygen, double final_solute_concentration_1, double final_solute_concentration_2,double final_solute_concentration_3,double final_solute_concentration_4, double final_solute_concentration_5)
{
  //function to set dirichlet nodes during simulation
  //TODO: Do I still need this? commented out for now, but would need restructuring
  ////Dirichlet node for all the voxels located outside of the ring
//   std::vector<double> dirichlet_solute_end_state( 3 , 0.0 );
// 	dirichlet_solute_end_state[0]=9.0;//oxygen;//oxygen
// 	dirichlet_solute_end_state[1]=8.0;//final_solute_concentration_1;//hm
// 	dirichlet_solute_end_state[2]=10.0;//final_solute_concentration_2;//eg
// 	//dirichlet_solute_end_state[3]=final_solute_concentration_3;//gly
//   //dirichlet_solute_end_state[4]=final_solute_concentration_4;//pbs
// 	//dirichlet_solute_end_state[5]=final_solute_concentration_5;//suc
// 	//dirichlet_solute[4]=0.0;//sucrose
// 	//Dirichlet nodes for all the voxels located inside of the ring
// /*	std::vector<double> dirichlet_solute_start_state( 3 , 0 );
// 	dirichlet_solute_start_state[0]=0.0;//annoying oxygen no idea how to get rid of it but its turned off
// 	dirichlet_solute_start_state[1]=initial_solute_concentration_1;//final_external_osmolarity;//eg
// 	dirichlet_solute_start_state[2]=initial_solute_concentration_2;//gly
// */
// /*
// 	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )
// 	{
// 		if(dist(microenvironment.voxels(i).center, {0.0,0.0,0.0})>=130)
// 		{
// 			microenvironment.update_dirichlet_node( i,dirichlet_solute_end_state );
// 		}
// 		else
// 		{
// 			microenvironment.remove_dirichlet_node(i);
// 		}
// 	}
// */
// 	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )//something here isn't right TODO: Check this
// 	{
// 		if(dist(microenvironment.voxels(i).center, {0.0,0.0,0.0})>=10)
// 		{
// 			// microenvironment.density_vector(i)=dirichlet_solute_end_state;//<<-- is this right?
// 		}
// 		else
// 		{
// 			//microenvironment.remove_dirichlet_node(i);
// 		}
// 	}
	return;
}
//TODO: update velocity to use Adams_Bashforth
void activate_nodes(double radius_of_activation)
{
  
	for( int i=0; i < microenvironment.number_of_voxels() ; i++ )//
	{
		if(dist(microenvironment.voxels(i).center, {0.0,0.0,0.0})>=radius_of_activation)
		{
		  microenvironment.add_dirichlet_node(i, default_microenvironment_options.Dirichlet_condition_vector);//<<-- is this right?
      microenvironment.density_vector(i)=default_microenvironment_options.Dirichlet_condition_vector;
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
  // exterior_voxel_index=get_exterior_voxels(pCell);
  Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index];
  if(norm(SPcell->previous_position-pCell->position)>default_microenvironment_options.dx && std::abs(SPcell->previous_radius-pCell->phenotype.geometry.radius)>default_microenvironment_options.dx)
  {
    get_intersecting_voxels(pCell,&(SPcell->uptake_voxels));
    // #pragma omp critical
    // {
    //   SPcell->uptake_voxels.insert(SPcell->uptake_voxels.end(),exterior_voxel_index.begin(),exterior_voxel_index.end() );
    // }

  } 
  
  if(SPcell->uptake_voxels.size()<2)//for cells smaller than or approx. equal to one voxel
  {
    average=microenvironment.nearest_density_vector( microenvironment.nearest_voxel_index(pCell->position))[solute_index];
  }
  else 
  {
    // #pragma omp reduction(+:sum)
    for (size_t i = 0; i < SPcell->uptake_voxels.size(); i++)
    {
       // std::cout<<" concentration "<<microenvironment.nearest_density_vector(exterior_voxel_index[i] )[solute_index]<<"\n";
        //std::cout<<" exterior voxel index"<< exterior_voxel_index[i]<<"\n";
        //std::cout<<" concentration for sum: "<<microenvironment.nearest_density_vector(exterior_voxel_index[i] )[solute_index]<<"\n";
      if(microenvironment.nearest_density_vector(SPcell->uptake_voxels[i])[solute_index]<1e-16)// to deal with numerical instability from tiny inirial diffusion values inside rasterized cell
      {
          sum+=0.0;
      }
      else {
        sum+=microenvironment.nearest_density_vector( SPcell->uptake_voxels[i])[solute_index];
      }
    }
    #pragma omp critical
    {
      average=sum/SPcell->uptake_voxels.size();
    }
  } 
  return average; 
}



void cells_in_me(Cell *pCell, std::vector<Cell*> *return_cells_in_me) // uses mechanics vectors to search for cells that are within or equal to pCell radius can also be used for bounding boxes
{
  //NOTE: current physicell uses the same cartesian mesh for mechanics and diffusion?? could use functions for diffusion bounding box
  std::vector<int> search_voxels={};
  diffusion_bounding_box(pCell,&search_voxels);
  std::vector<Cell *> agents_in_voxel={};
  // #pragma omp private(agents_in_voxel)
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
    return_cells_in_me->insert(return_cells_in_me->end(), agents_in_voxel.begin(), agents_in_voxel.end());
  }
      // std::cout<<"voxel " << temp_voxel.mesh_index<< " is inside my radius
      // and located at: "<< temp_voxel.center[0]<<",
      // "<<temp_voxel.center[1]<<", "<<temp_voxel.center[2]<<std::endl;
    
  // std::cout<<agents_in_me.size()<< " agents in me:
  // "<<agents_in_me[0]<<std::endl;
  return;
}
void find_basement_membrane_voxels(std::vector<double> center_of_sphere, std::vector<double> radius_of_sphere) {
  return;
}
void Hookes_law_force_vector(std::vector<double> &current_spring_vector, std::vector<double> &equilibrium_spring_vector,double &spring_constant, std::vector<double> *return_force) 
{
  std::vector<double> force(3,0.0);
  force = -1.0*spring_constant * (equilibrium_spring_vector-current_spring_vector);//f=-k*delta_X if delta_X is negative spring was stretched
  #pragma omp critical
  {
    return_force->assign(force.begin(),force.end());
  }
  return;
}
void Hookes_law_force_magnitude(std::vector<double> &my_position, std::vector<double> &neighbor_position, double &equilibrium_length, double &current_length, double &spring_constant, std::vector<double> *return_force)
{
  std::vector<double> force(3,0.0);
  std::vector<double> force_direction= neighbor_position-my_position;
  double delta_x=0;
  bool is_towardNeighbor=false;
  if(equilibrium_length<0)//both negative or both positive
  {
    if(current_length<0)
    {
      delta_x=std::abs(equilibrium_length)-std::abs(current_length);
      if(std::abs(delta_x)<1e-16){delta_x=0;return;}
      else if (delta_x<0){
        is_towardNeighbor=true;
      }
    }
    else {
      delta_x=std::abs(equilibrium_length-current_length);
      is_towardNeighbor=true;
    }
  }
  else if(equilibrium_length>=0){
    if(current_length>=0)
    {
      delta_x=std::abs(equilibrium_length)-std::abs(current_length);
      if(std::abs(delta_x)<1e-16){delta_x=0;return;}
      else if(delta_x>=0){
        is_towardNeighbor=true;
      }
    }
    else {
      delta_x=std::abs(equilibrium_length-current_length);
    }
  }
  // std::cout<<"delta X:"<< delta_x<<"\n";
  if(is_towardNeighbor==true)// force is along force_direction
  {
      
    (*return_force)[0]=1.0*(spring_constant *delta_x)* (normalize(force_direction))[0];
    (*return_force)[1]=1.0*(spring_constant *delta_x)* (normalize(force_direction))[1];
    (*return_force)[2]=1.0*(spring_constant *delta_x)* (normalize(force_direction))[2];
    
  }
  else { //force is along -force_direction
    (*return_force)[0]=-1.0*(spring_constant *delta_x)* (normalize(force_direction))[0];
    (*return_force)[1]=-1.0*(spring_constant *delta_x)* (normalize(force_direction))[1];
    (*return_force)[2]=-1.0*(spring_constant *delta_x)* (normalize(force_direction))[2];
  
    
  }

  // return_force->assign(force.begin(),force.end());
  return;
}
void update_all_forces(Cell* pCell, double dt, double spring_constant) {
  Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index];
  SPcell->previous_radius=pCell->phenotype.geometry.radius;
  // SPcell->previous_velocity=pCell->velocity;
  SPcell->previous_position=pCell->position;
  //these both update velocity
  // double test_radius=70;
  double basement_radius=95;//follicle_radius=oocyte radius+ 3* granulosa radius+ granulosa radius TODO: set this from xml
  std::vector<double> basement_center={0.0,0.0,0.0};
  custom_add_potentials_for_pCells(pCell,dt);
  non_connected_neighbor_pressure(pCell,dt,spring_constant);
  // function is designed so you could have a changing basement membrane but we set it static 
  basement_membrane_mechanics(SPcell,basement_radius,basement_center,dt);  
  return;
}
void non_connected_neighbor_pressure(Cell* pCell, double dt, double spring_constant) {
  // similar to the new function in PhysiCell 1.10 automatically pushes on cells
  // that press up against the cell; however this version doesnt require the
  // cell to be the same size as mechanics voxel uses same spring value as
  // connections, avoids double pushing on connected cells
  // probably much slower than it could be, optomize in future
//first find all possible neighbors
  Spring_Cell* SpCell=spring_cell_by_pCell_index[pCell->index];//access the spring cell for this pCell
  std::vector<Cell *> possible_neighbors ={}; 
  cells_in_me(pCell,&possible_neighbors); // vector containing cells that could be interacting
                                                               // with pCell that aren't in neighboors
  std::vector<Cell *> non_connected_neighbors={};
  //collect all the possible neighboors that aren't connected
  // #pragma omp private(possible_neighbors,found_neighbors)
  for (int j = 0; j < possible_neighbors.size(); j++) {
    for (int i = 0; i < spring_cell_by_pCell_index[pCell->index]->m_springs.size(); i++) {
      if(possible_neighbors.size()>0 && j<possible_neighbors.size()) //ugly way to deal with changing size should probably be a while loop instead
      {
        if (possible_neighbors[j] == spring_cell_by_pCell_index[pCell->index]->m_springs[i]->m_pNeighbor) {
          possible_neighbors.erase(possible_neighbors.begin()+j);
        }
      }
    }
  }
  #pragma omp critical
  {
    non_connected_neighbors.assign(possible_neighbors.begin(),possible_neighbors.end());
  }
  // std::cout<<"non connected neighbors list"<<
  // non_connected_neighbors.size()<< std::endl;
  double sum_x_acceleration =0.0;// pCell->velocity[0];
  double sum_y_acceleration =0.0; //pCell->velocity[1];
  double sum_z_acceleration =0.0;// pCell->velocity[2];
  // #pragma omp private(found_neighbors) reduction(+:sum_x_velocity,sum_y_velocity,sum_z_velocity)
  // go through all my possible non-connected neighboors and if they overlap me in space apply a force to me
  if (non_connected_neighbors.size()>0) {
    for (int i = 0; i < non_connected_neighbors.size(); i++) {
      Cell* neighbor=non_connected_neighbors[i];
      // std::cout<<neighbor<<std::endl;
      double neighbor_effective_radius=neighbor->phenotype.geometry.radius-neighbor->custom_data["allowed_overlap"];
      double pCell_effective_radius=pCell->phenotype.geometry.radius-pCell->custom_data["allowed_overlap"];
      double membrane_distance= norm(neighbor->position-pCell->position)-neighbor_effective_radius-pCell_effective_radius;
      // if statement assumes strict interpentration of neighboors membrane into
      // pCells membrane to exert spring pressure
      if (membrane_distance<0) 
      {
        // std::cout<<"NON-CONNECTED CALLED!"<<"\n";
        std::vector<double> force_on_pCell(3,0.0);
        double rest_length=0;
        Hookes_law_force_magnitude(pCell->position,neighbor->position,rest_length,membrane_distance,spring_constant, &force_on_pCell);
        double neighbor_mass = neighbor->phenotype.volume.total; // mass_of_cell(neighboor);
        //TODO: make mass of cell function
        double pCell_mass = pCell->phenotype.volume.total;
        //#pragma omp private(neighbor) for
        //test that the velocity component isn't going to cause a floating decimal issue if its super small
        sum_x_acceleration += (force_on_pCell[0] / pCell_mass);
        sum_y_acceleration += (force_on_pCell[1] / pCell_mass);
        sum_z_acceleration += (force_on_pCell[2] / pCell_mass);
      }
      // TODO: possibly, eventually use Adams_Bashforth_ODE_2nd_Order for velocity
      // velocity_x= previous_x+ 
    }
  }
#pragma omp critical //probably not needed but the function could be run with pragma omp parallel for
  {
    // std::cout<<"Sums: "<< sum_x_acceleration<<", "<<sum_y_acceleration<<", "<<sum_z_acceleration<<"\n";
    // velocity_x= previous_x+ 
    pCell->velocity[0] = SpCell->previous_velocity[0]+sum_x_acceleration*dt;
    pCell->velocity[1] = SpCell->previous_velocity[1]+sum_y_acceleration*dt;
    pCell->velocity[2] = SpCell->previous_velocity[2]+sum_z_acceleration*dt;
    SpCell->previous_velocity=pCell->velocity;
  }
  return;
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
std::vector<std::vector<double>> create_spheroid(double cell_radius, double sphere_radius) 
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
          /* output file of initial positions
          std::ofstream ofs;
          ofs.open ("temp_points_location.csv", std::ofstream::out |
          std::ofstream::app); ofs <<(sqrt(norm_squared(tempPoint)))<<"," <<"\n";
          ofs.close();
          */
          cells.push_back(tempPoint);
        }
      }
    }
  }
  return cells;
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
            /* output file of initial positions
            std::ofstream ofs;
            ofs.open ("temp_points_location.csv", std::ofstream::out |
            std::ofstream::app); ofs <<(sqrt(norm_squared(tempPoint)))<<"," <<"\n";
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

void break_TZPs(Spring_Cell* SP_Oocyte,double max_breakage_distance) 
{
  //TAG: possible issue here
  for (int i=0; i<SP_Oocyte->m_springs.size();i++) {//loop through the oocytes springs
    Spring* S=SP_Oocyte->m_springs[i];//S is the ith spring
    if(signed_distance_between_membranes(SP_Oocyte->m_my_pCell,S->m_pNeighbor)>max_breakage_distance)
    {
      Spring_Cell* SP_Neighbor=spring_cell_by_pCell_index[S->m_pNeighbor->index];
      SP_Oocyte->remove_spring(SP_Neighbor);
      SP_Neighbor->remove_spring(SP_Oocyte);
    }
  } 
  return;
}
void update_interior_concentrations(Spring_Cell* SPcell)
{
  double osmotically_active_volume=SPcell->water_volume+SPcell->solute_volume;
  // std::cout<<"osmotically_active_volume: "<< osmotically_active_volume<<"\n"; 
  // std::cout<<"should be the same volume: "<< SPcell->m_my_pCell->phenotype.volume.total-SPcell->solid_volume<<"\n";
  // std::cout<<"total volume: "<< SPcell->m_my_pCell->phenotype.volume.total<<"\n";
  for(size_t i=0; i<SPcell->solute_moles.size(); i++)// loop through vector of solute_moles indexed by solute 
  {
    SPcell->interior_molarity[i]=SPcell->solute_moles[i]/osmotically_active_volume;// calculate molarity
  }
  if(SPcell->simulation_selected==1)
    {// only EG selected
      SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
      SPcell->interior_component_molality[1]=molarity_to_molality(SPcell->interior_molarity[1], "EG");
      SPcell->interior_osmolality=ternary_virial(SPcell->interior_component_molality[0],SPcell->interior_component_molality[1],"NaCl","EG");
  }
  else if(SPcell->simulation_selected==2)
    {// only EG selected
      SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
      SPcell->interior_component_molality[1]=molarity_to_molality(SPcell->interior_molarity[1], "GLY");
      SPcell->interior_osmolality=ternary_virial(SPcell->interior_component_molality[0],SPcell->interior_component_molality[1],"NaCl","GLY");
  }
  else if(SPcell->simulation_selected==3)
    {// only EG selected
      SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
      SPcell->interior_osmolality=binary_virial(SPcell->interior_component_molality[0],"NaCl");
  }
  else if(SPcell->simulation_selected==4)
    {// only EG selected
      SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
      SPcell->interior_component_molality[1]=molarity_to_molality(SPcell->interior_molarity[1], "EG");
      SPcell->interior_component_molality[2]=molarity_to_molality(SPcell->interior_molarity[2], "GLY");
      SPcell->interior_osmolality=binary_virial(SPcell->interior_component_molality[0],"NaCl")+ternary_virial(SPcell->interior_component_molality[1],SPcell->interior_component_molality[2],"EG","GLY");
  }
  else if(SPcell->simulation_selected==5)
    {// karlsson EG selected
      SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
      SPcell->interior_component_molality[1]=molarity_to_molality(SPcell->interior_molarity[1], "EG");
      SPcell->interior_osmolality=ternary_virial(SPcell->interior_component_molality[0],SPcell->interior_component_molality[1],"NaCl","EG");
  }
  else if(SPcell->simulation_selected==6)
    {// karlsson dmso selected
      SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
      SPcell->interior_component_molality[1]=molarity_to_molality(SPcell->interior_molarity[1], "DMSO");
      SPcell->interior_osmolality=ternary_virial(SPcell->interior_component_molality[0],SPcell->interior_component_molality[1],"NaCl","DMSO");
  }
  else if(SPcell->simulation_selected==7)
    {// karlsson PROH selected
      SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
      SPcell->interior_component_molality[1]=molarity_to_molality(SPcell->interior_molarity[1], "PROH");
      SPcell->interior_osmolality=ternary_virial(SPcell->interior_component_molality[0],SPcell->interior_component_molality[1],"NaCl","PROH");
  }

  else{
    std::cout<<"BAD SIMULATION SELECTED CHECK THE XML FILE!"<<"\n";
  }
  return;
}
void two_parameter_single_step(Cell* pCell, Phenotype &phenotype, double dt) // two_p
{
  Spring_Cell* SPcell=spring_cell_by_pCell_index[pCell->index]; //find the spring_cell containing pCell by its index
  update_exterior_concentrations(SPcell);//update and get exterior concentrations
  // std::cout<<"External "<<SPcell->exterior_osmolality<<"\n";
  // std::cout<<"Internal "<<SPcell->interior_osmolality<<"\n";

  // std::cout<<"pre moles: "<<SPcell->solute_moles<<"\n";
  SPcell->dVw_Osmolality();// calculate dVw/dt
  SPcell->dN_molarity();// calculate dN/dt
  SPcell->two_p_forward_step(dt);// forward step Adams-bashforth 2nd order uses dN and dVw to get new moles and water volume
  SPcell->two_p_update_volume();// update pCell volume and Spring Cell water and solute volume
//   double osmotically_active_volume=SPcell->water_volume+SPcell->solute_volume;
//   // std::cout<<"osmotically_active_volume: "<< osmotically_active_volume<<"\n"; 
//   // std::cout<<"should be the same volume: "<< SPcell->m_my_pCell->phenotype.volume.total-SPcell->solid_volume<<"\n";
//   // std::cout<<"total volume: "<< SPcell->m_my_pCell->phenotype.volume.total<<"\n";
//   for(size_t i=0; i<SPcell->solute_moles.size(); i++)// loop through vector of solute_moles indexed by solute 
//   {
//
//     // std::cout<<"moles: "<< SPcell->solute_moles[i]<<"\n";
//     SPcell->interior_molarity[i]=SPcell->solute_moles[i]/osmotically_active_volume;// calculate molarity
//     // std::cout<<"molarities: "<< SPcell->interior_molarity[i]<<"\n";
//   }
//   if(SPcell->simulation_selected==1)
//     {// only EG selected
//       SPcell->interior_component_molality[0]=molarity_to_molality(SPcell->interior_molarity[0], "NaCl");
//       SPcell->interior_component_molality[1]=molarity_to_molality(SPcell->interior_molarity[1], "EG");
//   //     
//   //    // std::cout<<"interior_component_molality 0: "<< SPcell->interior_component_molality[0]<<"\n";   
//   //    // std::cout<<"interior_component_molality 1: "<< SPcell->interior_component_molality[1]<<"\n";   
//      // double vir=ternary_virial(SPcell->interior_component_molality[0],SPcell->interior_component_molality[1],"NaCl","EG"); 
//      // std::cout<<"Interior OSMOLALITY: "<< vir<<"\n";   
//       SPcell->interior_osmolality=ternary_virial(SPcell->interior_component_molality[0],SPcell->interior_component_molality[1],"NaCl","EG");
// //SPcell->interior_component_molality[0]*1.68+SPcell->interior_component_molality[1];  
//   }
/* code for testing function 
  int counting_test=0;
      // std::cout<<"AREA: "<<SPcell->surface_area<<"\n";
  if(osmotically_active_volume<551060&& SPcell->interior_osmolality>0.5|| SPcell->interior_osmolality<0.01)
  {
    if(counting_test%200==0)
    { 
      // std::cout<<"External: "<<SPcell->exterior_osmolality<<"\n";
      // std::cout<<"Internal: "<<SPcell->interior_osmolality<<"\n";
      // std::cout<<"Total Volume: "<<SPcell->m_my_pCell->phenotype.volume.total<<"\n";
      // std::cout<<"water volume: "<<SPcell->water_volume<<"\n";
      // std::cout<<"solid volume: "<<SPcell->solid_volume<<"\n";
      // std::cout<<"solute volume: "<<SPcell->solute_volume<<"\n";
      // std::cout<<"number of cores: "<<PhysiCell_settings.omp_num_threads<<"\n";
      // std::cout<<"AREA: "<<SPcell->surface_area<<"\n";
    } 
    counting_test++;
  }
*/
  update_interior_concentrations(SPcell);

  uptake(SPcell); 
  return;
}
void uptake(Spring_Cell* SPcell)
{
  double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
  double uptake_voxel_num=SPcell->uptake_voxels.size();
  double reduce_water_volume=SPcell->water_uptake/uptake_voxel_num;
  std::vector<double> specific_volumes;
  std::vector<double> solute_uptake_per_voxel;
  // std::cout<<"Number of voxels: "<< SPcell->uptake_voxels.size()<<"\n"; 
  // std::cout<<"Uptakes: "<<SPcell->solute_uptake<<"in :"<<uptake_voxel_num<<"voxels."<<"\n";
  // std::cout<<"reduce_water_volume: "<<reduce_water_volume<<"in :"<<uptake_voxel_num<<"voxels."<<"\n";
  solute_uptake_per_voxel.resize(SPcell->solute_uptake.size(),0.0);
  specific_volumes.resize(SPcell->solute_uptake.size(),0.0);
  for (size_t i =0; i<SPcell->solute_moles.size();i++) {
    solute_uptake_per_voxel[i]=SPcell->solute_uptake[i]/uptake_voxel_num;
    specific_volumes[i]=SPcell->solute_specific_volume[i];
  } 
  // std::cout<<"uptake per voxel: "<<solute_uptake_per_voxel<<"in :"<<uptake_voxel_num<<"voxels."<<"\n";
  //reduce voxel concentration by cell uptake
  for (size_t i =0; i<SPcell->uptake_voxels.size();i++) {
    uptake_in_one_voxel(SPcell->uptake_voxels[i],reduce_water_volume,solute_uptake_per_voxel, SPcell->solute_specific_volume);
  }
  return;
}
void output_voxel_uptakes(Spring_Cell* SPcell,std::vector<Voxel> output_voxels)
{
  Voxel temp_voxel;
  for(int i=0; i<output_voxels.size(); i++)
  {

    temp_voxel=output_voxels[i];
    std::ofstream ofs;
    ofs.open ("./output/an_uptake.csv", std::ofstream::out | std::ofstream::app);
    std::vector<double> concentration=microenvironment.nearest_density_vector(temp_voxel.center);

    ofs <<PhysiCell_globals.current_time<<", "<<temp_voxel.mesh_index<<", "<<temp_voxel.center[0]<<", "<<temp_voxel.center[1]<<", "<<temp_voxel.center[2]<<", "<<concentration[1]<<"\n";
    ofs.close();
  }
    return;
}

void uptake_in_one_voxel(int voxel, double water_uptake_per_voxel, std::vector<double> solute_uptake_per_voxel, std::vector<double> specific_volumes )
{

  double voxel_volume=default_microenvironment_options.dx*default_microenvironment_options.dy*default_microenvironment_options.dz;
  std::vector<double> temp_density_vec=microenvironment.density_vector(voxel);
  std::vector <double> moles_in_voxel;
  moles_in_voxel.resize(solute_uptake_per_voxel.size(),0.0);
  for(int i=0; i<solute_uptake_per_voxel.size(); i++)
  {
    moles_in_voxel[i]=(temp_density_vec[i]*voxel_volume);//fmole/um^3* um^3 
  }
  // std::ofstream ofs;
  // ofs.open("./output/uptake_in_one_voxel.csv", std::ofstream::out|std::ofstream::app);
  // ofs <<PhysiCell_globals.current_time<<", "<<microenvironment.voxels(voxel).mesh_index<<", "<<microenvironment.voxels(voxel).center[0]<<", "<<microenvironment.voxels(voxel).center[1]<<", "<<microenvironment.voxels(voxel).center[2]<<", "<<microenvironment.density_vector(voxel)[1]<<", ";
   
  std::vector <double> new_moles=moles_in_voxel-solute_uptake_per_voxel;
  double new_water_volume=voxel_volume-water_uptake_per_voxel;
  // ofs<< new_moles[1]<<", "<<moles_in_voxel[1]<<", "<<solute_uptake_per_voxel[1]<<", "<<new_water_volume<<", "<< voxel_volume<<", "<<water_uptake_per_voxel<<", "; 
  // std::cout<<"new water volume: "<< new_water_volume<<"\n";
  for(int i=0; i<solute_uptake_per_voxel.size();i++)
  {
    // std::cout<<"old old: "<<microenvironment.density_vector(voxel)[i]<<"\n";
    // std::cout<<"old: "<<temp_density_vec[i]<<"\n";
    microenvironment.density_vector(voxel)[i]=new_moles[i]/new_water_volume;
    // std::cout<<"new: "<<microenvironment.density_vector(voxel)[i]<<"\n";
  }
  // ofs<<microenvironment.density_vector(voxel)[1]<<"\n ";

   // ofs.close();

  return;
}
[[deprecated]] void rasterize_my_uptake(Cell* pCell, double solute_index)//old version
{
  std::cout<<"WARNING!!! OLD VERSION OF UPTAKING ARE YOU SURE YOU WANT TO USE THIS??"<<"\n";
  // double total_voxel_volume=0.0;
  // double uptake_per_voxel=0.0;
  // //test double
  // double custom_cell_uptake=2.22;//kg/um^3 dt already factored in at calculate 2p
  // std::vector<int> interior_voxels=get_interior_voxels(pCell);
  // //get my interior voxels
  // if(interior_voxels.size()==1)
  // {
  //   //smaller than voxel
  //   total_voxel_volume=pCell->get_container()->underlying_mesh.dV;
  //   uptake_per_voxel=custom_cell_uptake*total_voxel_volume/(pCell->phenotype.volume.total);
  // }
  // else
  // {
  //   //if multivoxel divide 
  //   total_voxel_volume=interior_voxels.size()*pCell->get_container()->underlying_mesh.dV;
  //   uptake_per_voxel=custom_cell_uptake*total_voxel_volume/(pCell->phenotype.volume.total);
  // }
  // //for my voxels uptake (secretion=-uptake)
  // #pragma omp critical
  // {
  //   for (size_t i = 0; i < interior_voxels.size(); i++)
  //   {
  //     microenvironment.density_vector(interior_voxels[i])[solute_index]+=uptake_per_voxel;
  //     
  //   }
  // }
  return;
  
}
//as of May 23 2023 mechanics voxels seem to equal microenvironment voxels
std::vector<int> basement_membrane_voxels={};//external list for checking if a cell is intersecting the basement membrane

void get_basement_membrane_intersection(std::vector <double> &center_point, double &radius, std::vector<int> *basement_voxels )
{
  double voxel_size=default_microenvironment_options.dx;
  //double inner_radius=radius-voxel_size;
  double outer_radius=radius+(voxel_size);//~2 voxels thick
  //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector <int> bounding_voxels={};
    spherical_bounding_box(center_point,outer_radius,&bounding_voxels);
    std::vector<int> intersecting_voxels={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=microenvironment.voxels(bounding_voxels[i]).center;
      // std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
      std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      get_voxel_corners(test_voxel_center,test_corners);
      int sum=0;
      // #pragma omp private(sum,exterior_voxels)
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
      basement_voxels->assign(intersecting_voxels.begin(), intersecting_voxels.end());
    }

    return;
}

void initialize_basement_membrane(std::vector <double> center_point, double radius)
{
  //calculate voxels that represent the BM
  get_basement_membrane_intersection(center_point,radius,&basement_membrane_voxels);
  return;
}
// void variable_moore_neighborhood(Cell* pCell, double radius)
// { //this version is not currently used
  //works just like bounding diffusion box to get mechnical voxels around pCell with any radius
  //returns voxel by index
  //NOTE!!! currently PhysiCells microenvironment voxels are used for both mechanics and diffusion
 /* 
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
 */ 
  // return;
// }
void max_interaction_variable_moore_neighborhood(Cell* pCell, double &maximum_interaction_distance, std::vector<int> *voxel_neighboorhood)
{
  double radius=pCell->phenotype.geometry.radius+maximum_interaction_distance;
  //works just like bounding diffusion box to get mechnical voxels around pCell with any radius
  //returns voxel by index
  //NOTE!!! currently PhysiCells microenvironment voxels are used for both mechanics and diffusion
  if(radius<=pCell->phenotype.geometry.radius)
  {
    std::cout<< "radius: "<< radius<<"\n";
    std::cout<< "WARNING!! Are you sure you want to use a moore neighborhood smaller than the cell?"<<"\n";
  }

  
  //bounding box of voxels in the microenvironment
 //for speed get center voxel figure out x,y,z offset and only check local voxels
    std::vector<double> center_voxels_center= microenvironment.nearest_voxel(pCell->position).center; // get center voxel
    std::vector<double> offset=center_voxels_center-pCell->position;
    //figure out the corners of my bounding box of voxels to search
    std::vector<double> lower_point(3,0.0);
    std::vector<double> upper_point(3,0.0);
    for (size_t i = 0; i < 3; i++)
    {
      lower_point[i]= center_voxels_center[i]-offset[i]-ceil(radius);
      upper_point[i]= center_voxels_center[i]+offset[i]+ceil(radius);
 
    }

  std::vector<int> bounding_box_by_index={};
  general_voxel_bounding_box(voxel_neighboorhood,lower_point,upper_point,default_microenvironment_options.dx,(pCell->get_container()->underlying_mesh));
  //std::cout<<"Bounding box"<<bounding_box_by_index<<"\n";
  
  return;
}
void general_voxel_bounding_box(std::vector<int> *return_bounding_box,std::vector<double> &starting_position, std::vector <double>&ending_position, double &voxel_length, BioFVM::Cartesian_Mesh &a_mesh)
{
  //function to find a box of voxels in a cartesian mesh given the lowest corner to the highest, should work in 2D and 3D
  std::vector<int> bounding_box_by_index={-1};
  int starting_voxel_index=a_mesh.nearest_voxel_index(starting_position);// should work with pCell-get_container or microenvironment
  int ending_voxel_index=a_mesh.nearest_voxel_index(ending_position);
  std::vector<double> starting_voxel_center=a_mesh.nearest_voxel(starting_position).center;
  std::vector<double> ending_voxel_center=a_mesh.nearest_voxel(ending_position).center;
  if(starting_voxel_index==ending_voxel_index)
  {
    //contained entirely in 1 voxel my exterior is my voxel
    bounding_box_by_index[0]=starting_voxel_index;
  }
  else{
   std::vector <int> num_voxels_across(3,0.0);
   //loop through voxel centers, get index and store it
   std::vector <int> dimensions={0,0,0};
   for (size_t i = 0; i < 3; i++)
   {
    num_voxels_across[i]=int((ceil(ending_voxel_center[i])-floor(starting_voxel_center[i]))/voxel_length); 
    if(floor(starting_voxel_center[i])==ceil(ending_voxel_center[i]))
    {
      dimensions[i]=1;//if a dimension is flat mark it with 1
    }
   }
   
    int upper_z=num_voxels_across[2]+dimensions[2];
    int upper_y=num_voxels_across[1]+dimensions[1];
    int upper_x=num_voxels_across[0]+dimensions[0];
  
    #pragma omp critical
    { //to avoid race conditions on resize if multiple cells are resized, unsure if this is an issue 
      bounding_box_by_index.resize(upper_x*upper_y*upper_z);
    }
    std::vector<double> position(3,0.0);
    int count=0;
    for(int zi=0; zi<upper_z; zi++)
    { 
      position[2]=double(starting_voxel_center[2]+(double(zi-dimensions[2])*voxel_length));
      for (int yi = 0; yi <upper_y ; yi++)
      {
        position[1]=double(starting_voxel_center[1]+(double(yi-dimensions[1])*voxel_length));
        for (int xi = 0; xi <upper_x ; xi++)
        {
            position[0]=double(starting_voxel_center[0]+(double(xi-dimensions[0])*voxel_length));
            bounding_box_by_index[count]=a_mesh.nearest_voxel_index(position);
            count++;
        }
      }
    }
  } 
  #pragma omp critical
  {
    return_bounding_box->assign(bounding_box_by_index.begin(),bounding_box_by_index.end());
  }
  return;
}
/* function just connects the spring cells doesn't check distances, this might be better places in springs.cpp     */
void connect_spring_cells(Spring_Cell* SpCell_1, Spring_Cell* SpCell_2)
{
  // if(SpCell_1->m_my_pCell->type==1 || SpCell_2->m_my_pCell->type==1){
    // std::cout<<"Connections made from "<< SpCell_1->m_my_pCell->type_name<<" to "<< SpCell_2->m_my_pCell->type_name<<"\n";
  // }
  double spring_length=signed_distance_between_membranes(SpCell_1->m_my_pCell,SpCell_2->m_my_pCell);
  SpCell_1->add_spring(SpCell_2,spring_length);
  SpCell_2->add_spring(SpCell_1,spring_length);
  return;
}
void cells_in_neighborhood(Cell* pCell, double &maximum_interaction_distance, std::vector<Cell*> *neighboors)
{
  std::vector <int> voxels_to_search_indicies={};
  max_interaction_variable_moore_neighborhood(pCell,maximum_interaction_distance,&voxels_to_search_indicies);
  std::vector <Cell*> cells_found_private{};
  std::vector <Cell*> cells_found{};
  // #pragma omp private(cells_found_private)
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
    neighboors->assign(cells_found.begin(),cells_found.end()); 
  }
  return;
}
/* function to add the forces from all the connected cells, this needs to be combined with non-connected cell forces in a thread safe way*/
void custom_add_potentials(Spring_Cell* SpCell, double dt)
{
  std::vector <double> sum_of_forces_private(3,0.0);
  std::vector <double> sum_of_forces(3,0.0);
  // #pragma omp for private(sum_of_forces_private)
  for (size_t i = 0; i < SpCell->m_springs.size(); i++)//loop through my springs
  {
    double spring_stretch=signed_distance_between_membranes(SpCell->m_my_pCell,SpCell->m_springs[i]->m_pNeighbor);
    //spring_stretch 
    std::vector <double> force(3,0.0);
    //std::vector<double> force_direction= (SpCell->m_springs[i]->m_pNeighbor->position)-(SpCell->m_my_pCell->position);
      Hookes_law_force_magnitude(SpCell->m_my_pCell->position,SpCell->m_springs[i]->m_pNeighbor->position,SpCell->m_springs[i]->m_length,spring_stretch,SpCell->m_my_pCell->custom_data["cell_k"], &force);
    //if any component of the force is less than 1.0*10^-16 set that equal to 0 to avoid floating decimal issues
    sum_of_forces_private+= force;
    if(SpCell->m_my_pCell->index==5)
    {
      // std::cout<<"spring stretch: "<<spring_stretch<<"\n";
       // std::cout<<"equilibrium_length: "<<SpCell->m_springs[i]->m_length<<"\n";
       // std::cout<<"difference: "<< SpCell->m_springs[i]->m_length-spring_stretch<<"\n";
       // std::cout<<"FORCE: "<<force<<"\n";
    //std::cout<<"force_direction "<<force_direction<<"\n";
    //std::cout<<"force_direction "<<force_direction<<"\n";
    // f/m=a=dv/dt
    }
  }
  
  #pragma omp critical
  {
    if(SpCell->m_my_pCell->index==5)
    {
      // std::cout<<"previous_velocity: "<<SpCell->previous_velocity<<std::endl;
     // std::cout<<"sum of forces: "<<sum_of_forces_private<<std::endl;
      // std::cout<<"size of neighbor: "<<SpCell->m_springs.size()<<"\n";
      // std::cout<<"first neighbor: "<<SpCell->m_springs[0]->m_pNeighbor<<"\n";
      // std::cout<<"equilibrium_length: "<<SpCell->m_springs[0]->m_length<<"\n";
    }
      // sum_of_forces.insert(sum_of_forces.end(),sum_of_forces_private.begin(),sum_of_forces_private.end());
    sum_of_forces[0]=(dt*sum_of_forces_private[0]/SpCell->m_my_pCell->phenotype.volume.total); //f/m=a=(dv/dt) in x
    sum_of_forces[1]=(dt*sum_of_forces_private[1]/SpCell->m_my_pCell->phenotype.volume.total); //f/m=a=(dv/dt) in y
    sum_of_forces[2]=(dt*sum_of_forces_private[2]/SpCell->m_my_pCell->phenotype.volume.total); //f/m=a=(dv/dt) in z 
    
    for(int i=0; i<3; i++)//if force components are tiny, round that component to 0 to avoid floating point issues
    {
      if(std::abs(sum_of_forces[i])<=1e-16)
      {
        sum_of_forces[i]=0;
      }
    }  
    SpCell->m_my_pCell->velocity=SpCell->previous_velocity+sum_of_forces;//Forward_Euler
    SpCell->previous_velocity=SpCell->m_my_pCell->velocity;
    if(SpCell->m_my_pCell->index==5)
    {
      // std::cout<<"acceleration: "<<sum_of_forces<<std::endl;
      // std::cout<<"position: "<<SpCell->m_my_pCell->position<<std::endl;
      // std::cout<<"velocity: "<<SpCell->m_my_pCell->velocity<<std::endl;
    }
  }
  // std::cout<< "some velocities: "<< SpCell->m_my_pCell->velocity<<"\n";
  return;
}
void custom_add_potentials_for_pCells(Cell* pCell, double dt)
{
  Spring_Cell* current_spring_cell=spring_cell_by_pCell_index[pCell->index];
  custom_add_potentials(current_spring_cell, dt);
  return;
}
void spring_cells_in_neighborhood(Spring_Cell* SpCell, double maximum_interaction_distance,std::vector <Spring_Cell*> *neighboors )
{
  std::vector <int> voxels_to_search_indicies={};
  max_interaction_variable_moore_neighborhood(SpCell->m_my_pCell,maximum_interaction_distance,&voxels_to_search_indicies);
  std::vector <Cell*> cells_found_private{};
  std::vector <Cell*> cells_found{};
  std::vector <Spring_Cell*> spring_cells_found{};
  // #pragma omp private(cells_found_private)
  for (size_t i = 0; i < voxels_to_search_indicies.size(); i++)
  {
    for (size_t j = 0; j < SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]].size(); j++) 
    {
      if(SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]][j]!=SpCell->m_my_pCell)
      {
        Cell* tempCell=(SpCell->m_my_pCell->get_container()->agent_grid[voxels_to_search_indicies[i]][j]);
        // if(norm(tempCell->position-SpCell->m_my_pCell->position)<=(tempCell->phenotype.geometry.radius+SpCell->m_my_pCell->phenotype.geometry.radius))
        // {
          cells_found_private.push_back(tempCell);
        // }
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
  #pragma omp critical
  {
    neighboors->assign(spring_cells_found.begin(),spring_cells_found.end());
  }
  return;
}
void initialize_spring_connections()//also initialize mechanics and 2p vectors which probably should be done elsewhere
{
  // std::cout<<"initializing connections!"<<"\n"; 
  //using modified moore neighborhood search for "spring cell" neighbors and connect them with spring objects
  for(size_t i=0;i<all_spring_cells.size();i++) //loop through all spring cells
  {
    //std::cout<<"Second!!"<<std::endl;
    Spring_Cell* SpCell=all_spring_cells[i];
    SpCell->set_2p_initial_conditions();
    SpCell->initialize_mechanics(); 
    //std::cout<<" spcell "<<SpCell<<std::endl;
    std::vector<Spring_Cell*> SpNeighbors={};
    double maximum_interaction_distance=SpCell->m_my_pCell->custom_data["neighborhood_radius"]+SpCell->m_my_pCell->phenotype.geometry.radius;
    spring_cells_in_neighborhood(SpCell,maximum_interaction_distance,&SpNeighbors);
   
    for (size_t j = 0; j < SpNeighbors.size(); j++)
    {
      //use center to center to establish connections to deal with oocyte-to-granulosa and granulosa-granulosa size difference
      if (norm(SpCell->m_my_pCell->position-SpNeighbors[j]->m_my_pCell->position)<=maximum_interaction_distance) {
        connect_spring_cells(SpCell,SpNeighbors[j]);
      }
    }
  }
  return;
}
void create_spring_cell_from_cell(Cell* pCell, Phenotype& phenotype)
{ //not used at the moment
  
  return;
}
/*convert from molarity to molality using the polynomial Ax^3+Bx^2+Cx+D, D is 0 where they intercept 
 * from python files in virial folder using CRC data located there for 20 degrees C --assumed ~= at 23 C
 * To convert from  molarity  to  molality  for  EG  the coefficients are:  [ 0.01277023 -0.02080002  1.13628309  0.        ]
 * To convert from  molarity  to  molality  for  GLY  the coefficients are:  [0.00665793 0.06926148 1.00285049 0.        ]
 * To convert from  molarity  to  molality  for  NaCl  the coefficients are:  [ 5.80771883e-05 -1.93492494e-02  1.00003795e+00  0.00000000e+00]
 * */
std::unordered_map<std::string,double> molal_conversion_coeff_A={{"EG",0.01277023},{"GLY",0.00665793},{"NaCl",5.80771883e-05}};
std::unordered_map<std::string,double> molal_conversion_coeff_B={{"EG",-0.02080002},{"GLY",0.06926148},{"NaCl",-1.93492494e-02}};
std::unordered_map<std::string,double> molal_conversion_coeff_C={{"EG",1.13628309},{"GLY",1.00285049},{"NaCl",1.00003795}};
std::unordered_map<std::string,double> molal_conversion_coeff_D={{"EG",0.0},{"GLY",0.0},{"NaCl",0.0}};

double molarity_to_molality(double molarity, std::string component_name)
{
  double molality=0.0;
  molality= molal_conversion_coeff_A[component_name]*molarity*molarity*molarity+molal_conversion_coeff_B[component_name]*molarity*molarity+molal_conversion_coeff_C[component_name]*molarity+molal_conversion_coeff_D[component_name];
  return molality;
}
/*
 * The virial osmotic coefficients for the cubic virial osmotic equation and molar_mass
 * */
std::unordered_map<std::string,double> virial_coeff_B={{"EG",0.037},{"GLY",0.023},{"NaCl",0.044}};
std::unordered_map<std::string,double> virial_coeff_C={{"EG",-0.001},{"GLY",0.0},{"NaCl",0.0}};
std::unordered_map<std::string,double> kdiss={{"EG",1.0},{"GLY",1.0},{"NaCl",1.678}};
std::unordered_map<std::string,double> molar_mass={{"EG",62.07},{"GLY",92.09},{"NaCl",58.44}};
double binary_virial(double molality, std::string component_name)
{//osmolality=mi+B_i*m_i^2+C_i*m_i^3
  double osmolality=molality+virial_coeff_B[component_name]*molality*molality+virial_coeff_C[component_name]*molality*molality*molality;
  return osmolality;
}

/* The following equation computes the osmolality of a ternary solution with m ininitial molalilities using B and C coefficients,with kdiss for ions*/


double ternary_virial(double molality_1, double molality_2, std::string component_1, std::string component_2)
{
    //compute virial equation uses vectors for future expansion
  double osmolality=0.0;
  // #pragma omp reduction(+:osmolality) 
  // avoid thread issues by mearly computing all the terms
  std::vector <double> m={molality_1*kdiss[component_1],molality_2*kdiss[component_2]};
  std::vector <double> B={virial_coeff_B[component_1],virial_coeff_B[component_2]};
  std::vector <double> C={virial_coeff_C[component_1],virial_coeff_C[component_2]};
  osmolality=m[0]+m[1]+B[0]*m[0]*m[0]+B[1]*m[1]*m[1]+(B[0]+B[1])*m[0]*m[1]+C[0]*m[0]*m[0]*m[0]+C[1]*m[1]*m[1]*m[1]+3*(std::pow((C[0]*C[0]*C[1]),(1.0/3.0))*m[1]*m[1]*m[2])+3*(std::pow((C[0]*C[1]*C[1]),(1.0/3.0))*m[1]*m[2]*m[2]);
  // for (size_t i = 0; i < m.size(); i++)
  // {
  //   osmolality+=m[i];
  //   for (size_t j = 0; j < B.size(); j++)
  //   {
  //       osmolality+=((B[i]+B[j])/2)*m[i]*m[j];
  //       for (size_t k = 0; k < C.size(); k++)
  //       {
  //         std::cout<<"Cs"<< std::abs(C[i]*C[j]*C[k])<<"\n";
  //         std::cout<<"test pow "<< (std::pow(std::abs(C[i]*C[j]*C[k]),(0.333333)))<<"\n";
  //         std::cout<<"test ms: "<<m[i]*m[j]*m[k]<<"\n";
  //           osmolality+=(std::pow((C[i]*C[j]*C[k]),(0.333333)))*m[i]*m[j]*m[k];
  //       }
  //       
  //   }
  // }
  // std::cout<<"Osmole "<< osmolality<< "\n";
    return osmolality;
}
void update_exterior_concentrations(Spring_Cell* SPcell)
{
  // this function passes the exterior concentrations to the Spring_Cell and specifies what solutes are being simulated
  // std::cout<<"Simulation selected"<< SPcell->simulation_selected<<"\n";
  if(SPcell->simulation_selected==1)
  {
    // std::cout<<"Salt molarity: "<< concentration_at_boundary(SPcell->m_my_pCell,0)<<"\n";
    SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);// read in the averaged exterior molarity in the voxels along the cell boundary
    SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); //convert the molarity into molality using the best fit polynomial
    
    // std::cout<<"Eg molarity: "<< concentration_at_boundary(SPcell->m_my_pCell,1)<<"\n";
    SPcell->exterior_molarity[1]=concentration_at_boundary(SPcell->m_my_pCell,1);
    
    // std::cout<<"Test molal: "<< molarity_to_molality((SPcell->exterior_molarity[1]),"EG")<<"\n";
    SPcell->exterior_component_molality[1]=molarity_to_molality(SPcell->exterior_molarity[1],"EG");
    
   // std::cout<<"Test virial: "<< ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","EG")<<"\n";
    SPcell->exterior_osmolality=ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","EG");
    // SPcell->exterior_component_molality[0]*1.68+SPcell->exterior_component_molality[1];//ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","EG");
  }
  else if(SPcell->simulation_selected==2)
  { 
    SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);
    SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); 
    SPcell->exterior_molarity[1]=concentration_at_boundary(SPcell->m_my_pCell,1);
    SPcell->exterior_component_molality[1]=molarity_to_molality(SPcell->exterior_molarity[1],"GLY"); 
    SPcell->exterior_osmolality=ternary_virial(SPcell->exterior_component_molality[0],SPcell->exterior_component_molality[1],"NaCl","GLY");
  }
  else if(SPcell->simulation_selected==3)
  { 
    SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);
    SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); 
    SPcell->exterior_osmolality=binary_virial(SPcell->exterior_component_molality[0],"NaCl");
  }
  else if(SPcell->simulation_selected==4)
  { 
    SPcell->exterior_molarity[0]=concentration_at_boundary(SPcell->m_my_pCell,0);
    SPcell->exterior_component_molality[0]=molarity_to_molality(SPcell->exterior_molarity[0],"NaCl"); 
    SPcell->exterior_molarity[1]=concentration_at_boundary(SPcell->m_my_pCell,1);
    SPcell->exterior_component_molality[1]=molarity_to_molality(SPcell->exterior_molarity[1],"EG"); 
    SPcell->exterior_molarity[2]=concentration_at_boundary(SPcell->m_my_pCell,2);
    SPcell->exterior_component_molality[2]=molarity_to_molality(SPcell->exterior_molarity[2],"GLY");
    SPcell->exterior_osmolality=binary_virial(SPcell->exterior_component_molality[0],"NaCl")+ternary_virial(SPcell->exterior_component_molality[1],SPcell->exterior_component_molality[2],"EG","GLY");
  }
  return;
}
void initialize_spring_cells()//make all cells spring cells
{
  // std::cout<<"Initializing Spring Cells"<<"\n";
  //encapsulate all cells in the super class spring cell
  //find initial spring lengths and connect neighbors
  //set up basement membrane and connected exterior cells
  std::vector <double> basement_membrane_center={0.0, 0.0, 0.0};
  // double test_radius=70;
  double basement_membrane_radius=95;
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
  int count_cells=0;
  std::cout<<"basement voxels size: "<< basement_membrane_voxels.size()<<"\n";
  //TODO: should be thread safe but verify
  for (int i=0; i< basement_membrane_voxels.size();i++)
  {
    for( int j=0; j<(all_spring_cells).size();j++)
    {
      Spring_Cell* SpCell=(all_spring_cells)[j];
      if(is_in_voxel(SpCell->m_my_pCell,basement_membrane_voxels[i]))
      {

        double length=basement_membrane_radius-(norm(SpCell->m_my_pCell->position-basement_membrane_center)+SpCell->m_my_pCell->phenotype.geometry.radius);
        SpCell->is_basement_connected=true;
        SpCell->basement_length=length;
        count_cells++;
      }
    }
  }
  std::cout<<"Cells connected to basement: "<< count_cells<<"\n";
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
    if (pCell->get_container()->agent_grid[voxel_index][i] == pCell) 
    {

      // std::cout<<"connecting cells to basement membrane: "<< i <<"\n";
      value=true;
    }
  }
  return value;
}
void basement_membrane_mechanics(Spring_Cell* SpCell, double basement_membrane_radius, std::vector <double> basement_membrane_center, double dt)
{
  // function is designed so you could have a changing basement membrane but we set it static
  //basement membrane radius is the center of the membrane which is 2-3 voxels thick
  std::vector<double>force_direction=normalize(SpCell->m_my_pCell->position);//force along normal
  double current_length=basement_membrane_radius-(norm(SpCell->m_my_pCell->position-basement_membrane_center)+SpCell->m_my_pCell->phenotype.geometry.radius);
  if(SpCell->is_basement_connected==true)//normal case don't need to calculated edge cases if only this
  {
      std::vector<double> force(3,0.0);

      // std::abs(current_length);
      Hookes_law_force_magnitude(basement_membrane_center,SpCell->m_my_pCell->position,SpCell->basement_length,current_length,SpCell->m_my_pCell->custom_data["basement_k"],&force);
      #pragma omp critical
      {
        
        std::vector<double> acceleration(3,0.0);
        acceleration[0]=force[0]/SpCell->m_my_pCell->phenotype.volume.total;//a_x=F_x/m
        acceleration[1]=force[1]/SpCell->m_my_pCell->phenotype.volume.total;//a_y=F_y/m
        acceleration[2]=force[2]/SpCell->m_my_pCell->phenotype.volume.total;//a_z=F_z/m
        // if(SpCell->m_my_pCell->index==5)
      // { std::cout<<"position: "<<SpCell->m_my_pCell->position<<"acceleration: "<< acceleration<<"\n";
        // std::cout<<"basement length: "<<SpCell->basement_length<<"\n";
        // std::cout<<"current length: "<<current_length<<"\n";
      // }
        SpCell->m_my_pCell->velocity=SpCell->previous_velocity+(dt*acceleration);//overloaded double*vector double
        SpCell->previous_velocity=SpCell->m_my_pCell->velocity;
      }
      
  }
  //edge case where cells that don't initially neighbor the BM are pushed into the BM
  else if(SpCell->is_basement_connected!=true && current_length<0) 
  {
    // std::cout<<"Basement: "<<current_length<<"\n";
    if(std::abs(current_length)>(basement_membrane_radius+default_microenvironment_options.dx))//cell membrane to BM center is greater than a whole voxel
    {
     //updated to use outer_radius from basement membrane function
      std::cout<< "WARNING!! CELL PASSED THROUGH BASEMENT MEMBRANE!"<<std::endl;
      SpCell->is_outside=true; 
    }
    //basement membrane acts on all cells that get within the basement_membrane_voxels
    for (int i=0; i< basement_membrane_voxels.size();i++)
    {
      if(is_in_voxel(SpCell->m_my_pCell,basement_membrane_voxels[i]))
      {
        std::vector<double> force(3,0.0);
        //current_length is from edge of cell to center of BM
        double equilibrium_length=0.0;
        Hookes_law_force_magnitude(basement_membrane_center,SpCell->m_my_pCell->position,equilibrium_length,current_length,SpCell->m_my_pCell->custom_data["basement_k"],&force);
        std::vector<double> acceleration(3,0.0);
        acceleration[0]=force[0]/SpCell->m_my_pCell->phenotype.volume.total;//a_x=F_x/m
        acceleration[1]=force[1]/SpCell->m_my_pCell->phenotype.volume.total;//a_y=F_y/m
        acceleration[2]=force[2]/SpCell->m_my_pCell->phenotype.volume.total;//a_z=F_z/m
        // std::cout<<"membrane acceleration: "<< acceleration<<"\n";
        SpCell->m_my_pCell->velocity=SpCell->previous_velocity+(dt*acceleration);//overloaded double*vector double
        SpCell->previous_velocity=SpCell->m_my_pCell->velocity;

      } 
    }
    //
  }
    return;
} 
void basement_membrane_mechanics(Cell* pCell, double basement_membrane_radius) //for future versions in follicle distance test will be faster currently not working
{
  // currently not used/ use Spring Cell version above
  // std::vector <int>  test_voxel_indicies{};
  // //mechanics and diffusion are using the same mesh as of may 23 2023
  // std::vector <int> interior_indicies=get_interior_voxels(pCell);
  // std::vector <int> exterior_indicies=get_exterior_voxels(pCell);
  // if(exterior_indicies.size()!=1)//allows for testing if multivoxel cell intersects the basement membrane
  // {
  //   #pragma omp critical
  //   {
  //     test_voxel_indicies.resize(interior_indicies.size()+exterior_indicies.size());
  //   }
  //   int exterior_count=0;
  //   for(int i=0;i<test_voxel_indicies.size();i++)
  //   {
  //     if(i<interior_indicies.size())
  //     {
  //       test_voxel_indicies[i]=interior_indicies[i];
  //     }
  //     else
  //     {
  //       test_voxel_indicies[i]=exterior_indicies[exterior_count];
  //       exterior_count++;
  //     }
  //   }
  // }
  // for(int j=0; j<test_voxel_indicies.size();j++)
  // {
  //   for(int k=0; k<basement_membrane_voxels.size();k++)
  //   {
  //     if(test_voxel_indicies[j]==basement_membrane_voxels[k])
  //     {
  //       
  //     }
  //   }
  // }
  return;
}

void spherical_bounding_box(std::vector <double> &center_point, double &radius,std::vector <int> *return_bounding_box)
{
  //gets voxels of a bounding box around an arbitrary sphere
  //NOTE!!! currently PhysiCells microenvironment voxels are both mechanics and diffusion
  std::vector<int> bounding_box_by_index={};
  //bounding box of voxels in the microenvironment
  //for speed get center voxel figure out x,y,z offset and only check local voxels
  std::vector<double> center_voxels_center=microenvironment.nearest_voxel(center_point).center;
  std::vector<double> offset=center_voxels_center-center_point;
  //figure out the corners of my bounding box of voxels to search
  std::vector<double> lower_point(3,0.0);
  std::vector<double> upper_point(3,0.0);
  for (size_t i = 0; i < 3; i++)
  {
    lower_point[i]= center_voxels_center[i]-offset[i]-(ceil(radius));
    upper_point[i]= center_voxels_center[i]+offset[i]+ceil(radius);
  }
  std::vector<double> starting_voxel_center=microenvironment.nearest_voxel(lower_point).center;
  std::vector<double> ending_voxel_center=microenvironment.nearest_voxel(upper_point).center;
  general_voxel_bounding_box(return_bounding_box,starting_voxel_center,ending_voxel_center,default_microenvironment_options.dx,microenvironment.mesh);
  return;
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
{ //not currently used
  // int vox=microenvironment.find_density_index(solute_name);
  double moles=0.0;
  // //get water volume of voxel
  // std::vector <double> solute_specific_volumes={0.0,1.1,2.2,3.3,4.4};//grams/mole//TODO put in real values and move to config
  //
  // double total_volume=microenvironment.mesh.voxels[voxel_id].volume;
  // double sum=0.0;//get volume taken up by all solutes in a voxel
  // #pragma omp private(sum)
  // for (size_t i = 0; i < microenvironment.density_vector(vox).size(); i++)
  // {
  //   if(solute_specific_volumes.size()!=microenvironment.density_vector(vox).size())
  //   {
  //     std::cout<<"MISSING/EXTRA SPECIFIC VOLUME!"<<"\n";
  //     return moles;
  //   }
  //   //sum+=microenvironment.density_vector(vox)[i]*solute_molar_masses[vox];//mole/kg/m^3/kg*10^18um^3/m^3
  // }
  // 

  return moles;
}
void output_all_voxels_concentrations()
{
  Voxel temp_voxel;
  for(int i=0; i<microenvironment.number_of_voxels(); i++)
  {
    temp_voxel=microenvironment.mesh.voxels[i];
    std::vector <unsigned int> output_indicies=microenvironment.cartesian_indices(i);
    std::ofstream ofs;
    ofs.open ("./output/concentrations.csv", std::ofstream::out | std::ofstream::app);
    std::vector<double> concentration=microenvironment.nearest_density_vector(temp_voxel.center);

    ofs <<PhysiCell_globals.current_time<<", "<<temp_voxel.mesh_index<<", "<<temp_voxel.center[0]<<", "<<temp_voxel.center[1]<<", "<<temp_voxel.center[2]<<", "<<concentration[1]<<"\n";
    ofs.close();
  }
    return;
}
void output_all_voxels_as_cartesian_index()
{
  for(int i=0; i<microenvironment.number_of_voxels(); i++)
  {
    std::vector <unsigned int> output_indicies=microenvironment.cartesian_indices(i);
    std::ofstream ofs;
    ofs.open ("./output/voxel_indicies.csv", std::ofstream::out | std::ofstream::app);

    ofs <<i<<", "<<output_indicies[0]<<", "<<output_indicies[1]<<", "<<output_indicies[2]<<"\n";
    ofs.close();
  }
    return;
}
void output_cell_and_voxels(std::vector <int> voxel_list, Cell* pCell)
{
  #pragma omp critical
  {
  std::vector <std::vector <std::vector<bool>>> matplotlib_values;
  matplotlib_values.resize(40,std::vector<std::vector<bool> >(40,std::vector<bool>(40)));
  // matplotlib_values.resize(40,std::vector<bool> (40));
  for(int i=0; i<40; i++)//x
  {
    for(int j=0; j<40; j++)//y
    {
      for(int k=0; k<40; k++)
      {
        matplotlib_values[i][j].push_back(false);
      }
    }
  }
  std::ofstream ofs;
  ofs.open ("./output/pcell_and_indicies.csv", std::ofstream::out | std::ofstream::app);
  ofs <<pCell<<", "<<pCell->phenotype.geometry.radius<<", "<<pCell->position<<"\n"<<"\n";
  ofs.close();
  for(int i=0; i<40; i++)//x
  {
    for(int j=0; j<40; j++)//y
    {
      for(int k=0; k<40;k++)//z
      {
        for(int m=0; m<voxel_list.size();m++)
        {

          if(microenvironment.voxel_index(i,j,k)==voxel_list[m])
          {
              matplotlib_values[i][j][k]=true;
          }
        }
      }
    }
  }

  ofs.open ("./output/pcell_and_indicies.csv", std::ofstream::out | std::ofstream::app);

  ofs<<"vox_array=np.array([";
  for(int i=0; i<40; i++)//x
  {
    ofs<<"[";
    for(int j=0; j<40; j++)//y
    {
      ofs<<"[";
      for(int k=0; k<40;k++)//z
      {
          ofs<<matplotlib_values[i][j][k]<<", ";
      }
      ofs<<"],";
    }
    ofs<<"],";
  }
  ofs<<"])";
  ofs.close();
  }
  return;
}
