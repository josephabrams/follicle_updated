#include "./multivoxel_utilities.h"

/*Function gets a voxel bounding box, with offset based on the cells' position relative to the center of a voxel; not gauranteed to be the least bounding box, but will contain the cell and it's immediate neighborhood - there is currently no check to see if all voxels are valid.
 --
 Parameters: 
  - return_bounding_box - std::vector<int>* that holds the voxels contained in the bounding box by their index
  - starting_position - vector point that corresponds to the "smallest x,y,z, the "bottom, leftmost, frontmost" corner of the cube
  - ending_position - position of top, rightmost, backmost corner of the bounding box
  - voxel length is the length of one dimension of the regular cartesian cube forming the mesh
  - a_mesh can be either the microenvironment or mechanics voxel mesh, or some custom mesh of the BioFVM::Cartesian_Mesh type
  
  Note: it may be possible to switch this function to use std::set in the future for a considerable speed up
 */

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
    //bounding_box is 1 voxel
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
    { //to avoid race conditions on resize if multiple cells are resized, unsure if this is a real issue but better safe than sorry for now
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
  #pragma omp critical // can be parrallelized 
  {
    return_bounding_box->assign(bounding_box_by_index.begin(),bounding_box_by_index.end());
  }
  return;
}
/*
  Function that gets the bounding box of microevironment voxels for a cell and returns it in bounding_box_by_index, a vector of voxels by their index. Will for sure contain the cell, but will be slightly larger than the cell  * */
void diffusion_bounding_box(Cell* pCell, std::vector<int>* bounding_box_by_index)
{
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
/*
 - function to calculate the corner positions of a given voxel, give the voxel center, uses the microenviroment mesh to determine voxel size,
 fills in std::vector<std::vector<double>> &return_corners as return_corners[corner_index][{x,y,z}]  
 */

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
/*
 * Function fills return_intersection_voxel_indicies with those voxels that the cell membrane passes through,
 * an intersection voxel is one in which at least one voxel corner is contained within the cell, but 
 * not all 8 corners are inside the cell radius. Voxel center may be inside or outside the cell.
 * */
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

/*
 * Function fills return_exterior_voxel_indicies with those voxels that some but not all the coners of the voxel are within the cell and the voxel center is outside the cell
 * */
void get_exterior_voxels(Cell* pCell, std::vector<double>* return_exterior_voxel_indicies)
{
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
/* Function returns return_interior_voxel_indicies with the index of all the voxels for whom their center falls within the cell radius*/
void get_interior_voxels(Cell* pCell, std::vector<int>* return_interior_voxel_indicies) //voxel center is inside the cell
{
   //if voxel is edge voxel count as exterior
  //a voxel is exterior if some of its corners fall within the sphere but not all 
  //take tight bounding box and remove voxels where all or none of the corners are <Radius
    std::vector<int> bounding_voxels={};
    diffusion_bounding_box(pCell,&bounding_voxels);
    
    std::vector<int> interior_voxels={};

    for (size_t i = 0; i < bounding_voxels.size(); i++)
    {
      std::vector<double> test_voxel_center=pCell->get_container()->underlying_mesh.voxels[bounding_voxels[i]].center;
      // std::vector <std::vector <double>> test_corners=get_voxel_corners(test_voxel_center);
      // std::vector<std::vector <double>> test_corners(8,std::vector<double>(3,0.0));
      // get_voxel_corners(test_voxel_center,test_corners);
     
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


/*MECHANICS FOR MULTIVOXEL OBJECTS ***************************************************************************/

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
/* Diffusion reaction term for multivoxel cells*************************************************/
// This is up to the user and non-trivial.
