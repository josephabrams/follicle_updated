
#ifndef __MULTIVOXEL_UTILITIES_H__
#define __MULTIVOXEL_UTILITIES_H__

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

void general_voxel_bounding_box(std::vector<int> *return_bounding_box,std::vector<double> &starting_position, std::vector <double>&ending_position, double &voxel_length, BioFVM::Cartesian_Mesh &a_mesh);

void diffusion_bounding_box(Cell* pCell, std::vector<int>* bounding_box_by_index);

void get_voxel_corners(std::vector<double> &voxel_center, std::vector<std::vector<double>> &return_corners );

void get_intersecting_voxels(Cell* pCell,std::vector<int>* return_intersecting_voxel_indicies);
void get_exterior_voxels(Cell* pCell, std::vector<double>* return_exterior_voxel_indicies);

void get_interior_voxels(Cell* pCell, std::vector<int>* return_interior_voxel_indicies); //voxel center is inside the cell

void output_cell_and_voxels(std::vector <int> voxel_list, Cell* pCell);

/*MECHANICS FOR MULTIVOXEL OBJECTS ***************************************************************************/

void spherical_bounding_box(std::vector <double> &center_point, double &radius,std::vector <int> *return_bounding_box);

void max_interaction_variable_moore_neighborhood(Cell* pCell, double &maximum_interaction_distance, std::vector<int> *voxel_neighboorhood);


void cells_in_neighborhood(Cell* pCell, double &maximum_interaction_distance, std::vector<Cell*> *neighboors);
#endif //__MULTIVOXEL_UTILITIES_H__
