#ifndef __FOLLICLE_UTILITIES_H__
#define __FOLLICLE_UTILITIES_H__

#pragma once

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "./springs.h"
using namespace BioFVM;
using namespace PhysiCell;
using namespace Springs;

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius, double inner_radius);
std::vector <int> general_voxel_bounding_box(std::vector <double> starting_position, std::vector <double>ending_position, double voxel_length, BioFVM::Cartesian_Mesh a_mesh);
void create_spring_cell_from_cell(Cell* pCell, Phenotype& phenotype);
void initialize_neighboring_spring_connections();
void attach_neighboring_springs(Cell *pCell, double max_spring_length); 
std::vector<double> displacement_between_membranes(Cell* pCell_1, Cell* pCell_2);
double distance_between_membranes(Cell* pCell_1, Cell* pCell_2);
void non_connected_neighbor_pressure(Cell* pCell, double dt);
void two_parameter_single_step(Cell* pCell,Phenotype& phenotype,double dt);//calculates current volume based on previous volume and parameters
std::vector<int> get_intersecting_voxels(Cell* pCell);
std::vector<int> get_exterior_voxels(Cell * pCell);
std::vector<int> get_interior_voxels(Cell * pCell);
double concentration_at_boundary(Cell * pCell, int solute_index);

std::vector<std::vector<double> >
get_voxel_corners (std::vector<double> voxel_center);

std::vector<int> diffusion_bounding_box (Cell *pCell);

std::vector<Cell *> cells_in_me (Cell *pCell);

void find_basement_membrane_voxels (std::vector<double> center_of_sphere, std::vector<double> radius_of_sphere);

std::vector<double> Hookes_law_force (std::vector<double> direction, double rest_length, double current_length, double spring_constant);

double Adams_Bashforth_ODE_2nd_Order (double y_value, double prev_df_dt, double df_dt, double step_size);
std::vector<double> Adams_Bashforth_ODE_2nd_Order(std::vector<double> Y_values, std::vector<double> prev_df_dts, std::vector<double> df_dts, double step_size);
double Forward_Euler (double f_value, double step_size, double df_dt);
//experimental add springs with length to cells using Spring_Cell class

void solute_loading( double oxygen, double final_solute_concentration_1, double final_solute_concentration_2,double final_solute_concentration_3,double final_solute_concentration_4, double final_solute_concentration_5);

void rasterize_my_uptake(Cell* pCell, double solute_index);

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius,double inner_radius);
extern std::vector <int> basement_membrane_voxels;
std::vector<int> get_basement_membrane_intersection(std::vector <double> center_point, double radius);
void initialize_basement_membrane(std::vector <double> center_point, double radius);
void initialize_spring_cells();
void initialize_basement_membrane_connections( std::vector <double> basement_membrane_center, double basement_membrane_radius);
std::vector <int> variable_moore_neighborhood(Cell* pCell, double radius);
std::vector <int> max_interaction_variable_moore_neighborhood(Cell* pCell, double maximum_interaction_distance);
void connect_spring_cells(Spring_Cell* SpCell_1, Spring_Cell* SpCell_2);
std::vector <Cell*> cells_in_neighborhood(Cell* pCell, double maximum_interaction_distance);
void custom_add_potentials(Spring_Cell* SpCell);
std::vector <Spring_Cell*> spring_cells_in_neighborhood(Spring_Cell* SpCell, double maximum_interaction_distance);
void initialize_spring_connections();
void initialize_spring_cells();
void initialize_basement_membrane_connections( std::vector <double> basement_membrane_center, double basement_membrane_radius);
bool is_in_voxel(Cell* pCell, Voxel* pVoxel);
bool is_in_voxel(Cell* pCell, int voxel_index);
void basement_membrane_mechanics(Spring_Cell* SpCell, double basement_membrane_radius, std::vector <double> basement_membrane_center);
void basement_membrane_mechanics(Cell* pCell, double basement_membrane_radius);
std::vector <int> spherical_bounding_box(std::vector <double> center_point, double radius);
void print_voxels_for_quick_plotting(Cell* pCell,std::vector <int> bounding_voxels, std::vector <int> sub_section);
void custom_add_potentials_for_pCells(Cell* pCell);
void internal_concentration(Cell* pCell, std::vector<double)

#endif // __FOLLICLE_UTILITIES_H__