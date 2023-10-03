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
//std::vector <int> general_voxel_bounding_box(std::vector <double> starting_position, std::vector <double>ending_position, double voxel_length, BioFVM::Cartesian_Mesh a_mesh);

void general_voxel_bounding_box(std::vector<int> *return_bounding_box,std::vector<double> &starting_position, std::vector <double>&ending_position, double &voxel_length, BioFVM::Cartesian_Mesh &a_mesh);
void create_spring_cell_from_cell(Cell* pCell, Phenotype& phenotype);
void initialize_neighboring_spring_connections();
void attach_neighboring_springs(Cell *pCell, double max_spring_length); 
std::vector<double> displacement_between_membranes(Cell* pCell_1, Cell* pCell_2);
double distance_between_membranes(Cell* pCell_1, Cell* pCell_2);
void non_connected_neighbor_pressure(Cell* pCell, double dt, double spring_constant);
void two_parameter_single_step(Cell* pCell,Phenotype& phenotype,double dt);//calculates current volume based on previous volume and parameters
//std::vector<int> get_intersecting_voxels(Cell* pCell);
void get_intersecting_voxels(Cell* pCell,std::vector<int>* return_intersecting_voxel_indicies);
//std::vector<int> get_exterior_voxels(Cell * pCell);
//std::vector<int> get_interior_voxels(Cell * pCell);
void get_exterior_voxels(Cell *pCell, std::vector<int>* return_exterior_voxel_indicies);
void get_interior_voxels(Cell *pCell, std::vector<int>* return_interior_voxel_indicies);
double concentration_at_boundary(Cell * pCell, int solute_index);

void update_all_forces(Cell* pCell, double dt, double spring_constant);
void uptake_in_one_voxel(int voxel, double water_uptake_per_voxel, std::vector<double> solute_uptake_per_voxel, std::vector<double> specific_volumes );
double molarity_to_molality(double molarity, std::string component_name);
// std::vector<std::vector<double> >get_voxel_corners (std::vector<double> voxel_center);
void get_voxel_corners(std::vector<double> &voxel_center, std::vector<std::vector<double>> *return_corners );
double binary_virial(double molality, std::string component_name);
double ternary_virial(double molality_1, double molality_2, std::string component_1, std::string component_2);
// std::vector<int> diffusion_bounding_box (Cell *pCell);
void diffusion_bounding_box(Cell* pCell, std::vector<int>* bounding_box_by_index);

// std::vector<Cell *> cells_in_me (Cell *pCell);
void cells_in_me(Cell *pCell, std::vector<Cell*> *return_cells_in_me); // uses mechanics vectors to search for cells that are within or equal to pCell radius can also be used for bounding boxes

void find_basement_membrane_voxels (std::vector<double> center_of_sphere, std::vector<double> radius_of_sphere);

void Hookes_law_force (std::vector<double> &direction, double &rest_length, double &current_length, double &spring_constant, std::vector<double> *return_force);
void output_all_voxels_concentrations();
double Adams_Bashforth_ODE_2nd_Order (double y_value, double prev_df_dt, double df_dt, double step_size);
std::vector<double> Adams_Bashforth_ODE_2nd_Order(std::vector<double> Y_values, std::vector<double> prev_df_dts, std::vector<double> df_dts, double step_size);
double Forward_Euler (double f_value, double step_size, double df_dt);
//experimental add springs with length to cells using Spring_Cell class
//
void uptake(Spring_Cell* SPcell);
void solute_loading( double oxygen, double final_solute_concentration_1, double final_solute_concentration_2,double final_solute_concentration_3,double final_solute_concentration_4, double final_solute_concentration_5);
void activate_nodes(double radius_of_activation);
void activate_edges(double domain_size);

void set_dirchlet_nodes_cylinder(double radius, double height,std::vector<double> center, std::vector<double> dirchlet_condition);//cylinder specifies voxels by voxel center
void rasterize_my_uptake(Cell* pCell, double solute_index);
void update_exterior_concentrations(Spring_Cell* SpCell);
void update_interior_concentrations(Spring_Cell* SPcell);
std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius,double inner_radius);
extern std::vector <int> basement_membrane_voxels;
void get_basement_membrane_intersection(std::vector <double> &center_point, double &radius, std::vector<int> *basement_voxels );
void initialize_basement_membrane(std::vector <double> center_point, double radius);
void initialize_spring_cells();
void initialize_basement_membrane_connections( std::vector <double> basement_membrane_center, double basement_membrane_radius);
// void variable_moore_neighborhood(Cell* pCell, double radius, std::vector<int> *voxel_neighborhood);
void max_interaction_variable_moore_neighborhood(Cell* pCell, double &maximum_interaction_distance, std::vector<int> *voxel_neighborhood);
void connect_spring_cells(Spring_Cell* SpCell_1, Spring_Cell* SpCell_2);
void break_TZPs(Spring_Cell* Oocyte,double max_breakage_distance); 
void cells_in_neighborhood(Cell* pCell, double maximum_interaction_distance, std::vector <Cell*> *neighboors);
void custom_add_potentials(Spring_Cell* SpCell, double dt);
void spring_cells_in_neighborhood(Spring_Cell* SpCell, double maximum_interaction_distance, std::vector <Spring_Cell*> *neighboors);
void initialize_spring_connections();
void initialize_basement_membrane_connections( std::vector <double> basement_membrane_center, double basement_membrane_radius);
bool is_in_voxel(Cell* pCell, Voxel* pVoxel);
bool is_in_voxel(Cell* pCell, int voxel_index);
void basement_membrane_mechanics(Spring_Cell* SpCell, double basement_membrane_radius, std::vector <double> basement_membrane_center, double dt);
// void basement_membrane_mechanics(Cell* pCell, double basement_membrane_radius);
void spherical_bounding_box(std::vector <double> &center_point, double &radius,std::vector <int> *return_bounding_box);
void print_voxels_for_quick_plotting(Cell* pCell,std::vector <int> bounding_voxels, std::vector <int> sub_section);
void custom_add_potentials_for_pCells(Cell* pCell, double dt);
void internal_concentration(Cell* pCell, std::string solute_name);
void output_all_voxels_as_cartesian_index();
void output_cell_and_voxels(std::vector <int> voxel_list, Cell* pCell);
#endif // __FOLLICLE_UTILITIES_H__
