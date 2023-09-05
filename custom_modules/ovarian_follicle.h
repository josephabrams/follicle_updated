#ifndef __OVARIAN_FOLLICLE_H__
#define __OVARIAN_FOLLICLE_H__

#pragma once
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>

using namespace BioFVM;
using namespace PhysiCell;



extern Cell_Definition oocyte_cell;
extern Cell_Definition hepato_cell;
extern Cell_Definition granulosa_cell;
extern std::string output_filenames;
extern double k_granulosa;
extern double k_oocyte;
extern double k_basement;
std::vector<std::string> my_coloring_function( Cell* pCell );
std::vector<std::string> follicle_coloring_function( Cell* pCell );
////////////////////////////////Custom Cell Construction/////////////////////////////////////////////////////////////////////////

void spring_cell_functions(Cell* pCell, Phenotype& phenotype, double dt);

void create_cell_types( void );
void create_oocyte_cell_type(void);
void create_granulosa_cell_type( void );
void create_hepato_cell_type(void);
void hepato_cell_rule( Cell* pCell, Phenotype& phenotype, double dt );
void hepato_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt );
void hepato_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt );
void oocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double dt );
void oocyte_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt );
void oocyte_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt );

void granulosa_cell_rule( Cell* pCell, Phenotype& phenotype, double dt );
void granulosa_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt );
void granulosa_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt );
///////////////////////////////PhysiCell///////////////////////////////////////////////////////////////////
void setup_just_oocyte(void);
void setup_4_granulosa_test_case(void);
void setup_microenvironment( void );
void setup_secondary_stage_follicle(void);
void setup_toy_granulosa_model(void);

void setup_tissue( void );
struct points{
  double x;
  double y;
  double z;
  points(double i,double j, double k);
};

void setup_curve();
struct poly_corner_2D {
  std::vector<double> my_position;
  std::vector<double> connected_point_1;
  std::vector<double> connected_point_2;
  std::vector<double> segment_1;//{slope,b} 
  std::vector<double> segment_2;
  poly_corner_2D(std::vector<double> position);
  void calculate_connections(std::vector<std::vector<double>> &points);
  void calculate_segments();
  bool is_within(std::vector<double> polygon_center_point);
};
struct polygon{
  std::vector<poly_corner_2D*> corners;
  polygon(int number_of_sides);
  void populate(std::vector<poly_corner_2D*> corner_points);
};
struct bez_curve_3{
  int index;
  std::vector<std::vector<double>> end_points;
  std::vector<std::vector<double>> control_points;
  std::vector<double> ratios;
  void initalize(int number);
  void set_control_points(std::vector<double> p1,std::vector<double> p2);
  void set_end_points(std::vector<double> p1,std::vector<double> p2);
  void set_ratios(double r1,double r2, double r3, double r4);

  void get_point(double t, std::vector<double> *return_point);
  void points_on_curve(double &t_values_as_ratio,std::vector<std::vector<double>> &end_points, 
                       std::vector<std::vector<double>> &control_points, std::vector<double> &ratios, 
                       std::vector<std::vector<double>> *return_points);
  void points_on_me(double &t_values_as_ratio,std::vector<std::vector<double>> *return_points);
};
void cells_on_bez(std::vector<double> t_values,bez_curve_3 *bez, std::vector<std::vector<double>>* cell_positions, double radius);
void cells_on_path(std::vector<std::vector<double>> path_points, std::vector<std::vector<double>>* cell_positions, Cell* pCell);//approximate curvy path with line segments, if segment is smaller than a cell, place cells at points
void cells_on_line(std::vector<double>&line_start,std::vector<double>&line_end, std::vector<std::vector<double>>* cell_positions, Cell* pCell);
// void subtract_overlapping_cells_by_list(std::vector<std::vector<double>>* cell_list);//can almost certainly use voxels to speed up// if center-to_center with other cell < (radius-overlap) remove it from the cell_list 
void subtract_overlapping_cells();
void add_venule(std::vector<double> position, double radius);
void add_portals(std::vector<std::vector<double>> positions, double radius);
void add_sinusoids(double density);
std::vector<std::vector<double>> create_hepato_positions(double cell_radius,
                                                         double sphere_radius, 
                                                         std::vector<std::vector<double>> hole_centers, 
                                                         std::vector<double> hole_radii,
                                                         std::vector<Cell*> sinusoids
                                                         );
void setup_liver();

#endif // __OVARIAN_FOLLICLE_H__
