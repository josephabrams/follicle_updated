#include "./ovarian_follicle.h"

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius, double inner_radius);
void attach_neighboring_springs(Cell* pCell, double max_spring_length);
void initialize_neighboring_spring_connections();
std::vector<double> displacement_between_memebranes(Cell* pCell_1, Cell* pCell_2);
double distance_between_membranes(Cell* pCell_1, Cell* pCell_2);
void non_connected_neighbor_pressure(Cell* pCell, double dt);
void two_parameter_single_step(Cell* pCell,Phenotype& phenotype,double dt);//calculates current volume based on previous volume and parameters