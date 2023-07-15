#include "./ovarian_follicle.h"

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius, double inner_radius);
void attach_neighboring_springs(Cell* pCell, double max_spring_length);
void initialize_neighboring_spring_connections();
std::vector<double> displacement_between_membranes(Cell* pCell_1, Cell* pCell_2);
double distance_between_membranes(Cell* pCell_1, Cell* pCell_2);
void non_connected_neighbor_pressure(Cell* pCell, double dt);
void two_parameter_single_step(Cell* pCell,Phenotype& phenotype,double dt);//calculates current volume based on previous volume and parameters

std::vector<int> get_exterior_voxels(Cell * pCell);

double concentration_at_boundary(Cell * pCell, int solute_index);

std::vector<std::vector<double> >
get_voxel_corners (std::vector<double> voxel_center);

std::vector<int> diffusion_bounding_box (Cell *pCell);

std::vector<Cell *> cells_in_me (Cell *pCell);
