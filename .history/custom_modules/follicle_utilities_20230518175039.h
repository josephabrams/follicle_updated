#include "./ovarian_follicle.h"

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius, double inner_radius);
void attach_neighboring_springs(Cell* pCell, double max_spring_length);
void initialize_neighboring_spring_connections();
std::vector<double> displacement_between_membranes(Cell* pCell_1, Cell* pCell_2);
double distance_between_membranes(Cell* pCell_1, Cell* pCell_2);
void non_connected_neighbor_pressure(Cell* pCell, double dt);
void two_parameter_single_step(Cell* pCell,Phenotype& phenotype,double dt);//calculates current volume based on previous volume and parameters

std::vector<int> get_exterior_voxels(Cell * pCell);
std::vector<int> get_interior_voxels(Cell * pCell);
double concentration_at_boundary(Cell * pCell, int solute_index);

std::vector<std::vector<double> >
get_voxel_corners (std::vector<double> voxel_center);

std::vector<int> diffusion_bounding_box (Cell *pCell);

std::vector<Cell *> cells_in_me (Cell *pCell);

void find_basement_membrane_voxels (std::vector<double> center_of_sphere, std::vector<double> radius_of_sphere);

std::vector<double> Hookes_law_force (std::vector<double> direction,
                                      double rest_length,
                                      double current_length,
                                      double spring_constant);

double Adams_Bashforth_ODE_2nd_Order (double y_value, double prev_df_dt, double df_dt, double step_size);

double Forward_Euler (double f_value, double step_size, double df_dt);
//experimental add springs with length to cells using Spring_Cell class



class Spring
{
  private:
  public:
    double m_length{};
    Spring_Cell* m_neighbor{};
    Spring(Spring_Cell* neighbor, double length);
};
class Spring_Cell
{
//class for adding springs with length between membranes instead of center to center
//class does nothing but add linked data structure
private:

public:
  Cell* m_my_pCell{};
  std::vector<Spring*> m_springs{};
  Spring_Cell(Cell* my_pCell);
  //std::vector<double> spring_lengths={};
  //std::vector<double> neighbor={};
  void add_spring(Spring_Cell* other_pSCell, double spring_length);
  void remove_spring(Spring_Cell* other_pSCell);

};
std::vector <Spring_Cell*> all_spring_cells{};