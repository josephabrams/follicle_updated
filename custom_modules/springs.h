#ifndef __SPRINGS_H__
#define __SPRINGS_H__

#pragma once
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include<unordered_map>
using namespace BioFVM;
using namespace PhysiCell;

namespace Springs{
// class Cryo_State_Variables
// {
//   public:
    
//     static double temperature;
//     static double gas_constant;
//     double surface_area;
//     double next_water_volume;//for ABM 2nd order
//     double water_volume;
//     double dVw;//cell water flux
//     double previous_dVw;
//     std::vector<double> solute_moles;
//     std::vector<double> dN;//cell mole flux of solutes

//     double exterior_osmolality;//total exterior osmolality salt+CPA (mole/kg)
//     double interior_osmolality;// total internal osmolality salt+CPA
//     std::vector<double> internal_molarity;
//     std::vector<double> external_molarity;
//     std::vector<double> Ps;
//     double Lp;

//   //2p methods
//     Cryo_State_Variables(int solute_num);
//     void get_external_molarity(int component, Cell* pCell);
//     double compute_EG_molality(double EG_molarity);
//     double compute_NaCl_molality(double EG_molarity);
    


// };
class Spring //spring between cells that should be from membrane to membrane
{
  private:
  public:
    double m_length;
    Cell* m_pNeighbor;
    Spring(Cell* pNeighbor, double length);
};

class Spring_Cell
{
//class for adding springs with length between membranes instead of center to center
//class does nothing but add linked data structure
private:

public:
  bool is_basement_connected;
  double basement_length;
  Cell* m_my_pCell;
  //Phenotype m_phenotype;
  int index;
  std::vector<Spring*> m_springs;
  int simulation_selected;
  //constants available to every spring cell, future version could have changing temp
  static const double temperature;
  static const double gas_constant;
  static const std::vector <double> solute_specific_volume;
  static const std::vector <double> kdiss;//virial ionic constant
  static const std::vector <double> virial_B;
  static const std::vector <double> virial_C;
  static const std::vector <double> molar_mass;
  double solute_volume;
  double surface_area;
  double next_water_volume;//for ABM 2nd order
  double water_volume;
  double solid_volume;
  double dVw;//cell water flux
  double previous_dVw;
  double previous_radius;
  std::vector<double> solute_moles;
  std::vector<double> next_solute_moles;
  std::vector<double> previous_position;
  std::vector<double> previous_velocity;
  std::vector<double> previous_acceleration;
  std::vector<double> uptake;//molar uptake/secretion of solutes
  std::vector<double> dN;//cell mole flux of solutes
  std::vector<double> previous_dN;//previous mole flux of solutes
  bool use_virial;
  double exterior_osmolality;//total exterior osmolality salt+CPA (mole/kg)
  double interior_osmolality;// total internal osmolality salt+CPA
  std::vector<double> interior_molarity;
  std::vector<double> exterior_molarity;
  std::vector<double> interior_component_molality;
  std::vector<double> exterior_component_molality;
  std::vector<double> Ps;
  double Lp;
  double toxicity;
  std::vector <double> solute_uptake;
  double water_uptake;
  std::vector <int> uptake_voxels;
  int initial_number_of_connections;
//   //2p methods
//     Cryo_State_Variables(int solute_num);
//     void get_external_molarity(int component, Cell* pCell);//sample microenvironment
   //void compute_molality();
//     double compute_NaCl_molality(double EG_molarity);
    //     void compute_GLY_molality();
    //     void compute_solution_osmolality();//ternary virial
    //     void compute_component_osmolality();//binary virial
    //     void compute_mixed_solution_osmolality; //ternary virial+binary virial for HM
    //for mixture use virial for GLY and EG but add HM
  void dN_molarity(); 
  void dVw_Osmolality();
  void two_p_forward_step(double dt);
    // void two_parameter_update
    // void cell_volume_update
    // void compute_internal_component_molarity //moles/(V_total-V_b)
  // methods
  void initialize();
  Spring_Cell(Cell* pCell);//(Cell* pCell, Phenotype& phenotype);
  ~Spring_Cell();
  //std::vector<double> spring_lengths={};
  //std::vector<double> neighbor={};
   
  void two_p_update_volume();
  void set_2p_initial_conditions();
  void initialize_mechanics();
  void add_spring(Spring_Cell* other_pSCell, double spring_length);
  void remove_spring(Spring_Cell* other_pSCell);
 // 

};



extern std::vector <Spring_Cell* > all_spring_cells;
extern std::vector <Spring_Cell*> spring_cell_by_pCell_index;
Spring_Cell* create_spring_cell(Cell* pCell);//( Cell* pCell, Phenotype& phenotype);
void delete_spring_cell( int index);

};
#endif // __SPRINGS_H__
