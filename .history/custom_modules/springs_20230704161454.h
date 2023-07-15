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
class Cryo_State_Variables
{
  public:
    
    double temperature;
    static double gas_constant;
    double surface_area;
    double next_water_volume;//for ABM 2nd order
    double water_volume;
    double dVw;//cell water flux
    double previous_dVw;
    std::vector<double> solute_moles;
    std::vector<double> dN;//cell mole flux of solutes

    double exterior_osmolality;//total exterior osmolality salt+CPA (mole/kg)
    double interior_osmolality;// total internal osmolality salt+CPA
    std::vector<double> internal_molarity;
    std::vector<double> external_molarity;
    std::vector<double> Ps;
    double Lp;

  //2p methods
    Cryo_State_Variables(int solute_num);
    void get_external_molarity(int component, Cell* pCell);
    double compute_EG_molality(double EG_molarity);
    double compute_NaCl_molality(double EG_molarity);
    


};
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
  // methods
  void initialize();
  Spring_Cell(Cell* pCell);//(Cell* pCell, Phenotype& phenotype);
  ~Spring_Cell();
  //std::vector<double> spring_lengths={};
  //std::vector<double> neighbor={};
  void add_spring(Spring_Cell* other_pSCell, double spring_length);
  void remove_spring(Spring_Cell* other_pSCell);
  Cryo_State_Variables cryo;

};



extern std::vector <Spring_Cell* > all_spring_cells;
extern std::vector <Spring_Cell*> spring_cell_by_pCell_index;
Spring_Cell* create_spring_cell(Cell* pCell);//( Cell* pCell, Phenotype& phenotype);
void delete_spring_cell( int index);

};
#endif // __SPRINGS_H__