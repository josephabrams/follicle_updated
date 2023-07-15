
//This file and the associated header contain a custom class to enfold the standard cell class in physicell
//Author Joseph Abrams
// Note this was built for PhysiCell v1.10 but it should work with future versions
// TODO Someday: this class could be improved with the addition of some vector <int> for passing voxel indicies around


#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <memory>
#include <stdexcept>

#include "springs.h"
using namespace BioFVM;
using namespace PhysiCell;


namespace Springs{

std::vector <Spring_Cell*> all_spring_cells;
std::vector <Spring_Cell*> spring_cell_by_pCell_index;

/* eventually make its own class objects*/
// Cryo_State_Variables::Cryo_State_Variables(int solute_num)
// {
//     temperature=25;
//     gas_constant=8.314;
//     surface_area=0;
//     next_water_volume=0;//for ABM 2nd order
//     water_volume=0;
//     dVw=0;//cell water flux
//     previous_dVw=0;
//     solute_moles.resize(solute_num,0.0);
//     dN.resize(solute_num,0.0);//cell mole flux of solutes

//     exterior_osmolality=0.0;//total exterior osmolality salt+CPA (mole/kg)
//     interior_osmolality=0.0;// total internal osmolality salt+CPA
//     internal_molarity.resize(solute_num,0.0);
//     external_molarity.resize(solute_num,0.0);
//     Ps.resize(solute_num,0.0);
//     Lp=0.0;
//   return;
// }

void Spring_Cell::initialize()
{
  //initializing stuff after constructor is called
  return;
}

Spring_Cell::Spring_Cell(Cell* pCell)//(Cell* pCell,Phenotype& phenotype)
{
  
  is_basement_connected=false;
  basement_length=0.0;
  m_my_pCell=pCell;

//     temperature;
//     gas_constant;
//     surface_area;
//     next_water_volume;//for ABM 2nd order
//     water_volume;
//     dVw;//cell water flux
//     previous_dVw;
l//     std::vector<double> solute_moles;
l//     std::vector<double> dN;//cell mole flux of solutes

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
  index=0;
   return;
}
Spring_Cell::~Spring_Cell()
{
  delete_spring_cell(this->index);
  return;
}
void Spring_Cell::add_spring(Spring_Cell* other_pSCell, double spring_length)
{
  #pragma omp critical//could speed up with private merge
  {
    std::cout<<"new spring made"<<std::endl;
    Spring* new_spring;
    new_spring=new Spring(other_pSCell->m_my_pCell,spring_length);
    this->m_springs.push_back(new_spring);
  }
  return;
}

void Spring_Cell::remove_spring(Spring_Cell* other_pSCell)
{
  //copied from removing default spring attachments
  #pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < this->m_springs.size() )
		{
			
			if( this->m_springs[i]->m_pNeighbor == other_pSCell->m_my_pCell )
			{
        int n = this->m_springs.size();
        Spring* spring_ptr= this->m_springs[n-i];
				// copy last entry to current position 
				this->m_springs[i] = this->m_springs[n-1]; 
        
				// shrink by one 
				this->m_springs.pop_back(); 
        delete spring_ptr;
				found = true; 
			}
			i++; 
		}
	}
	return; 
}

Spring::Spring(Cell* pNeighbor, double length)
{
  
  m_pNeighbor=pNeighbor;
  m_length=length;
  std::cout<< "spring constructed"<<std::endl;
  return;
}

Spring_Cell* create_spring_cell(Cell* pCell)//( Cell* pCell,Phenotype& phenotype )
{
	//enfold already existing Cell into super class spring cell
	//auto pNewSpring_Cell=std::make_shared<Spring_Cell>;
  Spring_Cell* pNew; 
  #pragma omp critical //could be made multithread probably if needed but shouldn't be an issue
  { 
    
	  pNew = new Spring_Cell(pCell);//(pCell,phenotype);		
	  all_spring_cells.push_back( pNew ); 
	  pNew->index=all_spring_cells.size()-1;
    spring_cell_by_pCell_index[pNew->m_my_pCell->index]=pNew;
  	
  }
  return pNew; 
}
void delete_spring_cell( int index)
{
  #pragma omp critical
  {
	  Spring_Cell* pDeleteMe = all_spring_cells[index]; 
	  // move last item to index location  
	  all_spring_cells[ all_spring_cells.size()-1 ]->index=index;
	  all_spring_cells[index] = all_spring_cells[ all_spring_cells.size()-1 ];
	  // shrink the vector
	  all_spring_cells.pop_back();	
	  delete pDeleteMe; 
	  
  }
  return; 
}

};
//std::vector <Spring_Cell*> all_spring_cells{};