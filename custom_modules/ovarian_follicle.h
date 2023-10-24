#ifndef __OVARIAN_FOLLICLE_H__
#define __OVARIAN_FOLLICLE_H__

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
#endif // __OVARIAN_FOLLICLE_H__
