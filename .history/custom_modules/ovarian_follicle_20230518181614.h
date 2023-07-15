#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include "./f"
using namespace BioFVM;
using namespace PhysiCell;

extern Cell_Definition oocyte_cell;
extern Cell_Definition granulosa_cell;
extern std::string output_filenames;
std::vector<std::string> my_coloring_function( Cell* pCell );
std::vector<std::string> follicle_coloring_function( Cell* pCell );
////////////////////////////////Custom Cell Construction/////////////////////////////////////////////////////////////////////////
void create_cell_types( void );
void create_oocyte_cell_type(void);
void oocyte_cell_rule( Cell* pCell, Phenotype& phenotype, double dt );
void oocyte_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt );
void oocyte_velocity_rule( Cell* pCell, Phenotype& phenotype, double dt );
void create_granulosa_cell_type( void );
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