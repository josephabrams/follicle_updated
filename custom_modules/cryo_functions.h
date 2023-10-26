/* Author: Joseph E.S. Abrams
 * Please cite if you use this code
 * This is part of a library of functions and paramaters for simualting cells and tissue undergoing cryopreservation*/

#ifndef __CRYO_FUNCTIONS_H__
#define __CRYO_FUNCTIONS_H__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <vector>
#include "cryo_utilities.h"
#include "cryo_parameters.h"
#include "follicle_utilities.h"
#include "springs.h"
#include <unordered_map>
#include <algorithm> 

#define PI 3.14159265
#define R 0.08205 // granulosa (10^-3 J/mole*k)

using namespace BioFVM;
using namespace PhysiCell;
using namespace Springs;
namespace cryo_parameters{


void set_2p_initial_conditions(Cell* pCell, Cryo_Parameters &cryo_p);
};
#endif //__CRYO_FUNCTIONS_H__
