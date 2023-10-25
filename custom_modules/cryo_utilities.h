/*Author Joseph E.S. Abrams
 * Please cite if you use this code.
 * Virial Equations for converting between molarity and molality of some common cryoprotectants.
 * The polynomical coefficients were determined using the python code found in ./virial_eq/virial.py 
 * these use the data from literature that has been saved into the associated csv files, if you need
 * other CPAs you'll need to add the data and run that python code.
 * */
#ifndef __CRYO_UTILITIES_H__
#define __CRYO_UTILITIES_H__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <cstddef>
#include <unordered_map>
#include <memory>
#include <vector>
#include "./multivoxel_utilities.h"

using namespace BioFVM;
using namespace PhysiCell;
double molarity_to_molality(double molarity, std::string component_name);
double binary_virial(double molality, std::string component_name);
double ternary_virial(double molality_1, double molality_2, std::string component_1, std::string component_2);

#endif //__CRYO_UTILITIES_H__
