
#ifndef __MULTIVOXEL_H__
#define __MULTIVOXEL_H__
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
#include "./springs.h"
using namespace BioFVM;
using namespace PhysiCell;
namespace multiVoxel{

/*Function gets a voxel bounding box, with offset based on the cells' position relative to the center of a voxel; not gauranteed to be the least bounding box, but will contain the cell and it's immediate neighborhood*/



};


#endif // __MULTIVOXEL_H__
