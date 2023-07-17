# follicle_updated
- Pleae see wiki for current status and detailed overview (pending)
- Repo will eventually contain all the current code for the ovarian follicle project.

## This is an updated version that:
1. corrects all race conditions
2. is independant of the PhysiCell core code so serves as a true add on (tested up to PhysiCell 1.12 with adjustments made to makefile and config)
3. uses an improved approximation of solution concentrations
4. deals with membrane pressure independant of connections to avoid any failure to apply cell forces
5. uses dynamic voxel neighboorhoods to more efficiently handle mechanics and mass transport
6. number 5 also gaurantees that any abnormally large deviations in cell position, basement membrane mechanics or volume change are captured
7. Uses independant mechanics and mass transport linked to PhysiCell by encapsulating the Cell object in a Spring_Cell class
8. Initial conditions moved into the PhysiCell_settings.config xml file to avoid recompiling between simulations
9. Functions to make plotting and analysis in R/python easier (PhysiCell is relying more and more on matplotlib and this new version reflects that)
Note: this version is not yet compatible with the PhysiCell GUI
Also note: Many of the files still contain old versions of code and uncalled functions or objects. Eventually this repo will get cleaned up, but for now the goal is mearly to contain all the working code.
A final note: Because the encapsulation in Spring_Cell objects has been tested when Cells are destroyed or go out of domain, the code behavior is unknown, undefined and vulnerable to memory leaks (does not use smart pointers) so consider yourself warned if you allow this in simulation. Time permitting, I will update and address this. 
