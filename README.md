*This is a fork of PhysiCell* from https://github.com/MathCancer/PhysiCell/ 
The code was written originally using v1.7 however, it has been updated to v1.10 of PhysiCell

**Note: This is a work in progress. Some core code was minorly altered to allow svg outputs in seconds instead of minutes and to the operator override for streaming out vectors (this version, they cout as comma separated). All custom modules will function without these minor changes, and they were simply for my convenience. I am currently refactoring and cleaning up this repository, so expect frequent changes - note that branches beside the main branch are essentially my developing other features for my PhD project/side projects.**

**If you use this code in part or all, please cite both me and the original PhysiCell since this would all be impossible without them. (like and subscribe for academics?) There are citation files for myself and core PhysiCell located in the repository.**


(slightly outdated see above) The only base PhysiCell file modified is modules/PhysiCell_Pathology.cpp, SVG plots have been altered to display seconds, minutes and hours instead of minutes, hours and days. This is optional to the working of the code. 
Several experimental files in custom_modules are not utilized. See the main.cpp include statements for which files are being used. Eventually, time permitting, I will clean up this repo, refactor the code to just the essentials, and expand the documentation.
Until then, please cite this repository and PhysiCell if you use this code.


# follicle_updated
- Please see wiki for the current status and detailed overview (pending)
- Repo will eventually contain all the current code for the ovarian follicle project.

## This is an updated version that:
1. corrects all race conditions
2. is independent of the PhysiCell core code, so serves as a true add-on (tested up to PhysiCell 1.12 with adjustments made to makefile and config)
3. uses an improved approximation of solution concentrations
4. deals with membrane pressure independent of connections to avoid any failure to apply cell forces
5. uses dynamic voxel neighbourhoods to more efficiently handle mechanics and mass transport with differently sized cells
6. Number 5 also guarantees that any abnormally large deviations in cell position, basement membrane mechanics or volume change are captured (baring cells leaving the domain, destruction of out of scope cells will produce undefined behavior, but a fix is coming)
7. Uses independent mechanics and mass transport linked to PhysiCell by encapsulating the Cell object in a Spring_Cell class
8. Initial conditions moved into the PhysiCell_settings.config xml file to avoid recompiling between simulations
9. Functions to make plotting and analysis in python easier (PhysiCell is relying more and more on matplotlib and this new version reflects that)
Note: this version is not yet compatible with the PhysiCell GUI
Also note: Many of the files still contain old versions of code and uncalled functions or objects. Eventually, this repo will get cleaned up, but for now, the goal is merely to contain all the working code.
A final note: Because the encapsulation in Spring_Cell objects has not been tested when Cells are destroyed or go out of domain, the code behavior is unknown, undefined and vulnerable to memory leaks (does not use smart pointers) so consider yourself warned if you allow this in simulation. Time permitting, I will update and address this. 
