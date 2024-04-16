#ifndef __ADDON_H__
#define __ADDON_H__
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include "./quaternion.h"
#include <stdexcept>
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;
class Base_Addon_Class{
private:
public:
  Cell* m_pCell;
  bool pCell_is_safe;
  bool is_updated;
  bool divided;
  Cell* m_daughter;
  Base_Addon_Class(Cell* pCell);
  virtual ~Base_Addon_Class();
  Base_Addon_Class* get_instance();
  virtual void update_state()=0;
  virtual void on_division()=0;
};

class Addon_Factory{
  public:
    virtual Base_Addon_Class* Factory_Method(Cell* pCell)=0;
    Base_Addon_Class* Create_Addon_Instance(Cell* pCell);  
};

class Addon{
private:
  
public:
  Addon(Addon const&) = delete; //cannot copy
  Addon& operator=(Addon const&) = delete; // cannot assign
  Addon_Factory* m_addon_factory;
  std::unordered_map< int, Base_Addon_Class* > class_instances_by_pCell;//store instance
  std::unordered_map< int, Base_Addon_Class* > detached_instances_by_pCell;//store detached instance
  Addon(Addon_Factory* custom_class_type);
  //Addon_factory* class_one = new Class_One_Creator();
  ~Addon();
  void spawn_instance(Cell* pCell);//make an instance using the factory and then run the creator
    //m_addon_factory->Create_Addon_Instance
  Base_Addon_Class* get_instance(Cell* pCell);
  void copy_instance_to_daughter(Base_Addon_Class* instance);
  void detach_instance(Cell* pCell);
  void check_pCell_safety(Cell* pCell);// if pCell goes out of domain or is destructed deal with addons, maybe an observer is needed
  void update_custom_class(Cell* pCell); //run the class
     
};
#endif
