#include "./addon.h"

Base_Addon_Class::Base_Addon_Class(Cell* pCell): m_pCell{pCell},pCell_is_safe{true},is_updated{false},m_daughter{nullptr} {
  if(pCell->is_active && !pCell->is_out_of_domain)
  {
    m_pCell=pCell;
    pCell_is_safe=true;
  }
  else
  {
   throw std::invalid_argument("assigning Addon to invalid pCell"); 
  }
}
Base_Addon_Class* Addon_Factory::Create_Addon_Instance(Cell* pCell){
  Base_Addon_Class* ptr = this->Factory_Method(pCell);
  return ptr;
}
Addon::Addon(Addon_Factory* custom_class_type):m_addon_factory{custom_class_type}{
}

Addon::~Addon(){
  //delete all the instances located in map 
  // this could be very slow let the stack handle Addon instances at the end of simulations
  
  for ( auto it = class_instances_by_pCell.begin(); it != class_instances_by_pCell.end(); ++it )
  {  delete it->second;}
  for ( auto it = detached_instances_by_pCell.begin(); it != detached_instances_by_pCell.end(); ++it )
  {  delete it->second;}
  
}

void Addon::spawn_instance(Cell* pCell){
  Base_Addon_Class* ptr= m_addon_factory->Create_Addon_Instance(pCell);
  std::pair<int,Base_Addon_Class*> new_instance (pCell->index,ptr);
  class_instances_by_pCell.insert(new_instance);
}
Base_Addon_Class* Addon::get_instance(Cell* pCell){//could overload [] if code gets really cludgey
  return class_instances_by_pCell[pCell->index];
  
}
void Addon::copy_instance_to_daughter(Base_Addon_Class* instance){
  this->spawn_instance(instance->m_daughter);
  class_instances_by_pCell[instance->m_pCell->index]->on_division();
  class_instances_by_pCell[instance->m_daughter->index]->on_division();
  instance->m_daughter=nullptr; //reset 
}

void Addon::detach_instance(Cell* pCell){
  Base_Addon_Class* detach_ptr=class_instances_by_pCell[pCell->index];
  std::pair<int,Base_Addon_Class*> detached_instance (pCell->index,detach_ptr);
  detached_instances_by_pCell.insert(detached_instance);
  class_instances_by_pCell.erase(pCell->index);

}

void Addon::check_pCell_safety(Cell* pCell){
  Base_Addon_Class* instance=class_instances_by_pCell[pCell->index];
  if(!pCell->is_active || pCell->is_out_of_domain){
    
    instance->pCell_is_safe=false; 
    detach_instance(pCell);
  }
  else if(pCell->phenotype.death.dead){
    instance->pCell_is_safe=false; 
    detach_instance(pCell);
    //do dead stuff if you've overwritten what PhysiCell will do by default
  }
  else if(pCell->phenotype.flagged_for_removal){
    #pragma critical
    {
      instance->pCell_is_safe=false;
      detach_instance(pCell);
    }
  }
  else if(pCell->phenotype.flagged_for_division){
     #pragma critical
    {
      instance->pCell_is_safe=false;
      detach_instance(pCell);// maybe you want to detach, maybe not,
      //instance->m_daughter should be set here need some way to get child if you want to attach new addon to it
      if (instance->m_daughter!=nullptr)
      {
        copy_instance_to_daughter(instance);// be sure addon has on_division specified if you want this to work correctly 
      }
    }
  }  
  return;
}

void Addon::update_custom_class(Cell* pCell){
  Base_Addon_Class* instance=class_instances_by_pCell[pCell->index];
  if(!(instance->is_updated) || instance->pCell_is_safe){
    instance->update_state();
  }
  check_pCell_safety(instance->m_pCell);
  return;
}


