#include "./addon.h"
class Daughter_Tracker: public Base_Addon_Class{
  Cell* m_pCell;
  void initialize(){
    //when quantity is divided by 2 daughter and pcell will match
    int tracer= m_pCell->index*2;
    m_pCell->custom_data.add_variable("Daughter_Tracker",tracer);
  }
void update_state() override{
    is_updated=true;
  }

};
