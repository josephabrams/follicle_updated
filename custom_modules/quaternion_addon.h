#include "./addon.h"
class Quaternion_Orientation: public Base_Addon_Class{
  public:
  Cell* m_pCell;
  Quaternions::Quaternion quaternion_orientation;
  void sync_orientation();
  void sync_quaternion_orientation();
  Quaternion_Orientation(Cell* pCell);
  void update_state() override{
    is_updated=true;
  }
  void on_division() override{
    return;
  }
};
class Quaternion_Orientation_Creator: public Addon_Factory{
  Base_Addon_Class* Factory_Method(Cell* pCell) override{
    return new Quaternion_Orientation(pCell);
  }


};

