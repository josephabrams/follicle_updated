#include "./addon.h"
class Example_Addon_Class: public Base_Addon_Class{
  public:
    Example_Addon_Class(Cell* pCell);
    int example_int;
    void update_state() override{
    is_updated=true;
  }
    void on_division() override{
    this->example_int = this->example_int/2;
    return;
  }
};


class Example_Addon_Creator : public Addon_Factory{
  Base_Addon_Class* Factory_Method(Cell* pCell) override {
    return new Example_Addon_Class(pCell);
  }
};
