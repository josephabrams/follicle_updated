#include "./quaternion_test.h"
using namespace Quaternions;
int main(){
  std::cout<< "RUNNING!\n";
  std::cout<<test_construction();
  std::cout<<test_operations();
  return 0;
}
