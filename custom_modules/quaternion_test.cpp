#include "./quaternion.h"
#include "./quaternion_test.h"
using namespace Quaternions;

bool test_construction(){
  std::cout<< "\nTESTING CONSTRUCTORS: \n\n";
  
  Quaternion q1{0,1,1,1};
  std::cout<<"Constructed q1: "<< q1<<"\n";
  std::cout<<"is pure Quaternion: " <<q1.is_pure<<"\n"; 
  
  Quaternion q2{1,1,1};
  std::cout<<"Constructed q2: "<< q2<<"\n";
  std::cout<<"is pure Quaternion: " <<q2.is_pure<<"\n"; 
  
  std::vector<double> v1 {1.1,2.2,3.3};
  Quaternion q3(v1);
  std::cout<<"Constructed q3: "<< q3<<"\n";
  std::cout<<"is pure Quaternion: " <<q3.is_pure<<"\n"; 
  
  std::vector<double> v2 {1.1,2.2,3.3, 4.4};
  Quaternion q4(v2);
  std::cout<<"Constructed q4: "<< q4<<"\n";
  std::cout<<"is pure Quaternion: " <<q4.is_pure<<"\n"; 
  return true;
}

bool test_operations(){
  std::cout<< "\nTESTING OPERATIONS: \n\n";
  Quaternion Q1{100, 121, 225,1000};
  std::cout<< "Q1: "<< Q1<<"\n";
  double norm=Q1.norm();
  Quaternion Q1_normal= Q1.normalize();

  std::cout<< Q1_normal <<" Is the result of Q1 divided by "<< norm<<"\n";
  Q1.normalize_by_ref();
  std::cout<< Q1 << " normalized!\n";
  
  Quaternion Q2{8888888, 999999, 77777,66666};
  Quaternion Q_sum = Q1+Q2;
  std::cout<< "Q_sum: " << Q_sum<<"\n";

  Quaternion Q_subtract= Q1-Q2;
  std::cout<< "Q_subtract: " <<Q_subtract<< "\n";

  Quaternion Q_product= Q1*Q2;
  std::cout<< "Q_product: " <<Q_product<< "\n";

  Quaternion Q_scalar_product= Q2*2.5;
  std::cout<< "Q_scalar_product: " <<Q_scalar_product<< "\n";

  Quaternion Q_scalar_division= Q2/2.5;
  std::cout<< "Q_scalar_division: " <<Q_scalar_division<< "\n";
  
  Q1.set(Q2);
  std::cout<< "Q1.set() : " <<Q1<< "\n";

  Quaternion Q1_conj=Q1.conjugate();
  std::cout<< "Q1_conj: " <<Q1_conj<< "\n";

  Quaternion Q_multiply= Q2.multiply(Q1,Q2);
  std::cout<< "Q_multiply: " <<Q_multiply<< "\n";

  Quaternion z_axis{0,0,25};
  std::cout<< "z_axis: " <<z_axis<< "\n";
  
  Q2.rotate_by_ref(z_axis, (3.14159));
  std::cout<< "Q2 rotated by PI: " <<Q2<< "\n";

  Quaternion Q_rotation=Q1.rotation_by_radians(z_axis, (3.14159/2));
  std::cout<< "Q_rotation is Q1 rotated by PI/2: " <<Q_rotation<< "\n";

  return true;
}
