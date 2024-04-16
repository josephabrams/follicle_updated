#include "./quaternion.h"
// #include <cmath>
// #include <stdexcept>
// #include <cassert>
// #include <cmath>
namespace Quaternions{
Quaternion::Quaternion(double a, double b, double c, double d): m_a{a}, m_b{b}, m_c{c}, m_d{d}, is_pure{false} {
  if(a-0.0==0.0){this->is_pure=true;}
}
Quaternion::Quaternion(double b, double c, double d): m_a{0}, m_b{b}, m_c{c}, m_d{d}, is_pure{true}{}
Quaternion::Quaternion(): m_a{0}, m_b{0}, m_c{0}, m_d{0}, is_pure{false} {}//private_constructor is_pure is ambiguous here set to false by default
Quaternion::Quaternion(std::vector <double> Vec) {

  if (Vec.size()==3 || (Vec.size()==4 && (Vec[0]-0.0==0.0))) {
    this->m_a=0; 
    this->m_b=Vec[0]; 
    this->m_c=Vec[1]; 
    this->m_d=Vec[2];
    is_pure=true;
  }
  else if (Vec.size()==4) {
    this->m_a=Vec[0]; 
    this->m_b=Vec[1]; 
    this->m_c=Vec[2];
    this->m_d=Vec[3];
    is_pure=false;
  }
  else{
    throw std::invalid_argument("Quaternion constructed with bad vector size!");
  }
}

Quaternion Quaternion::operator+(const Quaternion& Q2) const{
  Quaternion result;
  result.m_a= (this->m_a + Q2.m_a);
  result.m_b= (this->m_b + Q2.m_b);
  result.m_c= (this->m_c + Q2.m_c);
  result.m_d= (this->m_d + Q2.m_d);
  return result;
}
Quaternion Quaternion::operator-(const Quaternion& Q2) const{
  Quaternion result;
  result.m_a= (this->m_a - Q2.m_a);
  result.m_b= (this->m_b - Q2.m_b);
  result.m_c= (this->m_c - Q2.m_c);
  result.m_d= (this->m_d - Q2.m_d);
  return result;
}
Quaternion Quaternion::operator*(const Quaternion& Q2) const{
  Quaternion result;
  result.m_a= (this->m_a * Q2.m_a)-(this->m_b * Q2.m_b)-(this->m_c * Q2.m_c)-(this->m_d * Q2.m_d);
  result.m_b= (this->m_a * Q2.m_b)+(this->m_b * Q2.m_a)+(this->m_c * Q2.m_d)-(this->m_d * Q2.m_c);//i
  result.m_c= (this->m_a * Q2.m_c)+(this->m_c * Q2.m_a)+(this->m_d * Q2.m_b)-(this->m_b * Q2.m_d);//j
  result.m_d= (this->m_a * Q2.m_d)+(this->m_d * Q2.m_a)+(this->m_b * Q2.m_c)-(this->m_c * Q2.m_b);//k
  return result;
}

Quaternion Quaternion::operator*(const double& D) const{
  Quaternion result;
  result.m_a= (this->m_a * D);
  result.m_b= (this->m_b * D);
  result.m_c= (this->m_c * D);
  result.m_d= (this->m_d * D);
  return result;
}
Quaternion Quaternion::operator/(const double& D) const{
  Quaternion result;
  result.m_a= (this->m_a / D);
  result.m_b= (this->m_b / D);
  result.m_c= (this->m_c / D);
  result.m_d= (this->m_d / D);
  return result;
}
void Quaternion::set(Quaternion& Q1){
   this->m_a =Q1.m_a;
   this->m_b =Q1.m_b;
   this->m_c =Q1.m_c;
   this->m_d =Q1.m_d;
}
const double Quaternion::norm(){
  
  double result=(this->m_a * this->m_a)+ (this->m_b * this->m_b) + (this->m_c*this->m_c) + (this->m_d*this->m_d);
  result=sqrt(result);
  return result;
}
const Quaternion Quaternion::normalize(){
  Quaternion normal_vec= *this/ this->norm();
  return normal_vec;
}
void Quaternion::normalize_by_ref(){
  *this = *this/this->norm();
}
const Quaternion Quaternion::conjugate(){
  Quaternion result;
  result.m_a= this->m_a;
  result.m_b= (this->m_b * -1.0);
  result.m_c= (this->m_c * -1.0);
  result.m_d= (this->m_d * -1.0);
  return result;
}
const Quaternion Quaternion::multiply(Quaternion &Q1, Quaternion &Q2){
  Quaternion result= Q1*Q2;
  return result;
}
const Quaternion Quaternion::rotation_by_radians( Quaternion vector_axis, double angle)
{
  if(!vector_axis.is_pure){
    throw std::invalid_argument(" Rotation about vector axis must be a pure Quaternion, i.e. a vector! "); 
  }
  double real_part= std::cos(angle*0.5);
  Quaternion normal_vec= vector_axis/vector_axis.norm();
  normal_vec= normal_vec*(std::sin(angle*0.5));
  Quaternion q(real_part,normal_vec.m_b,normal_vec.m_c,normal_vec.m_d);
  Quaternion* p=this;
  Quaternion result=*p*q*q.conjugate();
  return result; 
}

void Quaternion::rotate_by_ref( Quaternion &vector_axis, double angle)
{
  Quaternion temp=rotation_by_radians(vector_axis, angle);
  this->set(temp); 
}

std::ostream& operator<<(std::ostream& os, const Quaternion& Q )
{
 
  os << Q.m_a << " + ("<< Q.m_b <<")i + ("<< Q.m_c<<")j + ("<< Q.m_d<<")k"; 
 return os; 
}




}; //namespace ends 
