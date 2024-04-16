#ifndef __QUATERNIONS_H__
#define __QUATERNIONS_H__
#include <ostream>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>
namespace Quaternions{
struct Quaternion{
  private:
  Quaternion();

  public:
  //a+bi+cj+dk
  double m_a;// real
  double m_b;// i
  double m_c;// j
  double m_d;// k
  bool is_pure;
  
  Quaternion(double a, double b, double c, double d);
  Quaternion(double b, double c, double d);
  Quaternion(std::vector <double> Vec);
  
  Quaternion operator+(const Quaternion& Q2) const;
  Quaternion operator-(const Quaternion& Q2) const;
  Quaternion operator*(const Quaternion& Q2) const;
  Quaternion operator*(const double& D) const;
  Quaternion operator/(const double& D) const;
  void set(Quaternion& Q1); //safely change values of a Quaternion
  
  const double norm();
  const Quaternion normalize();
  void normalize_by_ref();
  const Quaternion conjugate();
  const Quaternion multiply(Quaternion &Q1, Quaternion &Q2);
  const Quaternion rotation_by_radians(Quaternion vector_axis, double angle);
  void rotate_by_ref(Quaternion &vector_axis, double angle);
  friend std::ostream& operator<<(std::ostream& os, const Quaternion& Q);
  //Quaternion inverse();
  //double dot();
  //convert_euclid_angles();
  //
};

};

#endif // !QUATERNIONS_H_

