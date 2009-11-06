#include "Utils/Point.hh"

namespace CoupledField
{

Point & Point::operator=(const Point& t) {
  for (UInt i=0; i<3; i++)
    data[i]=t.data[i];
  return *this;
}

Point & Point::operator+=(const Point& t) {
  UInt i;
  for (i=0; i<3; i++)
    data[i]+=t.data[i];
  return *this;
}

Point Point::operator+(const Point &t) {
  return Point(data[0]+t.data[0], data[1]+t.data[1], data[2] + t.data[2]);
}

Point Point::operator+(const Point &t) const {
  return Point(data[0]+t.data[0], data[1]+t.data[1], data[2] + t.data[2]);
}

Point & Point::operator-=(const Point & t) {
  for (UInt i=0; i<3; i++)
    data[i]-=t.data[i];
  return *this;
}

Point Point::operator-(const Point &t) {
  return Point(data[0]-t.data[0], data[1]-t.data[1], data[2]-t.data[2]);
}

Point Point::operator-(const Point &t) const {
  return Point(data[0]-t.data[0], data[1]-t.data[1], data[2]-t.data[2]);
}

Point& Point::operator*=(const Double factor) {
  for (UInt i=0; i<3; i++)
    data[i] *= factor;
  return *this;
}

Point Point::operator*(const Double factor) {
  return Point(data[0]*factor, data[1]*factor, data[2]*factor);
}

Point Point::operator*(const Double factor) const {
  return Point(data[0]*factor, data[1]*factor, data[2]*factor);
}

Point& Point::operator/=(const Double factor) {
  for (UInt i=0; i<3; i++)
    data[i] /= factor;
  return *this;
}

Point Point::operator/(const Double factor) {
  return Point(data[0]/factor, data[1]/factor, data[2]/factor);
}

Point Point::operator/(const Double factor) const {
  return Point(data[0]/factor, data[1]/factor, data[2]/factor);
}

std::string Point::ToString() const
{
   std::ostringstream os;
   os << "(" << data[0] << ";" << data[1] << ";" << data[2] << ")";
   return os.str();
}


}
