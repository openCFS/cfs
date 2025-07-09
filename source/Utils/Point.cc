#include <stddef.h>
#include <ostream>

#include "General/Exception.hh"
#include "MatVec/Vector.hh"
#include "Utils/Point.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

Point& Point::operator+=(const Point& t) {
  assert(data.GetSize() == t.data.GetSize());
  for (unsigned int i=0; i<data.GetSize(); i++)
    data[i]+=t.data[i];
  return *this;
}

Point& Point::operator-=(const Point& t) {
  assert(data.GetSize() == t.data.GetSize());
  for (UInt i=0; i<data.GetSize(); i++)
    data[i]-=t.data[i];
  return *this;
}

Point& Point::operator*=(double factor) {
  for (UInt i=0; i<data.GetSize(); i++)
    data[i] *= factor;
  return *this;
}

Point Point::operator*(double factor) const {
  return Point(data[0]*factor, data[1]*factor, data[2]*factor);
}

double Point::Dot(const Point& other) const
{
  double sum = 0.0;
  sum  = data[0] * other.data[0];
  sum += data[1] * other.data[1];
  sum += data[2] * other.data[2];
  return sum;
}

std::string Point::ToString() const
{
   std::ostringstream os;
   os << "(";
   if(data.GetSize() > 0) 
     os << data[0] << "," << data[1] << "," << data[2];
   os << ")";
   return os.str();
}

int Point::GetCartesianOrientation(const Vector<double>* vec)
{
  assert(vec != NULL);
  return GetCartesianOrientation(vec->GetPointer(), vec->GetSize());
}


int Point::GetCartesianOrientation(const double* vec, int length)
{
  assert(vec != NULL);
  assert(length == 2 || length == 3);
  assert(length == 3 || ((vec[0] != 0.0 && vec[1] == 0.0)
                       ||(vec[0] == 0.0 && vec[1] != 0.0)));

  assert(length == 2 || ((vec[0] != 0.0 && vec[1] == 0.0 && vec[2] == 0.0)
                       ||(vec[0] == 0.0 && vec[1] != 0.0 && vec[2] == 0.0)
                       ||(vec[0] == 0.0 && vec[1] == 0.0 && vec[2] != 0.0)));
  for(int i = 0; i < length; i++)
    if(vec[(unsigned int) i] != 0)
      return i;

  throw Exception("vector is not Cartesian: " + StdVector<double>::ToString(length, vec));
}

}
