#ifndef POINT_HH_
#define POINT_HH_

#include "General/environment.hh" // for the CFS types Double and UInt

namespace CoupledField
{

/** Useful 3D coordinate class.
* One wonders if this doesn't exist in boost ?? */
class Point
{
public:
  //! constructor
  Point() {
    SetZero();
  }

  Point(const double val)
  {
    data[0] = val;
    data[1] = val;
    data[2] = val;
  }

  Point(const Double x, const Double y, const Double z) {
    data[0]=x;
    data[1]=y;
    data[2]=z;
  }


  //!destructor
  ~Point() {}

  /** resets the values */
  void SetZero() {
    for(UInt i=0; i<3; i++)
      data[i]=0.0;
  }

  //!
  Point & operator=(const Point & t);

  //!
  Point & operator+=(const Point & t);
  Point operator+(const Point &t);
  Point operator+(const Point &t) const;

  //!
  Point & operator-=(const Point &t);
  Point operator-(const Point &t);
  Point operator-(const Point &t) const;

  /** scale the point */
  Point& operator*=(const Double factor);
  Point operator*(const Double factor);
  Point operator*(const Double factor) const;

  /** scale the point */
  Point&  operator/=(const Double factor);
  Point operator/(const Double factor);
  Point operator/(const Double factor) const;

  //! return coordinate number i
  Double &operator[](UInt i) {
    assert(i < 3);
    return data[i];
  }

  //! return coordinate number i
  Double operator[](UInt i) const {
    assert(i < 3);
    return data[i];
  }

  /** Interpret point as vector and calc length
   * Is not cached but always calculated! */
  inline Double CalcLength() const
  {
    return sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
  }

  //! calculate distance between two points
  inline static Double Dist(const Point& a, const Point& b) {
    Double preSqrt = 0.0;
    for(UInt i=0; i<3; i++)
      preSqrt+= (a.data[i]-b.data[i]) * (a.data[i]-b.data[i]);
    return sqrt(preSqrt);
  }

  //! calculate distance to another point
  Double Dist(const Point& other) const {
    return Point::Dist(*this, other);
  }

  /** Lists the content
  * @return the form "(0.3;4.3;0.0)" but no digit control */
  std::string ToString() const;

  Double data[3];
};


}
#endif /* POINT_HH_ */
