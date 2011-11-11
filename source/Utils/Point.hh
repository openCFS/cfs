#ifndef POINT_HH_
#define POINT_HH_

#include "General/Environment.hh" // for the CFS types Double and UInt

namespace CoupledField
{
template<class TYPE> class Vector;

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

  /** Assumes a Cartesian orientation and gives the direction, 0-based!*/
  int GetCartesianOrientation() const {
    return GetCartesianOrientation(&data[0]);
  }

  /** Assumes a Cartesian orientation and gives the direction, 0-based!
   * @param vec assumed to be of size 3 - pointer because of circular inclusion :(
   * @return 0-based index of (firt) non-zero index. */
  static int GetCartesianOrientation(const Vector<double>* vec);

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

  /** Difference to another point but w/o return value to save memory in loops
   * dist = this - other */
  void Dist(const Point& other, Point& dist) const
  {
    dist[0] = data[0] - other.data[0];
    dist[1] = data[1] - other.data[1];
    dist[2] = data[2] - other.data[2];
  }

  /** Lists the content
  * @return the form "(0.3;4.3;0.0)" but no digit control */
  std::string ToString() const;

  Double data[3];

private:
  /** common implementation
   * @param assume to be of size 3! */
  static int GetCartesianOrientation(const double* vec);

};


}
#endif /* POINT_HH_ */
