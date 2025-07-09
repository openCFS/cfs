#ifndef POINT_HH_
#define POINT_HH_

#include <assert.h>
#include <cmath>
#include <string>

#include "General/defs.hh"
#include "MatVec/Vector.hh"

namespace CoupledField
{
template<class TYPE> class Vector;

/** Useful 3D coordinate class. Initialized by default to 0.
 * 
* One wonders if this doesn't exist in boost ?? */
class Point
{
public:
  //! constructor


  Point(const double val = 0.0) {
    data.Resize(3, val);
  }

  Point(double x, double y, double z = 0.0) {
    data.Resize(3);
    data[0]=x;
    data[1]=y;
    data[2]=z;
  }

  //!destructor
  ~Point() {}

  /** resets the values */
  void SetZero() {
    data.Init(0.0);
  }

  void Set(double x, double y, double z = 0.0) {
    data[0] = x;
    data[1] = y;
    data[2] = z;  
  } 


  /** Assumes a Cartesian orientation and gives the direction, 0-based!*/
  int GetCartesianOrientation() const {
    return GetCartesianOrientation(&data[0]);
  }

  /** Assumes a Cartesian orientation and gives the direction, 0-based!
   * @param vec assumed to be of size 3 - pointer because of circular inclusion :(
   * @return 0-based index of (first) non-zero index. */
  static int GetCartesianOrientation(const Vector<double>* vec);

  /** Returns data vector */
  inline const Vector<double>& GetCoordVector() { return data;}

  // the default assignment operator is sufficient: inline Point& operator=(const Point& rhs)
  inline bool operator==(const Point& t) const {
    assert(data.GetSize() == t.data.GetSize());
    return data[0] == t.data[0] && data[1] == t.data[1] && data[2] == t.data[2];
  }

  inline Point& operator+=(const Point& t);
  inline Point operator+(const Point& t) const {
    assert(data.GetSize() == t.data.GetSize());
    return Point(data[0]+t.data[0], data[1]+t.data[1], data[2] + t.data[2]);
  }  

  inline Point& operator-=(const Point& t);
  Point operator-(const Point& t) const {
    assert(data.GetSize() == t.data.GetSize());
    return Point(data[0]-t.data[0], data[1]-t.data[1], data[2]-t.data[2]);
  }

  inline Point& operator*=(double factor);
  inline Point operator*(double factor) const;

  /** scale the point */
  inline Point&  operator/=(double factor) {
    data[0] /= factor;
    data[1] /= factor;
    data[2] /= factor;
    return *this; 
  }
  
  inline Point operator/(double factor) const {
    return Point(data[0]/factor, data[1]/factor, data[2]/factor);
  }

  //! return coordinate number i
  inline Double &operator[](UInt i) {
    assert(i < data.GetSize());
    return data[i];
  }

  //! return coordinate number i
  inline Double operator[](UInt i) const {
    assert(i < data.GetSize());
    return data[i];
  }

  /** Interpret point as vector and calc length
   * Is not cached but always calculated! */
  inline Double CalcLength() const
  {
    assert(data.GetSize() == 3);
    return std::sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
  }

  //! calculate distance between two points. Needs to be inline to avoid linking issues
  inline static Double Dist(const Point& a, const Point& b) {
    double preSqrt = 0.0;
    assert(a.data.GetSize() == b.data.GetSize());
    for(UInt i=0; i<a.data.GetSize(); i++)
      preSqrt+= (a.data[i]-b.data[i]) * (a.data[i]-b.data[i]);
    return std::sqrt(preSqrt);
  }

  /** calculate the distance of two 2D points (x,y) and (a,b)*/
  inline static double Dist(double x, double y, double a, double b) {
    return std::sqrt((a-x)*(a-x)+(b-y)*(b-y));
  }

  double Dot(const Point& other) const;

  //! calculate distance to another point
  double Dist(const Point& other) const {
    return Point::Dist(*this, other);
  }

  /** Difference to another point but w/o return value to save memory in loops
   * dist = this - other */
  void Dist(const Point& other, Point& dist) const
  {
    assert(data.GetSize() == other.data.GetSize());
    assert(data.GetSize() == dist.data.GetSize());
    dist[0] = data[0] - other.data[0];
    dist[1] = data[1] - other.data[1];
    dist[2] = data[2] - other.data[2];
  }
  
  size_t GetHash() const {
    size_t hash = std::hash<double>()(data[0]) - std::hash<double>()(data[1]); 
    hash += std::hash<double>()(data[2]);
    return hash;
  }

  /** Lists the content
  * @return the form "(0.3;4.3;0.0)" but no digit control */
  std::string ToString() const;

  Vector<double> data;

private:
  /** common implementation
   * @param assume to be of size 3! */
  static int GetCartesianOrientation(const double* vec, int length = 3);
};

}

namespace std {

template<>
struct hash<CoupledField::Point>
{
  size_t operator()(const CoupledField::Point & obj) const { return obj.GetHash(); }
};
}
#endif /* POINT_HH_ */
