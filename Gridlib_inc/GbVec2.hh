/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

// ----------------------------------------------------------------------------
//
//  KLASSE  Vec2
//
//  Ein einfacher Punkt im 2D mit Zugriffsfunktionen wie bei einem Vektor   
//  Abwandlung des Vis3D aus Vision
//
//-----------------------------------------------------------------------------

#ifndef  GBVEC2_HH
#define  GBVEC2_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template<class T>
class GbVec2
{
public:
  //! Various constructors: this class is intended for floating point
  //! numbers. Note the special treatment of integers
  GbVec2 ();
  explicit GbVec2 (T s);
  GbVec2 (T x, T y);
  GbVec2 (const T p[2]);
  GbVec2 (const GbVec2<T>& p);

  //! Destructor
  ~GbVec2();

  // comparison operators (support fuzzy arithmetic if EPSILON > 0)
  INLINE GbBool operator== (const GbVec2<T>& p) const;
  INLINE GbBool operator!= (const GbVec2<T>& p) const;
  INLINE GbBool operator<  (const GbVec2<T>& p) const;
  INLINE GbBool operator<= (const GbVec2<T>& p) const;
  INLINE GbBool operator>  (const GbVec2<T>& p) const;
  INLINE GbBool operator>= (const GbVec2<T>& p) const;

  // arithmetic operators
  INLINE GbVec2<T>& operator+= (const GbVec2<T>& rkVector);
  INLINE GbVec2<T>& operator+= (const T& fScalar);
  INLINE GbVec2<T>& operator-= (const GbVec2<T>& rkVector);
  INLINE GbVec2<T>& operator-= (const T& fScalar);
  INLINE GbVec2<T>& operator*= (const int& s);
  INLINE GbVec2<T>& operator*= (const T& fScalar);
  INLINE GbVec2<T>& operator/= (const int& s);
  INLINE GbVec2<T>& operator/= (const T& fScalar);

  // unary operators
  INLINE GbVec2<T> operator- () const;

  // binary operators
  INLINE GbVec2<T> operator+ (const GbVec2<T>& rkVector) const;
  INLINE GbVec2<T> operator+ (const T& fScalar) const;
  INLINE GbVec2<T> operator- (const GbVec2<T>& rkVector) const;
  INLINE GbVec2<T> operator- (const T& fScalar) const;
  INLINE GbVec2<T> operator* (const int& s) const;
  INLINE GbVec2<T> operator* (const T& fScalar) const;
  INLINE GbVec2<T> operator/ (const int& s) const;
  INLINE GbVec2<T> operator/ (const T& fScalar) const;

  friend GbVec2<T> operator* (T fScalar, const GbVec2<T>& rkVector);

  // norm
  INLINE T normalize(T fTolerance = std::numeric_limits<T>::epsilon());
  INLINE T getNorm() const;
  INLINE T getSquareNorm() const;
  INLINE GbVec2<T> getNormalized(T fTolerance = std::numeric_limits<T>::epsilon()) const;

  // scalar product
  INLINE T dot (const GbVec2<T>& rkVector) const;
  INLINE T operator|(const GbVec2<T>& rkVector) const;

  // vector product
  INLINE GbVec2<T> cross () const;  // returns (y,-x)
  INLINE GbVec2<T> unitCross () const;  // returns (y,-x)/sqrt(x*x+y*y)
  INLINE T operator^(const GbVec2<T>& rkVector) const;  // return the z component as if the vectors were embedded in 3D

  // Access operators
  INLINE T &        operator[] (int i);
  INLINE const T &  operator[] (int i) const;

  //! Read methods
  INLINE const T * getVec2() const;
  INLINE void      getVec2(T &x, T &y) const;
  
  //! Direct manipulation
  INLINE void setX(T x);
  INLINE void setY(T y);
  INLINE void setVec2(const T f[]);
  INLINE void setVec2(T x,T y);
//  INLINE static GbVec2<T> from3D(GbVec3<T> &v);

  // Gram-Schmidt orthonormalization.
  INLINE static void orthonormalize (GbVec2<T> v[/*2*/]);

  // I/O operators
  friend std::ostream& operator<<(std::ostream& s, const GbVec2<T>& v);
  friend std::istream& operator>>(std::istream& s, GbVec2<T>& v);

  // special points
  static const GbVec2<T> ZERO;
  static const GbVec2<T> UNIT_X;
  static const GbVec2<T> UNIT_Y;

  // fuzzy arithmetic (set to > 0 to enable)
  static T EPSILON;

private:
  T xy[2];
};

template <class T>
GbVec2<T> 
operator* (T fScalar, const GbVec2<T>& rkVector)
{
  return GbVec2<T>(fScalar*rkVector[0],fScalar*rkVector[1]);
}

// output vector in human readable form
template<class T>
std::ostream&
operator<<(std::ostream& s, const GbVec2<T>& v)
{
  s<<"["<<v[0]<<", "<<v[1]<<"]";
  return s;
}

// read vector in same format as output
template<class T>
std::istream&
operator>>(std::istream& s, GbVec2<T>& v)
{
  char c;
  char dummy[3];
  T x,y;

  s>>c>>x>>dummy>>y>>c;
  v.setVec2(x,y);
  return s;
}

#ifdef __GNUG__

#include "GbVec2.in"
#include "GbVec2.T"

#else

#pragma instantiate GbVec2<float>
#pragma instantiate GbVec2<double>
#pragma instantiate GbVec2<float> operator*(float s, const GbVec2<float>& a)
#pragma instantiate GbVec2<double> operator*(double s, const GbVec2<double>& a)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbVec2<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbVec2<double>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbVec2<float>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbVec2<double>&)

#ifndef OUTLINE
#include "GbVec2.in"
#endif

#endif  // g++

#endif  // GBVEC2_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.3  2001/09/12 09:17:38  prkipfer
| added comparision operators to be able to use it with STL containers
|
| Revision 1.2  2001/08/16 16:53:21  prkipfer
| improved type safety for template parameter
|
| Revision 1.1  2001/06/20 16:14:39  prkipfer
| introduced 2D Vector class
|
|
+---------------------------------------------------------------------*/
