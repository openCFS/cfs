/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

// ----------------------------------------------------------------------------
//
//  KLASSE  Vec3
//
//  Ein einfacher Punkt im 3D mit Zugriffsfunktionen wie bei einem Vektor   
//  Original Source: MoPoint.hh (C) 1996-1998 Swen Campagna
//
//-----------------------------------------------------------------------------

#ifndef  GBVEC3_HH
#define  GBVEC3_HH

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

template <class T>
class GbVec3
{
public:
  //! Various constructors: this class is intended for floating point
  //! numbers. Note the special treatment of integers
  GbVec3();
  explicit GbVec3(T s);
  GbVec3(T x, T y, T z);
  GbVec3(const T p[3]);
  GbVec3(const GbVec3<T> &p);

  //! Destructor
  ~GbVec3();

  //! Comparision operators (support fuzzy arithmetic when EPSILON > 0)
  INLINE GbBool operator==(const GbVec3<T>& p) const;
  INLINE GbBool operator!=(const GbVec3<T>& p) const;
  INLINE GbBool operator< (const GbVec3<T>& p) const;
  INLINE GbBool operator<=(const GbVec3<T>& p) const;
  INLINE GbBool operator> (const GbVec3<T>& p) const;
  INLINE GbBool operator>=(const GbVec3<T>& p) const;
  
  //! Arithmetic operations
  INLINE GbVec3<T>& operator+=(const GbVec3<T>& p);
  INLINE GbVec3<T>& operator+=(const T& s);
  INLINE GbVec3<T>& operator-=(const GbVec3<T>& p);
  INLINE GbVec3<T>& operator-=(const T& s);
  INLINE GbVec3<T>& operator*=(const GbVec3<T>& p);
  INLINE GbVec3<T>& operator*=(const int& s);
  INLINE GbVec3<T>& operator*=(const T& s);
  INLINE GbVec3<T>& operator/=(const GbVec3<T>& p);
  INLINE GbVec3<T>& operator/=(const int& s);
  INLINE GbVec3<T>& operator/=(const T& s);

  //! Binary operators
  INLINE GbVec3<T>  operator+(const GbVec3<T>& p) const;
  INLINE GbVec3<T>  operator+(const T& p) const;
  INLINE GbVec3<T>  operator-(const GbVec3<T>& p) const;
  INLINE GbVec3<T>  operator-(const T& p) const;
  INLINE GbVec3<T>  operator*(const GbVec3<T>& p) const;
  INLINE GbVec3<T>  operator*(const int& s) const;
  INLINE GbVec3<T>  operator*(const T& s) const;
  INLINE GbVec3<T>  operator/(const GbVec3<T>& p) const;
  INLINE GbVec3<T>  operator/(const int& s) const;
  INLINE GbVec3<T>  operator/(const T& s) const;

  friend GbVec3<T>  operator*(T s, const GbVec3<T>& a);

  //! Unary operators
  INLINE GbVec3<T>  operator-() const;

  //! Vector product
  INLINE GbVec3<T>  operator^(const GbVec3<T>& a) const;
  INLINE GbVec3<T>  cross(const GbVec3<T> &a) const;
  INLINE GbVec3<T>  unitCross(const GbVec3<T>& a) const;

  //! Scalar product
  INLINE        T operator|(const GbVec3<T> &a) const;
  INLINE        T dot(const GbVec3<T> &a) const;

  //! Access operators
  INLINE T &        operator[] (int i);
  INLINE const T &  operator[] (int i) const;

  //! Read methods
  INLINE const T * getVec3() const;
  INLINE void      getVec3(T &x, T &y, T &z) const;
  
  //! Direct manipulation
  INLINE void setX(T x);
  INLINE void setY(T y);
  INLINE void setZ(T z);
  INLINE void setVec3(const T f[]);
  INLINE void setVec3(T x,T y,T z);

  //! Normal vector handling
  INLINE T    getNorm() const;
  INLINE T    getSquareNorm() const;
  INLINE T    normalize(T tolerance = std::numeric_limits<T>::epsilon());
  INLINE GbVec3<T> getNormalized(T tolerance = std::numeric_limits<T>::epsilon()) const;

  //! Vector computation methods
  INLINE static T det(const GbVec3<T> &a, const GbVec3<T> &b, const GbVec3<T> &c);
  INLINE static T angle(const GbVec3<T> &a, const GbVec3<T> &b);

  //! Orthonormal stuff
  INLINE static void orthonormalize(GbVec3<T> v[/*3*/]);
  INLINE static void generateOrthonormalBasis(GbVec3<T>& u, GbVec3<T>& v, GbVec3<T>& w, GbBool unitLenW = true);

  //! Projection normal to a vector
  INLINE GbVec3<T>	  getOrthogonalVector() const;
  INLINE const GbVec3<T>& projectNormalTo(const GbVec3<T> &v);
  INLINE GbVec3<T>	  getReflectedAt(const GbVec3<T>& n) const;
  INLINE GbVec3<T>	  getRefractedAt(const GbVec3<T>& n, T index, GbBool& total_reflection) const;
  
  //! This operator displays the coordinate values of the vector
  friend std::ostream& operator<<(std::ostream&, const GbVec3<T>&);
  friend std::istream& operator>>(std::istream&, GbVec3<T>&);

  // Change dimensions
//  INLINE static GbVec3<T> from3DProj(GbVec3Proj<T> &v ) ;

  //! special vectors
  static const GbVec3<T> ZERO;
  static const GbVec3<T> UNIT_X;
  static const GbVec3<T> UNIT_Y;
  static const GbVec3<T> UNIT_Z;

  //! to control exactness of comparisions set this to > 0
  static T EPSILON;
    
private:
  //! Private data storage
  T xyz[3];
};

template<class T>
GbVec3<T>
operator* (T s, const GbVec3<T>& a)
{
  T x = s*a[0];
  T y = s*a[1];
  T z = s*a[2];
  return GbVec3<T>(x,y,z);
}

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbVec3<T>& v)
{
  s<<"["<<v[0]<<", "<<v[1]<<", "<<v[2]<<"]";
  return s;
}

template<class T>
std::istream&
operator>>(std::istream& s, GbVec3<T>& v)
{
  char c;
  char dummy[3];
  T x,y,z;

  s>>c>>x>>dummy>>y>>dummy>>z>>c;
  v.setVec3(x,y,z);
  return s;
}

#ifdef __GNUG__

#include "GbVec3.in"
#include "GbVec3.T"

#else

#pragma instantiate GbVec3<float>
#pragma instantiate GbVec3<double>
#pragma instantiate GbVec3<float> operator*(float s, const GbVec3<float>& a)
#pragma instantiate GbVec3<double> operator*(double s, const GbVec3<double>& a)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbVec3<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbVec3<double>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbVec3<float>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbVec3<double>&)

#ifndef OUTLINE
#include "GbVec3.in"
#endif

#endif  // g++

#endif  // GBVEC3_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.10  2001/09/12 09:17:39  prkipfer
| added comparision operators to be able to use it with STL containers
|
| Revision 1.9  2001/08/16 16:53:22  prkipfer
| improved type safety for template parameter
|
| Revision 1.8  2001/06/18 10:48:21  prkipfer
| beautified inline syntax
|
| Revision 1.7  2001/02/13 10:58:39  prkipfer
| exchanged fixed numeric tolerances by machine dependent ones
|
| Revision 1.6  2001/01/02 14:46:10  prkipfer
| extended functionality to work with new classes
|
| Revision 1.5  2000/07/28 15:24:26  uflabsik
| friend bug fixed
|
| Revision 1.4  2000/07/21 10:40:44  prkipfer
| dropped g++ support
|
| Revision 1.3  2000/07/20 13:32:38  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.2  2000/06/14 15:39:13  prkipfer
| improved base classes and added funcstruct processing for mesh
|
| Revision 1.1.1.1  2000/06/08 16:24:44  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
