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
  GbVec3(T x, T y, T z);
  GbVec3(const T p[3]);
  GbVec3(const GbVec3<T> &p);

  //! Destructor
  ~GbVec3();

  //! Operators
  INLINE GbVec3<T>& operator=(const GbVec3<T>& p);
  INLINE GbBool     operator==(const GbVec3<T>& p) const;
  INLINE GbBool     operator!=(const GbVec3<T>& p) const;
  INLINE GbVec3<T>& operator+=(const GbVec3<T>& p);
  INLINE GbVec3<T>& operator-=(const GbVec3<T>& p);
  INLINE GbVec3<T>& operator*=(int s);
  INLINE GbVec3<T>& operator*=(T s);
  INLINE GbVec3<T>& operator/=(int s);
  INLINE GbVec3<T>& operator/=(T s);
  INLINE GbVec3<T>  operator+(const GbVec3<T>& p) const;
  INLINE GbVec3<T>  operator-(const GbVec3<T>& p) const;
  INLINE GbVec3<T>  operator-() const;
  INLINE GbVec3<T>  operator*(int s) const;
  INLINE GbVec3<T>  operator*(T s) const;
  INLINE GbVec3<T>  operator/(int s) const;
  INLINE GbVec3<T>  operator/(T s) const;
  friend GbVec3<T>  operator*(T s, const GbVec3<T>& a);

  //! Cross product
  INLINE GbVec3<T>  operator%(const GbVec3<T>& a) const;
  INLINE GbVec3<T>  cross(const GbVec3<T> &a) const;
  INLINE GbVec3<T>  unitCross(const GbVec3<T>& a) const;

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
  // tolerance max. 1e-06
  INLINE T    normalize(T tolerance = std::numeric_limits<T>::epsilon());

  //! Vector computation methods
  INLINE        T operator*(const GbVec3<T> &a) const;
  INLINE        T dot(const GbVec3<T> &a) const;
  INLINE static T det(const GbVec3<T> &a, const GbVec3<T> &b, const GbVec3<T> &c);
  INLINE static T angle(const GbVec3<T> &a, const GbVec3<T> &b);

  //! Orthonormal stuff
  INLINE static void orthonormalize(GbVec3<T> v[3]);
  INLINE static void generateOrthonormalBasis(GbVec3<T>& u, GbVec3<T>& v, GbVec3<T>& w, GbBool unitLenW = true);

  //! This operator displays the coordinate values of the vector
  friend std::ostream& operator<<(std::ostream&, const GbVec3<T>&);

  //! special vectors
  static const GbVec3<T> ZERO;
  static const GbVec3<T> UNIT_X;
  static const GbVec3<T> UNIT_Y;
  static const GbVec3<T> UNIT_Z;
    
private:
  //! Private data storage
  T xyz[3];
};

template<class T>
GbVec3<T>
operator* (T s, const GbVec3<T>& a)
{
  GbVec3<T> kProd;
  kProd[0] = s*a[0];
  kProd[1] = s*a[1];
  kProd[2] = s*a[2];
  return kProd;
}

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbVec3<T>& v)
{
  s<<typeid(v).name()<<": p("<<v[0]<<","<<v[1]<<","<<v[2]<<")"<<std::endl;
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

#ifndef OUTLINE
#include "GbVec3.in"
#endif

#endif  // g++

#endif  // GBVEC3_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
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
