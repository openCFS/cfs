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
//  Abwandlung des GbVec3
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
  GbVec2 (T x, T y);
  GbVec2 (const T p[2]);
  GbVec2 (const GbVec2<T>& p);

  //! Destructor
  ~GbVec2();

  // assignment and comparison
  INLINE GbVec2<T>& operator= (const GbVec2<T>& p);
  INLINE GbBool     operator== (const GbVec2<T>& p) const;
  INLINE GbBool     operator!= (const GbVec2<T>& p) const;
  INLINE GbVec2<T>& operator+= (const GbVec2<T>& rkVector);
  INLINE GbVec2<T>& operator-= (const GbVec2<T>& rkVector);
  INLINE GbVec2<T>& operator*= (int s);
  INLINE GbVec2<T>& operator*= (T fScalar);
  INLINE GbVec2<T>& operator/= (int s);
  INLINE GbVec2<T>& operator/= (T fScalar);
  INLINE GbVec2<T> operator+ (const GbVec2<T>& rkVector) const;
  INLINE GbVec2<T> operator- (const GbVec2<T>& rkVector) const;
  INLINE GbVec2<T> operator- () const;
  INLINE GbVec2<T> operator* (int s) const;
  INLINE GbVec2<T> operator* (T fScalar) const;
  INLINE GbVec2<T> operator/ (int s) const;
  INLINE GbVec2<T> operator/ (T fScalar) const;
  friend GbVec2<T> operator* (T fScalar, const GbVec2<T>& rkVector);

  // vector operations
  INLINE T getNorm() const;
  INLINE T getSquareNorm() const;
  INLINE T dot (const GbVec2<T>& rkVector) const;
  // tolerance max. 1e-06
  INLINE T normalize(T fTolerance = std::numeric_limits<T>::epsilon());
  INLINE GbVec2<T> cross () const;  // returns (y,-x)
  INLINE GbVec2<T> unitCross () const;  // returns (y,-x)/sqrt(x*x+y*y)

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

  // Gram-Schmidt orthonormalization.
  INLINE static void orthonormalize (GbVec2<T> v[2]);

  friend std::ostream& operator<<(std::ostream& s, const GbVec2<T>& v);

  // special points
  static const GbVec2<T> ZERO;
  static const GbVec2<T> UNIT_X;
  static const GbVec2<T> UNIT_Y;

private:
  T xy[2];
};

template <class T>
GbVec2<T> 
operator* (T fScalar, const GbVec2<T>& rkVector)
{
  return GbVec2<T>(fScalar*rkVector[0],fScalar*rkVector[1]);
}

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbVec2<T>& v)
{
  s<<typeid(v).name()<<": p("<<v[0]<<","<<v[1]<<")"<<std::endl;
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

#ifndef OUTLINE
#include "GbVec2.in"
#endif

#endif  // g++

#endif  // GBVEC2_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/06/20 16:14:39  prkipfer
| introduced 2D Vector class
|
|
+---------------------------------------------------------------------*/
