/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBQUATERNION_HH
#define  GBQUATERNION_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbMatrix3.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbQuaternion
{
public:
  //! Construction and destruction
  GbQuaternion();
  GbQuaternion(const T q[4]);
  GbQuaternion(T fW, T fX, T fY, T fZ);
  GbQuaternion(const GbQuaternion<T>& rkQ);

  //! Access
  INLINE T W() const;
  INLINE GbVec3<T> XYZ() const;
  INLINE void set(T w, T x, T y, T z);

  //! Conversion between quaternions, matrices, and angle-axes
  void fromRotationMatrix (const GbMatrix3<T>& kRot);
  void toRotationMatrix (GbMatrix3<T>& kRot) const;
  void fromAngleAxis (const T& rfAngle, const GbVec3<T>& rkAxis);
  void toAngleAxis (T& rfAngle, GbVec3<T>& rkAxis) const;
  void fromAxes (const GbVec3<T>* akAxis);
  void toAxes (GbVec3<T>* akAxis) const;

  //! Arithmetic operations
  INLINE GbQuaternion<T>& operator= (const GbQuaternion<T>& rkQ);
  INLINE GbQuaternion<T> operator+ (const GbQuaternion<T>& rkQ) const;
  INLINE GbQuaternion<T> operator- (const GbQuaternion<T>& rkQ) const;
  INLINE GbQuaternion<T> operator* (const GbQuaternion<T>& rkQ) const;
  INLINE GbQuaternion<T> operator* (T fScalar) const;
  INLINE GbQuaternion<T> operator- () const;
  friend GbQuaternion<T> operator*(T fScalar, const GbQuaternion<T>& rkQ);

  //! Functions of a quaternion
  INLINE T dot (const GbQuaternion<T>& rkQ) const;  // dot product
  INLINE T norm () const;  // squared-length
  INLINE GbQuaternion<T> inverse () const;  // apply to non-zero quaternion
  INLINE GbQuaternion<T> unitInverse () const;  // apply to unit-length quaternion
  GbQuaternion<T> exp () const;
  GbQuaternion<T> log () const;

  //! Rotation of a vector by a quaternion
  GbVec3<T> operator* (const GbVec3<T>& rkVector) const;

  //! Spherical linear interpolation
  static GbQuaternion<T> slerp (T fT, const GbQuaternion<T>& rkP, const GbQuaternion<T>& rkQ);

  static GbQuaternion<T> slerpExtraSpins (T fT, const GbQuaternion<T>& rkP, 
					  const GbQuaternion<T>& rkQ, int iExtraSpins);

  //! Setup for spherical quadratic interpolation
  static void intermediate (const GbQuaternion<T>& rkQ0, const GbQuaternion<T>& rkQ1, const GbQuaternion<T>& rkQ2,
			    GbQuaternion<T>& rka, GbQuaternion<T>& rkB);

  //! Spherical quadratic interpolation
  static GbQuaternion<T> squad (T fT, const GbQuaternion<T>& rkP, const GbQuaternion<T>& rkA, 
				const GbQuaternion<T>& rkB, const GbQuaternion<T>& rkQ);

  //! Special values
  static const GbQuaternion<T> ZERO;
  static const GbQuaternion<T> IDENTITY;

  //! Cutoff for sine near zero	
  static const T EPSILON;

  //! This operator displays the values
  friend std::ostream& operator<<(std::ostream&, const GbQuaternion<T>&);
  friend std::istream& operator>>(std::istream&, GbQuaternion<T>&);

private:
  T w_, x_, y_, z_;
};

template <class T>
GbQuaternion<T> 
operator* (T fScalar, const GbQuaternion<T>& rkQ)
{
  GbVec3<T> v = rkQ.XYZ();
  return GbQuaternion<T>(fScalar*rkQ.W(),fScalar*v[0],fScalar*v[1],fScalar*v[2]);
}

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbQuaternion<T>& q)
{
  GbVec3<T> v = q.XYZ();
  s<<"("<<q.W()<<", "<<v[0]<<", "<<v[1]<<", "<<v[2]<<")";
  return s;
}

template<class T>
std::istream&
operator>>(std::istream& s, GbQuaternion<T>& q)
{
  char c;
  char dummy[3];
  T w,x,y,z;
  
  s>>c>>w>>dummy>>x>>dummy>>y>>dummy>>z>>c;
  q.set(w,x,y,z);
  return s;
}

#ifdef __GNUG__

#include "GbQuaternion.in"
#include "GbQuaternion.T"

#else

#pragma instantiate GbQuaternion<float>
#pragma instantiate GbQuaternion<double>
#pragma instantiate GbQuaternion<float>  operator* (float fScalar, const GbQuaternion<float>& rkQ)
#pragma instantiate GbQuaternion<double> operator* (double fScalar, const GbQuaternion<double>& rkQ)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbQuaternion<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbQuaternion<double>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbQuaternion<float>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbQuaternion<double>&)

#ifndef OUTLINE
#include "GbQuaternion.in"
#endif

#endif  // g++

#endif  // GBQUATERNION_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2001/08/16 16:53:20  prkipfer
| improved type safety for template parameter
|
| Revision 1.2  2001/06/15 08:38:10  prkipfer
| removed useless qualifiers
|
| Revision 1.1  2001/01/02 14:49:47  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
