/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBPLANE_HH
#define  GBPLANE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbVec3.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbPlane
{
public:
  // construction
  INLINE GbPlane ();
  INLINE GbPlane (const GbVec3<T>& rkNormal, T fConstant);
  INLINE GbPlane (const GbVec3<T>& rkNormal, const GbVec3<T>& rkPoint);
  GbPlane (const GbVec3<T>& rkPoint0, const GbVec3<T>& rkPoint1, const GbVec3<T>& rkPoint2);

  INLINE void setNormal(GbVec3<T>& n);
  INLINE const GbVec3<T> getNormal () const;

  INLINE void setConstant(T c);
  INLINE const T getConstant () const;

  //! access plane P as P[0] = N.x, P[1] = N.y, P[2] = N.z, P[3] = c
  INLINE T operator[] (int i) const;

  /*! The "positive side" of the plane is the half space to which the plane
      normal points.  The "negative side" is the other half space.  The flag
      "no side" indicates the plane itself.
  */
  typedef enum _side {
    NO_SIDE,
    POSITIVE_SIDE,
    NEGATIVE_SIDE
  } Side;

  Side whichSide (const GbVec3<T>& rkPoint) const;

  /*! This is a pseudodistance.  The sign of the return value is positive if
      the point is on the positive side of the plane, negative if the point
      is on the negative side, and zero if the point is on the plane.  The
      absolute value of the return value is the true distance only when the
      plane normal is a unit length vector.
  */
  INLINE T distanceTo (const GbVec3<T>& rkPoint) const;

protected:
  // plane is Dot(m_normal,(x,y,z)) = constant_
  GbVec3<T> normal_;
  T constant_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbPlane<T>& v)
{
  s<<typeid(v).name()<<": normal("<<v.getNormal()[0]<<","<<v.getNormal()[1]<<","<<v.getNormal()[2]<<") ";
  s<<"const("<<v.getConstant()<<")"<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GbPlane.in"
#include "GbPlane.T"

#else

#pragma instantiate GbPlane<float>
#pragma instantiate GbPlane<double>
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbPlane<float>& v)
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbPlane<double>& v)

#ifndef OUTLINE
#include "GbPlane.in"
#endif

#endif  // g++

#endif  // GBPLANE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/01/09 11:12:29  prkipfer
| fixed reference copy bug
|
| Revision 1.1  2001/01/02 14:55:04  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
