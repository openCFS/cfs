/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBFRUSTUM_HH
#define  GBFRUSTUM_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbVec3.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  Orthogonal frustum.  Let E be the origin, L be the left vector, U be
  the up vector, and D be the direction vector.  Let l > 0 and u > 0 be
  the extents in the L and U directions, respectively.  Let n and f be
  the extents in the D direction with 0 < n < f.  The four corners of the
  frustum in the near plane are E + s0*l*L + s1*u*U + n*D where |s0| =
  |s1| = 1 (four choices).  The four corners of the frustum in the far
  plane are E + (f/n)*(s0*l*L + s1*u*U) where |s0| = |s1| = 1 (four
  choices).
*/

template <class T>
class GbFrustum
{
public:
  GbFrustum ();

  INLINE GbVec3<T>& origin ();
  INLINE const GbVec3<T>& origin () const;

  INLINE GbVec3<T>& lVector ();
  INLINE const GbVec3<T>& lVector () const;
  INLINE GbVec3<T>& uVector ();
  INLINE const GbVec3<T>& uVector () const;
  INLINE GbVec3<T>& dVector ();
  INLINE const GbVec3<T>& dVector () const;

  INLINE T& lBound ();
  INLINE const T& lBound () const;
  INLINE T& uBound ();
  INLINE const T& uBound () const;
  INLINE T& dMin ();
  INLINE const T& dMin () const;
  INLINE T& dMax ();
  INLINE const T& dMax () const;

  void update ();
  INLINE T getDRatio () const;
  INLINE T getMTwoLF () const;
  INLINE T getMTwoUF () const;

  void computeVertices (GbVec3<T> akVertex[8]) const;

protected:
  GbVec3<T> origin_, lVector_, uVector_, dVector_;
  T lBound_, uBound_, dMin_, dMax_;

  // derived quantities
  T dRatio_;
  T mTwoLF_, mTwoUF_;
};

#ifdef __GNUG__

#include "GbFrustum.in"
#include "GbFrustum.T"

#else

#pragma instantiate GbFrustum<float>
#pragma instantiate GbFrustum<double>
// #pragma instantiate std::ostream& operator<<(std::ostream&, const GbFrustum<float>&)
// #pragma instantiate std::ostream& operator<<(std::ostream&, const GbFrustum<double>&)

#ifndef OUTLINE
#include "GbFrustum.in"
#endif

#endif  // g++

#endif  // GBFRUSTUM_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/06/18 10:43:04  prkipfer
| introduced Frustum object
|
|
+---------------------------------------------------------------------*/
