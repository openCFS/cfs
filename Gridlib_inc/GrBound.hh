/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRBOUND_HH
#define GRBOUND_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GbVec3.hh"
#include "GbPlane.hh"
#include "GbMatrix3.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GrBound
{
public:
  INLINE GrBound ();

  INLINE void setCenter(GbVec3<T>& c);
  INLINE void setRadius(T& r);
  INLINE const GbVec3<T> getCenter () const;
  INLINE const T getRadius () const;

  // access bound B as B[0] = C.x, B[1] = C.y, B[2] = C.z, B[3] = r
  INLINE T operator[] (int i) const;

  INLINE GrBound<T>& operator= (const GrBound<T>& rkBound);
  GrBound<T>& operator+= (const GrBound<T>& rkBound);

  void computeFromData (int uiQuantity, const GbVec3<T>* akData);

  GrBound<T> transformBy (const GbMatrix3<T>& rkRotate,
			  const GbVec3<T>& rkTranslate, T fScale) const;

  typename GbPlane<T>::Side whichSide (const GbPlane<T>& rkPlane) const;

  static const T TOLERANCE;

protected:
  GbVec3<T> center_;
  T radius_;
};

#ifdef __GNUG__

#include "GrBound.in"
#include "GrBound.T"

#else

#pragma instantiate GrBound<float>
#pragma instantiate GrBound<double>

#ifndef OUTLINE
#include "GrBound.in"
#endif

#endif  // g++


#endif // GRBOUND_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 09:02:29  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:26  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
