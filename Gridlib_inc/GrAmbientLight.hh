/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRAMBIENTLIGHT_HH
#define GRAMBIENTLIGHT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbMatrix3.hh"
#include "GrColor.hh"
#include "GrLight.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GrAmbientLight 
  : public GrLight<T>
{
public:
  INLINE GrAmbientLight ();

  INLINE virtual typename GrLight<T>::Type getType () const;

  INLINE virtual void computeDiffuse (const GbMatrix3<T>& rkWorldRotate,
				      const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				      const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				      int uiQuantity, const GbBool* abVisible,
				      GrColor<T>* akDiffuse);

  INLINE virtual void computeSpecular (const GbMatrix3<T>& rkWorldRotate,
				       const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				       const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				       int uiQuantity, const GbBool* abVisible,
				       const GbVec3<T>& rkCameraModelLocation, GrColor<T>* akSpecular);
};

#ifdef __GNUG__

#include "GrAmbientLight.in"
#include "GrAmbientLight.T"

#else

#pragma instantiate GrAmbientLight<float>
#pragma instantiate GrAmbientLight<double>

#ifndef OUTLINE
#include "GrAmbientLight.in"
#endif

#endif  // g++

#endif // GRAMBIENTLIGHT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 09:02:29  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:25  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
