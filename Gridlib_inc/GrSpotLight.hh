/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRSPOTLIGHT_HH
#define GRSPOTLIGHT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbMatrix3.hh"
#include "GrColor.hh"
#include "GrPointLight.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GrSpotLight 
  : public GrPointLight<T>
{
public:
  GrSpotLight ();

  INLINE virtual typename GrLight<T>::Type getType () const;

  INLINE GbVec3<T> direction ();
  INLINE T exponent ();
  void setAngle (T fAngle);
  INLINE T getAngle () const;

  virtual void computeDiffuse (const GbMatrix3<T>& rkWorldRotate,
			       const GbVec3<T>& rkWorldTranslate, T fWorldScale,
			       const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
			       int uiQuantity, const GbBool* abVisible,
			       GrColor<T>* akDiffuse);

  virtual void computeSpecular (const GbMatrix3<T>& rkWorldRotate,
				const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				int uiQuantity, const GbBool* abVisible,
				const GbVec3<T>& rkCameraModelLocation, GrColor<T>* akSpecular);

protected:
  GbVec3<T> direction_;
  T angle_, cosSqr_, sinSqr_;
  T exponent_;
};

#ifdef __GNUG__

#include "GrSpotLight.in"
#include "GrSpotLight.T"

#else

#pragma instantiate GrSpotLight<float>
#pragma instantiate GrSpotLight<double>

#ifndef OUTLINE
#include "GrSpotLight.in"
#endif

#endif  // g++

#endif // GRSPOTLIGHT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:58  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:04:44  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:36  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
