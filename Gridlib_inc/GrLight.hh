/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRLIGHT_HH
#define GRLIGHT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbMatrix3.hh"
#include "GbVec3.hh"
#include "GrColor.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GrLight 
{
public:
  INLINE virtual ~GrLight ();

  typedef enum _type {
    LT_AMBIENT,
    LT_DIRECTIONAL,
    LT_POINT,
    LT_SPOT,
    LT_QUANTITY
  } Type;

  virtual Type getType () const = 0;

  INLINE GrColor<T> ambient ();
  INLINE GrColor<T> diffuse ();
  INLINE GrColor<T> specular ();

  INLINE GbBool attenuate ();
  INLINE T constant ();
  INLINE T linear ();
  INLINE T quadratic ();

  INLINE GbBool isOn ();
  INLINE int getId();
  INLINE T intensity ();

  virtual void computeDiffuse (const GbMatrix3<T>& rkWorldRotate,
			       const GbVec3<T>& rkWorldTranslate, T fWorldScale,
			       const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
			       int uiQuantity, const GbBool* abVisible,
			       GrColor<T>* akDiffuse) = 0;

  virtual void computeSpecular (const GbMatrix3<T>& rkWorldRotate,
				const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				int uiQuantity, const GbBool* abVisible,
				const GbVec3<T>& rkCameraModelLocation, GrColor<T>* akSpecular) = 0;

protected:
  // abstract base class
  GrLight ();

  GrColor<T> ambient_;
  GrColor<T> diffuse_;
  GrColor<T> specular_;
  T constant_;
  T linear_;
  T quadratic_;
  T intensity_;
  GbBool attenuate_;
  GbBool on_;
  int id_;
  static int nextId_;
};

#ifdef __GNUG__

#include "GrLight.in"
#include "GrLight.T"

#else

#pragma instantiate GrLight<float>
#pragma instantiate GrLight<double>

#ifndef OUTLINE
#include "GrLight.in"
#endif

#endif  // g++

#endif // GRLIGHT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:04:43  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:33  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
