/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRPOINTLIGHT_HH
#define GRPOINTLIGHT_HH

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
class GrPointLight 
  : public GrLight<T>
{
public:
  GrPointLight ();

  INLINE virtual typename GrLight<T>::Type getType () const;

  INLINE GbVec3<T> location ();

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
  GbVec3<T> location_;
};

#ifdef __GNUG__

#include "GrPointLight.in"
#include "GrPointLight.T"

#else

#pragma instantiate GrPointLight<float>
#pragma instantiate GrPointLight<double>

#ifndef OUTLINE
#include "GrPointLight.in"
#endif

#endif  // g++

#endif // GRPOINTLIGHT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:04:44  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:35  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
