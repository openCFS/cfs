/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRDIRECTIONALLIGHT_HH
#define GRDIRECTIONALLIGHT_HH

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
class GrDirectionalLight 
  : public GrLight<T>
{
public:
  GrDirectionalLight ();

  INLINE virtual typename GrLight<T>::Type getType () const;

  INLINE GbVec3<T> direction ();

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
};

#ifdef __GNUG__

#include "GrDirectionalLight.in"
#include "GrDirectionalLight.T"

#else

#pragma instantiate GrDirectionalLight<float>
#pragma instantiate GrDirectionalLight<double>

#ifndef OUTLINE
#include "GrDirectionalLight.in"
#endif

#endif  // g++

#endif // GRDIRECTIONALLIGHT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/15 09:02:30  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:30  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
