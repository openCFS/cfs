/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRGEOPRIMITIVECAMERA_HH
#define GRGEOPRIMITIVECAMERA_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrCamera.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrGeoPrimitiveCamera 
  : public GrCamera
{
public:
  // construction and destruction
  GrGeoPrimitiveCamera (float fWidth, float fHeight);

  INLINE float getFarDivideNear () const;
  INLINE float getFarDivideFarMinusNear () const;

  void computeModelToViewTransform (const GbMatrix3<float>& rkWRot,
				    const GbVec3<float>& rkWTrn, float fWScale,
				    GbMatrix3<float>& rkModelToViewMat, GbVec3<float>& rkModelToViewTrn);

  virtual void onFrustumChange ();
  virtual void onViewPortChange ();
  virtual void onFrameChange ();

protected:
  GrGeoPrimitiveCamera ();

  float floatWidth_, floatHeight_;

  // frustum-related values
  float inverseNear_;           // 1/N
  float farDivideNear_;         // F/N
  float farDivideFarMinusNear_; // F/(F-N)
  float twoDivideRightMinusLeft_;             // 2/(R-L)
  float rightPlusLeftDivideNormalTimesRightMinusLeft_;            // (R+L)/(N*(R-L))
  float twoDivideTopMinusBottom_;             // 2/(T-B)
  float topPlusBottomDivideNormalTimesTopMinusBottom_;            // (T+B)/(N*(T-B))

  // viewport (left, bottom, width, height)
  int intViewportLeft_, intViewportBottom_, intViewportWidth_, intViewportHeight_;

  // projection-times-rotation matrix
  void updateProjRotMatrix ();
  GbMatrix3<float> projectionTimesRotationMatrix_;
};

#ifndef OUTLINE
#include "GrGeoPrimitiveCamera.in"
#endif

#endif // GRGEOPRIMITIVECAMERA_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/02/13 13:29:31  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
