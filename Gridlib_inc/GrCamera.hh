/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRCAMERA_HH
#define GRCAMERA_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"

#include "GbVec3.hh"
#include "GbMatrix3.hh"
#include "GbPlane.hh"
#include "GrBound.hh"

class GrSpatial;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


class GrCamera 
{
public:
  // construction
  GrCamera ();

  // view frustum
  void setFrustum (float fNear, float fFar, float fLeft,
		   float fRight, float fTop, float fBottom);
  void setFrustumNear (float fNear);
  void setFrustumFar (float fFar);
  void getFrustum (float& rfNear, float& rfFar, float& rfLeft,
		   float& rfRight, float& rfTop, float& rfBottom) const;
  INLINE float getFrustumNear () const;
  INLINE float getFrustumFar () const;
  INLINE float getFrustumLeft () const;
  INLINE float getFrustumRight () const;
  INLINE float getFrustumTop () const;
  INLINE float getFrustumBottom () const;
  float getMaxCosSqrFrustumAngle () const;

  // view port (contained in [0,1]^2)
  void setViewPort (float fLeft, float fRight, float fTop, float fBottom);
  void getViewPort (float& rfLeft, float& rfRight, float& rfTop, float& rfBottom);
  INLINE float getViewPortLeft () const;
  INLINE float getViewPortRight () const;
  INLINE float getViewPortTop () const;
  INLINE float getViewPortBottom () const;

  // camera frame (world coordinates)
  //   default location = (0,0,0)
  //   default left = (1,0,0)
  //   default up = (0,1,0)
  //   default direction = (0,0,1)
  void setFrame (const GbVec3<float>& rkLocation, const GbVec3<float>& rkLeft,
		 const GbVec3<float>& rkUp, const GbVec3<float>& rkDirection);
  void setFrame (const GbVec3<float>& rkLocation, const GbMatrix3<float>& rkAxes);
  void setLocation (const GbVec3<float>& rkLocation);
  void setAxes (const GbVec3<float>& rkLeft, const GbVec3<float>& rkUp, const GbVec3<float>& rkDirection);
  void setAxes (const GbMatrix3<float>& rkAxes);
  INLINE const GbVec3<float>& getLocation () const;
  INLINE const GbVec3<float>& getLeft () const;
  INLINE const GbVec3<float>& getUp () const;
  INLINE const GbVec3<float>& getDirection () const;

  // update frustum, viewport, and frame
  void update ();

  // access to stack of world culling planes
  INLINE int getPlaneQuantity () const;
  INLINE const GbPlane<float>* getPlanes () const;
  void pushPlane (const GbPlane<float>& rkPlane);
  void popPlane ();

protected:
  // update callbacks
  friend class GrRenderer;
  virtual void onResize(int iWidth, int iHeight);
  virtual void onFrustumChange ();
  virtual void onViewPortChange ();
  virtual void onFrameChange ();

  friend class GrSpatial;
  INLINE void setPlaneState (unsigned int uiPlaneState);
  INLINE unsigned int getPlaneState () const;
  GbBool culled (const GrBound<float>& kWorldBound);

  // culling support in GrPortal::Draw
  friend class GrPortal;
  GbBool culled (int iVertexQuantity, const GbVec3<float>* akVertex,
		 GbBool bIgnoreNearPlane);

  // view frustum
  float frustumNear_, frustumFar_;
  float frustumLeft_, frustumRight_, frustumTop_, frustumBottom_;

  // view port
  float viewportLeft_, viewportRight_, viewportTop_, viewportBottom_;

  // camera frame (world location and coordinate axes)
  GbVec3<float> lookFrom_, vecLeft_, vecUp_, lookTo_;

  // Temporary values computed in OnFrustumChange that are needed if a
  // call is made to OnFrameChange.
  float coeffLeft_[2], coeffRight_[2], coeffBottom_[2], coeffTop_[2];

  // Bit flag to store whether or not a plane is active in the culling
  // system.  A bit of 1 means the plane is active, otherwise the plane is
  // inactive.  The first 6 planes are the view frustum planes.  Indices
  // are:  0 = near, 1 = far, 2 = left, 3 = right, 4 = bottom, 5 = top.
  // The system currently supports at most 32 culling planes.
  unsigned int planeState_;

  // world planes
  enum {
    CAM_LEFT_PLANE = 0,
    CAM_RIGHT_PLANE = 1,
    CAM_BOTTOM_PLANE = 2,
    CAM_TOP_PLANE = 3,
    CAM_FAR_PLANE = 4,
    CAM_NEAR_PLANE = 5,
    CAM_FRUSTUM_PLANES = 6,
    CAM_MAX_WORLD_PLANES = 32
  };
  int nbrPlanes_;
  GbPlane<float> worldPlane_[CAM_MAX_WORLD_PLANES];
};

#ifndef OUTLINE
#include "GrCamera.in"
#endif

#endif // GRCAMERA_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.4  2001/06/19 16:30:20  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.3  2001/06/18 11:12:24  prkipfer
| enabled world coord culling
|
| Revision 1.2  2001/06/15 09:02:30  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:27  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
