/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRRENDERER_HH
#define GRRENDERER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GrCamera.hh"
#include "GrColor.hh"
#include "GoTriangleElement.hh"
#include "GoTetrahedronElement.hh"
#include "GrAlphaState.hh"
#include "GrDitherState.hh"
#include "GrFogState.hh"
#include "GrLightState.hh"
#include "GrMaterialState.hh"
#include "GrShadeState.hh"
#include "GrTextureState.hh"
#include "GrVertexColorState.hh"
#include "GrWireframeState.hh"
#include "GrZBufferState.hh"

class GrGeoPrimitiveMesh;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrRenderer 
{
public:
  virtual ~GrRenderer ();

  INLINE int getWidth () const;
  INLINE int getHeight () const;
  virtual void resize (int iWidth, int iHeight);

  GrCamera* setCamera (GrCamera* pkCamera);
  INLINE GrCamera* getCamera ();

  void setState (GrRenderState* aspkState[]);

  virtual INLINE void setBackgroundColor (const GrColor<float>& rkColor);
  INLINE GrColor<float> getBackgroundColor () const;
  virtual void clearBackBuffer () = 0;
  virtual void displayBackBuffer () = 0;
  virtual void clearZBuffer () = 0;
  virtual void clearBuffers () = 0;

  virtual void draw (const GoTriangleElement& tri) {}
  virtual void draw (const GoTetrahedronElement& tri) {}
  virtual void draw (const GrGeoPrimitiveMesh& rkMesh) {}
  virtual void draw (int iX, int iY, const char* acText) {}

protected:
  // abstract base class
  GrRenderer (int iWidth, int iHeight);
  GrRenderer ();

  // state management
  virtual void initializeState ();
  virtual void setState (GrAlphaState* pkState) = 0;
  virtual void setState (GrDitherState* pkState) = 0;
  virtual void setState (GrFogState* pkState) = 0;
  virtual void setState (GrLightState* pkState) = 0;
  virtual void setState (GrMaterialState* pkState) = 0;
  virtual void setState (GrShadeState* pkState) = 0;
  virtual void setState (GrTextureState* pkState) = 0;
  virtual void setState (GrVertexColorState* pkState) = 0;
  virtual void setState (GrWireframeState* pkState) = 0;
  virtual void setState (GrZBufferState* pkState) = 0;

  // window parameters
  int width_, height_, wXh_;
  GrColor<float> backgroundColor_;

  // camera for establishing view frustum
  GrCamera* camera_;
};

#ifndef OUTLINE
#include "GrRenderer.in"
#endif

#endif // GRRENDERER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/15 09:19:34  prkipfer
| introduced resizing
|
| Revision 1.1  2001/02/13 14:06:47  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
