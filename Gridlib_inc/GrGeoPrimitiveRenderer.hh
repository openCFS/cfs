/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRGEOPRIMITIVERENDERER_HH
#define GRGEOPRIMITIVERENDERER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrRenderer.hh"
#include "GrGeoPrimitiveCamera.hh"
#include "GrGeoPrimitiveRendererMesh.hh"
#include "GoTriangleElement.hh"
#include "GoTetrahedronElement.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrGeoPrimitiveRenderer 
  : public GrRenderer
{
public:
  virtual ~GrGeoPrimitiveRenderer ();

  virtual void draw (const GrGeoPrimitiveMesh& rkTriMesh);
  virtual void draw (const GoTriangleElement& tri);
  virtual void draw (const GoTetrahedronElement& tri);

protected:
  // abstract base class for software renderers
  GrGeoPrimitiveRenderer (int iWidth, int iHeight);
  GrGeoPrimitiveRenderer () {}

  // clipping
  GbBool clip ();
  GbBool clipTriMesh ();
  GbBool clipNear ();
  GbBool clipFar ();
  GbBool clipLeft ();
  GbBool clipRight ();
  GbBool clipTop ();
  GbBool clipBottom ();

  // projection
  void projectVertices ();

  // rasterization
  virtual void rasterizeEdges () = 0;
  virtual void rasterizeTriangles () = 0;

  void computeEdgeBuffer (int iX0, int iY0, float* afLinAttr0, float* afPerAttr0,
			  int iX1, int iY1, float* afLinAttr1, float* afPerAttr1,
			  int* aiXBuffer, float** aafLinAttrBuffer,
			  float** aafPerAttrBuffer);

  GbBool computeEdgeBuffers ();

  // renderer state management
  virtual void setState (GrAlphaState* pkState);
  virtual void setState (GrDitherState* pkState);
  virtual void setState (GrFogState* pkState);
  virtual void setState (GrLightState* pkState);
  virtual void setState (GrMaterialState* pkState);
  virtual void setState (GrShadeState* pkState);
  virtual void setState (GrTextureState* pkState);
  virtual void setState (GrVertexColorState* pkState);
  virtual void setState (GrWireframeState* pkState);
  virtual void setState (GrZBufferState* pkState);

  INLINE GbBool alphaBlendActive ();
  INLINE GbBool zBufferActive ();

  // window parameters
  float floatWidth_, floatHeight_, floatXScale_, floatYScale_;

  // upcast GrCamera to GrGeoPrimitiveCamera
  INLINE GrGeoPrimitiveCamera* getCamera ();

  // renderer's local storage for trimesh (grows to maximum size needed)
  GrGeoPrimitiveRendererMesh mesh_;

  // edge buffering for scan line setup
  int edgeYMin_, edgeYMax_;
  int edgeX_[3], edgeY_[3];
  int* edgeXMin_;
  int* edgeXMax_;
  enum { MAX_ATTRIBUTES = 6 };  // u, v, w, r, g, b
  int nbrScanlineAttributes_;  // in [0..MAX_ATTRIBUTES-1]
  float** scanlineAttributeMin_;
  float** scanlineAttributeMax_;
  float* scanlineAttribute_[3];  // attribute packing
  int nbrPerspectiveAttributes_;  // in [0..MAX_ATTRIBUTES-1]
  float** perspectiveAttributeMin_;
  float** perspectiveAttributeMax_;
  float* perspectiveAttribute_[3];  // attribute packing

  // current render state
  GrAlphaState* currentAlphaState_;
  GrDitherState* currentDitherState_;
  GrFogState* currentFogState_;
  GrLightState* currentLightState_;
  GrMaterialState* currentMaterialState_;
  GrShadeState* currentShadeState_;
  GrTextureState* currentTextureState_;
  GrVertexColorState* currentVertexColorState_;
  GrWireframeState* currentWireframeState_;
  GrZBufferState* currentZBufferState_;
};

#ifndef OUTLINE
#include "GrGeoPrimitiveRenderer.in"
#endif

#endif // GRGEOPRIMITIVERENDERER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/06/18 11:29:42  prkipfer
| adapted to GrGeoPrimitiveRendererMesh
|
| Revision 1.2  2001/06/15 09:07:39  prkipfer
| introduced draw method for tets, switched to STL containers, added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 14:06:43  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
