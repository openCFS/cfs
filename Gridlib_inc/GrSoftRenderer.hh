/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRSOFTRENDERER_HH
#define GRSOFTRENDERER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrRenderer.hh"
#include "GrGeoPrimitiveRenderer.hh"
//#include "GrSoftColor.hh"
#include "GrColor.hh"
#include "GrSoftImage.hh"
#include "GrImage.hh"
#include "GrRayCaster.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


class GrSoftRenderer 
  : public GrGeoPrimitiveRenderer
{
public:
  GrSoftRenderer (int iWidth, int iHeight);
  virtual ~GrSoftRenderer ();

  virtual void setBackgroundColor (const GrColor<float>& rkColor);
  virtual void clearBackBuffer ();
  virtual void clearZBuffer ();
  virtual void clearBuffers ();

  virtual void draw(const GrRayCaster& rc);

protected:
  GrSoftRenderer () {}

  // rasterizer selection
  int getRasterizerIndex ();

  // edge rasterizing
  typedef void (GrSoftRenderer::*ERasterizer)(unsigned int, unsigned int);
  ERasterizer edgeRasterizer_;
  static ERasterizer edgeRasterizers_[32];
  virtual void rasterizeEdges ();

  void drawEdgeNULL (unsigned, unsigned);
  void drawEdgeC (unsigned uiV0, unsigned uiV1);
  void drawEdgeT (unsigned uiV0, unsigned uiV1);
  void drawEdgeTC (unsigned uiV0, unsigned uiV1);
  void drawEdgePT (unsigned uiV0, unsigned uiV1);
  void drawEdgePTC (unsigned uiV0, unsigned uiV1);
  void drawEdgeZC (unsigned uiV0, unsigned uiV1);
  void drawEdgeZT (unsigned uiV0, unsigned uiV1);
  void drawEdgeZTC (unsigned uiV0, unsigned uiV1);
  void drawEdgeZPT (unsigned uiV0, unsigned uiV1);
  void drawEdgeZPTC (unsigned uiV0, unsigned uiV1);
  void drawEdgeAC (unsigned uiV0, unsigned uiV1);
  void drawEdgeAT (unsigned uiV0, unsigned uiV1);
  void drawEdgeATC (unsigned uiV0, unsigned uiV1);
  void drawEdgeAPT (unsigned uiV0, unsigned uiV1);
  void drawEdgeAPTC (unsigned uiV0, unsigned uiV1);
  void drawEdgeAZC (unsigned uiV0, unsigned uiV1);
  void drawEdgeAZT (unsigned uiV0, unsigned uiV1);
  void drawEdgeAZTC (unsigned uiV0, unsigned uiV1);
  void drawEdgeAZPT (unsigned uiV0, unsigned uiV1);
  void drawEdgeAZPTC (unsigned uiV0, unsigned uiV1);

  // triangle rasterizers
  typedef void (GrSoftRenderer::*TRasterizer)(unsigned int, unsigned int, unsigned int);
  TRasterizer triangleRasterizer_;
  static TRasterizer triangleRasterizers_[32];
  virtual void rasterizeTriangles ();

  void drawTriNULL (unsigned, unsigned, unsigned);
  void drawTriC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriTC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriPT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriPTC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriZC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriZT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriZTC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriZPT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriZPTC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriATC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAPT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAPTC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAZC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAZT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAZTC (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAZPT (unsigned uiV0, unsigned uiV1, unsigned uiV2);
  void drawTriAZPTC (unsigned uiV0, unsigned uiV1, unsigned uiV2);

  GrColor<float>* frameBuffer_;
//  GrSoftColor* frameBuffer_;
//  GrSoftColor backgroundSoftColor_;

  // current texture image
  GrSoftImage image_;

  // depth buffer
  unsigned short* zBuffer_;
};

//#ifndef OUTLINE
//#include "GrSoftRenderer.in"
//#endif

#endif // GRSOFTRENDERER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:58  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.4  2001/06/19 16:30:00  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.3  2001/06/15 09:20:38  prkipfer
| minor changes to better overload interface
|
| Revision 1.2  2001/04/23 10:05:53  prkipfer
| moved raycaster draw method to parent class
|
| Revision 1.1  2001/03/23 17:21:39  prkipfer
| added software rasterizer support
|
|
+---------------------------------------------------------------------*/
