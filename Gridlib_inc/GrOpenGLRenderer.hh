/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GROPENGLRENDERER_HH
#define GROPENGLRENDERER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <GL/gl.h>
#include <GL/glu.h>

#include "GbTypes.hh"
#include "GrRenderer.hh"
#include "GrTexture.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrOpenGLRenderer 
  : public GrRenderer
{
public:
  // construction
  GrOpenGLRenderer (int iWidth, int iHeight);

  void getDriverInformation (GbString& aucVendor,GbString& aucRenderer, 
			     GbString& aucVersion,GbString& aucGluVersion,
			     GbString& aucExtensions);

  virtual void setBackgroundColor (const GrColor<float>& rkColor);
  virtual void clearBackBuffer ();
  virtual void displayBackBuffer () = 0;
  virtual void clearZBuffer ();
  virtual void clearBuffers ();

  virtual void draw (const GoTriangleElement& tri);
  virtual void draw (const GrGeoPrimitiveMesh& rkMesh);
  virtual void draw (int iX, int iY, const char* acText) = 0;

protected:
  // state management
  virtual void initializeState ();
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

  // maps from Magic enums to OpenGL enums
  static GLenum alphaSrcBlend_[GrAlphaState::SBF_QUANTITY];
  static GLenum alphaDstBlend_[GrAlphaState::DBF_QUANTITY];
  static GLenum alphaTest_[GrAlphaState::TF_QUANTITY];
  static GLenum fogDensity_[GrFogState::DF_QUANTITY];
  static GLenum fogApply_[GrFogState::AF_QUANTITY];
  static GLenum shade_[GrShadeState::SM_QUANTITY];
  static GLenum textureCorrection_[GrTexture::CM_QUANTITY];
  static GLenum textureApply_[GrTexture::AM_QUANTITY];
  static GLenum textureFilter_[GrTexture::FM_QUANTITY];
  static GLenum textureMipmap_[GrTexture::MM_QUANTITY];
  static GLenum imageComponents_[GrImage::IT_QUANTITY];
  static GLenum imageFormats_[GrImage::IT_QUANTITY];
  static GLenum zBufferCompare_[GrZBufferState::CF_QUANTITY];

#ifdef DEBUG
  static char* ms_acMessage;
  static GbBool OperationOkay ();
#endif
};

//#ifndef OUTLINE
//#include "GrOpenGLRenderer.in"
//#endif

#endif // GROPENGLRENDERER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:11:43  prkipfer
| enabled texturing, added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 14:06:46  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
