/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GROPENGLCAMERA_HH
#define GROPENGLCAMERA_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <GL/gl.h>
#include <GL/glu.h>

#include "GbTypes.hh"
#include "GrCamera.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrOpenGLCamera 
  : public GrCamera
{
public:
  GrOpenGLCamera (float fWidth, float fHeight);

  virtual void onResize(int iWidth, int iHeight);
  virtual void onFrustumChange ();
  virtual void onViewPortChange ();
  virtual void onFrameChange ();

protected:
  GrOpenGLCamera ();

  float width_, height_;

#ifdef DEBUG
  static char* ms_acMessage;
  static GbBool OperationOkay ();
#endif
};

//#ifndef OUTLINE
//#include "GrOpenGLCamera.in"
//#endif

#endif // GROPENGLCAMERA_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/15 09:04:43  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:34  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
