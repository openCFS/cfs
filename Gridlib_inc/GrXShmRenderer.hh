/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRXSHMRENDERER_HH
#define GRXSHMRENDERER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#undef USE_SHM
#ifdef USE_SHM
#include <sys/ipc.h>
#include <sys/shm.h>
#include <X11/extensions/XShm.h>
#endif

#include "GbTypes.hh"
#include "GrSoftRenderer.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrXShmRenderer 
  : public GrSoftRenderer
{
public:
  GrXShmRenderer (int iXPos, int iYPos, int iWidth, int iHeight);
  virtual ~GrXShmRenderer ();

  virtual void displayBackBuffer ();

protected:
  GrXShmRenderer () {}

  Display *display_;
  XVisualInfo vi_;
  int screen_,depth_;
  Window window_;
  XImage *xImage_;
  GC gc_; // ready for use in copyarea/putimage
#ifdef USE_SHM
  XShmSegmentInfo SHMInfo_;
#endif
  unsigned char* buffer_;
  GbBool doublebuffer_;
  
};

//#ifndef OUTLINE
//#include "GrXShmRenderer.in"
//#endif

#endif // GRXSHMRENDERER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/04/23 10:05:53  prkipfer
| moved raycaster draw method to parent class
|
| Revision 1.2  2001/04/17 12:34:31  prkipfer
| added redraw method for RayCaster
|
| Revision 1.1  2001/03/23 17:17:52  prkipfer
| added software rasterizer support
|
|
+---------------------------------------------------------------------*/
