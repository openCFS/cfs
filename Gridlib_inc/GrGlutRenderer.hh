/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRGLUTRENDERER_HH
#define GRGLUTRENDERER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <GL/glut.h>

#include "GrOpenGLRenderer.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrGlutRenderer 
  : public GrOpenGLRenderer
{
public:
  // construction and destruction
  GrGlutRenderer (int iWindowID, int iWidth, int iHeight);

  void activate ();
  int getWindowID () const;

  virtual void displayBackBuffer ();
  virtual void draw (int iX, int iY, const char* acText);

protected:
  // GLUT identifier for the window of the renderer
  int windowID_;
};

//#ifndef OUTLINE
//#include "GrGlutRenderer.in"
//#endif

#endif // GRGLUTRENDERER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/02/13 14:06:44  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
