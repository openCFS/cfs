/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRIMAGERENDERER_HH
#define GRIMAGERENDERER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrSoftRenderer.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrImageRenderer
  : public GrSoftRenderer
{
public:
  GrImageRenderer (int iWidth, int iHeight);
  virtual ~GrImageRenderer ();

  virtual void displayBackBuffer ();

  virtual GbBool save(const GbString& filename);

protected:
  GrImageRenderer () {}

  unsigned char* buffer_;
  GbBool doublebuffer_;
  
};

//#ifndef OUTLINE
//#include "GrImageRenderer.in"
//#endif

#endif // GRIMAGERENDERER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2001/04/23 10:05:52  prkipfer
| moved raycaster draw method to parent class
|
| Revision 1.2  2001/04/17 12:33:19  prkipfer
| added save method
|
| Revision 1.1  2001/03/23 17:14:26  prkipfer
| added software rasterizer support
|
|
+---------------------------------------------------------------------*/
