/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRRAYCASTER_HH
#define GRRAYCASTER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbImages.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


class GrRayCaster
{
public:
  GrRayCaster (GbImageUChar3D* pkImage, float fGamma, float fScale=1.0);
  ~GrRayCaster ();

  GbBool moveTrackBall (float fX0, float fY0, float fX1, float fY1, float fScale=1.0);
  void trace (int iSpacing);
  void drawWireFrame ();
  void correction (float fGamma);

  INLINE unsigned short getRendered (int uiX, int uiY);
  INLINE unsigned short getRendered (int uiIndex);
  INLINE GbImageUShort2D* getRenderResult () const;
  INLINE float& frame (int uiY, int uiX);

private:
  // storage for scaled volume data
  GbImageFloat3D* density_;

  // accumulator and render images
  int bound_, boundM1_;
  int boundH_;
  GbImageFloat2D* accu_;
  GbImageUShort2D* rendering_;

  // center point of image
  float centerX_, centerY_, centerZ_;

  // for gamma correction of rendered image values
  float gamma_;

  // frame field for eyepoint:  u = column 0, v = column 1, w = column 2
  float frameField_[3][3];
  float scaling_;

  GbBool clipped (float fP, float fQ, float& rfU1, float& rfU2);
  GbBool clip3D (float& rfX1, float& rfY1, float& rfZ1, 
		 float& rfX2, float& rfY2, float& rfZ2);
  void line3D (int iJ0, int iJ1, 
	       int iX0, int iY0, int iZ0, 
	       int iX1, int iY1, int iZ1);
  void line2D (GbBool bVisible, int iX0, int iY0, int iX1, int iY1);
};

#ifndef OUTLINE
#include "GrRayCaster.in"
#endif

#endif // GRRAYCASTER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/06/15 09:12:38  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.2  2001/04/23 10:05:17  prkipfer
| added scaling support
|
| Revision 1.1  2001/04/17 12:36:40  prkipfer
| added RayCaster
|
|
+---------------------------------------------------------------------*/
