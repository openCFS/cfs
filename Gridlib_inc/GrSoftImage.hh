/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRSOFTIMAGE_HH
#define GRSOFTIMAGE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrSoftColor.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


class GrSoftImage
{
public:
  GrSoftImage ();

  void set (int iXSize, int iYSize, GrSoftColor* akColor, 
	    GbBool bXClamp, GbBool bYClamp);

  INLINE int getXSize () const;
  INLINE int getYSize () const;
  INLINE const GrSoftColor* getColors () const;

  INLINE int getX (float fX) const;
  INLINE int getY (float fY) const;

  GrSoftColor getColor (int iX, int iY);
  GrSoftColor getColor (float fX, float fY);

protected:
  int xSize_, ySize_;
  int xSizeM1_, ySizeM1_;
  int logXSize_, logYSize_;
  float floatXSizeM1_, floatYSizeM1_;
  GbBool xClamp_, yClamp_;
  GrSoftColor* color_;
};

#ifndef OUTLINE
#include "GrSoftImage.in"
#endif

#endif // GRSOFTIMAGE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:58  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/03/23 16:55:22  prkipfer
| added software rasterizer support
|
|
+---------------------------------------------------------------------*/
