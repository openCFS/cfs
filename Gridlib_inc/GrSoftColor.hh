/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRSOFTCOLOR_HH
#define GRSOFTCOLOR_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <typeinfo>
#include <iostream>

#include "GbTypes.hh"
#include "GrColor.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrSoftColor
{
public:
  INLINE GrSoftColor (unsigned short usColor = 0);
  GrSoftColor (float fR, float fG, float fB);
  GrSoftColor (const GrColor<float>& rkColor);

  // (iR,iG,iB) must be 8-bit-per-channel color
  INLINE GrSoftColor (int iR, int iG, int iB);

  // returned values are 8-bit
  INLINE int getR () const;
  INLINE int getG () const;
  INLINE int getB () const;

  INLINE operator unsigned short ();
  operator GrColor<float> ();

  GrSoftColor& operator= (const GrColor<float>& rkColor);
  GrSoftColor& operator= (const GrSoftColor& rkColor);

  static const GrSoftColor WHITE;
  static const GrSoftColor BLACK;
  static const float SCALE;

protected:
  // format is 1555 (1-bit unused, 5-bits R, 5-bits G, 5-bits B)
  unsigned short color_;

  static const float inverse31_;
};

#ifndef OUTLINE
#include "GrSoftColor.in"
#endif

#endif // GRSOFTCOLOR_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/03/23 16:55:21  prkipfer
| added software rasterizer support
|
|
+---------------------------------------------------------------------*/
