/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBIMAGES_HH
#define GBIMAGES_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbImageElement.hh"
#include "GbTImage.hh"
#include "GbTImage2D.hh"
#include "GbTImage3D.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


typedef GbTImage<GbEchar> GbImageChar;
typedef GbTImage<GbEuchar> GbImageUChar;
typedef GbTImage<GbEshort> GbImageShort;
typedef GbTImage<GbEushort> GbImageUShort;
typedef GbTImage<GbEint> GbImageInt;
typedef GbTImage<GbEuint> GbImageUInt;
typedef GbTImage<GbElong> GbImageLong;
typedef GbTImage<GbEulong> GbImageULong;
typedef GbTImage<GbEfloat> GbImageFloat;
typedef GbTImage<GbEdouble> GbImageDouble;
typedef GbTImage<GbErgb5> GbImageRGB5;
typedef GbTImage<GbErgb8> GbImageRGB8;

typedef GbTImage2D<GbEchar> GbImageChar2D;
typedef GbTImage2D<GbEuchar> GbImageUChar2D;
typedef GbTImage2D<GbEshort> GbImageShort2D;
typedef GbTImage2D<GbEushort> GbImageUShort2D;
typedef GbTImage2D<GbEint> GbImageInt2D;
typedef GbTImage2D<GbEuint> GbImageUInt2D;
typedef GbTImage2D<GbElong> GbImageLong2D;
typedef GbTImage2D<GbEulong> GbImageULong2D;
typedef GbTImage2D<GbEfloat> GbImageFloat2D;
typedef GbTImage2D<GbEdouble> GbImageDouble2D;
typedef GbTImage2D<GbErgb5> GbImageRGB52D;
typedef GbTImage2D<GbErgb8> GbImageRGB82D;

typedef GbTImage3D<GbEchar> GbImageChar3D;
typedef GbTImage3D<GbEuchar> GbImageUChar3D;
typedef GbTImage3D<GbEshort> GbImageShort3D;
typedef GbTImage3D<GbEushort> GbImageUShort3D;
typedef GbTImage3D<GbEint> GbImageInt3D;
typedef GbTImage3D<GbEuint> GbImageUInt3D;
typedef GbTImage3D<GbElong> GbImageLong3D;
typedef GbTImage3D<GbEulong> GbImageULong3D;
typedef GbTImage3D<GbEfloat> GbImageFloat3D;
typedef GbTImage3D<GbEdouble> GbImageDouble3D;
typedef GbTImage3D<GbErgb5> GbImageRGB53D;
typedef GbTImage3D<GbErgb8> GbImageRGB83D;

//#ifndef OUTLINE
//#include "GbImages.in"
//#endif


#endif // GBIMAGES_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/03/20 09:50:16  prkipfer
| introduced image handling tool
|
|
+---------------------------------------------------------------------*/
