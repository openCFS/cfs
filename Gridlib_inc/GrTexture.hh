/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRTEXTURE_HH
#define GRTEXTURE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrColor.hh"
#include "GrImage.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrTexture 
{
public:
  GrTexture ();
  virtual ~GrTexture ();

  INLINE void setImage (GrImage* pkImage);
  INLINE GrImage* getImage ();


  typedef enum _correctionMode {
    CM_AFFINE,
    CM_PERSPECTIVE,
    CM_QUANTITY
  } CorrectionMode;

  INLINE CorrectionMode& correction ();


  typedef enum _applyMode {
    AM_REPLACE,
    AM_DECAL,
    AM_MODULATE,
    AM_BLEND,
    AM_QUANTITY
  } ApplyMode;

  INLINE ApplyMode& apply ();
  INLINE GrColor<float> blendColor ();


  typedef enum _wrapMode {
    WM_CLAMP_S_CLAMP_T,
    WM_CLAMP_S_WRAP_T,
    WM_WRAP_S_CLAMP_T,
    WM_WRAP_S_WRAP_T,
    WM_QUANTITY
  } WrapMode;

  INLINE WrapMode& wrap ();


  typedef enum _filterMode {
    FM_NEAREST,
    FM_LINEAR,
    FM_QUANTITY
  } FilterMode;

  INLINE FilterMode& filter ();


  typedef enum _mipmapMode {
    MM_NONE,
    MM_NEAREST,
    MM_LINEAR,
    MM_NEAREST_NEAREST,
    MM_NEAREST_LINEAR,
    MM_LINEAR_NEAREST,
    MM_LINEAR_LINEAR,
    MM_QUANTITY
  } MipmapMode;

  INLINE MipmapMode& mipmap ();

  INLINE float priority ();
  INLINE unsigned int userData ();
  INLINE void setUserData(unsigned int);

protected:
  GrImage* image_;
  CorrectionMode correction_;
  ApplyMode apply_;
  GrColor<float> blendColor_;
  WrapMode wrap_;
  FilterMode filter_;
  MipmapMode mipmap_;
  float priority_;
  unsigned int userData_;
};

#ifndef OUTLINE
#include "GrTexture.in"
#endif

#endif // GRTEXTURE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/19 16:30:22  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.1  2001/02/13 13:29:37  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
