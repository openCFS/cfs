/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRTEXTURESTATE_HH
#define GRTEXTURESTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GrTexture.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrTextureState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrTextureState);

public:
  GrTextureState ();
  virtual ~GrTextureState ();

  virtual Type getType () const;

  enum { MAX_TEXTURES = 4 };

  void setTexture (int uiIndex, GrTexture* pkTexture);
  GrTexture* getTexture (int uiIndex);
  GrTexture* removeTexture (int uiIndex);
  void removeAll ();

  int getNumTextures() const;

protected:
  GrTexture* textures_[MAX_TEXTURES];
};

//#ifndef OUTLINE
//#include "GrTextureState.in"
//#endif

#endif // GRTEXTURESTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/15 09:21:24  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/02/13 14:06:48  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
