/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRMATERIALSTATE_HH
#define GRMATERIALSTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GrColor.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrMaterialState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrMaterialState);

public:
  GrMaterialState ();
  virtual Type getType () const;

  GrColor<float> emissive ();
  GrColor<float> ambient ();
  GrColor<float> diffuse ();
  GrColor<float> specular ();
  float shininess ();
  float alpha ();

protected:
  GrColor<float> emissive_;
  GrColor<float> ambient_;
  GrColor<float> diffuse_;
  GrColor<float> specular_;
  float shininess_;
  float alpha_;
};

#ifndef OUTLINE
#include "GrMaterialState.in"
#endif

#endif // GRMATERIALSTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/02/13 14:06:45  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
