/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRLIGHTSTATE_HH
#define GRLIGHTSTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GrRenderState.hh"
#include "GrLight.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrLightState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrLightState);

public:
  GrLightState ();
  virtual ~GrLightState ();

  virtual Type getType () const;

  enum { MAX_LIGHTS = 8 };

  int attach (GrLight<float>* pkLight);
  int detach (GrLight<float>* pkLight);
  GrLight<float>* detach (int uiIndex);
  void detachAll ();

  int getNumLights() const;
  GrLight<float>* getLight (int uiIndex);

protected:
  GrLight<float>* lights_[MAX_LIGHTS];
};

//#ifndef OUTLINE
//#include "GrLightState.in"
//#endif

#endif // GRLIGHTSTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 09:10:24  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 14:06:45  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
