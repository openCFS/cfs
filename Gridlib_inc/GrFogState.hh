/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRFOGSTATE_HH
#define GRFOGSTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrColor.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrFogState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrFogState);

public:
  GrFogState ();
  virtual Type getType () const;

  GbBool enabled ();
  float start ();
  float end ();
  float density ();
  GrColor<float> color ();

  typedef enum _densityFunction {
    DF_LINEAR,
    DF_EXP,
    DF_EXPSQR,
    DF_QUANTITY
  } DensityFunction;

  DensityFunction& dFunction ();

  typedef enum _applyFunction {
    AF_PER_VERTEX,
    AF_PER_PIXEL,
    AF_QUANTITY
  } ApplyFunction;

  ApplyFunction& aFunction ();

protected:
  GbBool enabled_;
  float start_, end_, density_;
  GrColor<float> color_;
  DensityFunction dFunction_;
  ApplyFunction aFunction_;
};

#ifndef OUTLINE
#include "GrFogState.in"
#endif

#endif // GRFOGSTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/02/13 14:06:41  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
