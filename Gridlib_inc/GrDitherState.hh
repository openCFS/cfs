/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRDITHERSTATE_HH
#define GRDITHERSTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrDitherState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrDitherState);

public:
  GrDitherState ();
  virtual Type getType () const;

  GbBool enabled ();

protected:
  GbBool enabled_;
};

#ifndef OUTLINE
#include "GrDitherState.in"
#endif

#endif // GRDITHERSTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/02/13 14:06:40  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
