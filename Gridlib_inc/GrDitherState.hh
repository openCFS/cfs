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
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/02/13 14:06:40  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
