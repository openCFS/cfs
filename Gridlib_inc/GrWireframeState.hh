/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRWIREFRAMESTATE_HH
#define GRWIREFRAMESTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrWireframeState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrWireframeState);

public:
  GrWireframeState ();
  virtual Type getType () const;

  INLINE GbBool enabled ();
  INLINE void enable(GbBool e);

protected:
  GbBool enabled_;
};

#ifndef OUTLINE
#include "GrWireframeState.in"
#endif

#endif // GRWIREFRAMESTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/02/13 14:06:50  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
