/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRSHADESTATE_HH
#define GRSHADESTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrShadeState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrShadeState);

public:
  GrShadeState ();
  virtual Type getType () const;

  typedef enum _shadeMode {
    SM_FLAT,
    SM_SMOOTH,
    SM_QUANTITY
  } ShadeMode;

  INLINE ShadeMode& shade ();

protected:
  ShadeMode eShade_;
};

#ifndef OUTLINE
#include "GrShadeState.in"
#endif

#endif // GRSHADESTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/02/13 14:06:48  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
