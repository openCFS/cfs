/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRZBUFFERSTATE_HH
#define GRZBUFFERSTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrZBufferState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrZBufferState);

public:
  GrZBufferState ();
  virtual Type getType () const;

  INLINE GbBool enabled ();
  INLINE void enable(GbBool e);
  INLINE GbBool writeable ();
  INLINE void write(GbBool e);

  typedef enum _compareFunction {
    CF_NEVER,
    CF_LESS,
    CF_EQUAL,
    CF_LEQUAL,
    CF_GREATER,
    CF_NOTEQUAL,
    CF_GEQUAL,
    CF_ALWAYS,
    CF_QUANTITY
  } CompareFunction;

  INLINE CompareFunction& compare ();
  INLINE void setCompare(CompareFunction c);

protected:
  GbBool enabled_, writeable_;
  CompareFunction compare_;
};

#ifndef OUTLINE
#include "GrZBufferState.in"
#endif

#endif // GRZBUFFERSTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/26 12:25:45  prkipfer
| minor fixes for clean Linux compile
|
| Revision 1.1  2001/02/13 14:06:51  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
