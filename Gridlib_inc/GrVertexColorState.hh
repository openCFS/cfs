/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRVERTEXCOLORSTATE_HH
#define GRVERTEXCOLORSTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GrColor.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrVertexColorState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrVertexColorState);

public:
  GrVertexColorState ();
  virtual Type getType () const;

  typedef enum _sourceMode {
    SM_IGNORE,
    SM_EMISSIVE,
    SM_DIFFUSE,  // includes ambient + diffuse
    SM_QUANTITY
  } SourceMode;

  typedef enum _lightingMode {
    LM_EMISSIVE,
    LM_DIFFUSE,  // includes emissive + ambient + diffuse
    LM_QUANTITY
  } LightingMode;

  SourceMode& source ();
  LightingMode& lighting ();

protected:
  SourceMode source_;
  LightingMode lighting_;
};

#ifndef OUTLINE
#include "GrVertexColorState.in"
#endif

#endif // GRVERTEXCOLORSTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/02/13 14:06:49  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
