/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRALPHASTATE_HH
#define GRALPHASTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrRenderState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrAlphaState 
  : public GrRenderState
{
  GR_DECLARE_DEFAULT_STATE(GrAlphaState);

public:
  GrAlphaState ();
  virtual Type getType () const;

  // blending
  typedef enum _srcBlendFunction {
    SBF_ZERO,
    SBF_ONE,
    SBF_DST_COLOR,
    SBF_ONE_MINUS_DST_COLOR,
    SBF_SRC_ALPHA,
    SBF_ONE_MINUS_SRC_ALPHA,
    SBF_DST_ALPHA,
    SBF_ONE_MINUS_DST_ALPHA,
    SBF_SRC_ALPHA_SATURATE,
    SBF_QUANTITY
  } SrcBlendFunction;

  typedef enum _dstBlendFunction {
    DBF_ZERO,
    DBF_ONE,
    DBF_SRC_COLOR,
    DBF_ONE_MINUS_SRC_COLOR,
    DBF_SRC_ALPHA,
    DBF_ONE_MINUS_SRC_ALPHA,
    DBF_DST_ALPHA,
    DBF_ONE_MINUS_DST_ALPHA,
    DBF_QUANTITY
  } DstBlendFunction;

  INLINE GbBool blendEnabled ();
  INLINE void blendEnable(GbBool b);
  SrcBlendFunction& srcBlend ();
  DstBlendFunction& dstBlend ();

  // testing
  typedef enum _testFunction {
    TF_NEVER,
    TF_LESS,
    TF_EQUAL,
    TF_LEQUAL,
    TF_GREATER,
    TF_NOTEQUAL,
    TF_GEQUAL,
    TF_ALWAYS,
    TF_QUANTITY
  } TestFunction;

  INLINE GbBool testEnabled ();
  INLINE void testEnable(GbBool b);
  TestFunction& test ();
  float reference ();

protected:
  GbBool blendEnabled_;
  SrcBlendFunction srcBlend_;
  DstBlendFunction dstBlend_;

  GbBool testEnabled_;
  TestFunction test_;
  float reference_;  // in [0,1]
};

#ifndef OUTLINE
#include "GrAlphaState.in"
#endif

#endif // GRALPHASTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/02/13 14:06:39  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
