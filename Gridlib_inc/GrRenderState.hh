/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRRENDERSTATE_HH
#define GRRENDERSTATE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrRenderState
{
public:
  // supported render states
  typedef enum _type {
    RS_ALPHA,
    RS_DITHER,
    RS_FOG,
    RS_LIGHT,
    RS_MATERIAL,
    RS_SHADE,
    RS_TEXTURE,
    RS_VERTEXCOLOR,
    RS_WIREFRAME,
    RS_ZBUFFER,
    RS_MAX_STATE
  } Type;

  enum { RS_QUANTITY = RS_MAX_STATE*sizeof(GrRenderState*) };

  virtual Type getType () const = 0;
  virtual void push (GrRenderState*& rspkState, GrRenderState*& rpkDefaultState);
  virtual void pop (GrRenderState*& rspkState, GrRenderState*& rpkDefaultState);
  virtual void copyTo (GrRenderState*& rspkState);

  class List 
  {
  public:
    GrRenderState* state_;
    List* next_;
  };

  static GrRenderState** getDefaultStates () { return defaultStates_; }
  static void push (GrRenderState*& rspkState);
  static void pop (GrRenderState*& rspkState);
  static void copyTo (GrRenderState* aspkState[]);

protected:
  // abstract base class
  GrRenderState ();
  virtual ~GrRenderState ();

  // default states
  static GrRenderState* defaultStates_[RS_MAX_STATE];
};

//#ifndef OUTLINE
//#include "GrRenderState.in"
//#endif

#endif // GRRENDERSTATE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/18 11:31:50  prkipfer
| introduced render state stack
|
| Revision 1.1  2001/02/13 14:06:46  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
