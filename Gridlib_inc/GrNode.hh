/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRNODE_HH
#define GRNODE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrSpatial.hh"
#include "GrRenderer.hh"
#include "GrBound.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrNode 
  : public GrSpatial
{
public:
  // construction and destruction
  GrNode (int iQuantity = 1, int iGrowBy = 1);
  virtual ~GrNode ();

  // children
  INLINE int getQuantity () const;
  INLINE int getUsed () const;
  int attachChild (GrSpatial* pkChild);
  int detachChild (GrSpatial* pkChild);
  GrSpatial* detachChildAt (int i);
  GrSpatial* setChild (int i, GrSpatial* pkChild);
  GrSpatial* getChild (int i);

  // support for searching by name
//  virtual Object* getObjectByName (const char* acName);

protected:
  // allow internal calls to update functions and drawing
  friend class GrSpatial;
  friend class GrRenderer;

  // geometric updates
  virtual void updateWorldData (float fAppTime);
  virtual void updateWorldBound ();

  // render state update
  virtual void updateRenderState ();

  // drawing
  virtual void draw (GrRenderer& rkRenderer);


  // children
  std::vector<GrSpatial*> child_;
  int growBy_, used_;
};

#ifndef OUTLINE
#include "GrNode.in"
#endif

#endif // GRNODE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/19 16:30:21  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.1  2001/06/18 11:23:14  prkipfer
| introduced new classes for spatial partitioning and hierarchical culling
|
|
+---------------------------------------------------------------------*/
