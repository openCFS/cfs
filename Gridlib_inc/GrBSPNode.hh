/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRBSPNODE_HH
#define GRBSPNODE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrNode.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrBSPNode 
  : public GrNode
{
public:
  // Construction.  The base class GrNode has *three* children and is not
  // allowed to grow.  The first and last children (indices 0 and 2) are
  // the left and right children of the binary tree.  The middle child slot
  // is where a derived class can attach additional geometry such as the
  // triangles that are coplanar with the splitting plane.
  GrBSPNode ();

  // These methods should be used instead of the attach/detach methods in
  // the GrNode base class.  The left child corresponds to the positive
  // side of the separating plane.  The right child corresponds to the
  // negative side of the plane.
  GrSpatial* attachLeftChild (GrSpatial* pkChild);
  GrSpatial* attachRightChild (GrSpatial* pkChild);
  GrSpatial* detachLeftChild ();
  GrSpatial* detachRightChild ();
  GrSpatial* getLeftChild ();
  GrSpatial* getRightChild ();

  // plane access
  INLINE GbPlane<float>& modelPlane ();
  INLINE const GbPlane<float>& modelPlane () const;
  INLINE const GbPlane<float>& getWorldPlane () const;

  // determine BSP node whose represented region contains the point
  GrBSPNode* getContainingNode (const GbVec3<float>& rkPoint);

  // An application can attach relevant data to each BSP node.  The
  // callback is executed in the Draw pass when it is determined that the
  // node is potentially visible.  The application is responsible for
  // memory management associated with the data.
  typedef void (*CallbackFunction)(void);
  INLINE CallbackFunction& callback ();
  INLINE void*& data ();

protected:
  // geometric updates
  virtual void updateWorldData (float fAppTime);

  // drawing
  virtual void draw (GrRenderer& rkRenderer);

  GbPlane<float> modelPlane_;
  GbPlane<float> worldPlane_;
  CallbackFunction callback_;
  void* data_;
};

#ifndef OUTLINE
#include "GrBSPNode.in"
#endif

#endif // GRBSPNODE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/19 16:30:19  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.1  2001/06/18 11:23:09  prkipfer
| introduced new classes for spatial partitioning and hierarchical culling
|
|
+---------------------------------------------------------------------*/
