/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRCONVEXREGIONMANAGER_HH
#define GRCONVEXREGIONMANAGER_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrBSPNode.hh"

class GrCamera;
class GrConvexRegion;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


class GrConvexRegionManager 
  : public GrBSPNode
{
public:
  // Construction.  The BSP tree should be built so that the leaf nodes are
  // where the GrConvexRegion objects are located.
  GrConvexRegionManager (GbBool bUseEyePlusNear = false);

  // The middle child of GrConvexRegionManager is where the representation
  // of the outside of the set of regions is stored.  This can be an
  // arbitrary subgraph, not just drawable geometry.
  GrSpatial* attachOutside (GrSpatial* pkOutside);
  GrSpatial* detachOutside ();
  GrSpatial* getOutside ();

  // Determine region that contains the point.  If the point is outside
  // the set of regions, the return values is null.
  GrConvexRegion* getContainingRegion (const GbVec3<float>& rkPoint);

  // getContainingRegion is called in Draw to determine in which
  // room the eye point is.  When transitions are allowed between inside
  // and outside, there is the problem when the eye point is inside
  // (outside) and the near face of the view frustum is outside (inside).
  // In this case the wrong contents are drawn.  To avoid this, set
  // the eye_plus_near Boolean value to 'true', in which case the point
  // E+n*D is used instead of E (E = eye point, n = near distance, D =
  // camera direction).
  GbBool& useEyePlusNear ();

protected:
  // drawing
  virtual void draw (GrRenderer& rkRenderer);

  GbBool useEyePlusNear_;
};

#ifndef OUTLINE
#include "GrConvexRegionManager.in"
#endif

#endif // GRCONVEXREGIONMANAGER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/06/18 11:23:12  prkipfer
| introduced new classes for spatial partitioning and hierarchical culling
|
|
+---------------------------------------------------------------------*/
