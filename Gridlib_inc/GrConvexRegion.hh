/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRCONVEXREGION_HH
#define GRCONVEXREGION_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GrBSPNode.hh"

class GrPortal;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrConvexRegion
  : public GrBSPNode
{
public:
  // Construction and destruction.  ConvexRegion accepts responsibility
  // for deleting the input array.
  GrConvexRegion (int iPortalQuantity, GrPortal** apkPortal,
		  GrSpatial* pkRepresentation);

  virtual ~GrConvexRegion ();

  // The left and right children of GrConvexRegion must be null.  The
  // middle child is where the representation of the region is stored.
  // This can be an arbitrary subgraph, not just drawable geometry.
  GrSpatial* attachRepresentation (GrSpatial* pkRepresentation);
  GrSpatial* detachRepresentation ();
  GrSpatial* getRepresentation ();

  // portal access
  INLINE int getPortalQuantity () const;
  INLINE GrPortal* getPortal (int i) const;

protected:
  // streaming
  GrConvexRegion ();

  // geometric updates
  virtual void updateWorldData (float fAppTime);

  // Drawing.  GrConvexRegionManager starts the region graph traversal
  // with the region containing the eye point.  GrPortal continues the
  // traversal.
  friend class GrConvexRegionManager;
  friend class GrPortal;
  virtual void draw (GrRenderer& rkRenderer);

  // portals of the region (these are not set up to be shared)
  int nbrPortals_;
  GrPortal** portal_;

  // for region graph traversal
  GbBool visited_;
};

#ifndef OUTLINE
#include "GrConvexRegion.in"
#endif

#endif // GRCONVEXREGION_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/19 16:30:20  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.1  2001/06/18 11:23:11  prkipfer
| introduced new classes for spatial partitioning and hierarchical culling
|
|
+---------------------------------------------------------------------*/
