/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRPORTAL_HH
#define GRPORTAL_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbMatrix3.hh"
#include "GrConvexRegion.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*! 
  The portal is a unidirectional connector between two regions.  The
  vertices of the portal must satisfy the following constraints:
  1. They must form a planar convex polygon (quantity >= 3 is implied).
  2. They must be counterclockwise ordered when looking through the
     portal to the adjacent region.
  3. They must be in the model space coordinates for the region that
     contains the portal.
*/

class GrPortal
{
public:
  // Construction and destruction.  GrPortal accepts responsibility for
  // deleting the input array.
  GrPortal (int iVertexQuantity, GbVec3<float>* akModelVertex,
	    GrConvexRegion* pkAdjacentRegion, GbBool bOpen);

  ~GrPortal ();

  // member access
  INLINE GrConvexRegion*& adjacentRegion();
  INLINE GbBool& open();

protected:
  // streaming
  GrPortal();

  // updates
  friend class GrConvexRegion;
  void updateWorldData (const GbMatrix3<float>& rkRot, const GbVec3<float>& rkTrn,
			float fScale);

  // drawing
  void draw (GrRenderer& rkRenderer);

  // portal vertices
  int nbrVertices_;
  GbVec3<float>* modelVertex_;
  GbVec3<float>* worldVertex_;

  // region to which portal leads
  GrConvexRegion* adjacentRegion_;

  // portals can be open or closed
  GbBool isOpen_;
};

#ifndef OUTLINE
#include "GrPortal.in"
#endif

#endif // GRPORTAL_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/06/18 11:23:15  prkipfer
| introduced new classes for spatial partitioning and hierarchical culling
|
|
+---------------------------------------------------------------------*/
