/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRBSPTREE_HH
#define GRBSPTREE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbVec2.hh"
#include "GrColor.hh"
#include "GrBSPNode.hh"
#include <list>

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


// Creation of a BSP tree from a list of triangles.  The internal nodes of
// the tree are GrBSPNode objects and store the coplanar triangles (if any)
// of the splitting planes.


class GrBSPTree 
  : public GrBSPNode
{
public:
  class Triangle {
  public:
    Triangle ();

    GbVec3<float> vertex_[3];
    GbVec3<float>* normal_[3];
    GrColor<float>* color_[3];
    GbVec2<float>* texture_[3];
  };

  typedef std::list<Triangle*> TriangleList;

  // Construction.  The input list is consumed by the construction, so the
  // application should keep a copy of the list if it needs to be used
  // elsewhere.
  GrBSPTree (TriangleList& rkList);

protected:
  GrBSPTree ();

  void createTree (TriangleList& rkList);

  void splitTriangle (Triangle* pkTri, TriangleList& rkPositive,
		      TriangleList& rkNegative, TriangleList& rkCoincident);

  void clipTriangle (int i0, int i1, int i2, Triangle* pkTri,
		     float afDistance[3], TriangleList& rkPositive,
		     TriangleList& rkNegative);

  void addTriangle (TriangleList& rkList, const GbVec3<float>& rkV0,
		    const GbVec3<float>& rkV1, const GbVec3<float>& rkV2, GbBool bHasNormals,
		    const GbVec3<float>* pkN0, const GbVec3<float>* pkN1, const GbVec3<float>* pkN2,
		    GbBool bHasColors, const GrColor<float>* pkC0, const GrColor<float>* pkC1,
		    const GrColor<float>* pkC2, GbBool bHasTextures, const GbVec2<float>* pkT0,
		    const GbVec2<float>* pkT1, const GbVec2<float>* pkT2);
};

// #ifndef OUTLINE
// #include "GrBSPTree.in"
// #endif

#endif // GRBSPTREE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/06/26 12:25:25  prkipfer
| minor fixes for clean Linux compile
|
| Revision 1.2  2001/06/19 16:30:19  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.1  2001/06/18 11:23:10  prkipfer
| introduced new classes for spatial partitioning and hierarchical culling
|
|
+---------------------------------------------------------------------*/
