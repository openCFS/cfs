/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRGEOPRIMITIVEMESH_HH
#define GRGEOPRIMITIVEMESH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbVec2.hh"
#include "GrColor.hh"
#include "GrRenderer.hh"
#include "GrBound.hh"
#include "GrSpatial.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrGeoPrimitiveMesh
  : public GrSpatial
{
public:
  // Construction and destruction.  GrGeoPrimitiveMesh accepts responsibility for
  // deleting the input arrays.
  GrGeoPrimitiveMesh (int uiVertexQuantity, GbVec3<float>* akVertex,
		      GbVec3<float>* akNormal, GrColor<float>* akColor, 
		      GbVec2<float>* akTexture);

  GrGeoPrimitiveMesh ();

  virtual ~GrGeoPrimitiveMesh ();

  // vertex access
  INLINE int getVerticesUsed () const;
  INLINE GbVec3<float>& vertex (int uiV);
  INLINE GbVec3<float>& normal (int uiV);
  INLINE GrColor<float>& color (int uiV);
  INLINE GbVec2<float>& texture (int uiV);
  INLINE std::vector<GbVec3<float> >& getVertices ();
  INLINE std::vector<GbVec3<float> >& getNormals ();
  INLINE std::vector<GrColor<float> >& getColors ();
  INLINE std::vector<GbVec2<float> >& getTextures ();
  INLINE GbBool hasNormals () const;
  INLINE GbBool hasColors () const;
  INLINE GbBool hasTextures () const;

  // flag for invalid vertices
  enum { INVALID_INDEX = -1 };

  class EdgeRecord {
  public:
    INLINE EdgeRecord ();

    INLINE int& vertex (int i);
    INLINE int& triangle (int i);
    INLINE int getVertexIndex (int uiV);
    INLINE int getTriangleIndex (int uiT);

  protected:
    int vertex_[2];
    int triangle_[2];
  };

  class TriangleRecord {
  public:
    INLINE TriangleRecord ();

    INLINE int& vertex (int i);
    INLINE int& edge (int i);
    INLINE int getVertexIndex (int uiV);
    INLINE int getEdgeIndex (int uiE);

  protected:
    int vertex_[3];
    int edge_[3];
  };

  // edge access
  INLINE int getEdgesUsed () const;
  INLINE EdgeRecord& edge (int uiE);
  INLINE std::vector<EdgeRecord>& getEdges();

  // triangle access
  INLINE int getTrianglesUsed () const;
  INLINE TriangleRecord& triangle (int uiT);
  INLINE std::vector<TriangleRecord>& getTriangles ();
  void getTriangle (int uiT, int& uiV0, int& uiV1, int& uiV2);

  // facet planes
  INLINE const GbPlane<float>& getFacetPlane (int uiT) const;
  void updateFacetPlanes () { /* TO DO */ }

  // mesh construction
  void insertTriangle (int uiV0, int uiV1, int uiV2);
  void removeTriangle (int uiT);

  // bounds
  virtual void updateModelBound ();

protected:
  // mesh construction
  int addEdge (int uiV0, int uiV1, int uiT);

  // drawing
  virtual void draw (GrRenderer& rkRenderer);

  // vertices
  std::vector<GbVec3<float> > vertex_;
  std::vector<GbVec3<float> > normal_;
  std::vector<GrColor<float> > color_;
  std::vector<GbVec2<float> > texture_;
  GbBool hasNormals_, hasColors_, hasTextures_;

  // edges
  std::vector<EdgeRecord> edge_;

  // triangles
  std::vector<TriangleRecord> triangle_;
  std::vector<GbPlane<float> > facetPlane_;

  // local stuff
  GrBound<float> bound_;
};

#ifndef OUTLINE
#include "GrGeoPrimitiveMesh.in"
#endif

#endif // GRGEOPRIMITIVEMESH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.4  2001/06/26 12:25:43  prkipfer
| minor fixes for clean Linux compile
|
| Revision 1.3  2001/06/19 16:29:58  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.2  2001/06/18 11:28:27  prkipfer
| now using STL, removed unnecessary signedness and rerooted mesh to be spatial object
|
| Revision 1.1  2001/02/13 14:06:42  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
