/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRGEOPRIMITIVERENDERERMESH_HH
#define GRGEOPRIMITIVERENDERERMESH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbTStorage.hh"
#include "GbTClassStorage.hh"
#include "GrGeoPrimitiveMesh.hh"
#include "GrVertexColorState.hh"
#include "GrLightState.hh"
#include "GrMaterialState.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrGeoPrimitiveRendererMesh
  : public GrGeoPrimitiveMesh
{
public:
  // construction and destruction
  GrGeoPrimitiveRendererMesh ();
  virtual ~GrGeoPrimitiveRendererMesh ();

  // model vertices and lighting
  INLINE GbVec3<float>& modelVertex (int uiV);
  INLINE std::vector<GbVec3<float> >& getModelVertices ();
  INLINE GbBool hasModelVertices () const;
  INLINE GbBool hasLights () const;
  void getLightingStatus (GrLightState* pkLState,
			  GrMaterialState* pkMState, 
			  GrVertexColorState* pkVState,
			  GrColor<float>* akVertexColors);

  // visibility
  INLINE GbBool getVertexVisible (int uiV) const;
  INLINE GbBool getEdgeVisible (int uiE) const;
  INLINE GbBool getTriangleVisible (int uiT) const;
  INLINE float& depth (int uiV);

protected:
  friend class GrGeoPrimitiveRenderer;

  // renderer utilities
  void copy (GrGeoPrimitiveMesh& rkMesh);

  // extension of mesh from clipped data
  int addVertex ();
  int addEdge ();
  int addTriangle ();

  // mesh growth
  virtual void growVertices (int uiNewQuantity);
  virtual void growEdges (int uiNewQuantity);
  virtual void growTriangles (int uiNewQuantity);

  // vertices
  std::vector<GbBool> vertexVisible_;
  std::vector<GbBool> vertexNeedsLighting_;
  std::vector<float> distance_;
  std::vector<int> oldEdge_;
  std::vector<int> newEdge_;
  std::vector<float> depth_;

  // edges
  std::vector<GbBool> edgeVisible_;
  std::vector<int> clipVertex_;
  std::vector<float> clipParam_;

  // triangles
  std::vector<GbBool> triangleVisible_;

  // interpolation and lighting of attributes
  void interpolate (int uiV0, int uiV1, float fParam, int uiClipV);
  void interpolateVertexAttributes (int uiV);
  void markVertexForLighting (int uiOldVerticesUsed, int uiV);
  void lightVertices (int uiOldVerticesUsed);
  void computeVertexAttributes (int uiOldVerticesUsed);

  // lighters
  typedef void (GrGeoPrimitiveRendererMesh::*LightFunction)(int);
  static LightFunction lightFx_[GrVertexColorState::SM_QUANTITY][GrVertexColorState::LM_QUANTITY];

  void sourceIgnoreLightingEmissive (int uiQuantity);
  void sourceIgnoreLightingDiffuse (int uiQuantity);
  void sourceEmissiveLightingEmissive (int uiQuantity);
  void sourceEmissiveLightingDiffuse (int uiQuantity);
  void sourceDiffuseLightingEmissive (int uiQuantity);
  void sourceDiffuseLightingDiffuse (int uiQuantity);

  void computeDiffuseSpecular (int uiQuantity, GrColor<float>& rkAmbient,
			       GbBool& rbHasSpecular);

  // current lighting state and members for lighting calculations
  GrLightState* lightState_;
  GrMaterialState* materialState_;
  GrVertexColorState* vertexColorState_;
  std::vector<GbVec3<float> > modelVertex_;
  std::vector<GrColor<float> > diffuse_;
  std::vector<GrColor<float> > specular_;
  GbVec3<float> cameraModelLocation_;
  GbBool hasModelVertices_;
  GbBool hasLights_;
};

#ifndef OUTLINE
#include "GrGeoPrimitiveRendererMesh.in"
#endif

#endif // GRGEOPRIMITIVERENDERERMESH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.4  2001/06/19 16:29:58  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.3  2001/06/18 11:31:11  prkipfer
| adapted to GrGeoPrimitiveRenderer, using STL and culling support
|
| Revision 1.2  2001/06/15 09:08:43  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 14:06:43  prkipfer
| introduced rendering subsystem
|
|
+---------------------------------------------------------------------*/
