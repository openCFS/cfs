/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOVOLUMEMESH_HH
#define GOVOLUMEMESH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GoMesh.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

//! class GoVolumeMesh
/*!
  Implementation of a hierarchical hybrid 3d mesh consisting of
  tetrahedra, hexahedra, octahedra, prisms and pyramids.
 */
class GoVolumeMesh 
  : public GoMesh
{
public:
  // Constructor
  GoVolumeMesh();
  virtual ~GoVolumeMesh();

  // Query topology
  virtual GoVertex<float> *         nextVertex(GoVertex<float> *v);
  virtual GoGeometryElement<float> *nextElement(GoVertex<float> *v);
  virtual GoEdge<float> *           nextEdge(GoVertex<float> *v);

  // Actions on the mesh
  virtual void switchNormals();
  virtual void updateNeighbourPositionsOnLevel(int l);

  // collapse operations
  virtual std::vector<GoGeometryElement<float> *> *halfEdgeCollapse(GoGeometryElement<float> *, int, int, int&, std::vector<GoVertex<float> *> *);
  virtual bool faceCollapse(GoGeometryElement<float> *elem, int);
  virtual bool elementCollapse(GoGeometryElement<float> *);

  // split operations
  virtual bool vertexSplit(GoVertex<float> *start, std::vector<GoVertex<float> *> *allNeighbours, float x, float y, float z );

  // This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoVolumeMesh&);

protected:
  virtual void updateFaceCache(int level=0);
};


//  #ifndef OUTLINE
//  #include "GoVolumeMesh.in"
//  #endif

#endif // GOVOLUMEMESH_HH

