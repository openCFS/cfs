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

#include <vector>
#include <iostream>
#include <typeinfo>
#include <algorithm>

#include "GoVertex.hh"
#include "GoGeometryElement.hh"
#include "GoMesh.hh"
#include "GbMeshFunctions.hh"

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
  virtual GoVertex<float> *nextVertex(GoVertex<float> *v);
  virtual GoGeometryElement<float>   *nextElement(GoVertex<float> *v);

  // Actions on the mesh
  virtual void setupNeighbours();
  virtual void switchNormals();

  // collapse operations
  virtual std::vector<GoGeometryElement<float> *> *halfEdgeCollapse(GoGeometryElement<float> *, int, int, int&, std::vector<GoVertex<float> *> *);
  virtual bool faceCollapse(GoGeometryElement<float> *elem, int);
  virtual bool elementCollapse(GoGeometryElement<float> *);

  // split operations
  virtual bool vertexSplit(GoVertex<float> *start, std::vector<GoVertex<float> *> *allNeighbours, float x, float y, float z );

  
  virtual std::vector<GoGeometryElement<float> *> * getAllNeighbourElements(GoVertex<float>* vertex);

  // This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoVolumeMesh&);

private:
  struct TableEntry {
    int i,j,k;
    GoGeometryElement<float> *f;
    TableEntry *next;
  };
  
  void searchTable(int, int, int, GoGeometryElement<float> *, TableEntry **);
  void searchTable(int, int, int, int, GoGeometryElement<float> *, TableEntry **);
};


//  #ifndef OUTLINE
//  #include "GoVolumeMesh.in"
//  #endif

#endif // GOVOLUMEMESH_HH

