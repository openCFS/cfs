/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOTETRAHEDRONMESH_HH
#define GOTETRAHEDRONMESH_HH

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
#include "GoTetrahedronElement.hh"
#include "GoMesh.hh"
#include "GbMeshFunctions.hh"
#include "GbEdgeCollapse.hh"
/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

//! class GoTetrahedronMesh
/*!

 */
class GoTetrahedronMesh 
  : public GoMesh
{
public:
  //! Constructor
  GoTetrahedronMesh();
  virtual ~GoTetrahedronMesh();

  // Query topology
  virtual GoVertex<float> *nextVertex(GoVertex<float> *v);
  virtual GoGeometryElement<float>   *nextElement(GoVertex<float> *v);

  // Actions on the mesh
  virtual void setupNeighbours();
  virtual void switchNormals();
  virtual int consistency();

  // Gradient
  float computeGradient(GoVertex<float> *v);
  float computeGradient(GoGeometryElement<float> *e);

  // collapse operations
  virtual bool mayCollapseEdge(GoGeometryElement<float> *, int, int);
  virtual std::vector<GoGeometryElement<float> *> *halfEdgeCollapse(GoGeometryElement<float> *, int, int, int &, std::vector<GoVertex<float> *> *);
  virtual void halfEdgeCollapse(int startId, int endId);

  // vertex split operations
  virtual int vertexSplit( GbHalfEdgeCollapse &); //int id, const std::vector<int> &n, const GbVec3<float> &pos, int bound, int);


  virtual std::vector<GoGeometryElement<float> *> * getAllNeighbourElements(GoVertex<float>* vertex);

  virtual std::vector<GoVertex<float> *> * getAllNeighbourVertices(GoVertex<float>* vertex);

  virtual std::vector<GoGeometryElement<float> *> * getBoundaryNeighbourElements(GoVertex<float>* vertex);

  virtual std::vector<GoVertex<float> *> * getBoundaryNeighbourVertices(GoVertex<float>* vertex);
  virtual void setUndeletedElement(GoVertex<float>* vertex);


  //! This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoTetrahedronMesh&);

private:
  struct TableEntry {
    int i,j;
    GoGeometryElement<float> *f;
    TableEntry *next;
  };
  
  void searchTable(int, int, int, GoGeometryElement<float> *, TableEntry **);
};


//  #ifndef OUTLINE
//  #include "GoTetrahedronMesh.in"
//  #endif

#endif // GOTETRAHEDRONMESH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.9  2001/05/03 08:21:42  uflabsik
| *** empty log message ***
|
| Revision 1.8  2001/05/03 08:05:51  uflabsik
| functionality added
|
| Revision 1.7  2000/12/12 10:03:50  prkipfer
| moved consistency check to GoMesh
|
| Revision 1.6  2000/11/07 17:27:55  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.5  2000/10/31 16:29:07  uflabsik
| new version
|
| Revision 1.4  2000/10/23 09:35:20  prkipfer
| moved bounding box calculation to GoMesh
|
| Revision 1.3  2000/07/20 13:32:41  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.2  2000/07/11 07:55:38  uflabsik
| setupNeighbours implemented
|
| Revision 1.1  2000/07/03 16:15:13  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
