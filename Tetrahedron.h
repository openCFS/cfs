// - C++ -
/***************************************************************************
    File        : Tetrahedron.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H


// Grid
#include "Pool.h"
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Jacobian.h"
// include topology information
#include "TetrahedronTopology.h"

namespace grd {

class GridLevel;

class Tetrahedron : public Element {
public:
  // Constructors
  Tetrahedron();

  // Destructor
  virtual ~Tetrahedron() {}

  // Element refinement
  virtual void refine(GridLevel& gridLevel);
  virtual void coarsen(GridLevel& gridLevel);

  virtual void close(ConformingClosure& closure);
  virtual void close(int refP,ConformingClosure& closure);

  virtual bool queryChildrenRefinement();
  virtual bool queryEdgeRefinement();
  virtual bool queryChildrenEdgeRefinement();

  virtual bool queryChildrenMarks();

  virtual int  getCommonFace(Element* theElem);
  virtual void getRefElements(int face,Vertex** theVertex,Edge** theEdge,
                                     Element** theChild,int* theChildFace);

  virtual int getNoOfChildren() { return 5; }

  // Topology
  virtual int getNoOfVertices() { return 4; }
  virtual int getNoOfEdges() { return 6; }
  virtual int getNoOfFaces() { return 4; }

  virtual void setVertex(int i,Vertex* vt);
  virtual void setEdge(int i,Edge* edg);
  virtual void setNeighbor(int face,Element* theNeighbor);

  virtual Vertex*  getVertex(int vertex);
  virtual Edge*    getEdge(int edge);
  virtual Element* getNeighbor(int face);

  virtual void setVertices(Vertex** vt);
  virtual void getVertices(Vertex** vt);

  virtual void setEdges(Edge** edg);
  virtual void getEdges(Edge** edg);

  // Topology
  virtual int type();
  virtual int noOfVertexAtFace(int fc);
  virtual Vertex* getVertexAtFace(int fc,int vt);
  virtual Vertex* getVertexAtEdge(int ed,int vt);
  virtual Edge*   getEdgeAtFace(int fc,int ed);

//  // Geometry
//  virtual void getGlobalCoords(const double* lc,double* coord);
//  virtual void getJacobian(Jacobian& Jac);

  // overload new and delete operators
  void * operator new(size_t) { return pool.alloc(); }
  void   operator delete(void *p) { pool.free(p); }

private:
  // geometry and topology
  Vertex*  v[4];
  Edge*    e[6];
  Element* n[4];

  // private function for closure computation
  int getRefinementPattern()
  {
    int refMask = 1;
    int refPattern = 0;
    for (int i = 0; i < 6; i++) {
      if (e[i]->isRefined()) {
        refPattern = (refPattern | refMask);
      }
      refMask = (refMask<<1);
    }
    return refPattern;
  }

  // memory management
  static Pool<Tetrahedron> pool;
};


} // namespace grd

#endif // TETRAHEDRON_H
