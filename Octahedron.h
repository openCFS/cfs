// - C++ -
/***************************************************************************
    File        : Octahedron.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef OCTAHEDRON_H
#define OCTAHEDRON_H

#include <vector>

// Grid
#include "Pool.h"
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "OctahedronTopology.h"

namespace grd {

class GridLevel;

class Octahedron : public Element {
public:
  // Constructors
  Octahedron();

  // Destructor
  virtual ~Octahedron() {}

  // Elment refinement
  virtual void refine(GridLevel& gridLevel);
  virtual void coarsen(GridLevel& gridLevel);

  virtual void close(ConformingClosure& closure);
  virtual void close(int refP,ConformingClosure& closure);

  virtual bool queryChildrenRefinement();
  virtual bool queryEdgeRefinement();
  virtual bool queryChildrenEdgeRefinement();

  virtual bool queryChildrenMarks();

  virtual int getCommonFace(Element* theElem);
  virtual void getRefElements(int face,Vertex** theVertex,
                                     Edge** theEdge,Element** theChild,
                                     int* theChildFace);

  virtual int getNoOfChildren() { return 14; }

  // Geometry
  virtual int getNoOfVertices() { return 7; }
  virtual int getNoOfEdges() { return 12; }
  virtual int getNoOfFaces() { return 8; }

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

  // Regular splitting
  void getTetras(std::vector<Element*>& tetra);

  // overload new and delete operators
  void * operator new(size_t) { return pool.alloc(); }
  void   operator delete(void *p) { pool.free(p); }

private:
  // geometry and topology; include barycenter
  Vertex*  v[7];
  Edge*    e[12];
  Element* n[8];

  // closure computation
  inline int getRefinementPattern();

  // memory management
  static Pool<Octahedron> pool;
};


// Description
//
inline int
Octahedron::getRefinementPattern()
{
  int i;
  int refMask = 1;
  int refPattern = 0;

  for (i = 0; i < 12; i++) {
    if (e[i]->isRefined()) {
      refPattern = (refPattern | refMask);
    }

    refMask = (refMask<<1);
  }

  return refPattern;
}
} // namespace grd

#endif //
