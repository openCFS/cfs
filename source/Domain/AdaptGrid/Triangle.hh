// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// - C++ -
/***************************************************************************
    File        : Triangle.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef TRIANGLE_H
#define TRIANGLE_H


// Grid
#include "Pool.hh"
#include "Vertex.hh"
#include "Edge.hh"
#include "Element.hh"
#include "TriangleTopology.hh"


namespace grd {


class GridLevel;

class Triangle : public Element {
public:
  // Constructors
  Triangle();

  // Destructor
  virtual ~Triangle() {}

  // Element refinement
  virtual void refine(GridLevel& gridLevel);
  virtual void coarsen(GridLevel& gridLevel);
  virtual void close(GridLevel& gridLevel);

  // private function for closure computation
  virtual int getRefinementPattern() const;
  
  virtual bool queryEdgeRefinement()         const;
  virtual bool queryChildrenEdgeRefinement() const;

  virtual int  getCommonFace(Element* theElem) const;
  virtual void getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const;

  // Topology
  virtual int getNoOfVertices() const { return 3; }
  virtual int getNoOfEdges()    const { return 3; }
  virtual int getNoOfFaces()    const { return 1; }

  virtual void setVertex(int i,Vertex* vt);
  virtual void setEdge(int i,Edge* edg);
  virtual void setNeighbor(int face,Element* theNeighbor);

  virtual Vertex*  getVertex(int vertex) const;
  virtual Edge*    getEdge(int edge) const;
  virtual Element* getNeighbor(int face) const;

  virtual void setVertices(Vertex* const* vt);
  virtual void getVertices(Vertex** vt) const;

  virtual void setEdges(Edge** edg);
  virtual void getEdges(Edge** edg) const;

  // Topology
  virtual int     type() const;
  virtual int     noOfVerticesAtFace(int fc) const;
  virtual int     topologicalVertexAtFace(int fc,int vt) const;
  virtual Vertex* getVertexAtFace(int fc,int vt) const;
  virtual int     noOfEdgesAtFace(int fc) const;
  virtual Vertex* getVertexAtEdge(int ed,int vt) const;
  virtual Edge*   getEdgeAtFace(int fc,int ed) const;

  // overload new and delete operators
  //void * operator new(size_t) { return pool.alloc(); }
  //void   operator delete(void *p) { pool.free(p); }

private:
  // geometry and topology
  Vertex*  v[3];
  Edge*    e[3];

  // connection to volume
  Element* n;
  
  // connect
  void connect(GridLevel& gridLevel);

  // memory management
  //static Pool<Triangle> pool;
};


} // namespace grd

#endif // TRIANGLE_H
