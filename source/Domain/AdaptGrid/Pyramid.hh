// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// - C++ -
/***************************************************************************
    File        : Pyramid.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Nov 29 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef PYRAMID_H
#define PYRAMID_H

#include <vector>

// Grid
#include "Pool.hh"
#include "Vertex.hh"
#include "Edge.hh"
#include "Element.hh"
#include "PyramidTopology.hh"

namespace grd {

class GridLevel;

class Pyramid : public Element {
public:
  // Constructors
  Pyramid();

  // Destructor
  virtual ~Pyramid() {}

  // Elment refinement
  virtual void refine(GridLevel& gridLevel);
  virtual void coarsen(GridLevel& gridLevel);
  virtual void close(GridLevel& gridLevel);

  // closure computation
  virtual int getRefinementPattern() const;

  virtual bool queryEdgeRefinement() const;
  virtual bool queryChildrenEdgeRefinement() const;

  virtual int getCommonFace(Element* theElem) const;
  virtual void getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const;

  // Geometry
  virtual int getNoOfVertices() const { return 5;  }
  virtual int getNoOfEdges()    const { return 8; }
  virtual int getNoOfFaces()    const { return 5;  }

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
  virtual int type() const;
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
  Vertex*  v[5];
  Edge*    e[8];
  Element* n[5];

  // memory management
  //static Pool<Pyramid> pool;
};


} // namespace grd
#endif // PYRAMID_H
