// - C++ -
/***************************************************************************
    File        : Hexahedron.h
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Feb 22 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef HEXAHEDRON_H
#define HEXAHEDRON_H

#include "Vertex.hh"
#include "Edge.hh"
#include "Element.hh"
#include "HexahedronTopology.hh"

namespace grd {


class Hexahedron : public Element {
public:
  // Constructors
  Hexahedron();

  // Destructor
  virtual ~Hexahedron() { }

  // Element refinement
  virtual void refine(GridLevel& gridLevel);
  virtual void coarsen(GridLevel& gridLevel);
  virtual void close(GridLevel& gridLevel);

  // closure computation
  virtual int getRefinementPattern() const;

  virtual bool queryEdgeRefinement() const;
  virtual bool queryChildrenEdgeRefinement() const;

  virtual int getCommonFace(Element* theElem) const;
  virtual void getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const;

  // Closure hexas
  inline Vertex* getFaceVertex(int face)       const;
  inline Vertex* getFaceVertex(Element* neigh) const;

  // Geometry
  virtual int getNoOfVertices() const { return 8; }
  virtual int getNoOfEdges()    const { return 12; }
  virtual int getNoOfFaces()    const { return 6; }

  virtual void setVertex(int i,Vertex* vt);
  virtual void setEdge(int i,Edge* edg);
  virtual void setNeighbor(int i,Element* elem);

  virtual Vertex*  getVertex(int i) const;
  virtual Edge*    getEdge(int i) const;
  virtual Element* getNeighbor(int i) const;

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

  // overload new and delete
  //void * operator new(size_t) { return pool.alloc(); }
  //void   operator delete(void *p) { pool.free(p); }

private:
  // gerometry and topology
  Vertex*  v[8];
  Edge*    e[12];
  Element* n[6];

  // conforming closure
  inline int getRefinementPattern(Vertex** vt) const;

  // memory management
  //static Pool<Hexahedron> pool;
};


/////////////////////////////////////////////////////////////
//
//  inline functions
//
/////////////////////////////////////////////////////////////

// Description
//
inline int
Hexahedron::getRefinementPattern(Vertex** vt) const
{
  int i;
  int refMask = 1;
  int refPattern = 0;

  for (i = 0; i < 12; i++) {
    if (e[i]->isRefined()) {
      refPattern = (refPattern | refMask);
      vt[i] = e[i]->getMidPoint();
    }
    refMask = (refMask<<1);
  }

  return refPattern;
}


// Description
//
inline Vertex*
Hexahedron::getFaceVertex(int face) const
{
  int childIndex;
  Vertex* theVertex;

  childIndex = hexChildOfParentFace[face][0];
  theVertex = child[childIndex]->getVertex(hexChildVertexOfParentFace[face]);
  return theVertex;
}



// Description
//
inline Vertex*
Hexahedron::getFaceVertex(Element* theNeighbor) const
{
  int face;
  int childIndex;
  Vertex* theVertex;

  for (face = 0; face < 6; face++) {
    if (n[face] == theNeighbor) {
      childIndex = hexChildOfParentFace[face][0];
      theVertex = child[childIndex]->getVertex(hexChildVertexOfParentFace[face]);
      return theVertex;
    }
  }

  return 0;
}


} // namespace grd
#endif // HEXAHEDRON_H
