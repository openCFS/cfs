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

#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
// Include element topology
#include "HexahedronTopology.h"

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

  virtual void close(ConformingClosure& closure);
  virtual void close(int refP,ConformingClosure& closure);
  // Only for hexs
  void close(GridLevel& gridLevel);

  virtual bool queryChildrenRefinement();
  virtual bool queryEdgeRefinement();
  virtual bool queryChildrenEdgeRefinement();

  virtual bool queryChildrenMarks();

  virtual int getCommonFace(Element* theElem);
  virtual void getRefElements(int face,Vertex** theVertex,
                                     Edge** theEdge,Element** theChild,
                                     int* theChildFace);

  virtual int getNoOfChildren() { return 8; }

  // close hexas
  //int close(apTessellation& tessellation);
  inline Vertex* getFaceVertex(int face);
  inline Vertex* getFaceVertex(Element* neigh);

  // Geometry
  virtual int getNoOfVertices() { return 8; }
  virtual int getNoOfEdges() { return 12; }
  virtual int getNoOfFaces() { return 6; }

  virtual void setVertex(int i,Vertex* vt);
  virtual void setEdge(int i,Edge* edg);
  virtual void setNeighbor(int i,Element* elem);

  virtual Vertex*  getVertex(int i);
  virtual Edge*    getEdge(int i);
  virtual Element* getNeighbor(int i);

  virtual void setVertices(Vertex** vt);
  virtual void getVertices(Vertex** vt);

  virtual void setEdges(Edge** edg);
  virtual void getEdges(Edge** edg);

  //inline void getGlobalCoords(const double* lc,double* coord);

  Vertex* getBarycenter() { return c; }
  void resetBarycenter() { c = 0; }
  void getNoOfIrregElem(int noElem[2]);

  // Topology
  virtual int type();
  virtual int noOfVertexAtFace(int fc);
  virtual Vertex* getVertexAtFace(int fc,int vt);
  virtual Vertex* getVertexAtEdge(int ed,int vt);
  virtual Edge*   getEdgeAtFace(int fc,int ed);

  // overload new and delete
  void * operator new(size_t) { return pool.alloc(); }
  void   operator delete(void *p) { pool.free(p); }

private:
  // gerometry and topology
  Vertex*  v[8];
  Edge*    e[12];
  Element* n[6];

  // for irreg elems, c means closure
  Vertex* c;

  // closure computation
  inline int getRefinementPattern();
  inline int getRefinementPattern(Vertex** vt);

  // memory management
  static Pool<Hexahedron> pool;
};


/////////////////////////////////////////////////////////////
//
//  inline functions
//
/////////////////////////////////////////////////////////////

//
// Description
//
inline int
Hexahedron::getRefinementPattern()
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


// Description
//
inline int
Hexahedron::getRefinementPattern(Vertex** vt)
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
Hexahedron::getFaceVertex(int face)
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
Hexahedron::getFaceVertex(Element* theNeighbor)
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
