// - C++ -
/***************************************************************************
    File        : Quadrangle.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Thu Feb 7 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef QUADRANGLE_H
#define QUADRANGLE_H



// Grid
#include "Pool.h"
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "QuadrangleTopology.h"


namespace grd {


class GridLevel;

class Quadrangle : public Element {
public:
  // Constructors
  Quadrangle();

  // Destructor
  virtual ~Quadrangle() {}

  // Element refinement
  virtual void refine(GridLevel& gridLevel);
  virtual void coarsen(GridLevel& gridLevel);

  virtual void close(ConformingClosure& closure);
  virtual void close(int refPattern,ConformingClosure& closure);

  virtual bool queryChildrenRefinement();
  virtual bool queryEdgeRefinement();
  virtual bool queryChildrenEdgeRefinement();

  virtual bool queryChildrenMarks();

  virtual int  getCommonFace(Element* theElem);
  virtual void getRefElements(int face,Vertex** theVertex,
                                     Edge** theEdge,Element** theChild,
                                     int* theChildFace);

  virtual int getNoOfChildren() { return 4; }

  // Topology
  virtual int getNoOfVertices() { return 4; }
  virtual int getNoOfEdges() { return 4; }
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
  Edge*    e[4];

  // connection to volume
  Element* n;

  // closure computation
  int getRefinementPattern()
  {
    int refMask = 1;
    int refPattern = 0;
    for (int i = 0; i < 4; i++)
    {
      if (e[i]->isRefined()) {
        refPattern = (refPattern | refMask);
      }
      refMask = (refMask<<1);
    }
    return refPattern;
  }
  // closure computation
  int getRefinementPattern(Vertex** vertex)
  {
    int refMask = 1;
    int refPattern = 0;
    for (int i = 0; i < 4; i++)
    {
      vertex[i] = v[i];
      if (e[i]->isRefined()) {
        refPattern = (refPattern | refMask);
        vertex[4+i] = e[i]->getMidPoint();
      }
      refMask = (refMask<<1);
    }
    return refPattern;
  }

  // memory management
  static Pool<Quadrangle> pool;
};


} // namespace grd

#endif // QUADRANGLE_H
