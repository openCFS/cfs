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
#include "Pool.h"
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "TriangleTopology.h"
#include "TriSkeleton.h"


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

  virtual void close(ConformingClosure& closure);
  virtual void close(int refP,ConformingClosure& closure);

  virtual bool queryChildrenRefinement();
  virtual bool queryEdgeRefinement();
  virtual bool queryChildrenEdgeRefinement();

  virtual bool queryChildrenMarks();

  virtual int  getCommonFace(Element* theElem);
  virtual void getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace);

  virtual int getNoOfChildren() { return 4; }

  // Topology
  virtual int getNoOfVertices() { return 3; }
  virtual int getNoOfEdges() { return 3; }
  virtual int getNoOfFaces() { return 3; }

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

  // Surface skeleton
  void setSkeleton(TriSkeleton* skl) { skeleton = skl; }
  TriSkeleton* getSkeleton() { return skeleton; }

  // overload new and delete operators
  void * operator new(size_t) { return pool.alloc(); }
  void   operator delete(void *p) { pool.free(p); }

private:
  // geometry and topology
  Vertex*  v[3];
  Edge*    e[3];

  // connection to volume
  Element* n;

  // private function for closure computation
  inline int getRefinementPattern()
  {
    int refMask = 1;
    int refPattern = 0;
    for (int i = 0; i < 3; i++)
    {
      if (e[i]->isRefined())
      { refPattern = (refPattern | refMask); }
      refMask = (refMask<<1);
    }
    return refPattern;
  }

  // Adapt to boundary surface
  TriSkeleton *skeleton;
  void setVertexPosition();

  // memory management
  static Pool<Triangle> pool;
};


} // namespace grd

#endif // TRIANGLE_H
