// - C++ -
/***************************************************************************
    File        : Edge.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef EDGE_H
#define EDGE_H

// Grid
#include "Pool.h"
#include "Utility.h"
#include "Vertex.h"


namespace grd {


class Edge {
public:
  // Constructors & Destructor
  Edge() {
    v[0] = v[1] = 0;
    child[0] = child[1] = 0;
    tag = false;
    iso = false;
  }
  Edge(Vertex* v0,Vertex* v1) {
    v[0] = v0;
    v[1] = v1;
    child[0] = child[1] = 0;
    tag = false;
    iso = false;
  }
  virtual ~Edge() { }

  // Refinement
  inline bool isRefined();
  inline bool isRegular();
  inline bool queryChildrenRefinement();

  inline void refine();
  inline void resetDescendant();

  // Geometry
  inline void setVertex(int i,Vertex* vt);
  inline void setVertices(Vertex* v0,Vertex* v1);
  inline void getVertices(Vertex** v0,Vertex** v1);
  inline Vertex* getVertex(int i);

  inline Vertex* getMidPoint();
  inline void    getChildren(Vertex* v0,Edge** ch0,Edge** ch1);
  inline void    getChildren(Edge** ch0,Edge** ch1);

  inline double length();
  inline void barycenter(double* bc);

  // reference count for database
  inline void ref();
  inline void unref();
  inline bool queryRef();


  // Iso-surface
  void setIsoMark()   { iso = true;  }
  void unsetIsoMark() { iso = false; }
  bool isIsoMarked()  { return iso;  }

  // memory management
  void * operator new(size_t) { return pool.alloc(); }
  void   operator delete(void *p) { pool.free(p); }

private:
  // edge geometry and topology
  Vertex* v[2];

  // edge tree, don't need parent
  Edge* child[2];

  // tag for vertex numbering
  bool tag;

  // tag for iso-sruface
  bool iso;

  // memory management
  static Pool<Edge> pool;
};



// ____________________________________________________________
//
//	
//	         inline functions
//
// ____________________________________________________________
//


// Description
//
inline void
Edge::setVertex(int i,Vertex* vt)
{
  v[i] = vt;
}


// Description
//
inline void
Edge::setVertices(Vertex* v0,Vertex* v1)
{
  v[0] = v0;
  v[1] = v1;
}


// Description
//
inline void
Edge::getVertices(Vertex** v0,Vertex** v1)
{
  *v0 = v[0];
  *v1 = v[1];
}


// Description
//
inline Vertex*
Edge::getVertex(int i)
{
  return v[i];
}


// Description
//
inline bool
Edge::isRefined()
{
  return (child[0] != 0);
}


// Description
//
inline bool
Edge::queryChildrenRefinement()
{
  if (child[0]) {
    if (child[0]->isRefined())
      return true;
    if (child[1]->isRefined())
      return true;
  }

  return false;
}


// Description
//
inline void
Edge::refine()
{
  Vertex* midPoint = new Vertex;
  for (int i = 0; i < 3; i++) {
    (*midPoint)[i] = 0.5*((*v[0])[i] + (*v[1])[i]);
  }

  child[0] = new Edge(v[0],midPoint);
  child[1] = new Edge(midPoint,v[1]);
}


// Description
//
void
Edge::resetDescendant()
{
  child[0] = 0;
  child[1] = 0;
}


// Description
//
inline Vertex*
Edge::getMidPoint()
{
  if (child[0])
    return (child[0]->v[1]);
  else {
    return 0;
  }
}



// Description
//
inline void
Edge::getChildren(Vertex* v0,Edge** ch0,Edge** ch1)
{
  if (child[0]) {
    if (v[0] == v0) {
      *ch0 = child[0];
      *ch1 = child[1];
    }
    else {
      *ch0 = child[1];
      *ch1 = child[0];
    }
  }
  else {
    *ch0 = 0;
    *ch1 = 0;
  }
}


// Description
//
inline void
Edge::getChildren(Edge** ch0,Edge** ch1)
{
  *ch0 = child[0];
  *ch1 = child[1];
}


// Description
//
inline double
Edge::length()
{
  double* pos1 = v[0]->getPosition();
  double* pos2 = v[1]->getPosition();
  double len = sqrt(square(pos1[0]-pos2[0]) + square(pos1[1]-pos2[1]) + square(pos1[2]-pos2[2]));
  return len;
}

// Description
//
inline void
Edge::barycenter(double* bc)
{
  double* pos1 = v[0]->getPosition();
  double* pos2 = v[1]->getPosition();
  bc[0] = 0.5*(pos1[0] + pos2[0]);
  bc[1] = 0.5*(pos1[1] + pos2[1]);
  bc[2] = 0.5*(pos1[2] + pos2[2]);
}


// Description
//
inline void
Edge::ref()
{
  tag = true;
}


// Description
//
inline void
Edge::unref()
{
  tag = false;
}


// Description
//
inline bool
Edge::queryRef()
{
  return tag;
}



} // namespace grd

#endif // EDGE_H
