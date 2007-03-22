// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
#include "Pool.hh"
#include "Utility.hh"
#include "Vertex.hh"


namespace grd {


class Edge {
public:
  // Constructors & Destructor
  Edge() {
    v[0] = v[1] = 0;
    child = 0;
    tag = 0;
  }
  Edge(Vertex* v0,Vertex* v1) {
    v[0] = v0;
    v[1] = v1;
    child = 0;
    tag = 0;
  }
  virtual ~Edge() { if (child) delete [] child; }

  // Refinement
  inline bool isRefined() const;
  inline bool isRegular() const;
  inline bool queryChildrenRefinement() const;

  inline void refine();
  inline void resetDescendant();

  // Geometry
  inline void setVertex(int i,Vertex* vt);
  inline void setVertices(Vertex* v0,Vertex* v1);
  inline void getVertices(Vertex** v0,Vertex** v1);
  inline Vertex* getVertex(int i);

  inline Vertex* getMidPoint();
  inline Edge*   getChild(int ch);
  inline void    getChildren(Vertex* v0,Edge** ch0,Edge** ch1);
  inline void    getChildren(Edge** ch0,Edge** ch1);

  inline double length();
  inline void barycenter(double* bc);

  // reference count for database
  inline void ref();
  inline void unref();
  inline bool queryRef() const;

  // Edge mark for general purpose
  inline void setMark();
  inline void unsetMark();
  inline bool isMarked() const;

  // memory management
  //void * operator new(size_t) { return pool.alloc(); }
  //void   operator delete(void *p) { pool.free(p); }

private:
  // edge geometry and topology
  Vertex* v[2];

  // edge tree, don't need parent
  Edge** child;

  // general flag
  unsigned char tag;

  // memory management
  //static Pool<Edge> pool;
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
Edge::isRefined() const
{
  return (child != 0);
}

// Description
// check if element is not refined, i.e. regular
inline bool
Edge::isRegular() const
{
  return (child == 0);
}

// Description
//
inline bool
Edge::queryChildrenRefinement() const
{
  if (child) {
    if (child[0]->isRefined())
      return true;
    else if (child[1]->isRefined())
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

  child = new Edge*[2];
  child[0] = new Edge(v[0],midPoint);
  child[1] = new Edge(midPoint,v[1]);
}


// Description
//
void
Edge::resetDescendant()
{
  if (child)
  {
    child[0] = 0;
    child[1] = 0;
    delete [] child;
    child = 0;
  }
}


// Description
//
inline Vertex*
Edge::getMidPoint()
{
  if (child)
    return (child[0]->v[1]);
  else 
    return 0;
}

// Description
//
inline Edge*   
Edge::getChild(int ch)
{
  return child[ch];
}

// Description
//
inline void
Edge::getChildren(Vertex* v0,Edge** ch0,Edge** ch1)
{
  if (child) {
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
  if (child)
  {
    *ch0 = child[0];
    *ch1 = child[1];
  }
  else
  {
    *ch0 = 0;
    *ch1 = 0;
  }
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


//
//  Flags
//
// Note: edge states and marks
// -------------------
// the variable tag is used bit wise
// Marks:
//   edge can be marked, e.g. for adaptive traversal is refined
//   and an edge can be referenced, e.g. for removing from the
//   edge list
// Edge States:
//  (i)   an edge is referenced or unreferenced
//  (ii)  an edge is marked or not marked
//
// Bit correspondence: 00 00 00 00
//   first bit for reference
//     00 00 00 00  not referenced
//     00 00 00 01  referenced
//   second bit for marked
//     00 00 00 xx  not marked
//     00 00 1x xx  marked
#define EDGE_MARK_REFERENCE    1
#define EDGE_MARK_GENERAL      2

// Description
//
inline void
Edge::ref()
{
  tag |= EDGE_MARK_REFERENCE;
}


// Description
//
inline void
Edge::unref()
{
  if (tag&EDGE_MARK_REFERENCE)
    tag ^= EDGE_MARK_REFERENCE;
}


// Description
//
inline bool
Edge::queryRef() const
{
  return (tag & EDGE_MARK_REFERENCE);
}



// Description
//
inline void
Edge::setMark()
{
  tag |= EDGE_MARK_GENERAL;
}

// Description
//
inline void
Edge::unsetMark()
{
  if (tag&EDGE_MARK_GENERAL)
    tag ^= EDGE_MARK_GENERAL;
}

// Description
//
inline bool
Edge::isMarked() const
{
  return (tag & EDGE_MARK_GENERAL);
}

} // namespace grd

#endif // EDGE_H
