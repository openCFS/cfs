/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOQUADELEMENTBASE_HH
#define GOQUADELEMENTBASE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GbVec3.hh"
#include "GbMemPool.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"

template <class T> class GoGeometryElement;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoQuadElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoQuadElementBase();
  ~GoQuadElementBase();

  //! Operations on vertices of the quad
  INLINE void setVertex(int i, GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex(int i) const;

  //! Operations on edges of the quad
  INLINE void setEdge(int i, GoEdge<T> *e);
  INLINE GoEdge<T> *getEdge(int i) const;

  //! The neighboring quads
  INLINE void setNeighbour(int i, GoGeometryElement<T> *face);
  INLINE GoGeometryElement<T> *getNeighbour(int i) const;

  //! The parent and the child triangles
  INLINE GoGeometryElement<T> *getChild(int i) const;
  INLINE void setChild(int i, GoGeometryElement<T> *face);

  //! Print memory pool statistics
  void poolStatistic() const;
  void shrink();

  //! The memory handling: we get all memory from a pool
  //! Note that these operators are automatic-inline
  static void *operator new(size_t size) { return pool_.alloc(size); }
  static void operator delete(void *p, size_t size) { pool_.free(p, size); }

  //! This operator pretty prints the vertex values
  friend std::ostream& operator<<(std::ostream&, const GoQuadElementBase<T>&);

private:
  //! Inhibit assignment to a quad
  GoQuadElementBase<T>& operator=(const GoQuadElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[4];
    GoEdge<T>* edges[4];
    GoGeometryElement<T>* neighbours[4];
    GoGeometryElement<T>* children[4];
  } QuadMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<QuadMem> pool_;
  QuadMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoQuadElementBase<T>& v)
{
//  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
//  s<<"\tn("<<v.getNormal()[0]<<","<<v.getNormal()[1]<<","<<v.getNormal()[2]<<") "<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoQuadElementBase.in"
#include "GoQuadElementBase.T"

#else

#pragma instantiate GoQuadElementBase<float>
#pragma instantiate GoQuadElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoQuadElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoQuadElementBase<double>&)
#pragma instantiate GbMemPool<GoQuadElementBase<float>::QuadMem>
#pragma instantiate GbMemPool<GoQuadElementBase<double>::QuadMem>

#ifndef OUTLINE
#include "GoQuadElementBase.in"
#endif

#endif // g++

#endif // GOQUADELEMENTBASE_HH

