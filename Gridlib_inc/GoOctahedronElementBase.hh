/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOOCTAHEDRONELEMENTBASE_HH
#define GOOCTAHEDRONELEMENTBASE_HH

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
class GoOctahedronElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoOctahedronElementBase();
  ~GoOctahedronElementBase();

  //! Operations on vertices of the triangle
  INLINE void setVertex(int i, GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex(int i) const;

  //! Operations on edges of the quad
  INLINE void setEdge(int i, GoEdge<T> *e);
  INLINE GoEdge<T> *getEdge(int i) const;

  //! The neighboring elements
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
  friend std::ostream& operator<<(std::ostream&, const GoOctahedronElementBase<T>&);

private:
  //! Inhibit assignment to a triangle
  GoOctahedronElementBase<T>& operator=(const GoOctahedronElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[6];
    GoEdge<T>* edges[12];
    GoGeometryElement<T>* neighbours[8];
    GoGeometryElement<T>* children[10];
  } OctahedronMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<OctahedronMem> pool_;
  OctahedronMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoOctahedronElementBase<T>& v)
{
//  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoOctahedronElementBase.in"
#include "GoOctahedronElementBase.T"

#else

#pragma instantiate GoOctahedronElementBase<float>
#pragma instantiate GoOctahedronElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoOctahedronElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoOctahedronElementBase<double>&)
#pragma instantiate GbMemPool<GoOctahedronElementBase<float>::OctahedronMem>
#pragma instantiate GbMemPool<GoOctahedronElementBase<double>::OctahedronMem>

#ifndef OUTLINE
#include "GoOctahedronElementBase.in"
#endif

#endif // g++

#endif // GOOCTAHEDRONELEMENTBASE_HH

