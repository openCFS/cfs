/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOHEXAHEDRONELEMENTBASE_HH
#define GOHEXAHEDRONELEMENTBASE_HH

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
class GoHexahedronElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoHexahedronElementBase();
  ~GoHexahedronElementBase();

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
  friend std::ostream& operator<<(std::ostream&, const GoHexahedronElementBase<T>&);

private:
  //! Inhibit assignment to a triangle
  GoHexahedronElementBase<T>& operator=(const GoHexahedronElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[8];
    GoEdge<T>* edges[12];
    GoGeometryElement<T>* neighbours[6];
    GoGeometryElement<T>* children[8];
  } HexahedronMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<HexahedronMem> pool_;
  HexahedronMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoHexahedronElementBase<T>& v)
{
//  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoHexahedronElementBase.in"
#include "GoHexahedronElementBase.T"

#else

#pragma instantiate GoHexahedronElementBase<float>
#pragma instantiate GoHexahedronElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoHexahedronElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoHexahedronElementBase<double>&)
#pragma instantiate GbMemPool<GoHexahedronElementBase<float>::HexahedronMem>
#pragma instantiate GbMemPool<GoHexahedronElementBase<double>::HexahedronMem>

#ifndef OUTLINE
#include "GoHexahedronElementBase.in"
#endif

#endif // g++

#endif // GOHEXAHEDRONELEMENTBASE_HH

