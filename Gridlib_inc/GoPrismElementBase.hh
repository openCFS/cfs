/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOPRISMELEMENTBASE_HH
#define GOPRISMELEMENTBASE_HH

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
class GoPrismElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoPrismElementBase();
  ~GoPrismElementBase();

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
  friend std::ostream& operator<<(std::ostream&, const GoPrismElementBase<T>&);

private:
  //! Inhibit assignment to a triangle
  GoPrismElementBase<T>& operator=(const GoPrismElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[6];
    GoEdge<T>* edges[9];
    GoGeometryElement<T>* neighbours[5];
    GoGeometryElement<T>* children[2];
  } PrismMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<PrismMem> pool_;
  PrismMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoPrismElementBase<T>& v)
{
//  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoPrismElementBase.in"
#include "GoPrismElementBase.T"

#else

#pragma instantiate GoPrismElementBase<float>
#pragma instantiate GoPrismElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoPrismElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoPrismElementBase<double>&)
#pragma instantiate GbMemPool<GoPrismElementBase<float>::PrismMem>
#pragma instantiate GbMemPool<GoPrismElementBase<double>::PrismMem>

#ifndef OUTLINE
#include "GoPrismElementBase.in"
#endif

#endif // g++

#endif // GOPRISMELEMENTBASE_HH

