/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOPYRAMIDELEMENTBASE_HH
#define GOPYRAMIDELEMENTBASE_HH

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
class GoPyramidElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoPyramidElementBase();
  ~GoPyramidElementBase();

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
  friend std::ostream& operator<<(std::ostream&, const GoPyramidElementBase<T>&);

private:
  //! Inhibit assignment to a triangle
  GoPyramidElementBase<T>& operator=(const GoPyramidElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[5];
    GoEdge<T>* edges[8];
    GoGeometryElement<T>* neighbours[5];
    GoGeometryElement<T>* children[5];
  } PyramidMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<PyramidMem> pool_;
  PyramidMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoPyramidElementBase<T>& v)
{
//  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoPyramidElementBase.in"
#include "GoPyramidElementBase.T"

#else

#pragma instantiate GoPyramidElementBase<float>
#pragma instantiate GoPyramidElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoPyramidElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoPyramidElementBase<double>&)
#pragma instantiate GbMemPool<GoPyramidElementBase<float>::PyramidMem>
#pragma instantiate GbMemPool<GoPyramidElementBase<double>::PyramidMem>

#ifndef OUTLINE
#include "GoPyramidElementBase.in"
#endif

#endif // g++

#endif // GOPYRAMIDELEMENTBASE_HH

