/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOBOUNDARYVERTEXBASE_HH
#define GOBOUNDARYVERTEXBASE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>
#include <vector>

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbMemPool.hh"

template <class T> class GoGeometryElement;
template <class T> class GoVertex;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoBoundaryVertexBase
{
public: 

  //! Constructor: all values are initialized to zero, if not specified
  GoBoundaryVertexBase(GoVertex<T> *vel = NULL);
  ~GoBoundaryVertexBase();

  //! The corresponding vertex in the volume mesh
  INLINE void setVertex(GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex() const;

  //! Print memory pool statistics
  void poolStatistic() const;
  void shrink();

  //! The memory handling: we get all memory from a pool
  //! Note that these operators are automatic-inline
  static void *operator new(size_t size) { return pool_.alloc(size); }
  static void operator delete(void *p, size_t size) { pool_.free(p, size); }

  //! This operator pretty prints the vertex values
  friend std::ostream& operator<<(std::ostream&, const GoBoundaryVertexBase<T>&);

private:
  //! Inhibit assignment to a vertex
  GoBoundaryVertexBase<T>& operator=(const GoBoundaryVertexBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertex;
  } BoundaryVertexMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<BoundaryVertexMem> pool_;
  BoundaryVertexMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoBoundaryVertexBase<T>& v)
{
  s<<typeid(v).name()<<" boundary link to vertex: "<<v.getVertex()->getId()<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoBoundaryVertexBase.in"
#include "GoBoundaryVertexBase.T"

#else

#pragma instantiate GoBoundaryVertexBase<float>
#pragma instantiate GoBoundaryVertexBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoBoundaryVertexBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoBoundaryVertexBase<double>&)
#pragma instantiate GbMemPool<GoBoundaryVertexBase<float>::BoundaryVertexMem>
#pragma instantiate GbMemPool<GoBoundaryVertexBase<double>::BoundaryVertexMem>

#ifndef OUTLINE
#include "GoBoundaryVertexBase.in"
#endif

#endif // g++

#endif // GOBOUNDARYVERTEXBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2002/03/18 10:00:23  prkipfer
| refactored element structure
|
| Revision 1.1  2001/02/13 11:00:39  prkipfer
| added GoBoundaryVertexBase class
|
|
+---------------------------------------------------------------------*/
