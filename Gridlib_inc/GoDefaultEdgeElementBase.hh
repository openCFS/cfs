/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GODEFAULTEDGEELEMENTBASE_HH
#define GODEFAULTEDGEELEMENTBASE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbMemPool.hh"

template <class T> class GoVertex;
template <class T> class GoGeometryElement;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoDefaultEdgeElementBase
{
public:

  //! Constructor: all values initialized to zero, if not spefified
  GoDefaultEdgeElementBase(int i = -1);
  ~GoDefaultEdgeElementBase();

  //! Operations on vertices of the edge
  INLINE void setVertex(int i, GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex(int i) const;
  int findVertex(GoVertex<T> *v) const;

  //! Status flags of this triangle indicating semantic defined by the
  //! object using this class
  INLINE void setFlag(GbEdgeStatusFlag f);
  INLINE void delFlag(GbEdgeStatusFlag f);
  INLINE GbBool testFlag(GbEdgeStatusFlag f) const;

  //! Integer to identify the edge
  //! Has no meaning to this class's implementation
  INLINE void setId(int i);
  INLINE int getId() const;

  //! The face object, this edge belongs to
  INLINE void setElement(GoGeometryElement<T> *f);
  INLINE GoGeometryElement<T> *getElement() const;

  //! The neighbor valence of the edge
  INLINE void setValence(int v);
  INLINE int getValence() const;

  //! Print memory pool statistics
  void poolStatistic() const;
  void shrink();

  //! The memory handling: we get all memory from a pool
  //! Note that these operators are automatic-inline
  static void *operator new(size_t size) { return pool_.alloc(size); }
  static void operator delete(void *p, size_t size) { pool_.free(p, size); }

  //! This operator pretty prints the vertex values
  friend std::ostream& operator<<(std::ostream&, const GoDefaultEdgeElementBase<T>&);

private:
  //! Inhibit assignment to edge
  GoDefaultEdgeElementBase<T>& operator=(const GoDefaultEdgeElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[2];
    GoGeometryElement<T>* element;
    int valence;
    GbEdgeStatusFlag flag;
    int id;
  } DefaultEdgeMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<DefaultEdgeMem> pool_;
  DefaultEdgeMem data_;
};


template<class T>
std::ostream&
operator<<(std::ostream& s, const GoDefaultEdgeElementBase<T>& v)
{
  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
  return s;
}


#ifdef __GNUG__

#include "GoDefaultEdgeElementBase.in"
#include "GoDefaultEdgeElementBase.T"

#else

#pragma instantiate GoDefaultEdgeElementBase<float>
#pragma instantiate GoDefaultEdgeElementBase<double>
#pragma instantiate GbMemPool<GoDefaultEdgeElementBase<float>::DefaultEdgeMem>
#pragma instantiate GbMemPool<GoDefaultEdgeElementBase<double>::DefaultEdgeMem>

#endif // g++

#endif // GODEFAULTEDGEELEMENTBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 14:47:44  prkipfer
| introduced edge elements
|
|
+---------------------------------------------------------------------*/
