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

template <class T> class GoGeometryElement;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoQuadElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoQuadElementBase(int i = -1, T x = 0, T y = 0, T z = 0);
  ~GoQuadElementBase();

  //! Operations on vertices of the quad
  INLINE void setVertex(int i, GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex(int i) const;
  int findVertex(GoVertex<T> *v) const;
  GoVertex<T> *otherVertex(GoGeometryElement<T> *f) const;

  //! Other faces related to this quad
  GoGeometryElement<T> *otherFace(GoVertex<T> *v) const;

  //! The face normal of the quad
  void computeNormal();
  INLINE void getNormal(T& xx, T& yy, T& zz) const;
  INLINE GbVec3<T> getNormal() const;
  INLINE void setNormal(const GbVec3<T>&);

  //! Status flags of this quad indicating semantic defined by the
  //! object using this class
  INLINE void setFlag(GbGeoStatusFlag f);
  INLINE void delFlag(GbGeoStatusFlag f);
  INLINE GbBool testFlag(GbGeoStatusFlag f) const;

  //! Integer to identify the quad
  //! Has no meaning to this class's implementation
  INLINE void setId(int i);
  INLINE int getId() const;

  //! modification operations on the object
  INLINE GbBool subdivide();

   //! The scene part this object belongs to in the partition
  INLINE void setPartition(int i);
  INLINE int getPartition() const;

  //! The neighboring quads
  INLINE void setNeighbour(int i, GoGeometryElement<T> *face);
  void setNeighbour(GoGeometryElement<T> *face);
  INLINE GoGeometryElement<T> *getNeighbour(int i) const;
  int findNeighbour(GoGeometryElement<T> *v) const;

  //! The parent and the child triangles
  INLINE GoGeometryElement<T> *getChild(int i) const;
  INLINE void setChild(int i, GoGeometryElement<T> *face);
  INLINE GoGeometryElement<T> *getParent() const;
  INLINE void setParent(GoGeometryElement<T> *face);

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
    GoGeometryElement<T>* neighbours[4];
    GoGeometryElement<T>* parent;
    GoGeometryElement<T>* children[4];
    T normal[3];
    GbGeoStatusFlag flag;
    int id;
    int partition;
  } QuadMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<QuadMem> pool_;
  QuadMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoQuadElementBase<T>& v)
{
  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
  s<<"\tn("<<v.getNormal()[0]<<","<<v.getNormal()[1]<<","<<v.getNormal()[2]<<") "<<std::endl;
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

