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

template <class T> class GoGeometryElement;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoPrismElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoPrismElementBase(int i = -1, T x = 0, T y = 0, T z = 0);
  ~GoPrismElementBase();

  //! Operations on vertices of the triangle
  INLINE void setVertex(int i, GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex(int i) const;
  int findVertex(GoVertex<T> *v) const;
  GoVertex<T> *otherVertex(GoGeometryElement<T> *f) const;

  //! Other faces related to this triangle
  GoGeometryElement<T> *otherFace(GoVertex<T> *v) const;

  //! The face normal of the triangle
  void computeNormal();
  INLINE void getNormal(T& xx, T& yy, T& zz) const;
  INLINE GbVec3<T> getNormal() const;
  INLINE void setNormal(const GbVec3<T>&);

  //! Status flags of this triangle indicating semantic defined by the
  //! object using this class
  INLINE void setFlag(GbGeoStatusFlag f);
  INLINE void delFlag(GbGeoStatusFlag f);
  INLINE GbBool testFlag(GbGeoStatusFlag f) const;

  //! Integer to identify the triangle
  //! Has no meaning to this class's implementation
  INLINE void setId(int i);
  INLINE int getId() const;

  //! modification operations on the object
  INLINE GbBool subdivide();

  //! The scene part this object belongs to in the partition
  INLINE void setPartition(int i);
  INLINE int getPartition() const;

  //! The neighboring elements
  INLINE void setNeighbour(int i, GoGeometryElement<T> *face);
  void setNeighbour(GoGeometryElement<T> *face);
  INLINE GoGeometryElement<T> *getNeighbour(int i) const;
  int findNeighbour(GoGeometryElement<T> *face) const;

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
    GoGeometryElement<T>* neighbours[5];
    GbGeoStatusFlag flag;
    int id;
    int partition;
  } PrismMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<PrismMem> pool_;
  PrismMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoPrismElementBase<T>& v)
{
  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
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

