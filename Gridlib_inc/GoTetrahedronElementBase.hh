/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOTETRAHEDRONELEMENTBASE_HH
#define GOTETRAHEDRONELEMENTBASE_HH

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
class GoTetrahedronElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoTetrahedronElementBase(int i = -1, T x = 0, T y = 0, T z = 0);
  ~GoTetrahedronElementBase();

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
  friend std::ostream& operator<<(std::ostream&, const GoTetrahedronElementBase<T>&);

private:
  //! Inhibit assignment to a triangle
  GoTetrahedronElementBase<T>& operator=(const GoTetrahedronElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[4];
    GoGeometryElement<T>* neighbours[4];
    T normal[3];
    GbGeoStatusFlag flag;
    int id;
    int partition;
  } TetraMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<TetraMem> pool_;
  TetraMem data_;

};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoTetrahedronElementBase<T>& v)
{
  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
  s<<"\tn("<<v.getNormal()[0]<<","<<v.getNormal()[1]<<","<<v.getNormal()[2]<<") "<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoTetrahedronElementBase.in"
#include "GoTetrahedronElementBase.T"

#else

#pragma instantiate GoTetrahedronElementBase<float>
#pragma instantiate GoTetrahedronElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoTetrahedronElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoTetrahedronElementBase<double>&)
#pragma instantiate GbMemPool<GoTetrahedronElementBase<float>::TetraMem>
#pragma instantiate GbMemPool<GoTetrahedronElementBase<double>::TetraMem>

#ifndef OUTLINE
#include "GoTetrahedronElementBase.in"
#endif

#endif // g++

#endif // GOTETRAHEDRONELEMENTBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.7  2001/01/02 15:16:47  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.6  2000/11/07 17:27:39  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.5  2000/09/28 14:15:13  uflabsik
| errors corrected
|
| Revision 1.4  2000/09/03 14:05:42  uflabsik
| functionality added
|
| Revision 1.3  2000/08/28 09:25:22  prkipfer
| added localized subdivision mechanism
|
| Revision 1.2  2000/07/21 10:40:45  prkipfer
| dropped g++ support
|
| Revision 1.1  2000/07/03 16:15:03  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
