/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GODEFAULTVERTEXBASE_HH
#define GODEFAULTVERTEXBASE_HH

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


/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoDefaultVertexBase
{
public: 

  //! Constructor: all values are initialized to zero, if not specified
  GoDefaultVertexBase(int i = -1, T x = 0, T y = 0, T z = 0, T nx = 0, T ny = 0, T nz = 0);
  ~GoDefaultVertexBase();

  //! The position in 3D space
  INLINE void setPosition(T x, T y, T z);
  INLINE void setPosition(const GbVec3<T>& v);
  INLINE void getPosition(T& xx, T& yy, T& zz) const;
  INLINE GbVec3<T> getPosition() const;

  //! A normal at this vertex in 3D space
  INLINE void setNormal(T x, T y, T z);
  INLINE void setNormal(const GbVec3<T>& v);
  INLINE void getNormal(T& xx, T& yy, T& zz) const;
  INLINE GbVec3<T> getNormal() const;

  //! Integer to identify the vertex
  //! Has no meaning to this class's implementation
  INLINE void setId(int i);
  INLINE int getId() const;

  //! Status flags of this vertex indicating semantic defined by the
  //! object using this class
  INLINE void setFlag(GbVertexFlag f);
  INLINE void delFlag(GbVertexFlag f);
  INLINE GbBool testFlag(GbVertexFlag f) const;

  //! The face object, this vertex belongs to
  INLINE void setElement(GoGeometryElement<T> *f);
  INLINE GoGeometryElement<T> *getElement() const;

  //! The neighbor valence of the vertex
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
  friend std::ostream& operator<<(std::ostream&, const GoDefaultVertexBase<T>&);

private:
  //! Inhibit assignment to a vertex
  GoDefaultVertexBase<T>& operator=(const GoDefaultVertexBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    T pos[3];
    T normal[3];
    int id;
    GbVertexFlag flag;
    GoGeometryElement<T> *face;
    int valence;
  } DefaultVertexMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<DefaultVertexMem> pool_;
  DefaultVertexMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoDefaultVertexBase<T>& v)
{
  s<<typeid(v).name()<<" id: "<<v.getId()<<" valence: "<<v.getValence()<<std::endl;
  s<<"\tv("<<v.getPosition()[0]<<","<<v.getPosition()[1]<<","<<v.getPosition()[2]<<") ";
  s<<"n("<<v.getNormal()[0]<<","<<v.getNormal()[1]<<","<<v.getNormal()[2]<<") "<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoDefaultVertexBase.in"
#include "GoDefaultVertexBase.T"

#else

#pragma instantiate GoDefaultVertexBase<float>
#pragma instantiate GoDefaultVertexBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoDefaultVertexBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoDefaultVertexBase<double>&)
#pragma instantiate GbMemPool<GoDefaultVertexBase<float>::DefaultVertexMem>
#pragma instantiate GbMemPool<GoDefaultVertexBase<double>::DefaultVertexMem>

#ifndef OUTLINE
#include "GoDefaultVertexBase.in"
#endif

#endif // g++

#endif // GODEFAULTVERTEXBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.3  2002/03/18 10:00:25  prkipfer
| refactored element structure
|
| Revision 1.2  2001/01/02 15:16:45  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.1  2000/11/07 17:27:30  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
