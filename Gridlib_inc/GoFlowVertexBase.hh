/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOFLOWVERTEXBASE_HH
#define GOFLOWVERTEXBASE_HH

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

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoFlowVertexBase
{
public: 

  //! Constructor: all values are initialized to zero, if not specified
  GoFlowVertexBase(T d = 0, T mx = 0, T my = 0, T mz = 0, T e = 0);
  ~GoFlowVertexBase();

  //! The scalar values at this vertex
  INLINE void setDensity(T v);
  INLINE T getDensity() const;
  INLINE void setEnergy(T v);
  INLINE T getEnergy() const;
  INLINE void setMomentum(const GbVec3<T>& v);
  INLINE GbVec3<T> getMomentum() const;

  //! Print memory pool statistics
  void poolStatistic() const;
  void shrink();

  //! The memory handling: we get all memory from a pool
  //! Note that these operators are automatic-inline
  static void *operator new(size_t size) { return pool_.alloc(size); }
  static void operator delete(void *p, size_t size) { pool_.free(p, size); }

  //! This operator pretty prints the vertex values
  friend std::ostream& operator<<(std::ostream&, const GoFlowVertexBase<T>&);

private:
  //! Inhibit assignment to a vertex
  GoFlowVertexBase<T>& operator=(const GoFlowVertexBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    T density;
    T momentum[3];
    T energy;
  } FlowVertexMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<FlowVertexMem> pool_;
  FlowVertexMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoFlowVertexBase<T>& v)
{
  s<<typeid(v).name()<<" density: "<<v.getDensity()<<std::endl;
  s<<"momentum: ("<<v.getMomentum()[0]<<","<<v.getMomentum()[1]<<","<<v.getMomentum()[2]<<")"<<std::endl;
  s<<"energy: "<<v.getEnergy()<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoFlowVertexBase.in"
#include "GoFlowVertexBase.T"

#else

#pragma instantiate GoFlowVertexBase<float>
#pragma instantiate GoFlowVertexBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoFlowVertexBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoFlowVertexBase<double>&)
#pragma instantiate GbMemPool<GoFlowVertexBase<float>::FlowVertexMem>
#pragma instantiate GbMemPool<GoFlowVertexBase<double>::FlowVertexMem>

#ifndef OUTLINE
#include "GoFlowVertexBase.in"
#endif

#endif // g++

#endif // GOFLOWVERTEXBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/01/02 15:16:46  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.1  2000/11/07 17:27:32  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
