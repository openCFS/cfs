/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOSCALARVERTEXBASE_HH
#define GOSCALARVERTEXBASE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>
#include <vector>

#include "GbTypes.hh"
#include "GbMemPool.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoScalarVertexBase
{
public: 

  //! Constructor: all values are initialized to zero, if not specified
  GoScalarVertexBase(T v = 0);
  ~GoScalarVertexBase();

  //! The scalar value at this vertex
  INLINE void setValue(T v);
  INLINE T getValue() const;

  //! Print memory pool statistics
  void poolStatistic() const;
  void shrink();

  //! The memory handling: we get all memory from a pool
  //! Note that these operators are automatic-inline
  static void *operator new(size_t size) { return pool_.alloc(size); }
  static void operator delete(void *p, size_t size) { pool_.free(p, size); }

  //! This operator pretty prints the vertex values
  friend std::ostream& operator<<(std::ostream&, const GoScalarVertexBase<T>&);

private:
  //! Inhibit assignment to a vertex
  GoScalarVertexBase<T>& operator=(const GoScalarVertexBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    T value;
  } ScalarVertexMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<ScalarVertexMem> pool_;
  ScalarVertexMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoScalarVertexBase<T>& v)
{
  s<<typeid(v).name()<<" value: "<<v.getValue()<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoScalarVertexBase.in"
#include "GoScalarVertexBase.T"

#else

#pragma instantiate GoScalarVertexBase<float>
#pragma instantiate GoScalarVertexBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoScalarVertexBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoScalarVertexBase<double>&)
#pragma instantiate GbMemPool<GoScalarVertexBase<float>::ScalarVertexMem>
#pragma instantiate GbMemPool<GoScalarVertexBase<double>::ScalarVertexMem>

#ifndef OUTLINE
#include "GoScalarVertexBase.in"
#endif

#endif // g++

#endif // GOSCALARVERTEXBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/01/02 15:16:47  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.1  2000/11/07 17:27:36  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
