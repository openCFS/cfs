/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBMEMPOOL_HH
#define GBMEMPOOL_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>
//  #include <sys/types.h>
#include <vector>

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template<class T>
class GbMemPool
{
public:
  //! Constructor: init nothing. All variables are class-static.
  //! Memory is allocated by alloc()
  GbMemPool();
  ~GbMemPool();

  //! Returns a pointer to enough memory for sizeof(one element).
  //! Try to allocate new memory block, if exhausted
  void *alloc(size_t n);

  //! Return the space to the pool but don't delete it from memory.
  //! Instead it can be reallocated by future alloc() calls.
  void free(void *p, size_t n);

  //! Query the number of elements currently in the pool
  static unsigned long getPoolSize() { return poolSize_; }
  void getStatistic();

  //! Query the increment of the pool, if the memory is exhausted:
  //! the pool will grow by this value * sizeof(one element)
  static int getPoolIncrement() { return poolIncrement_; }

  void shrink();

  //! This operator displays memory pool statistics
  friend std::ostream& operator<<(std::ostream&, const GbMemPool<T>&);

protected:
  //! Fold a pointer to unused space into the memory for one element.
  //! This implies, that the minimum size of one element is sizeof(pointer)
  union {
    GbMemPool<T> *next_;
    T data_;
  };

private:
  //! return memblocks allocated to system
  void freePoolMem();

  //! for pretty printing a number
  float formatNumber(unsigned long a, char* s);

  //! Number of elements the pool grows each time memory is allocated
  static const int poolIncrement_;
  //! Number of elements currently in the pool
  static unsigned long poolSize_;
  //! Pointer to the next free space
  static GbMemPool<T> *headOfFreeList_;
  //! Array of allocated memory blocks, needed for proper destruction
  static GbMemPool<T> *poolMem_[MEMPOOLBLOCKS];
  static unsigned int poolMemSize_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbMemPool<T>& p)
{
  unsigned long b = poolSize_*sizeof(T);
  unsigned long c = poolMemSize_*poolIncrement_*sizeof(T);

  s<<typeid(p).name()<<" statistics:"<<std::endl;
  s<<"items:  "<<poolSize_<<" elements X "<<sizeof(T)<<" bytes = "<<b<<" bytes"<<std::endl;
  s<<"memory: "<<poolMemSize_<<" blocks X "<<poolIncrement_<<" items = "<<c<<" bytes"<<std::endl;
  if(c>0)
  s<<"\t--> efficiency "<<(b*100)/c<<"%"<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GbMemPool.in"
#include "GbMemPool.T"

#else

#ifndef OUTLINE
#include "GbMemPool.in"
#endif

#endif // g++

#endif // GBMEMPOOL_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.7  2001/08/16 16:53:19  prkipfer
| improved type safety for template parameter
|
| Revision 1.6  2001/06/15 08:36:56  prkipfer
| fixed range overflow bug
|
| Revision 1.5  2000/07/21 10:40:44  prkipfer
| dropped g++ support
|
| Revision 1.4  2000/07/20 13:32:38  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.3  2000/07/17 10:48:36  prkipfer
| changed to static memory block handling because g++ doesn't like partially specialized templates
|
| Revision 1.2  2000/06/14 15:39:13  prkipfer
| improved base classes and added funcstruct processing for mesh
|
| Revision 1.1.1.1  2000/06/08 16:24:44  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
