/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GDDUMPABLEADAPTER_HH
#define GDDUMPABLEADAPTER_HH

#ifndef INLINE
#ifdef OUTLINE
#define INLINE
#else
#define INLINE inline
#endif
#endif

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

using namespace std;

#include "GdDumpable.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GdDumpableAdapter
  : public GdDumpable
{
public:
  GdDumpableAdapter(const T* t);

  //! Concrete dump method (simply delegates to the <dump> method
  //! of <class T>)
  virtual void dump() const;

  //! Delegate to methods in class T
  T *operator->();

private:
  //! Pointer to <this> of <class T>
  const T *this_;
};

template <class T>
void
dump(const T* t)
{
  t->dump();
}

//#ifndef OUTLINE
//#include "GdDumpableAdapter.in"
//#endif

#endif // GDDUMPABLEADAPTER_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1.1.1  2000/06/08 16:24:44  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
