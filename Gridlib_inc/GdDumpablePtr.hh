/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GDDUMPABLEPTR_HH
#define GDDUMPABLEPTR_HH

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

class GdDumpablePtr
{
public:
  GdDumpablePtr(const GdDumpable *dumper = 0);

  //! Smart pointer delegation method
  const GdDumpable *operator->() const;

  //! Assignment operator
  void operator=(const GdDumpable *dumper) const;

private:
  //! Points to the actual GdDumpable
  const GdDumpable *dumper_;
};

//#ifndef OUTLINE
//#include "GdDumpablePtr.in"
//#endif

#endif // GDDUMPABLEPTR_HH
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
