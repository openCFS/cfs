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
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1.1.1  2000/06/08 16:24:44  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
