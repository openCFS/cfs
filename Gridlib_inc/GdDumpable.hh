/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GDDUMPABLE_HH
#define GDDUMPABLE_HH

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

class GdDumpablePtr;
class GdObjectDB;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GdDumpable
{
  friend class GdObjectDB;
  friend class GdDumpablePtr;

public:
  GdDumpable(const void *);

  //! This pure virtual method must be filled in by a subclass
  virtual void dump() const = 0;

protected:
  virtual ~GdDumpable();

private:
  const void *this_;
};

//#ifndef OUTLINE
//#include "GdDumpable.in"
//#endif

#endif // GDDUMPABLE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1.1.1  2000/06/08 16:24:43  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
