/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBTSTORAGE_HH
#define  GBTSTORAGE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*! The class T is intended to be native data (int, float, char*, etc.).  The
    resize member function does memory copies knowing that T does not have
    side effects from construction or assignment.  The application is
    responsible for maintaining the quantity in storage.
*/
template <class T>
class GbTStorage
{
public:
  // construction and destruction
  INLINE GbTStorage ();
  INLINE GbTStorage (int uiQuantity);
  INLINE GbTStorage (T* atStorage);
  INLINE ~GbTStorage ();

  // element access
  INLINE T* getStorage () const;
  INLINE T& operator[] (int uiIndex) const;

  // reorganization (new storage is uninitialized)
  void resize (int uiOldQuantity, int uiNewQuantity, GbBool bCopy);

protected:
  T* storage_;
};

#ifdef __GNUG__

#include "GbTStorage.in"
#include "GbTStorage.T"

#else

#pragma instantiate GbTStorage<float>
#pragma instantiate GbTStorage<double>
#pragma instantiate GbTStorage<GbBool>
#pragma instantiate GbTStorage<unsigned int>
#pragma instantiate GbTStorage<int>

#ifndef OUTLINE
#include "GbTStorage.in"
#endif

#endif  // g++

#endif  // GBTSTORAGE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/08/16 17:06:23  prkipfer
| revised memory management
|
| Revision 1.1  2001/01/09 11:10:02  prkipfer
| new util: TStorage and TClassStorage
|
|
+---------------------------------------------------------------------*/
