/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBTCLASSSTORAGE_HH
#define  GBTCLASSSTORAGE_HH

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

/*! The class T is intended to be class data and must have the following
    member functions:
    T::T ()
    T& T::operator= (const T&)

    The application is responsible for maintaining the quantity in storage.
*/
template <class T>
class GbTClassStorage
{
public:
  // construction and destruction
  INLINE GbTClassStorage ();
  INLINE GbTClassStorage (int uiQuantity);
  INLINE GbTClassStorage (T* atStorage);
  INLINE ~GbTClassStorage ();

  // element access
  INLINE T* getStorage () const;
  INLINE T& operator[] (int uiIndex) const;

  // reorganization (new storage initialized by default constructor)
  void resize (int uiOldQuantity, int uiNewQuantity, GbBool bCopy);

protected:
  T* storage_;

private:
  // do not copy
  GbTClassStorage(const GbTClassStorage<T>& rhs);
};


#ifdef __GNUG__

#include "GbTClassStorage.in"
#include "GbTClassStorage.T"

#else

#ifndef OUTLINE
#include "GbTClassStorage.in"
#endif

#endif  // g++

#endif  // GBTCLASSSTORAGE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2001/08/16 17:06:23  prkipfer
| revised memory management
|
| Revision 1.2  2001/02/13 11:07:00  prkipfer
| added instantiations for rendering subsystem
|
| Revision 1.1  2001/01/09 11:10:01  prkipfer
| new util: TStorage and TClassStorage
|
|
+---------------------------------------------------------------------*/
