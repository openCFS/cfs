// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/08/16 17:06:23  prkipfer
| revised memory management
|
| Revision 1.1  2001/01/09 11:10:01  prkipfer
| new util: TStorage and TClassStorage
|
|
|
+---------------------------------------------------------------------*/

template <class T>
void 
GbTStorage<T>::resize (int uiOldQuantity, int uiNewQuantity, GbBool bCopy)
{
  if ( uiNewQuantity != uiOldQuantity ) {
    if ( uiNewQuantity > 0 ) {
      T* atNewStorage = new T[uiNewQuantity];

      if ( bCopy ) {
	if ( uiNewQuantity > uiOldQuantity ) {
	  if ( uiOldQuantity > 0 ) {
	    int uiCopySize = uiOldQuantity*sizeof(T);
	    memcpy(atNewStorage,storage_,uiCopySize);
	  }
	}
	else {
	  memcpy(atNewStorage,storage_,uiNewQuantity*sizeof(T));
	}
      }

      delete[] storage_;
      storage_ = atNewStorage;
    }
    else {
      delete[] storage_;
      storage_ = NULL;
    }
  }
}
