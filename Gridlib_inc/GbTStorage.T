// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
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
