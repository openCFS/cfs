// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/09 11:10:01  prkipfer
| new util: TStorage and TClassStorage
|
|
|
+---------------------------------------------------------------------*/

template <class T>
void 
GbTStorage<T>::resize (unsigned int uiOldQuantity,
		       unsigned int uiNewQuantity, GbBool bCopy)
{
  if ( uiNewQuantity != uiOldQuantity ) {
    if ( uiNewQuantity > 0 ) {
      T* atNewStorage = new T[uiNewQuantity];

      if ( bCopy ) {
	if ( uiNewQuantity > uiOldQuantity ) {
	  if ( uiOldQuantity > 0 ) {
	    unsigned int uiCopySize = uiOldQuantity*sizeof(T);
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
