// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/09 11:10:00  prkipfer
| new util: TStorage and TClassStorage
|
|
|
+---------------------------------------------------------------------*/

template <class T>
void 
GbTClassStorage<T>::resize (unsigned int uiOldQuantity,
			    unsigned int uiNewQuantity, GbBool bCopy)
{
  if ( uiNewQuantity != uiOldQuantity ) {
    if ( uiNewQuantity > 0 ) {
      T* atNewStorage = new T[uiNewQuantity];

      if ( bCopy ) {
	unsigned int uiIMax = ( uiNewQuantity > uiOldQuantity ?	uiOldQuantity : uiNewQuantity );

	for (unsigned int uiI = 0; uiI < uiIMax; uiI++)
	  atNewStorage[uiI] = storage_[uiI];
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
