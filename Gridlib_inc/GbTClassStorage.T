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
| Revision 1.1  2001/01/09 11:10:00  prkipfer
| new util: TStorage and TClassStorage
|
|
|
+---------------------------------------------------------------------*/

template <class T>
void 
GbTClassStorage<T>::resize (int uiOldQuantity, int uiNewQuantity, GbBool bCopy)
{
  if ( uiNewQuantity != uiOldQuantity ) {
    if ( uiNewQuantity > 0 ) {
      T* atNewStorage = new T[uiNewQuantity];

      if ( bCopy ) {
	int uiIMax = ( uiNewQuantity > uiOldQuantity ? uiOldQuantity : uiNewQuantity );

	for (int uiI = 0; uiI < uiIMax; uiI++)
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
