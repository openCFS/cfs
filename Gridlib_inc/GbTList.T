// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/04/23 10:02:23  prkipfer
| introduced morphological image operations and TList class
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbTList<T>::GbTList ()
  : quantity_(0)
  , front_(NULL)
  , iterator_(NULL)
{
}

template <class T>
GbTList<T>::~GbTList ()
{
  removeAll();
}
