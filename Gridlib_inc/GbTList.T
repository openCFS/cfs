// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
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
