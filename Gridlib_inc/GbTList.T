// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
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
