// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:27:30  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:19  prkipfer
| introduced image handling tool
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbTImage2D<T>::GbTImage2D (int uiXBound, int uiYBound, T* atData)
  : GbTImage<T>(2)
{
  int* auiBound = new int[2];
  auiBound[0] = uiXBound;
  auiBound[1] = uiYBound;
  setBounds(auiBound);
  setData(atData);
}

template <class T>
GbTImage2D<T>::GbTImage2D (const GbTImage2D& rkImage)
  : GbTImage<T>(rkImage)
{
}

template <class T>
GbTImage2D<T>::GbTImage2D (const char* acFilename)
  : GbTImage<T>(acFilename)
{
}
