// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
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
