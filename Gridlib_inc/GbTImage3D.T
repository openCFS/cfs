// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 09:27:31  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:21  prkipfer
| introduced image handling tool
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbTImage3D<T>::GbTImage3D (int uiXBound, int uiYBound, int uiZBound, T* atData)
  : GbTImage<T>(3)
{
  int* auiBound = new int[3];
  auiBound[0] = uiXBound;
  auiBound[1] = uiYBound;
  auiBound[2] = uiZBound;
  setBounds(auiBound);
  setData(atData);
}

template <class T>
GbTImage3D<T>::GbTImage3D (const GbTImage3D& rkImage)
  : GbTImage<T>(rkImage)
{
}

template <class T>
GbTImage3D<T>::GbTImage3D (const char* acFilename)
  : GbTImage<T>(acFilename)
{
}
