// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.6  2002/03/18 10:00:25  prkipfer
| refactored element structure
|
| Revision 1.5  2001/01/02 15:16:46  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.4  2000/12/12 10:05:47  prkipfer
| added new element and vertex flags
|
| Revision 1.3  2000/11/07 17:27:32  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.2  2000/09/28 14:13:47  uflabsik
| errors corrected
|
| Revision 1.1  2000/09/26 09:29:24  uflabsik
| initial version
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoHexahedronElementBase<T>::HexahedronMem> GoHexahedronElementBase<T>::pool_;

template <class T>
GoHexahedronElementBase<T>::GoHexahedronElementBase()
{
  memset(data_.vertices,0,8*sizeof(GoVertex<T>*));
  memset(data_.edges,0,12*sizeof(GoEdge<T>*));
  memset(data_.neighbours,0,6*sizeof(GoGeometryElement<T>*));
  memset(data_.children,0,8*sizeof(GoGeometryElement<T>*));
}

template <class T>
GoHexahedronElementBase<T>::~GoHexahedronElementBase()
{
}

template <class T>
void
GoHexahedronElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoHexahedronElementBase<T>::shrink()
{
  pool_.shrink();
}
