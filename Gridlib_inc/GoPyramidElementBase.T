// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.5  2002/03/18 10:00:26  prkipfer
| refactored element structure
|
| Revision 1.4  2000/12/12 10:05:47  prkipfer
| added new element and vertex flags
|
| Revision 1.3  2000/11/07 17:27:34  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.2  2000/09/28 14:15:13  uflabsik
| errors corrected
|
| Revision 1.1  2000/09/26 09:29:43  uflabsik
| initial version
|
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoPyramidElementBase<T>::PyramidMem> GoPyramidElementBase<T>::pool_;

template <class T>
GoPyramidElementBase<T>::GoPyramidElementBase()
{
  memset(data_.vertices,0,5*sizeof(GoVertex<T>*));
  memset(data_.edges,0,8*sizeof(GoEdge<T>*));
  memset(data_.neighbours,0,5*sizeof(GoGeometryElement<T>*));
  memset(data_.children,0,5*sizeof(GoGeometryElement<T>*));
}

template <class T>
GoPyramidElementBase<T>::~GoPyramidElementBase()
{
}

template <class T>
void
GoPyramidElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoPyramidElementBase<T>::shrink()
{
  pool_.shrink();
}
