// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.5  2002/03/18 10:00:25  prkipfer
| refactored element structure
|
| Revision 1.4  2000/12/12 10:05:47  prkipfer
| added new element and vertex flags
|
| Revision 1.3  2000/11/07 17:27:33  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.2  2000/09/28 14:15:13  uflabsik
| errors corrected
|
| Revision 1.1  2000/09/26 09:30:07  uflabsik
| initial version
|
| Revision 1.4  2000/09/03 14:05:42  uflabsik
| functionality added
|
| Revision 1.3  2000/08/28 09:25:22  prkipfer
| added localized subdivision mechanism
|
| Revision 1.2  2000/07/11 07:56:47  uflabsik
| *** empty log message ***
|
| Revision 1.1  2000/07/03 16:15:02  prkipfer
| extended IO subsystem, created class registry, added prism
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoPrismElementBase<T>::PrismMem> GoPrismElementBase<T>::pool_;

template <class T>
GoPrismElementBase<T>::GoPrismElementBase()
{
  memset(data_.vertices,0,6*sizeof(GoVertex<T>*));
  memset(data_.edges,0,9*sizeof(GoEdge<T>*));
  memset(data_.neighbours,0,5*sizeof(GoGeometryElement<T>*));
  memset(data_.children,0,2*sizeof(GoGeometryElement<T>*));
}

template <class T>
GoPrismElementBase<T>::~GoPrismElementBase()
{
}

template <class T>
void
GoPrismElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoPrismElementBase<T>::shrink()
{
  pool_.shrink();
}
