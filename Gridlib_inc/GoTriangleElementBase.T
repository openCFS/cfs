// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.8  2002/03/18 10:00:27  prkipfer
| refactored element structure
|
| Revision 1.7  2001/01/02 15:16:47  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.6  2000/12/12 10:05:48  prkipfer
| added new element and vertex flags
|
| Revision 1.5  2000/11/07 17:27:39  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.4  2000/09/03 14:05:43  uflabsik
| functionality added
|
| Revision 1.3  2000/08/28 09:25:22  prkipfer
| added localized subdivision mechanism
|
| Revision 1.2  2000/06/16 13:33:50  prkipfer
| further hierarchy fixes
|
| Revision 1.1  2000/06/14 15:39:14  prkipfer
| improved base classes and added funcstruct processing for mesh
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoTriangleElementBase<T>::TriangleMem> GoTriangleElementBase<T>::pool_;

template <class T>
GoTriangleElementBase<T>::GoTriangleElementBase()
{
  memset(data_.vertices,0,3*sizeof(GoVertex<T>*));
  memset(data_.edges,0,3*sizeof(GoEdge<T>*));
  memset(data_.neighbours,0,3*sizeof(GoGeometryElement<T>*));
  memset(data_.children,0,4*sizeof(GoGeometryElement<T>*));
}

template <class T>
GoTriangleElementBase<T>::~GoTriangleElementBase()
{
}

template <class T>
void
GoTriangleElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoTriangleElementBase<T>::shrink()
{
  pool_.shrink();
}
