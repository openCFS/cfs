// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
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
| Revision 1.2  2000/09/28 14:15:12  uflabsik
| errors corrected
|
| Revision 1.1  2000/09/26 09:29:35  uflabsik
| initial version
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoOctahedronElementBase<T>::OctahedronMem> GoOctahedronElementBase<T>::pool_;

template <class T>
GoOctahedronElementBase<T>::GoOctahedronElementBase()
{
  memset(data_.vertices,0,6*sizeof(GoVertex<T>*));
  memset(data_.edges,0,12*sizeof(GoEdge<T>*));
  memset(data_.neighbours,0,8*sizeof(GoGeometryElement<T>*));
  memset(data_.children,0,10*sizeof(GoGeometryElement<T>*));
}

template <class T>
GoOctahedronElementBase<T>::~GoOctahedronElementBase()
{
}

template <class T>
void
GoOctahedronElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoOctahedronElementBase<T>::shrink()
{
  pool_.shrink();
}
