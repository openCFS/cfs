// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.10  2002/03/18 10:00:26  prkipfer
| refactored element structure
|
| Revision 1.9  2001/09/12 09:20:04  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.8  2001/01/02 15:16:47  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.7  2000/12/12 10:05:47  prkipfer
| added new element and vertex flags
|
| Revision 1.6  2000/11/07 17:27:39  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.5  2000/09/28 14:15:13  uflabsik
| errors corrected
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
| extended IO subsystem, created class registry, added tetrahedron
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoTetrahedronElementBase<T>::TetraMem> GoTetrahedronElementBase<T>::pool_;

template <class T>
GoTetrahedronElementBase<T>::GoTetrahedronElementBase()
{
  memset(data_.vertices,0,4*sizeof(GoVertex<T>*));
  memset(data_.edges,0,6*sizeof(GoEdge<T>*));
  memset(data_.neighbours,0,4*sizeof(GoGeometryElement<T>*));
  memset(data_.children,0,8*sizeof(GoGeometryElement<T>*));
}

template <class T>
GoTetrahedronElementBase<T>::~GoTetrahedronElementBase()
{
  int i;

  for(i = 0; i < 4; i++) {
    if (data_.vertices[i])
      data_.vertices[i]->removeElement();
  }
  
  for(i = 0; i < 6; i++) {
    if (data_.edges[i])
      data_.edges[i]->removeElement();
  }
}

template <class T>
void
GoTetrahedronElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoTetrahedronElementBase<T>::shrink()
{
  pool_.shrink();
}
