// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
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
GoPrismElementBase<T>::GoPrismElementBase(int i, T x, T y, T z)
{
  int k;
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.id = i;
  data_.partition = 0;
  data_.flag = NORMAL_G;
  for (k=0; k<6; k++) data_.vertices[k] = NULL;
  for (k=0; k<5; k++) data_.neighbours[k] = NULL;
}

template <class T>
GoPrismElementBase<T>::~GoPrismElementBase()
{
}

template <class T>
void
GoPrismElementBase<T>::setNeighbour(GoGeometryElement<T>* f1) 
{
  if (f1 == NULL) return;

  int i=-20;

  if (f1->findVertex(data_.vertices[0]) == -1) {
    i+=10;
  }
  if (f1->findVertex(data_.vertices[1]) == -1) {
    i+=11;
  }
  if (f1->findVertex(data_.vertices[2]) == -1) {
    i+=12;
  }
  if (f1->findVertex(data_.vertices[3]) == -1) {
    i+=14;
  }
  if (f1->findVertex(data_.vertices[4]) == -1) {
    i+=18;
  }
  if (f1->findVertex(data_.vertices[5]) == -1) {
    i+=26;
  }

  if (i==38)
    setNeighbour(3 , f1);
  else if (i==13)
    setNeighbour(4 , f1);
  else if (i==4)
    setNeighbour(0, f1);
  else if (i==9)
    setNeighbour(1, f1);
  else if (i==18)
    setNeighbour(2, f1);
  else
    debugmsg("Prism neighbour setup failed");
    
}

//
//  finde Knoten in Dreieck
// 

template <class T>
int
GoPrismElementBase<T>::findVertex(GoVertex<T> *v) const
{
  if (data_.vertices[0] == v)
    return(0);
  if (data_.vertices[1] == v)
    return(1);
  if (data_.vertices[2] == v)
    return(2);
  if (data_.vertices[3] == v)
    return(3);
  if (data_.vertices[4] == v)
    return(4);
  if (data_.vertices[5] == v)
    return(5);
  return(-1);
}

//
//  finde Nachbarn in Dreieck
// 

template <class T>
int
GoPrismElementBase<T>::findNeighbour(GoGeometryElement<T> *f) const
{
  if (data_.neighbours[0] == f)
    return(0);
  if (data_.neighbours[1] == f)
    return(1);
  if (data_.neighbours[2] == f)
    return(2);
  if (data_.neighbours[3] == f)
    return(3);
  if (data_.neighbours[4] == f)
    return(4);
  return(-1);
}



////////////////////////////////////////////////////////////////
//  face3::nextFace()
//
//  gegen"uberliegendes Dreieck zu gegebenem Knoten
//

template <class T>
GoGeometryElement<T> *
GoPrismElementBase<T>::otherFace(GoVertex<T> *v) const
{
  return(NULL);
}

////////////////////////////////////////////////////////////////
//  face3::nextVertex()
//
//  gegen"uberliegendes Knoten zu gegebenem Dreieck
// 

template <class T>
GoVertex<T> *
GoPrismElementBase<T>::otherVertex(GoGeometryElement<T> *f) const
{
  return(NULL);
}

template <class T>
void
GoPrismElementBase<T>::computeNormal()
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
