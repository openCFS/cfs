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
GoPyramidElementBase<T>::GoPyramidElementBase(int i, T x, T y, T z)
{
  int k;
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.id = i;
  data_.partition = 0;
  data_.flag = NORMAL_G;
  for (k=0; k<5; k++) data_.vertices[k] = NULL;
  for (k=0; k<5; k++) data_.neighbours[k] = NULL;
}

template <class T>
GoPyramidElementBase<T>::~GoPyramidElementBase()
{
}

template <class T>
void
GoPyramidElementBase<T>::setNeighbour(GoGeometryElement<T>* f1) 
{
  debugmsg("pyra:setNeighbour");

  if (f1 == NULL) return;

  int i=0;

  if (f1->findVertex(data_.vertices[0]) == -1) {
    i+=1;
  }
  if (f1->findVertex(data_.vertices[1]) == -1) {
    i+=2;
  }
  if (f1->findVertex(data_.vertices[2]) == -1) {
    i+=4;
  }
  if (f1->findVertex(data_.vertices[3]) == -1) {
    i+=8;
  }
  if (f1->findVertex(data_.vertices[4]) == -1) {
    i+=16;
  }


  if (i==16)
    setNeighbour(4, f1);
  else if (i==12)
    setNeighbour(0, f1);
  else if (i==9)
    setNeighbour(1, f1);
  else if (i==3)
    setNeighbour(2, f1);
  else if (i==6)
    setNeighbour(3, f1);
  else
    debugmsg("Pyramid neighbour setup failed");

  std::cerr<< i <<std::endl;
    
}

//
//  finde Knoten in Dreieck
// 

template <class T>
int
GoPyramidElementBase<T>::findVertex(GoVertex<T> *v) const
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
  return(-1);
}

//
//  finde Nachbarn in Dreieck
// 

template <class T>
int
GoPyramidElementBase<T>::findNeighbour(GoGeometryElement<T> *f) const
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
GoPyramidElementBase<T>::otherFace(GoVertex<T> *v) const
{
  if (v == data_.vertices[0]) return(data_.neighbours[0]);
  if (v == data_.vertices[1]) return(data_.neighbours[1]);
  if (v == data_.vertices[2]) return(data_.neighbours[2]);
  return(data_.neighbours[3]);
}

////////////////////////////////////////////////////////////////
//  face3::nextVertex()
//
//  gegen"uberliegendes Knoten zu gegebenem Dreieck
// 

template <class T>
GoVertex<T> *
GoPyramidElementBase<T>::otherVertex(GoGeometryElement<T> *f) const
{
  if (f == data_.neighbours[0]) return(data_.vertices[0]);
  if (f == data_.neighbours[1]) return(data_.vertices[1]);
  if (f == data_.neighbours[2]) return(data_.vertices[2]);
  return(data_.vertices[3]);
}

template <class T>
void
GoPyramidElementBase<T>::computeNormal()
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
