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
GoOctahedronElementBase<T>::GoOctahedronElementBase(int i, T x, T y, T z)
{
  int k;
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.id = i;
  data_.partition = 0;
  data_.flag = NORMAL_G;
  for (k=0; k<6; k++) data_.vertices[k] = NULL;
  for (k=0; k<8; k++) data_.neighbours[k] = NULL;
}

template <class T>
GoOctahedronElementBase<T>::~GoOctahedronElementBase()
{
}

template <class T>
void
GoOctahedronElementBase<T>::setNeighbour(GoGeometryElement<T>* f1) 
{
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
  if (f1->findVertex(data_.vertices[5]) == -1) {
    i+=32;
  }

  if (i==44)  // 2,3,5
    setNeighbour(0 , f1);
  else if (i==41)    // 0,3,5
    setNeighbour(1 , f1);
  else if (i==35)    // 0,1,5
    setNeighbour(2, f1);
  else if (i==38)    // 1,2,5
    setNeighbour(3, f1);
  else if (i==28)  // 2,3,4
    setNeighbour(4 , f1);
  else if (i==25)    // 0,3,4
    setNeighbour(5 , f1);
  else if (i==19)    // 0,1,4
    setNeighbour(6, f1);
  else if (i==22)    // 1,2,4
    setNeighbour(7, f1);
  else
    debugmsg("Octahedron neighbour setup failed");
    
}

//
//  finde Knoten in Dreieck
// 

template <class T>
int
GoOctahedronElementBase<T>::findVertex(GoVertex<T> *v) const
{
  int i;
  for (i=0; i<6; i++)
    if (data_.vertices[i] == v)
      return(i); 

  return(-1);
}

//
//  finde Nachbarn in Dreieck
// 

template <class T>
int
GoOctahedronElementBase<T>::findNeighbour(GoGeometryElement<T> *f) const
{
  int i;
  for (i=0; i<8; i++)
    if (data_.neighbours[i] == f)
      return(i);

  return(-1);
}



////////////////////////////////////////////////////////////////
//  face3::nextFace()
//
//  gegen"uberliegendes Dreieck zu gegebenem Knoten
//

template <class T>
GoGeometryElement<T> *
GoOctahedronElementBase<T>::otherFace(GoVertex<T> *v) const
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
GoOctahedronElementBase<T>::otherVertex(GoGeometryElement<T> *f) const
{
  return(NULL);
}

template <class T>
void
GoOctahedronElementBase<T>::computeNormal()
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
