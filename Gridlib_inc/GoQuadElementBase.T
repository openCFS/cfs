// -*- C++ -*- automatisch in C++mode wechseln (emacs)


// static memory pool init
template <class T>
GbMemPool<typename GoQuadElementBase<T>::QuadMem> GoQuadElementBase<T>::pool_;

template <class T>
GoQuadElementBase<T>::GoQuadElementBase(int i, T x, T y, T z)
{
  int k;
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.id = i;
  data_.partition = 0;
  data_.flag = NORMAL_G;
  for (k=0; k<4; k++) data_.vertices[k] = NULL;
  for (k=0; k<4; k++) data_.neighbours[k] = NULL;
  data_.normal[0] = x;
  data_.normal[1] = y;
  data_.normal[2] = z;
}

template <class T>
GoQuadElementBase<T>::~GoQuadElementBase()
{
}

template <class T>
void
GoQuadElementBase<T>::setNeighbour(GoGeometryElement<T>* f1) 
{
  if (f1 == NULL) return;

  int i=-20;

  if (f1->findVertex(data_.vertices[0]) != -1) {
    i+=10;
  }
  if (f1->findVertex(data_.vertices[1]) != -1) {
    i+=11;
  }
  if (f1->findVertex(data_.vertices[2]) != -1) {
    i+=12;
  }
  if (f1->findVertex(data_.vertices[3]) != -1) {
    i+=14;
  }

  switch(i) {
    case 1:
      setNeighbour(0 , f1); break;
    case 3:
      setNeighbour(1 , f1); break;
    case 6:
      setNeighbour(2 , f1); break;
    case 4:
      setNeighbour(3 , f1); break;
    default:
      debugmsg("setNeighbour not correct");
  }
}

//
//  finde Knoten in Dreieck
// 

template <class T>
int
GoQuadElementBase<T>::findVertex(GoVertex<T> *v) const
{
  if (data_.vertices[0] == v)
    return(0);
  if (data_.vertices[1] == v)
    return(1);
  if (data_.vertices[2] == v)
    return(2);
  if (data_.vertices[3] == v)
    return(3);
  return(-1);
}

//
//  finde Nachbarn in Dreieck
// 

template <class T>
int
GoQuadElementBase<T>::findNeighbour(GoGeometryElement<T> *f) const
{
  if (data_.neighbours[0] == f)
    return(0);
  if (data_.neighbours[1] == f)
    return(1);
  if (data_.neighbours[2] == f)
    return(2);
  if (data_.neighbours[3] == f)
    return(3);
  return(-1);
}


////////////////////////////////////////////////////////////////
//  face3::nextFace()
//
//  gegen"uberliegendes Dreieck zu gegebenem Knoten
//

template <class T>
GoGeometryElement<T> *
GoQuadElementBase<T>::otherFace(GoVertex<T> *v) const
{
  if (v == data_.vertices[0]) return(data_.neighbours[0]);
  if (v == data_.vertices[1]) return(data_.neighbours[1]);
  if (v == data_.vertices[2]) return(data_.neighbours[3]);
  return(data_.neighbours[3]);
}

////////////////////////////////////////////////////////////////
//  face3::nextVertex()
//
//  gegen"uberliegendes Knoten zu gegebenem Dreieck
// 

template <class T>
GoVertex<T> *
GoQuadElementBase<T>::otherVertex(GoGeometryElement<T> *f) const
{
  if (f == data_.neighbours[0]) return(data_.vertices[0]);
  if (f == data_.neighbours[1]) return(data_.vertices[1]);
  if (f == data_.neighbours[2]) return(data_.vertices[2]);
  return(data_.vertices[3]);
}

template <class T>
void
GoQuadElementBase<T>::computeNormal()
{
  T x0,y0,z0, x1,y1, z1, x2,y2,z2, xn,yn,zn;
  T n;

//  debugmsg("computing face normal");
  data_.vertices[0]->getPosition(x0,y0,z0);
  data_.vertices[1]->getPosition(x1,y1,z1);
  data_.vertices[2]->getPosition(x2,y2,z2);
  
  xn = ((y1-y0)*(z2-z0) - (y2-y0)*(z1-z0));
  yn = ((z1-z0)*(x2-x0) - (z2-z0)*(x1-x0));
  zn = ((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
  
  n = GbMath<T>::Sqrt(xn*xn + yn*yn + zn*zn);
  
  data_.normal[0] = xn/n;
  data_.normal[1] = yn/n;
  data_.normal[2] = zn/n;

}

template <class T>
void
GoQuadElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoQuadElementBase<T>::shrink()
{
  pool_.shrink();
}
