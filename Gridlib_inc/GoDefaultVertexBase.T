// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/02/13 10:59:58  prkipfer
| more precise warning messages
|
| Revision 1.2  2000/12/12 10:05:46  prkipfer
| added new element and vertex flags
|
| Revision 1.1  2000/11/07 17:27:30  prkipfer
| introduced external polymorphism for vertices
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoDefaultVertexBase<T>::DefaultVertexMem> GoDefaultVertexBase<T>::pool_;


// constructor
template <class T>
GoDefaultVertexBase<T>::GoDefaultVertexBase(int i, T x, T y, T z)
{
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.flag = NORMAL_V;
  data_.id = i;
  data_.valence = 0;
  data_.face = NULL;
  data_.pos[0] = x; 
  data_.pos[1] = y;
  data_.pos[2] = z;
  data_.normal[0] = data_.normal[1] = data_.normal[2] = 0;
}

// destructor
template <class T>
GoDefaultVertexBase<T>::~GoDefaultVertexBase()
{
}

// compute normal at this vertex by averaging the normals
// of neighboring face objects
template <class T>
void
GoDefaultVertexBase<T>::computeNormal(std::vector<GoGeometryElement<T> *>& S)
{
  GbVec3<T> n(0.0, 0.0, 0.0);

  if (S.size()==0) {
    if (getElement()!=NULL) {
      setNormal(getElement()->getNormal());
      warningmsg("no neighbors ? possible mesh consistency problem or mesh not closed ?");
    }
    else {
      warningmsg("no face associated ? possible mesh setup problem ?");
    }
    return;
  }
  
  for (std::vector<GoGeometryElement<T> *>::iterator sIter = S.begin(); 
       sIter != S.end(); 
       ++sIter)

    n += (*sIter)->getNormal();

  n.normalize();
  setNormal(n);
}

template <class T>
void
GoDefaultVertexBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoDefaultVertexBase<T>::shrink()
{
  pool_.shrink();
}
