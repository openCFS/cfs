// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 14:47:45  prkipfer
| introduced edge elements
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoDefaultEdgeElementBase<T>::DefaultEdgeMem> GoDefaultEdgeElementBase<T>::pool_;


template <class T>
GoDefaultEdgeElementBase<T>::GoDefaultEdgeElementBase(int i)
{
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.id = i;
  data_.element = NULL;
  data_.valence = 0;
  data_.flag = NORMAL_E;
  data_.vertices[0] = data_.vertices[1] = NULL;
}

template <class T>
GoDefaultEdgeElementBase<T>::~GoDefaultEdgeElementBase()
{
}

template <class T>
int
GoDefaultEdgeElementBase<T>::findVertex(GoVertex<T> *v) const
{
  if (data_.vertices[0] == v)
    return(0);
  if (data_.vertices[1] == v)
    return(1);
  return(-1);
}


template <class T>
void
GoDefaultEdgeElementBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}


template <class T>
void
GoDefaultEdgeElementBase<T>::shrink()
{
  pool_.shrink();
}
