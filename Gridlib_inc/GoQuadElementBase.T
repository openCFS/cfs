// -*- C++ -*- automatisch in C++mode wechseln (emacs)


// static memory pool init
template <class T>
GbMemPool<typename GoQuadElementBase<T>::QuadMem> GoQuadElementBase<T>::pool_;

template <class T>
GoQuadElementBase<T>::GoQuadElementBase()
{
  memset(data_.vertices,0,4*sizeof(GoVertex<T>*));
  memset(data_.edges,0,4*sizeof(GoEdge<T>*));
  memset(data_.neighbours,0,4*sizeof(GoGeometryElement<T>*));
  memset(data_.children,0,4*sizeof(GoGeometryElement<T>*));
}

template <class T>
GoQuadElementBase<T>::~GoQuadElementBase()
{
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
