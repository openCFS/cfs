// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2000/11/07 17:27:31  prkipfer
| introduced external polymorphism for vertices
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoFlowVertexBase<T>::FlowVertexMem> GoFlowVertexBase<T>::pool_;


// constructor
template <class T>
GoFlowVertexBase<T>::GoFlowVertexBase(T d, T mx, T my, T mz, T e)
{
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.density = d;
  data_.momentum[0] = mx;
  data_.momentum[1] = my;
  data_.momentum[2] = mz;
  data_.energy = e;
}

// destructor
template <class T>
GoFlowVertexBase<T>::~GoFlowVertexBase()
{
}

template <class T>
void
GoFlowVertexBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoFlowVertexBase<T>::shrink()
{
  pool_.shrink();
}
