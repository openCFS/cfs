// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2000/11/07 17:27:35  prkipfer
| introduced external polymorphism for vertices
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoScalarVertexBase<T>::ScalarVertexMem> GoScalarVertexBase<T>::pool_;


// constructor
template <class T>
GoScalarVertexBase<T>::GoScalarVertexBase(T v)
{
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.value = v;
}

// destructor
template <class T>
GoScalarVertexBase<T>::~GoScalarVertexBase()
{
}

template <class T>
void
GoScalarVertexBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoScalarVertexBase<T>::shrink()
{
  pool_.shrink();
}
