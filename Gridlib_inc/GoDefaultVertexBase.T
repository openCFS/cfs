// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.5  2002/03/18 10:00:24  prkipfer
| refactored element structure
|
| Revision 1.4  2001/09/12 09:20:03  prkipfer
| introduced adaptive tet subdivision
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
GoDefaultVertexBase<T>::GoDefaultVertexBase(int i, T x, T y, T z, T nx, T ny, T nz)
{
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.flag = NORMAL_V;
  data_.id = i;
  data_.valence = 0;
  data_.face = NULL;
  data_.pos[0] = x; 
  data_.pos[1] = y;
  data_.pos[2] = z;
  data_.normal[0] = nx;
  data_.normal[1] = ny;
  data_.normal[2] = nz;
}

// destructor
template <class T>
GoDefaultVertexBase<T>::~GoDefaultVertexBase()
{
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
