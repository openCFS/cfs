// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2002/03/18 10:00:23  prkipfer
| refactored element structure
|
| Revision 1.1  2001/02/13 11:00:37  prkipfer
| added GoBoundaryVertexBase class
|
|
|
+---------------------------------------------------------------------*/

// static memory pool init
template <class T>
GbMemPool<typename GoBoundaryVertexBase<T>::BoundaryVertexMem> GoBoundaryVertexBase<T>::pool_;


// constructor
template <class T>
GoBoundaryVertexBase<T>::GoBoundaryVertexBase(GoVertex<T>* vel)
{
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
  data_.vertex = vel;
}

// destructor
template <class T>
GoBoundaryVertexBase<T>::~GoBoundaryVertexBase()
{
}

template <class T>
void
GoBoundaryVertexBase<T>::poolStatistic() const
{
#ifdef MINIMEM
  pool_.shrink();
#endif
  pool_.getStatistic();
}

template <class T>
void
GoBoundaryVertexBase<T>::shrink()
{
  pool_.shrink();
}
