// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
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
GoBoundaryVertexBase<T>::GoBoundaryVertexBase()
{
//  debugmsg("constructor size is "<<sizeof(this)<<" bytes");
//    data_.element = f;
}

// destructor
template <class T>
GoBoundaryVertexBase<T>::~GoBoundaryVertexBase()
{
}

// compute normal at this vertex by averaging the normals
// of neighboring boundary face objects
//  template <class T>
//  void
//  GoBoundaryVertexBase<T>::computeNormal(std::vector<GoGeometryElement<T> *>& S)
//  {
//    GbVec3<T> n(0.0, 0.0, 0.0);

//    if (S.size()==0) {
//      if (getElement()!=NULL) {
//        setNormal(getElement()->getNormal());
//        warningmsg("no neighbors ? possible mesh consistency problem ?");
//      }
//      else {
//        warningmsg("no face associated ? possible mesh setup problem ?");
//      }
//      return;
//    }
  
//    for (std::vector<GoGeometryElement<T> *>::iterator sIter = S.begin(); 
//         sIter != S.end(); 
//         ++sIter)

//      n += (*sIter)->getNormal();

//    n.normalize();
//    setNormal(n);
//  }

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
