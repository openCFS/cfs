// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/15 09:04:43  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:33  prkipfer
| introduced basic rendering tools
|
|
|
+---------------------------------------------------------------------*/

template <class T>
int GrLight<T>::nextId_ = 0;

template <class T>
GrLight<T>::GrLight()
  : ambient_(GrColor<T>::BLACK)
  , diffuse_(GrColor<T>::BLACK)
  , specular_(GrColor<T>::BLACK)
  , constant_(1.0)
  , linear_(0.0)
  , quadratic_(0.0)
  , intensity_(1.0)
  , attenuate_(false)
  , on_(true)
{
  id_ = nextId_;
  nextId_++;
  debugmsg("created light "<<id_);
}

