// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/18 11:13:37  prkipfer
| made operator/ inline
|
| Revision 1.1  2001/02/13 13:29:28  prkipfer
| introduced basic rendering tools
|
|
|
+---------------------------------------------------------------------*/

template <class T>
void 
GrColor<T>::clamp ()
{
  if ( r_ < 0.0 ) r_ = 0.0;
  else if ( r_ > 1.0 ) r_ = 1.0;

  if ( g_ < 0.0 ) g_ = 0.0;
  else if ( g_ > 1.0 ) g_ = 1.0;

  if ( b_ < 0.0 ) b_ = 0.0;
  else if ( b_ > 1.0 ) b_ = 1.0;
}

template <class T>
void 
GrColor<T>::scaleByMax ()
{
  T fMax = r_;
  if ( g_ > fMax ) fMax = g_;
  if ( b_ > fMax ) fMax = b_;

  if ( fMax > 1.0 ) {
    T fInvMax = 1.0/fMax;
    r_ *= fInvMax;
    g_ *= fInvMax;
    b_ *= fInvMax;
  }
}

