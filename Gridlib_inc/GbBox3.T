// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 14:59:03  prkipfer
| introduced new classes
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbBox3<T>::GbBox3()
  : center_()
{
  axis_[0]=GbVec3<T>::UNIT_X;
  axis_[1]=GbVec3<T>::UNIT_Y;
  axis_[2]=GbVec3<T>::UNIT_Z;
  extent_[0]=extent_[1]=extent_[2]=0;
}

template <class T>
void 
GbBox3<T>::computeVertices (GbVec3<T> akVertex[8]) const
{
  GbVec3<T> akEAxis[3] =
  {
    extent_[0]*axis_[0],
    extent_[1]*axis_[1],
    extent_[2]*axis_[2]
  };

  akVertex[0] = center_ - akEAxis[0] - akEAxis[1] - akEAxis[2];
  akVertex[1] = center_ + akEAxis[0] - akEAxis[1] - akEAxis[2];
  akVertex[2] = center_ + akEAxis[0] + akEAxis[1] - akEAxis[2];
  akVertex[3] = center_ - akEAxis[0] + akEAxis[1] - akEAxis[2];
  akVertex[4] = center_ - akEAxis[0] - akEAxis[1] + akEAxis[2];
  akVertex[5] = center_ + akEAxis[0] - akEAxis[1] + akEAxis[2];
  akVertex[6] = center_ + akEAxis[0] + akEAxis[1] + akEAxis[2];
  akVertex[7] = center_ - akEAxis[0] + akEAxis[1] + akEAxis[2];
}
