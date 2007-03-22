// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/***************************************************************************
    File        : Jacobian.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/


#include "Jacobian.hh"

namespace grd {

// Description:
//
Jacobian::Jacobian( )
{
  for (int i = 0; i < 3; i++)
  {
    jc[0][i] = 0.0;
    jc[1][i] = 0.0;
    jc[2][i] = 0.0;
  }
}


// Description:
//
bool
Jacobian::getInverse(Jacobian& sc)
{
  double det,idet;

  // det = getDeterminant();
  det = jc[0][0]*(jc[1][1]*jc[2][2] - jc[1][2]*jc[2][1]) -
        jc[0][1]*(jc[1][0]*jc[2][2] - jc[2][0]*jc[1][2]) +
        jc[0][2]*(jc[1][0]*jc[2][1] - jc[2][0]*jc[1][1]);

  if (fabs(det) > std::numeric_limits<double>::epsilon())
    idet = 1.0 / det;
  else {
      return false;
  }

  sc.jc[0][0] = idet*(jc[1][1]*jc[2][2] - jc[1][2]*jc[2][1]);
  sc.jc[0][1] = idet*(jc[0][2]*jc[2][1] - jc[0][1]*jc[2][2]);
  sc.jc[0][2] = idet*(jc[0][1]*jc[1][2] - jc[0][2]*jc[1][1]);
  sc.jc[1][0] = idet*(jc[1][2]*jc[2][0] - jc[1][0]*jc[2][2]);
  sc.jc[1][1] = idet*(jc[0][0]*jc[2][2] - jc[0][2]*jc[2][0]);
  sc.jc[1][2] = idet*(jc[0][2]*jc[1][0] - jc[0][0]*jc[1][2]);
  sc.jc[2][0] = idet*(jc[1][0]*jc[2][1] - jc[1][1]*jc[2][0]);
  sc.jc[2][1] = idet*(jc[0][1]*jc[2][0] - jc[0][0]*jc[2][1]);
  sc.jc[2][2] = idet*(jc[0][0]*jc[1][1] - jc[0][1]*jc[1][0]);

  return true;
}

} // namespace grd
