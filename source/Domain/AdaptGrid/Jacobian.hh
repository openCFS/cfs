// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// - C++ -
/***************************************************************************
    File        : Jacobian.h
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

//  -*- C++ -*-
/* ______________________________________________________________________
 *
 *  Class: Jacobian
 *
 *  File: Jacobian.h
 *
 *  Project: Adaptive Finite Element
 *
 *  Author(s): Roberto Grosso
 *
 *  Copyright (c) 1996 by University of Erlangen, IMMD IX.
 *  Am Weichselgarten 9, 91058 Erlangen, Germany.
 *  All rights reserved
 * ______________________________________________________________________
 */

#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <cmath>
#include <limits>

namespace grd {

class Jacobian {
public:
  // Constructors
  Jacobian();

  // Destructor
  virtual ~Jacobian() { }

  // Member Functions
  const double& operator() (int i,int j) const { return jc[i][j]; }
  double&       operator() (int i,int j)       { return jc[i][j]; }

  double getDeterminant()
  {
    double det = jc[0][0]*(jc[1][1]*jc[2][2] - jc[1][2]*jc[2][1]) -
                 jc[0][1]*(jc[1][0]*jc[2][2] - jc[2][0]*jc[1][2]) +
                 jc[0][2]*(jc[1][0]*jc[2][1] - jc[2][0]*jc[1][1]);

    return det;
  }

  bool getInverse(Jacobian& invJac);

private:
  double jc[3][3];
};


} // namespace grd

#endif // JACOBIAN_H
