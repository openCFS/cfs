// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// - C++ -
/***************************************************************************
    File        : TetrahedronGeometry.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Wed Sep 11 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef TETRAHEDRONGEOMETRY_H
#define TETRAHEDRONGEOMETRY_H

// Grid
#include "Pool.hh"
#include "Vertex.hh"
#include "Edge.hh"
#include "Element.hh"
#include "Jacobian.hh"

namespace grd {

class TetrahedronGeometry {
public:
  // Constructors
  TetrahedronGeometry() {}

  // Destrcutor
  virtual ~TetrahedronGeometry() {}

  // Methods
  inline void getGlobalCoords(const Element* t,const double* lc,double* coord);
  inline void getJacobian(const Element* t,Jacobian& Jac);
  inline double getVolume(const Element* t);
};

//____________________________________________________________
//
//         inline functions
//____________________________________________________________
//

// Description
//
inline void
TetrahedronGeometry::getGlobalCoords(const Element* t,const double* lc,double* coord)
{
  double* pos[4];

  double l0 = lc[0];
  double l1 = lc[1];
  double l2 = lc[2];
  double l3 = 1.0 - l0 - l1 - l2;


  pos[0] = (t->getVertex(0))->getPosition();
  pos[1] = (t->getVertex(1))->getPosition();
  pos[2] = (t->getVertex(2))->getPosition();
  pos[3] = (t->getVertex(3))->getPosition();

  coord[0] = l0*pos[0][0] + l1*pos[1][0] + l2*pos[2][0] + l3*pos[3][0];
  coord[1] = l0*pos[0][1] + l1*pos[1][1] + l2*pos[2][1] + l3*pos[3][1];
  coord[2] = l0*pos[0][2] + l1*pos[1][2] + l2*pos[2][2] + l3*pos[3][2];
}

inline void
TetrahedronGeometry::getJacobian(const Element* t,Jacobian& Jac)
{
  double* pos[4];
  pos[0] = (t->getVertex(0))->getPosition();
  pos[1] = (t->getVertex(1))->getPosition();
  pos[2] = (t->getVertex(2))->getPosition();
  pos[3] = (t->getVertex(3))->getPosition();

  Jac(0,0) = pos[0][0] - pos[3][0];
  Jac(0,1) = pos[0][1] - pos[3][1];
  Jac(0,2) = pos[0][2] - pos[3][2];

  Jac(1,0) = pos[1][0] - pos[3][0];
  Jac(1,1) = pos[1][1] - pos[3][1];
  Jac(1,2) = pos[1][2] - pos[3][2];

  Jac(2,0) = pos[2][0] - pos[3][0];
  Jac(2,1) = pos[2][1] - pos[3][1];
  Jac(2,2) = pos[2][2] - pos[3][2];
}

inline double
TetrahedronGeometry::getVolume(const Element* t)
{
  Jacobian jac;
  getJacobian(t,jac);
  double vol = jac.getDeterminant();
  vol *= 1.0/6.0;
  return vol;
}

} // namespace grd

#endif // TETRAHEDRONGEOMETRY_H
