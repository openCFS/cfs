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
#include "Pool.h"
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Jacobian.h"

namespace grd {

class TetrahedronGeometry {
public:
  // Constructors
  TetrahedronGeometry() {}

  // Destrcutor
  virtual ~TetrahedronGeometry() {}

  // Methods
  inline void getGlobalCoords(Element* t,const double* lc,double* coord);
  inline void getJacobian(Element* t,Jacobian& Jac);
};

//____________________________________________________________
//
//         inline functions
//____________________________________________________________
//

// Description
//
inline void
TetrahedronGeometry::getGlobalCoords(Element* t,const double* lc,double* coord)
{
  double* pos[4];

  double l1 = lc[0];
  double l2 = lc[1];
  double l3 = lc[2];
  double l0 = 1.0 - l1 - l2 - l3;


  pos[0] = (t->getVertex(0))->getPosition();
  pos[1] = (t->getVertex(1))->getPosition();
  pos[2] = (t->getVertex(2))->getPosition();
  pos[3] = (t->getVertex(3))->getPosition();

  coord[0] = l0*pos[0][0] + l1*pos[1][0] + l2*pos[2][0] + l3*pos[3][0];
  coord[1] = l0*pos[0][1] + l1*pos[1][1] + l2*pos[2][1] + l3*pos[3][1];
  coord[2] = l0*pos[0][2] + l1*pos[1][2] + l2*pos[2][2] + l3*pos[3][2];
}

inline void
TetrahedronGeometry::getJacobian(Element* t,Jacobian& Jac)
{
  Vertex* v[4];
  v[0] = t->getVertex(0);
  v[1] = t->getVertex(1);
  v[2] = t->getVertex(2);
  v[3] = t->getVertex(3);

  Jac(0,0) = (*v[1])[0] - (*v[0])[0];
  Jac(0,1) = (*v[2])[0] - (*v[0])[0];
  Jac(0,2) = (*v[3])[0] - (*v[0])[0];

  Jac(1,0) = (*v[1])[1] - (*v[0])[1];
  Jac(1,1) = (*v[2])[1] - (*v[0])[1];
  Jac(1,2) = (*v[3])[1] - (*v[0])[1];

  Jac(2,0) = (*v[1])[2] - (*v[0])[2];
  Jac(2,1) = (*v[2])[2] - (*v[0])[2];
  Jac(2,2) = (*v[3])[2] - (*v[0])[2];
}

} // namespace grd

#endif // TETRAHEDRONGEOMETRY_H