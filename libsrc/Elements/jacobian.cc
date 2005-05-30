#include <iostream>
#include <fstream>

#include "jacobian.hh"

namespace CoupledField
{

  template <UInt dim>
  Jacobian<dim>::Jacobian()
  {
    J.Resize(dim,dim);
    Jinv.Resize(dim,dim);
  }

  template<UInt dim>
  void Jacobian<dim>::GetJinvX(Vector<Double> & JinvX)
  {
    JinvX.Resize(dim);

    JinvX[0]=Jinv[0][0];
    JinvX[1]=Jinv[1][0];

    if (dim==3)
      JinvX[2]=Jinv[2][0];
  }

  template<UInt dim>
  void Jacobian<dim>::GetJinvY(Vector<Double> & JinvY)
  {
    JinvY.Resize(dim);

    JinvY[0]=Jinv[0][1];
    JinvY[1]=Jinv[1][1];

    if (dim==3)
      JinvY[2]=Jinv[2][1];
  }

  template<UInt dim>
  void Jacobian<dim>::GetJinvZ(Vector<Double> & JinvZ)
  {
    if (dim==2)
      Error("Function GetJinvZ is implemented only for 3D",__FILE__,__LINE__);
    else {
      JinvZ.Resize(3);

      JinvZ[0]=Jinv[0][2];
      JinvZ[1]=Jinv[1][2];
      JinvZ[2]=Jinv[2][2];
    }

  }
  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Jacobian<2>
#pragma instantiate Jacobian<3>
#endif

} // end of namespace
