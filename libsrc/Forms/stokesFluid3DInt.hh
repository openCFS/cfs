#ifndef FILE_STOKESFLUID3DINT
#define FILE_STOKESFLUID3DINT

#include "stokesFluidInt.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a 3D Stokes operator
class StokesFluid3DInt : public StokesFluidInt
{
public:
  /// Constructor
  StokesFluid3DInt(Double density, Double dynamicViscosity);

  /// 
  virtual ~StokesFluid3DInt();

  /// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);
 
};


}

#endif // FILE_LAPLACEINT
