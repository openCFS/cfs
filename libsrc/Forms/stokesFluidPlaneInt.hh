#ifndef FILE_STOKESFLUIDPLANEINT
#define FILE_STOKESFLUIDPLANEINT

#include "stokesFluidInt.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a Plane Stokes operator
class StokesFluidPlaneInt : public StokesFluidInt
{
public:
  /// Constructor
  StokesFluidPlaneInt(Double density, Double dynamicViscosity);

  /// 
  virtual ~StokesFluidPlaneInt();

  /// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);
 
};


}

#endif // FILE_LAPLACEINT
