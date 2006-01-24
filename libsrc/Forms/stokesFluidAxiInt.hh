#ifndef FILE_STOKESFLUIDAXIINT
#define FILE_STOKESFLUIDAXIINT

#include "stokesFluidInt.hh"

namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a Axi Stokes operator
class StokesFluidAxiInt : public StokesFluidInt
{
public:
  /// Constructor
  StokesFluidAxiInt(Double density, Double dynamicViscosity);

  /// 
  virtual ~StokesFluidAxiInt();

  /// Calculation of stiffmess matrix
  virtual void CalcElementMatrix(Matrix<Double> & ptCoord, Matrix<Double> & elemMat);
 
};


}

#endif // FILE_LAPLACEINT
