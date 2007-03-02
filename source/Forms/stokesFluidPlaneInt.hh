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
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
};


}

#endif // FILE_LAPLACEINT
