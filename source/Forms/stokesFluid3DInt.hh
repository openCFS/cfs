// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );

};


}

#endif // FILE_LAPLACEINT
