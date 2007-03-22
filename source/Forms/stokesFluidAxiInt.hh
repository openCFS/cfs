// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
  void CalcElementMatrix( Matrix<Double>& elemMat,
                          EntityIterator& ent1, 
                          EntityIterator& ent2 );
  
};


}

#endif // FILE_LAPLACEINT
