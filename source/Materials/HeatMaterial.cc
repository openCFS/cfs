// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "HeatMaterial.hh"

#include "Utils/LinInterpolate.hh"


namespace CoupledField
{

// ***********************
//   Default Constructor
// ***********************
HeatMaterial::HeatMaterial(MathParser* mp,
                           CoordSystem * defaultCoosy) 
: BaseMaterial(THERMIC, mp, defaultCoosy)
{
  //set the allowed material parameters
  isAllowed_.insert( DENSITY );
  isAllowed_.insert( HEAT_CONDUCTIVITY_SCALAR );
  isAllowed_.insert( HEAT_CONDUCTIVITY_TENSOR );
  isAllowed_.insert( HEAT_CONDUCTIVITY_1 );
  isAllowed_.insert( HEAT_CONDUCTIVITY_2 );
  isAllowed_.insert( HEAT_CONDUCTIVITY_3 );
  isAllowed_.insert( HEAT_CAPACITY );
  isAllowed_.insert( HEAT_REF_TEMPERATURE );
  isAllowed_.insert( NONLIN_DEPENDENCY );
}

HeatMaterial::~HeatMaterial() {
}

void HeatMaterial::Finalize() {
  MaterialType orthoCond[3] = {
      HEAT_CONDUCTIVITY_1,
      HEAT_CONDUCTIVITY_2,
      HEAT_CONDUCTIVITY_3
  };
  CalcFull3x3Tensor(HEAT_CONDUCTIVITY_SCALAR, orthoCond, HEAT_CONDUCTIVITY_TENSOR);
}

}
