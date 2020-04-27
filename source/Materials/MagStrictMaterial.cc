// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <Materials/MagStrictMaterial.hh>
#include "Domain/CoefFunction/CoefFunction.hh"

namespace CoupledField {

// ***********************
//   Default Constructor
// ***********************
MagStrictMaterial::MagStrictMaterial(MathParser* mp,
                               CoordSystem * defaultCoosy) 
: BaseMaterial(MAGNETOSTRICTIVE, mp, defaultCoosy)
{
  //set the allowed material parameters
  isAllowed_.insert( MAGNETOSTRICTION_TENSOR_h );
}

MagStrictMaterial::~MagStrictMaterial() {

}

}
