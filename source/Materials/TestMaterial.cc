// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "TestMaterial.hh"

namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  TestMaterial::TestMaterial(MathParser* mp,
                                     CoordSystem * defaultCoosy) 
  : BaseMaterial(TESTMAT, mp, defaultCoosy)
  {
    //set the allowed material parameters
    isAllowed_.insert( TEST_ALPHA );
    isAllowed_.insert( TEST_BETA );
  }

  TestMaterial::~TestMaterial() {
  }

}
