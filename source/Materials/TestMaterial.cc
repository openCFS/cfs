// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <map>
#include <set>
#include <string>
#include <utility>

#include "Materials/BaseMaterial.hh"
#include "TestMaterial.hh"
#include "limits.h"
#include "math.h"
#include "stdlib.h"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  TestMaterial::TestMaterial(MathParser* mp,
                                     CoordSystem * defaultCoosy) 
  : BaseMaterial(mp, defaultCoosy) {

    materialDatabaseName_ = "Test";

    //set the allowed material parameters
    isAllowed_.insert( TEST_ALPHA );
    isAllowed_.insert( TEST_BETA );
  }

  TestMaterial::~TestMaterial() {
  }

  void TestMaterial::SetScalar( Double param, MaterialType matType, 
      Global::ComplexPart dataType ) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val;
      if ( dataType == Global::REAL ) {
        val = Complex ( param, 0.0 );
      }
      else if (dataType == Global::IMAG ) {
        val = Complex ( 0.0, param );
        isComplex_.insert( matType );
      }
      else {
        std::string msg = "SetScalar-Double";
        dataTypeNotAllowed4SetGet ( dataType, msg );
      }

      scalarParams_[matType] = val;
    }
  }

  void TestMaterial::GetScalar( Double& param, MaterialType matType, 
				    Global::ComplexPart dataType ) const {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == Global::REAL ) {
	param = val.real();
      }
      else if ( dataType == Global::IMAG ) {
	param = val.imag();
      }
      else {
	std::string msg = "GetScalar-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }
}
