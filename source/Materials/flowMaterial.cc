// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "flowMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  FlowMaterial::FlowMaterial() : BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");

    materialDatabaseName_ = "FluidFlow";
    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( DYNAMIC_VISCOSITY );
  }

  FlowMaterial::~FlowMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }

  
  void FlowMaterial::SetScalar( Double param, MaterialType matType, 
				DataType dataType ) {

    ENTER_FCN( "FlowMaterial::SetScalar" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val;
      if ( dataType == REAL ) {
	val = Complex ( param, 0.0 );
      }
      else if (dataType == IMAG ) {
	val = Complex ( 0.0, param );
	isComplex_.insert( matType );
      }
      else {
	std::string msg = "SetScalar-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
      
      scalarParams_[matType] = val;
    }
  }


  void FlowMaterial::GetScalar( Double& param, MaterialType matType, 
				DataType dataType ) const {

    ENTER_FCN( "FlowMaterial::GetScalar" );

    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == REAL ) {
	param = val.real();
      }
      else if ( dataType == IMAG ) {
	param = val.imag();
      }
      else {
	std::string msg = "GetScalar-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }

}
