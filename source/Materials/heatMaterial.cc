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

#include "heatMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  HeatMaterial::HeatMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Heat";

    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( HEAT_CONDUCTIVITY );
    isAllowed_.insert( HEAT_CAPACITY );
  }

  HeatMaterial::~HeatMaterial() {

    materialDatabaseName_ = "Heat";

  }

  void HeatMaterial::SetScalar( Double param, MaterialType matType, 
				DataType dataType ) {


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
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
      
      scalarParams_[matType] = val;
    }

  }

  void HeatMaterial::SetScalar( Complex param, MaterialType matType, 
				DataType dataType ) {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val = pos->second;
      if ( dataType == REAL ) {
	//	param = val.real();
      }
      else if ( dataType == IMAG ) {
	//	param = val.imag();
	isComplex_.insert( matType );
      }
      else if ( dataType == COMPLEX ) {
	param = val;
	isComplex_.insert( matType );
      }
    }
  }



  void HeatMaterial::GetScalar( Double& param, MaterialType matType, 
				DataType dataType ) const {


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


  void HeatMaterial::GetScalar( Complex& param, MaterialType matType, 
					   DataType dataType ) const {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == REAL ) {
	Complex valReal = Complex (val.real(), 0.0);
	param = valReal;
      }
      else if ( dataType == IMAG ) {
	Complex valImag = Complex (0.0, val.imag());
	param = valImag;
      }
      else if ( dataType == COMPLEX ) {
	param = val;
      }
    }
  }

}
