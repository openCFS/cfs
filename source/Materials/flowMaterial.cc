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

  void FlowMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "FlowMaterial::SetScalar" );
    Error("SetScalar for 'String' not implemented",__FILE__,__LINE__);
  }


  void FlowMaterial::SetScalar( Double& param, const MaterialType& matType, 
				const DataType& dataType ) {

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


  void FlowMaterial::SetScalar( Complex& param, const MaterialType& matType, 
				const DataType& dataType ) {

    ENTER_FCN( "FlowMaterial::SetScalar" );

    dataTypeNotAllowed( dataType, matType );
  }


  void FlowMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
					   const DataType& dataType ) {
    
    ENTER_FCN( "FlowMaterial::SetTensor" );
    dataTypeNotAllowed( dataType, matType );
  }

  void FlowMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
					   const DataType& dataType ) {
    
    ENTER_FCN( "FlowMaterial::SetTensor" );
    dataTypeNotAllowed( dataType, matType );    
  }


  void FlowMaterial::GetScalar( Double& param, const MaterialType& matType, 
				const DataType& dataType ) const {

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

  void FlowMaterial::GetScalar( Complex& param, const MaterialType& matType, 
					   const DataType& dataType ) const {

    ENTER_FCN( "FlowMaterial::GetScalar" );
    dataTypeNotAllowed( dataType, matType );    
  }

}
