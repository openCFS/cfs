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

    ENTER_FCN("BaseMaterial::BaseMaterial");
    materialDatabaseName_ = "Heat";

    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( HEAT_CONDUCTIVITY );
    isAllowed_.insert( HEAT_CAPACITY );
  }

  HeatMaterial::~HeatMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");
    materialDatabaseName_ = "Heat";

  }

  void  HeatMaterial::SetScalar( Integer& param, const MaterialType& matType) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );
    Error("SetScalar for 'Integer' not implemented",__FILE__,__LINE__);
  }

  void  HeatMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );
    Error("SetScalar for 'String' not implemented",__FILE__,__LINE__);
  }


  void HeatMaterial::SetScalar( Double& param, const MaterialType& matType, 
				const DataType& dataType ) {

    ENTER_FCN( "HeatMaterial::SetScalar" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      Complex val;
      if ( dataType == REAL ) {
	val = Complex ( param, 0.0 );
      }
      else if (dataType == IMAG ) {
	val = Complex ( 0.0, param );
      }
      else {
	std::string msg = "SetScalar-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
      
      scalarParams_[matType] = val;
    }

  }

  void HeatMaterial::SetScalar( Complex& param, const MaterialType& matType, 
					   const DataType& dataType ) {

    ENTER_FCN( "HeatMaterial::SetScalar" );

    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == REAL ) {
	//	param = val.real();
      }
      else if ( dataType == IMAG ) {
	//	param = val.imag();
      }
      else if ( dataType == COMPLEX ) {
	param = val;
      }
    }
  }


  void HeatMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
					   const DataType& dataType ) {
    
    ENTER_FCN( "HeatMaterial::SetTensor" );
    std::string dim = "tensor";
    matTypeNotAllowed( matType, dim );
  }

  void HeatMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
					   const DataType& dataType ) {
    
    ENTER_FCN( "HeatMaterial::SetTensor" );
    std::string dim = "tensir";
    matTypeNotAllowed( matType, dim );
    
  }


  void HeatMaterial::GetScalar( Double& param, const MaterialType& matType, 
				const DataType& dataType ) const {

    ENTER_FCN( "HeatMaterial::GetScalar" );

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


  void HeatMaterial::GetScalar( Complex& param, const MaterialType& matType, 
					   const DataType& dataType ) const {

    ENTER_FCN( "HeatMaterial::GetScalar" );

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
