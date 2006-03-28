#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "acousticMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  AcousticMaterial::AcousticMaterial() : BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");
    materialDatabaseName_ = "Acoustic";

    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( ACOU_BULK_MODULUS );
    isAllowed_.insert( ACOU_SOUND_SPEED );
    isAllowed_.insert( ACOU_ALPHA );
    isAllowed_.insert( FRACTIONAL_EXPONENT );
    isAllowed_.insert( RAYLEIGH_ALPHA );
    isAllowed_.insert( RAYLEIGH_BETA );
    isAllowed_.insert( RAYLEIGH_FREQUENCY);
    isAllowed_.insert( LOSS_TANGENS_DELTA);
    isAllowed_.insert( BOVERA );
  }

  AcousticMaterial::~AcousticMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }

  void AcousticMaterial::SetScalar( Integer& param, const MaterialType& matType) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );
    Error("SetScalar for 'Integer' not implemented",__FILE__,__LINE__);
  }

  void AcousticMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );
    Error("SetScalar for 'String' not implemented",__FILE__,__LINE__);
  }


  void AcousticMaterial::SetScalar( Double& param, const MaterialType& matType, 
				    const DataType& dataType ) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );

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


  void AcousticMaterial::SetScalar( Complex& param, const MaterialType& matType, 
				    const DataType& dataType ) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );
    dataTypeNotAllowed( dataType, matType );
  }


  void AcousticMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
				    const DataType& dataType ) {
    
    ENTER_FCN( "AcousticMaterial::SetTensor" );
    dataTypeNotAllowed( dataType, matType );
  }

  void AcousticMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
				    const DataType& dataType ) {
    
    ENTER_FCN( "AcousticMaterial::SetTensor" );
    dataTypeNotAllowed( dataType, matType );    
  }


  void AcousticMaterial::GetScalar( Double& param, const MaterialType& matType, 
				    const DataType& dataType ) const {

    ENTER_FCN( "AcousticMaterial::GetScalar" );

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

  void AcousticMaterial::GetScalar( Complex& param, const MaterialType& matType, 
				    const DataType& dataType ) const {

    ENTER_FCN( "AcousticMaterial::GetScalar" );
    dataTypeNotAllowed( dataType, matType );    
  }

}
