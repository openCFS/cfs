#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "electroMagneticMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroMagneticMaterial::ElectroMagneticMaterial() : BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");
    materialDatabaseName_ = "Electromagnetics";

    //set the allowed material parameters
    isAllowed_.insert( MAG_PERMEABILITY );
    isAllowed_.insert( MAG_RELUCTIVITY );
    isAllowed_.insert( MAG_CONDUCTIVITY );
    isAllowed_.insert( ELEC_PERMITTIVITY );
    isAllowed_.insert( E_SATURATION );
    isAllowed_.insert( P_SATURATION );
    isAllowed_.insert( A_JILES );
    isAllowed_.insert( ALPHA_JILES );
    isAllowed_.insert( K_JILES );
    isAllowed_.insert( C_JILES );
    isAllowed_.insert( P_DIRECTION );
    isAllowed_.insert( HYST_MODEL );
  }

  ElectroMagneticMaterial::~ElectroMagneticMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }

  void ElectroMagneticMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "ElectroMagneticMaterial::SetScalar" );

    if ( matType == HYST_MODEL ) {
      hystType_ = param;
    }
    else {
      std::string dim = "string";
      matTypeNotAllowed( matType, dim );
    }
  }


  void ElectroMagneticMaterial::SetScalar( Integer& param, const MaterialType& matType) {

    ENTER_FCN( "ElectroMagneticMaterial::SetScalar" );

    if ( matType == P_DIRECTION ) {
      directionP_ = param;
    }
    else {
      std::string dim = "integer";
      matTypeNotAllowed( matType, dim );
    }
  }

  void ElectroMagneticMaterial::SetScalar( Double& param, const MaterialType& matType, 
					   const DataType& dataType ) {

    ENTER_FCN( "ElectroMagneticMaterial::SetScalar" );

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

    //check for permeability
    if ( matType == MAG_PERMEABILITY ) {
      if ( param < 1.255e-6 ) {
	Error("Mag. permeability cannot be smaller then the one of vacuum",
		__FILE__,__LINE__);
      }
      else {
	scalarParams_[MAG_RELUCTIVITY] = Complex( 1.0/param, 0.0 );
      }
    }

  }


  void ElectroMagneticMaterial::SetScalar( Complex& param, const MaterialType& matType, 
					   const DataType& dataType ) {

    ENTER_FCN( "ElectroMagneticMaterial::SetScalar" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      Complex val;
      if ( dataType == REAL ) {
	val = param.real();
      }
      else if (dataType == IMAG ) {
	val = param.imag();
      }
      else if ( dataType == COMPLEX ) {
	val = param;
      }
      
      scalarParams_[matType] = val;
    }

    //check for permeability
    if ( matType == MAG_PERMEABILITY ) {
      if ( param.real() < 1.255e-6 ) {
	Error("Mag. permeability cannot be smaller then the one of vacuum",
		__FILE__,__LINE__);
      }
      else {
	scalarParams_[MAG_RELUCTIVITY] = 1.0/param;
      }
    }
  
  }


  void ElectroMagneticMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
					   const DataType& dataType ) {
    
    ENTER_FCN( "ElectroMagneticMaterial::SetTensor" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
      matTypeNotAllowed( matType, dim );
    }
    else {
      Complex val;
      if ( dataType == REAL ) {
	//	val = param.real();
      }
      else if (dataType == IMAG ) {
	//	val = param.imag();
      }
      else {
	std::string msg = "SetScalar-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }

      //      tensorParams_[matType] = val;
    }

  }

  void ElectroMagneticMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
					   const DataType& dataType ) {
    
    ENTER_FCN( "ElectroMagneticMaterial::SetTensor" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      if ( dataType != COMPLEX ) {
	std::string msg = "SetTensor with Matrix<Complex>";
	setMakesNoSense( dataType, msg );
      }
      else {
	tensorParams_[matType] = param;
      }
    }

    if ( matType == ELEC_PERMITTIVITY ) {
      // currently we assume for scalar the 3-3-component
      scalarPermittivity_ = param[2][2];
    }
  }


  void ElectroMagneticMaterial::GetScalar( Integer& param, const MaterialType& matType) const {

    ENTER_FCN( "ElectroStaticMaterial::GetScalar" );

    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      if ( matType ==  P_DIRECTION ) {
	param = directionP_;
      }
      else {
      std::string dim = "integer";
      matTypeNotInDataBase( matType, dim );
      }
    }
  }

  void ElectroMagneticMaterial::GetScalar( Double& param, const MaterialType& matType, 
					   const DataType& dataType ) const {

    ENTER_FCN( "ElectroMagneticMaterial::GetScalar" );

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

  void ElectroMagneticMaterial::GetScalar( Complex& param, const MaterialType& matType, 
					   const DataType& dataType ) const {

    ENTER_FCN( "ElectroMagneticMaterial::GetScalar" );

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

  void ElectroMagneticMaterial::GetTensor( Matrix<Double>& param, 
					   const MaterialType& matType, 
					   const DataType& dataType,
					   const SubTensorType ) const {
    
    ENTER_FCN( "ElectroMagneticMaterial::GetTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
     Matrix<Complex> val = pos->second;
     if ( dataType == REAL ) {
       //	param = val.real();
      }
      else if ( dataType == IMAG ) {
	//	param = val.imag();
      }
      else {
	std::string msg = "GetTensor-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }

  void ElectroMagneticMaterial::GetTensor( Matrix<Complex>& param, 
					   const MaterialType& matType, 
					   const DataType& dataType,
					   const SubTensorType ) const {
    
    ENTER_FCN( "ElectroMagneticMaterial::GetTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Matrix<Complex> val = pos->second;
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
  
   void ElectroMagneticMaterial::Print(std::ostream & out) const {
    ENTER_FCN( "ElectroMagneticMaterial::Print" );
  }
  
}
