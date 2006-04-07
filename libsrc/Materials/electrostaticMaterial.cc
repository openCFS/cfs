#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "electrostaticMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroStaticMaterial::ElectroStaticMaterial() : BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");
    materialDatabaseName_ = "Electrostatic";

    //set the allowed material parameters
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

  ElectroStaticMaterial::~ElectroStaticMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }

  void ElectroStaticMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "ElectroStaticMaterial::SetScalar" );

    if ( matType == HYST_MODEL ) {
      hystType_ = param;
    }
    else {
      std::string dim = "string";
      matTypeNotAllowed( matType, dim );
    }
  }


  void ElectroStaticMaterial::SetScalar( Integer& param, const MaterialType& matType) {

    ENTER_FCN( "ElectroStaticMaterial::SetScalar" );

    if ( matType == P_DIRECTION ) {
      directionP_ = param;
    }
    else {
      std::string dim = "integer";
      matTypeNotAllowed( matType, dim );
    }
    if ( matType == HYST_MODEL ) {
      hystType_ = param;
    }
    else {
      std::string dim = "string";
      matTypeNotAllowed( matType, dim );
    }
  }

  void ElectroStaticMaterial::SetScalar( Double& param, const MaterialType& matType, 
					   const DataType& dataType ) {

    ENTER_FCN( "ElectroStaticMaterial::SetScalar" );

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


  void ElectroStaticMaterial::SetScalar( Complex& param, const MaterialType& matType, 
					 const DataType& dataType ) {

    ENTER_FCN( "ElectroStaticMaterial::SetScalar" );

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
  }


  void ElectroStaticMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
					 const DataType& dataType ) {
    
    ENTER_FCN( "ElectroStaticMaterial::SetTensor" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
      matTypeNotAllowed( matType, dim );
    }
    else {
      if ( dataType == REAL || dataType == IMAG ) {
	if ( tensorParams_[matType].GetSizeRow() == 0 ) {
	  tensorParams_[matType].Resize( param.GetSizeRow(), param.GetSizeCol() );
	}
	tensorParams_[matType].SetPart( dataType, param );

	// to be consistent to old structure
	if ( dataType == REAL ) {
	  scalarParams_[matType] = Complex( param[2][2], 0.0);
	}
	else {
	  scalarParams_[matType] = Complex( 0.0, param[2][2]);
	}
      }
      else {
	std::string msg = "SetTensor-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
  }

  void ElectroStaticMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
					 const DataType& dataType ) {
    
    ENTER_FCN( "ElectroStaticMaterial::SetTensor" );

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
	// to be consistent to old structure
	scalarParams_[matType] = param[2][2];
      }
    }    
  }


  void ElectroStaticMaterial::GetScalar( Integer& param, const MaterialType& matType, 
					   const DataType& dataType )  const {

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



  void ElectroStaticMaterial::GetScalar( Double& param, const MaterialType& matType, 
					   const DataType& dataType ) const {

    ENTER_FCN( "ElectroStaticMaterial::GetScalar" );

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


  void ElectroStaticMaterial::GetScalar( Complex& param, const MaterialType& matType, 
					 const DataType& dataType ) const {

    ENTER_FCN( "ElectroStaticMaterial::GetScalar" );


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

  void ElectroStaticMaterial::GetTensor( Matrix<Double>& param, 
					 const MaterialType& matType, 
					 const DataType& dataType,
					 const SubTensorType subTensor) const {

    ENTER_FCN( "ElectroStaticMaterial::GetTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Matrix<Complex> matTensor;
      if ( subTensor == FULL ) {
	matTensor = pos->second;
      }
      else {
	ComputeSubTensor(matTensor, matType, subTensor);
      }

      if ( dataType == REAL || dataType == IMAG) {
	param = matTensor.GetPart( dataType );
      }
      else {
	std::string msg = "GetTensor-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }

  void ElectroStaticMaterial::GetTensor( Matrix<Complex>& param, 
					 const MaterialType& matType, 
					 const DataType& dataType,
					 const SubTensorType subTensor) const {
    
    ENTER_FCN( "ElectroStaticMaterial::GetTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Matrix<Complex> matTensor;
      if ( subTensor == FULL ) {
	matTensor = pos->second;
      }
      else {
	ComputeSubTensor(matTensor, matType, subTensor);
      }

      if ( dataType == REAL || dataType == IMAG) {
	Matrix<Double> help; 
	help = matTensor.GetPart( dataType );
	param.Resize( matTensor.GetSizeRow(), matTensor.GetSizeCol() );
	param.SetPart( dataType, help );
      }
      else if ( dataType == COMPLEX ) {
	param = matTensor;
      }
    }
  }
  


  void ElectroStaticMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
					       const MaterialType& matType, 
					       const SubTensorType& subTensor) const {

    ENTER_FCN( "ElectroStaticMaterial::ComputeSubTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor axi or plane is the same
    matMatrix.Resize(2,2);
    pos->second.GetSubMatrix(matMatrix, 1, 1);
  }



}
