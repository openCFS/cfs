#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "piezoMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  PiezoMaterial::PiezoMaterial() : BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");

    //set the allowed material parameters
    isAllowed_.insert( PIEZO_TENSOR );
    isAllowed_.insert( E_SATURATION );
    isAllowed_.insert( P_SATURATION );
    isAllowed_.insert( A_JILES );
    isAllowed_.insert( ALPHA_JILES );
    isAllowed_.insert( K_JILES );
    isAllowed_.insert( C_JILES );
    isAllowed_.insert( P_DIRECTION );
    isAllowed_.insert( HYST_MODEL );
    isAllowed_.insert( NONLIN_COEFFICIENT );
    isAllowed_.insert( NONLIN_DEPENDENCY );
    isAllowed_.insert( NONLIN_APPROXIMATION_TYPE );
    isAllowed_.insert( NONLIN_DATA_NAME );


  }

  PiezoMaterial::~PiezoMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");
  }


  void PiezoMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );
    if ( matType == HYST_MODEL ) {
      stringParams_[matType] = param;
      isSet_.insert( matType );
    }
    else {
      if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
        std::string dim = "scalar";
        matTypeNotAllowed( matType, dim );
      }
      else {
        isSet_.insert( matType );
      }     
      stringParams_[matType] = param;
    }
  }


  void PiezoMaterial::SetScalar( Double& param, const MaterialType& matType, 
				 const DataType& dataType ) {

    ENTER_FCN( "PiezoMaterial::SetScalar" );

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


  void PiezoMaterial::SetScalar( Complex& param, const MaterialType& matType, 
					   const DataType& dataType ) {

    ENTER_FCN( "PiezoMaterial::SetScalar" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val;
      if ( dataType == REAL ) {
	val = param.real();
      }
      else if (dataType == IMAG ) {
	val = param.imag();
	isComplex_.insert( matType );
      }
      else if ( dataType == COMPLEX ) {
	val = param;
	isComplex_.insert( matType );
      }
      
      scalarParams_[matType] = val;
    }
  }


  void PiezoMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
				 const DataType& dataType ) {
    
    ENTER_FCN( "PiezoMaterial::SetTensor" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      if ( dataType == REAL || dataType == IMAG ) {
	if ( tensorParams_[matType].GetSizeRow() == 0 ) {
	  tensorParams_[matType].Resize( param.GetSizeRow(), param.GetSizeCol() );
          tensorParams_[matType].Init();
	}
	if ( tensorParamsOrig_[matType].GetSizeRow() == 0 ) {
	  tensorParamsOrig_[matType].Resize( param.GetSizeRow(), param.GetSizeCol() );
          tensorParamsOrig_[matType].Init();
	}

	tensorParams_[matType].SetPart( dataType, param );
	tensorParamsOrig_[matType].SetPart( dataType, param );
	if ( dataType == IMAG ) {
	  isComplex_.insert( matType );
	}
      }
      else {
	std::string msg = "SetTensor-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
  }

  void PiezoMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
				 const DataType& dataType ) {
    
    ENTER_FCN( "PiezoMaterial::SetTensor" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      if ( dataType != COMPLEX ) {
	std::string msg = "SetTensor with Matrix<Complex>";
	setMakesNoSense( dataType, msg );
      }
      else {
	tensorParams_[matType]     = param;
	tensorParamsOrig_[matType] = param;
	isComplex_.insert( matType );
      }
    }
  }

  void PiezoMaterial::GetScalar( std::string& param, const MaterialType& matType)  const {

    ENTER_FCN( "PiezoMaterial::GetScalar" );

    stringMap::const_iterator pos;
    pos = stringParams_.find( matType );
    std::string value;

    if ( pos == stringParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      param=pos->second;
    }
  }    
 
   void PiezoMaterial::GetScalar( Integer& param, const MaterialType& matType)  const {
    
     ENTER_FCN( "PiezoMaterial::GetScalar" );
    
     integerMap::const_iterator pos;
     pos = integerParams_.find( matType );
     std::string value;
    
     if ( pos == integerParams_.end() ) {
       std::string dim = "scalar";
       matTypeNotInDataBase( matType, dim );
     }
     else {
       param=pos->second;
     }
   }  


  void PiezoMaterial::GetScalar( Double& param, const MaterialType& matType, 
					   const DataType& dataType )  const {

    ENTER_FCN( "PiezoMaterial::GetScalar" );

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

  void PiezoMaterial::GetScalar( Complex& param, const MaterialType& matType, 
				 const DataType& dataType )  const {

    ENTER_FCN( "PiezoMaterial::GetScalar" );

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

  void PiezoMaterial::GetTensor( Matrix<Double>& param, 
				 const MaterialType& matType, 
				 const DataType& dataType,
				 const SubTensorType subTensor) const {

    ENTER_FCN( "PiezoMaterial::GetTensor" );

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

  void PiezoMaterial::GetTensor( Matrix<Complex>& param, 
				 const MaterialType& matType, 
				 const DataType& dataType,
				 const SubTensorType subTensor) const {
    
    ENTER_FCN( "PiezoMaterial::GetTensor" );

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
        param.Init();
	param.SetPart( dataType, help );
      }
      else if ( dataType == COMPLEX ) {
	param = matTensor;
      }
    }
  }


  void PiezoMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
				       const MaterialType& matType, 
				       const SubTensorType& subTensor) const {

    ENTER_FCN( "PiezoMaterial::ComputeSubTensors" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    Matrix<Complex> const &mat = pos->second;

    if ( subTensor == AXI ) {
      matMatrix.Resize(2,4);

      matMatrix[0][0] = mat[0][0];
      matMatrix[0][1] = mat[0][1];
      matMatrix[0][2] = mat[0][5];
      matMatrix[0][3] = mat[0][2];
      matMatrix[1][0] = mat[1][0];
      matMatrix[1][1] = mat[1][1];
      matMatrix[1][2] = mat[1][5];
      matMatrix[1][3] = mat[1][2];
    }
    else if ( subTensor == PLANE_STRAIN ||
              subTensor == PLANE_STRESS ) {
      matMatrix.Resize(2,3);

      matMatrix[0][0] = mat[0][0];
      matMatrix[0][1] = mat[0][1];
      matMatrix[0][2] = mat[0][5];
      matMatrix[1][0] = mat[1][0];
      matMatrix[1][1] = mat[1][1];
      matMatrix[1][2] = mat[1][5];

    } else {
      subTensorNotAvailable( matType, subTensor );
    }
  }
  
}
