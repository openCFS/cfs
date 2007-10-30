#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "pyroelectricMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  PyroelectricMaterial::PyroelectricMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Pyroelectric";

    //set the allowed material parameters
	isAllowed_.insert( PYROCOEFFICIENT_TENSOR );
  }

  PyroelectricMaterial::~PyroelectricMaterial() {

    materialDatabaseName_ = "Pyroelectric";

  }

  void PyroelectricMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
					 const DataType& dataType ) {
    

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

	// to be consistent to old structure
	if ( dataType == REAL ) {
	  scalarParams_[matType] = Complex( param[2][2], 0.0);
	}
	else {
	  scalarParams_[matType] = Complex( 0.0, param[2][2]);
	  isComplex_.insert( matType );
	}
      }
      else {
	std::string msg = "SetTensor-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
  }

  void PyroelectricMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
					 const DataType& dataType ) {
    

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
	// to be consistent to old structure
	scalarParams_[matType] = param[2][2];
	isComplex_.insert( matType );
      }
    }    
  }

  void PyroelectricMaterial::GetTensor( Matrix<Double>& param, 
					 const MaterialType& matType, 
					 const DataType& dataType,
					 const SubTensorType subTensor) const {


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

  void PyroelectricMaterial::GetTensor( Matrix<Complex>& param, 
					 const MaterialType& matType, 
					 const DataType& dataType,
					 const SubTensorType subTensor) const {
    

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


  void PyroelectricMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
					       const MaterialType& matType, 
					       const SubTensorType& subTensor) const {


    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor axi or plane is the same
    matMatrix.Resize(2,2);
    matMatrix.Init();
    pos->second.GetSubMatrix(matMatrix, 0, 0);
  }

  void PyroelectricMaterial::SetScalar(Double param, MaterialType matType, 
                                       DataType dataType ) {
    Error("SetScalar for 'Double' not implemented",__FILE__,__LINE__);
  }

  void PyroelectricMaterial::GetScalar(Double& param, MaterialType matType, 
                                       DataType dataType ) const{
    Error("GetScalar for 'Double' not implemented",__FILE__,__LINE__);
  }

  void PyroelectricMaterial::GetScalar(Complex& param, MaterialType matType, 
                                       DataType dataType ) const{
    Error("GetScalar for 'Complex' not implemented",__FILE__,__LINE__);
  }


}
