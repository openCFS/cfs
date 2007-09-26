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

#include "piezoMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  PiezoMaterial::PiezoMaterial() : BaseMaterial() {


    //set the allowed material parameters
    isAllowed_.insert( PIEZO_TENSOR );
    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
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

  }

  void PiezoMaterial::SetScalar(const std::string& param, MaterialType matType) {

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


  void PiezoMaterial::SetScalar( Double param, MaterialType matType, DataType dataType ) {


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


  void PiezoMaterial::SetScalar( Complex param, MaterialType matType, DataType dataType ) {


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


  void PiezoMaterial::SetTensor(const Matrix<Double>& param, MaterialType matType, 
                                DataType dataType ) {
    

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

  void PiezoMaterial::SetTensor(const Matrix<Complex>& param, MaterialType matType, 
				                         DataType dataType ) {
    

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

  void PiezoMaterial::GetScalar( std::string& param, MaterialType matType)  const {


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
 
   void PiezoMaterial::GetScalar( Integer& param, MaterialType matType)  const {
    
    
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

                      
  void PiezoMaterial::GetScalar( Double& param, MaterialType matType, DataType dataType )  const {


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

  void PiezoMaterial::GetScalar( Complex& param, MaterialType matType, DataType dataType )  const {


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

  void PiezoMaterial::GetTensor( Matrix<Double>& param, MaterialType matType, 
				                         DataType dataType, SubTensorType subTensor) const {


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

  void PiezoMaterial::GetTensor( Matrix<Complex>& param, MaterialType matType, 
				                         DataType dataType, SubTensorType subTensor) const {
    

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
				                               MaterialType matType, SubTensorType subTensor) const {


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
