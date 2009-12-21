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

    // micro-piezoelectric model
    isAllowed_.insert( MEAN_TEMPERATURE );
    isAllowed_.insert( SPON_POLARIZATION );
    isAllowed_.insert( SPON_STRAIN );
    isAllowed_.insert( EFIELD0 );
    isAllowed_.insert( STRESS0 );
    isAllowed_.insert( DCOUPLE0 );
    isAllowed_.insert( DRIVING_FORCE_90 );
    isAllowed_.insert( DRIVING_FORCE_180 );
    isAllowed_.insert( RATE_CONSTANT );
    isAllowed_.insert( VISCO_PLASTIC_INDEX );
    isAllowed_.insert( SATURATION_INDEX );
    isAllowed_.insert( SCALE_FORCE_ELEC );
    isAllowed_.insert( SCALE_FORCE_MECH );
    isAllowed_.insert( SCALE_FORCE_COUPLE );
    isAllowed_.insert( VOLUME_FRAC_INIT );
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


  void PiezoMaterial::SetScalar( Double param, MaterialType matType, Global::ComplexPart dataType ) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val;
      if ( dataType == Global::REAL ) {
	val = Complex ( param, 0.0 );
      }
      else if (dataType == Global::IMAG ) {
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


  void PiezoMaterial::SetScalar( Complex param, MaterialType matType, Global::ComplexPart dataType ) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val;
      if ( dataType == Global::REAL ) {
	val = param.real();
      }
      else if (dataType == Global::IMAG ) {
	val = param.imag();
	isComplex_.insert( matType );
      }
      else if ( dataType == Global::COMPLEX ) {
	val = param;
	isComplex_.insert( matType );
      }
      
      scalarParams_[matType] = val;
    }
  }


  void PiezoMaterial::SetTensor(const Matrix<Double>& param, MaterialType matType, 
                                Global::ComplexPart dataType ) {
    

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      if ( dataType == Global::REAL || dataType == Global::IMAG ) {
	if ( tensorParams_[matType].GetNumRows() == 0 ) {
	  tensorParams_[matType].Resize( param.GetNumRows(), param.GetNumCols() );
          tensorParams_[matType].Init();
	}
	if ( tensorParamsOrig_[matType].GetNumRows() == 0 ) {
	  tensorParamsOrig_[matType].Resize( param.GetNumRows(), param.GetNumCols() );
          tensorParamsOrig_[matType].Init();
	}

	tensorParams_[matType].SetPart( dataType, param );
	tensorParamsOrig_[matType].SetPart( dataType, param );
	if ( dataType == Global::IMAG ) {
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
				                         Global::ComplexPart dataType ) {
    

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      if ( dataType != Global::COMPLEX ) {
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

                      
  void PiezoMaterial::GetScalar( Double& param, MaterialType matType, Global::ComplexPart dataType )  const {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == Global::REAL ) {
	param = val.real();
      }
      else if ( dataType == Global::IMAG ) {
	param = val.imag();
      }
      else {
	std::string msg = "GetScalar-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }    
  }

  void PiezoMaterial::GetScalar( Complex& param, MaterialType matType, Global::ComplexPart dataType )  const {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Complex val = pos->second;
      if ( dataType == Global::REAL ) {
	Complex valReal = Complex (val.real(), 0.0);
	param = valReal;
      }
      else if ( dataType == Global::IMAG ) {
	Complex valImag = Complex (0.0, val.imag());
	param = valImag;
      }
      else if ( dataType == Global::COMPLEX ) {
	param = val;
      }
    }    
  }

  void PiezoMaterial::GetTensor( Matrix<Double>& param, MaterialType matType, 
				                         Global::ComplexPart dataType, SubTensorType subTensor) const {


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

      if ( dataType == Global::REAL || dataType == Global::IMAG) {
	param = matTensor.GetPart( dataType );
      }
      else {
	std::string msg = "GetTensor-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }

  void PiezoMaterial::GetTensor( Matrix<Complex>& param, MaterialType matType, 
				                         Global::ComplexPart dataType, SubTensorType subTensor) const {
    

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

      if ( dataType == Global::REAL || dataType == Global::IMAG) {
	Matrix<Double> help; 
	help = matTensor.GetPart( dataType );
	param.Resize( matTensor.GetNumRows(), matTensor.GetNumCols() );
        param.Init();
	param.SetPart( dataType, help );
      }
      else if ( dataType == Global::COMPLEX ) {
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
