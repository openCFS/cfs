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

#include <Utils/SmoothSpline.hh>
#include "electroMagneticMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroMagneticMaterial::ElectroMagneticMaterial() : BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");
    materialDatabaseName_ = "Electromagnetics";
    nlinFncBH_ = NULL;

    //set the allowed material parameters
    isAllowed_.insert( MAG_PERMEABILITY );
    isAllowed_.insert( MAG_PERMEABILITY_1 );
    isAllowed_.insert( MAG_PERMEABILITY_2 );
    isAllowed_.insert( MAG_PERMEABILITY_3 );
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
    isAllowed_.insert( DATA_ACCURACY );
    isAllowed_.insert( MAX_APPROX_VAL );
  }

  ElectroMagneticMaterial::~ElectroMagneticMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }

  void ElectroMagneticMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "ElectroMagneticMaterial::SetScalar" );

    if ( matType == HYST_MODEL ) {
      stringParams_[matType] = param;
      isSet_.insert( matType );
    }
    else {
      std::string dim = "string";
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

    //check for permeability
    if ( matType == MAG_PERMEABILITY ) {
      if ( param < 1.255e-6 ) {
	EXCEPTION("Mag. permeability cannot be smaller then the one of vacuum" );
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

    //check for permeability
    if ( matType == MAG_PERMEABILITY ) {
      if ( param.real() < 1.255e-6 ) {
	EXCEPTION("Mag. permeability cannot be smaller then the one of vacuum" );
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

  void ElectroMagneticMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
					   const DataType& dataType ) {
    
    ENTER_FCN( "ElectroMagneticMaterial::SetTensor" );

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
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

    if ( matType == ELEC_PERMITTIVITY ) {
      // to be consistent to old structure
      scalarParams_[matType] = param[2][2];
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
					   const SubTensorType subTensor) const {
    
    ENTER_FCN( "ElectroMagneticMaterial::GetTensor" );


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

  void ElectroMagneticMaterial::GetTensor( Matrix<Complex>& param, 
					   const MaterialType& matType, 
					   const DataType& dataType,
					   const SubTensorType subTensor) const {
    
    ENTER_FCN( "ElectroMagneticMaterial::GetTensor" );

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
  
  
  void ElectroMagneticMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
						 const MaterialType& matType, 
						 const SubTensorType& subTensor) const {

    ENTER_FCN( "ElectroMagnetic::ComputeSubTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor axi or plane is the same
    matMatrix.Resize(2,2);
    matMatrix.Init();
    pos->second.GetSubMatrix(matMatrix, 1, 1);
  }


  void ElectroMagneticMaterial::InitApproxCurves() {
    ENTER_FCN( "ElectroMagneticMaterial::InitApproxCurves" );

    // check, if we need to approx BH curve
    if (  needApproxMatCurves_.find( magBH ) != needApproxMatCurves_.end() ) {
      std::string nlfnc = GetNonlinFileName();
      nlinFncBH_ = new SmoothSpline(nlfnc, BH);

      //get accuracy of approximation
      Double dataAccuracy;
      GetScalar( dataAccuracy, DATA_ACCURACY, REAL );

      //get maximal approximation value
      Double maxApproxVal;
      GetScalar( maxApproxVal, MAX_APPROX_VAL, REAL );

      nlinFncBH_->SetAccuracy( dataAccuracy );
      nlinFncBH_->SetMaxY( maxApproxVal );   //maximal value of magnetic induction B
      nlinFncBH_->CalcBestParameter();
      nlinFncBH_->CalcApproximation();
      nlinFncBH_->Print(); 
    }
  }

}
