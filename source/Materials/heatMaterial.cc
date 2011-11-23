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

#include "heatMaterial.hh"
#include <Utils/LinInterpolate.hh>


namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  HeatMaterial::HeatMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Heat";

    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( HEAT_CONDUCTIVITY );
    isAllowed_.insert( HEAT_CONDUCTIVITY_TENSOR );
    isAllowed_.insert( HEAT_CAPACITY );
    isAllowed_.insert( NONLIN_APPROXIMATION_TYPE );
    isAllowed_.insert( DATA_ACCURACY );
    isAllowed_.insert( MAX_APPROX_VAL );
  }

  HeatMaterial::~HeatMaterial() {

    materialDatabaseName_ = "Heat";

  }

  void HeatMaterial::SetScalar( Double param, MaterialType matType, 
				Global::ComplexPart dataType ) {


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

  void HeatMaterial::SetScalar( Complex param, MaterialType matType, 
				Global::ComplexPart dataType ) {


    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      isSet_.insert( matType );

      Complex val = pos->second;
      if ( dataType == Global::REAL ) {
	//	param = val.real();
      }
      else if ( dataType == Global::IMAG ) {
	//	param = val.imag();
	isComplex_.insert( matType );
      }
      else if ( dataType == Global::COMPLEX ) {
	param = val;
	isComplex_.insert( matType );
      }
    }
  }


  void HeatMaterial::GetScalar( Double& param, MaterialType matType, 
				Global::ComplexPart dataType ) const {


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


  void HeatMaterial::GetScalar( Complex& param, MaterialType matType, 
					   Global::ComplexPart dataType ) const {


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
  
  void HeatMaterial::SetTensor(const Matrix<Double>& param, 
                               MaterialType matType, 
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

	// to be consistent to old structure
	if ( dataType == Global::REAL ) {
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
  
  void HeatMaterial::GetTensor( Matrix<Double>& param, 
					 MaterialType matType, 
					 Global::ComplexPart dataType,
					 SubTensorType subTensor) const {


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
  
  void HeatMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
					       const MaterialType& matType, 
					       const SubTensorType& subTensor) const {

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor axi or plane is the same
    matMatrix.Resize(2,2);
    matMatrix.Init();
    pos->second.GetSubMatrix(matMatrix, 0, 0);
    

  }
  

  void HeatMaterial::InitApproxCurves() {

    // check, if we need to approx curve
    if ( needApproxMatCurves_.find( HEAT_CONDUCTIVITY ) != needApproxMatCurves_.end() ) {
      std::string nlfncName = nonLinMatInfo_[HEAT_CONDUCTIVITY].fileName;
      if ( nlfncName != "" ) 
        std::cout << "Cond: FileName: " << nlfncName << std::endl;
        nlinFncConductivity_ = new LinInterpolate( nlfncName, HEAT_CONDUCTIVITY );
    }

    if ( needApproxMatCurves_.find( HEAT_CAPACITY ) != needApproxMatCurves_.end() ) {
      std::string nlfncName = nonLinMatInfo_[HEAT_CAPACITY].fileName;
      if ( nlfncName != "" ) 
        std::cout << "Capacity: FileName: " << nlfncName << std::endl;
        nlinFncCapacity_ = new LinInterpolate( nlfncName, HEAT_CAPACITY );
    }

  }   
}
