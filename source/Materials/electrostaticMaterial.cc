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

#include "Utils/preisach.hh"
#include "electrostaticMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroStaticMaterial::ElectroStaticMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Electrostatic";

    //set the allowed material parameters
    isAllowed_.insert( ELEC_PERMITTIVITY );
    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
    isAllowed_.insert( Y_REMANENCE );
    isAllowed_.insert( PREISACH_WEIGHTS );
    isAllowed_.insert( A_JILES );
    isAllowed_.insert( ALPHA_JILES );
    isAllowed_.insert( K_JILES );
    isAllowed_.insert( C_JILES );
    isAllowed_.insert( P_DIRECTION );
    isAllowed_.insert( NONLIN_COEFFICIENT );
    isAllowed_.insert( NONLIN_DEPENDENCY );
    isAllowed_.insert( NONLIN_APPROXIMATION_TYPE );
    isAllowed_.insert( NONLIN_DATA_NAME );
    isAllowed_.insert( MAG_PERMEABILITY ); //only for maxwell homogenization


  }

  ElectroStaticMaterial::~ElectroStaticMaterial() {


  }

  void ElectroStaticMaterial::SetScalar(const std::string& param, MaterialType matType) {


    if ( matType == HYST_MODEL || matType == P_DIRECTION ) {
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


  void ElectroStaticMaterial::SetScalar( Double param, MaterialType matType, 
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


  void ElectroStaticMaterial::SetScalar( Complex param, MaterialType matType, 
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


  void ElectroStaticMaterial::SetTensor(const Matrix<Double>& param, MaterialType matType, 
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

  void ElectroStaticMaterial::SetTensor(const Matrix<Complex>& param, MaterialType matType, 
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
        // to be consistent to old structure
        scalarParams_[matType] = param[2][2];
        isComplex_.insert( matType );
      }
    }    
  }


  void ElectroStaticMaterial::GetScalar( std::string& param, 
                                         MaterialType matType) const {


    stringMap::const_iterator pos;
    pos = stringParams_.find( matType );

    if ( pos == stringParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      param = pos->second;
    }    
  }


  void ElectroStaticMaterial::GetScalar( Double& param, MaterialType matType, 
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


  void ElectroStaticMaterial::GetScalar( Complex& param, MaterialType matType, 
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



  void ElectroStaticMaterial::GetScalar( Integer& param, MaterialType matType)  const {


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

  void ElectroStaticMaterial::GetTensor( Matrix<Double>& param, 
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

  void ElectroStaticMaterial::GetTensor( Matrix<Complex>& param, 
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
        Matrix<Double> help; 
        help = matTensor.GetPart( dataType );
        param.Resize( matTensor.GetNumRows(), matTensor.GetNumCols() );
        param.SetPart( dataType, help );
      }
      else if ( dataType == Global::COMPLEX ) {
        param = matTensor;
      }
    }
  }



  void ElectroStaticMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
                                               MaterialType matType, 
                                               SubTensorType subTensor) const {


    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor axi or plane is the same
    matMatrix.Resize(2,2);
    matMatrix.Init();
    pos->second.GetSubMatrix(matMatrix, 0, 0);


  }


  //============================== Hysteresis =====================================

  Double ElectroStaticMaterial::ComputeScalarDiffVal( UInt nrElem, Double Xval ) {

    Double matDiff, eps;

    UInt idx = globalElem2Local_[nrElem];
    Double Ycurrent = hyst_->computeValueAndUpdate(Xval, idx);

    //    std::cout << "epsDiff: " << " Xval=" << Xval << "  Yval=" << Ycurrent << std::endl;

    //compute differential material parameter
    Double dX = Xval - Xprevious_[idx];
    Double dY = Ycurrent -Yprevious_[idx];

    if ( (abs(dY) < 1e-12) || (abs(dX) < 1e-10) ) {
      GetScalar(eps,ELEC_PERMITTIVITY,Global::REAL);
      matDiff = eps;
    }
    else {
      matDiff = dY / dX;
    }

    return matDiff;
  }


  void ElectroStaticMaterial::SetPreviousHystVal( UInt nrElem, Double Xval ) {

    UInt idx = globalElem2Local_[nrElem];

    Xprevious_[idx] = Xval;
    Yprevious_[idx] = hyst_->computeValueAndUpdate( Xval, idx );
  }


  Double ElectroStaticMaterial::ComputeScalarHystVal( UInt nrElem, Double Xval ) {

    UInt idx    = globalElem2Local_[nrElem];
    Double Yval = hyst_->computeValueAndUpdate( Xval, idx );

    return Yval;
  }


}
