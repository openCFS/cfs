#include "ElectroStaticMaterial.hh"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <limits.h>
#include <string>

#include "Materials/Models/Preisach.hh"

namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroStaticMaterial::ElectroStaticMaterial(MathParser* mp,
                                               CoordSystem * defaultCoosy) 
  : BaseMaterial(mp, defaultCoosy) {

    materialDatabaseName_ = "Electrostatic";

    //set the allowed material parameters
    isAllowed_.insert( ELEC_PERMITTIVITY );
    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
    isAllowed_.insert( Y_REMANENCE );

    isAllowed_.insert( PREISACH_WEIGHTS );
    isAllowed_.insert( PREISACH_WEIGHTS_DIM );
    isAllowed_.insert( PREISACH_WEIGHTS_CONSTVALUE );
    isAllowed_.insert( PREISACH_WEIGHTS_TYPE );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_A );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H2 );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA2 );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_ETA );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE );

    isAllowed_.insert( PREISACH_WEIGHTS_TENSOR );
    isAllowed_.insert( PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_ONLY );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_CINATAN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_A );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_B );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_C );
    isAllowed_.insert( PREISACH_MAYERGOYZ_NUM_DIR );
    isAllowed_.insert( PREISACH_MAYERGOYZ_ISOTROPIC );
    isAllowed_.insert( PREISACH_MAYERGOYZ_CLIPOUTPUT );

    isAllowed_.insert( MAYERGOYZ_STARTAXIS_X );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Y );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Z );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_X_STRAIN );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Y_STRAIN );
    isAllowed_.insert( MAYERGOYZ_STARTAXIS_Z_STRAIN );

    isAllowed_.insert( PREISACH_PRESCRIBEOUTPUT );
    isAllowed_.insert( PREISACH_SCALEINITIALSTATE );
    isAllowed_.insert( SCALETOSAT );
    isAllowed_.insert( SCALETOSAT_STRAIN );
    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
    isAllowed_.insert( Y_REMANENCE );
    isAllowed_.insert( A_JILES );
    isAllowed_.insert( ALPHA_JILES );
    isAllowed_.insert( K_JILES );
    isAllowed_.insert( C_JILES );
    isAllowed_.insert( JILES_TEST );
    isAllowed_.insert( P_DIRECTION );
    isAllowed_.insert( P_DIRECTION_X );
    isAllowed_.insert( P_DIRECTION_Y );
    isAllowed_.insert( P_DIRECTION_Z );
    isAllowed_.insert( EVAL_VERSION );
    isAllowed_.insert( PRINT_PREISACH );
    isAllowed_.insert( PRINT_PREISACH_RESOLUTION );
    isAllowed_.insert( IS_TESTING );
    isAllowed_.insert( ANG_DISTANCE );
    isAllowed_.insert( ANG_CLIPPING );
    isAllowed_.insert( ANG_RESOLUTION );
    isAllowed_.insert( AMP_RESOLUTION );
		isAllowed_.insert( INITIAL_STATE );
    isAllowed_.insert( INITIAL_STATE_X );
    isAllowed_.insert( INITIAL_STATE_Y );
    isAllowed_.insert( INITIAL_STATE_Z );
    isAllowed_.insert( HYSTERESIS_DIM );
    isAllowed_.insert( HYST_STRAIN_FORM );
    isAllowed_.insert( HYST_IRRSTRAINS );
    isAllowed_.insert( HYST_IRRSTRAIN_C1 );
    isAllowed_.insert( HYST_IRRSTRAIN_C2 );
    isAllowed_.insert( HYST_IRRSTRAIN_C3 );
    isAllowed_.insert( HYST_IRRSTRAIN_CI );
    isAllowed_.insert( HYST_IRRSTRAIN_CI_SIZE );
    isAllowed_.insert( HYST_IRRSTRAIN_D0 );
    isAllowed_.insert( HYST_IRRSTRAIN_D1 );
    isAllowed_.insert( HYST_IRRSTRAIN_SCALETOSAT );
    isAllowed_.insert( HYST_IRRSTRAIN_PARAMSFORHALFRANGE );
    isAllowed_.insert( ROT_RESISTANCE );
    isAllowed_.insert( HYST_MODEL );
    isAllowed_.insert( HYST_COUPLING_DEFINED );
    isAllowed_.insert( X_SATURATION_STRAIN );
    isAllowed_.insert( Y_SATURATION_STRAIN );
    isAllowed_.insert( S_SATURATION );
    isAllowed_.insert( PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_TENSOR_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_A_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_H2_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_SIGMA2_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_ETA_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_TYPE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_CONSTVALUE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_DIM_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_STRAIN );
    isAllowed_.insert( HYST_TYPE_IS_PREISACH_STRAIN );
    isAllowed_.insert( HYST_TYPE_IS_PREISACH );
    isAllowed_.insert( ANG_DISTANCE_STRAIN );
    isAllowed_.insert( EVAL_VERSION_STRAIN );
    isAllowed_.insert( ROT_RESISTANCE_STRAIN );
    isAllowed_.insert( HYSTERESIS_DIM_STRAIN );
    isAllowed_.insert( HYST_MODEL_STRAIN );
    isAllowed_.insert( S_DIRECTION_Z );
    isAllowed_.insert( S_DIRECTION_Y );
    isAllowed_.insert( S_DIRECTION_X );
    isAllowed_.insert( S_DIRECTION );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_CINATAN_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_ONLY_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_C_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_B_STRAIN );
    isAllowed_.insert( PREISACH_WEIGHTS_ANHYST_A_STRAIN );
    isAllowed_.insert( ANG_CLIPPING_STRAIN );
    isAllowed_.insert( ANG_RESOLUTION_STRAIN );
    isAllowed_.insert( AMP_RESOLUTION_STRAIN );
    isAllowed_.insert( PREISACH_MAYERGOYZ_NUM_DIR_STRAIN );
    isAllowed_.insert( PREISACH_MAYERGOYZ_ISOTROPIC_STRAIN );
    isAllowed_.insert( PREISACH_MAYERGOYZ_CLIPOUTPUT_STRAIN );
    isAllowed_.insert( IRRSTRAIN_REUSE_P );
    isAllowed_.insert( NONLIN_COEFFICIENT );
    isAllowed_.insert( NONLIN_DEPENDENCY );
    isAllowed_.insert( NONLIN_APPROXIMATION_TYPE );
    isAllowed_.insert( NONLIN_DATA_NAME );
  }

  ElectroStaticMaterial::~ElectroStaticMaterial() {


  }
  
  void ElectroStaticMaterial::Finalize() {
      ComputeFullEpsTensor();
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

        // added this check to avoid seg-faults for tensors of size Nx1
        // ( normal material parameter (like permittivity) do not need this check
        //   but for some hysteresis parameter, Nx1 arrays are needed)
        if(param.GetNumRows() >= 2 && param.GetNumCols() >= 2){
          // to be consistent to old structure
          if ( dataType == Global::REAL ) {
             scalarParams_[matType] = Complex( param[2][2], 0.0); 
          }
          else {
            scalarParams_[matType] = Complex( 0.0, param[2][2]);
            isComplex_.insert( matType );
          }
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
  
  void ElectroStaticMaterial::ComputeFullEpsTensor() {
    Matrix<Complex> epsTensor(3,3);
    
    // check, if tensor was already set
    if( tensorParams_.find( ELEC_PERMITTIVITY ) == tensorParams_.end() 
        && tensorCoef_.find(ELEC_PERMITTIVITY ) == tensorCoef_.end() ) { 
      Complex eps;
      GetScalar(eps, ELEC_PERMITTIVITY, Global::COMPLEX);
      epsTensor[0][0] = eps;
      epsTensor[1][1] = eps;
      epsTensor[2][2] = eps;
      SetTensor( epsTensor, ELEC_PERMITTIVITY, Global::COMPLEX );
    }
  }


  //============================== Hysteresis =====================================
//
//  Double ElectroStaticMaterial::ComputeScalarDiffVal( UInt nrElem, Double Xval ) {
//
//    Double matDiff, eps;
//
//    UInt idx = globalElem2Local_[nrElem];
//    Double Ycurrent = hyst_->computeValueAndUpdate(Xval, idx);
//
//    //    std::cout << "epsDiff: " << " Xval=" << Xval << "  Yval=" << Ycurrent << std::endl;
//
//    //compute differential material parameter
//    Double dX = Xval - Xprevious_[idx];
//    Double dY = Ycurrent -Yprevious_[idx];
//
//    if ( (abs(dY) < 1e-12) || (abs(dX) < 1e-10) ) {
//      GetScalar(eps,ELEC_PERMITTIVITY,Global::REAL);
//      matDiff = eps;
//    }
//    else {
//      matDiff = dY / dX;
//    }
//
//    return matDiff;
//  }
//
//
//  void ElectroStaticMaterial::SetPreviousHystVal( UInt nrElem, Double Xval ) {
//
//    UInt idx = globalElem2Local_[nrElem];
//
//    Xprevious_[idx] = Xval;
//    Yprevious_[idx] = hyst_->computeValueAndUpdate( Xval, idx );
//  }
//
//
//  Double ElectroStaticMaterial::ComputeScalarHystVal( UInt nrElem, Double Xval ) {
//
//    UInt idx    = globalElem2Local_[nrElem];
//    Double Yval = hyst_->computeValueAndUpdate( Xval, idx );
//
//    return Yval;
//  }
}
