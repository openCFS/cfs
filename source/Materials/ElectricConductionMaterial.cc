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

#include "ElectricConductionMaterial.hh"
#include <Utils/LinInterpolate.hh>
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
#include "Utils/BiLinInterpolate.hh"


namespace CoupledField
{

// ***********************
//   Default Constructor
// ***********************
ElectricConductionMaterial::ElectricConductionMaterial(MathParser* mp,
                           CoordSystem * defaultCoosy) 
: BaseMaterial(mp, defaultCoosy) {

  materialDatabaseName_ = "ElectricConduction";

  //set the allowed material parameters
  isAllowed_.insert( ELEC_CONDUCTIVITY );
  isAllowed_.insert( ELEC_CONDUCTIVITY_TENSOR );
  isAllowed_.insert( NONLIN_APPROXIMATION_TYPE );
  isAllowed_.insert( DATA_ACCURACY );
  isAllowed_.insert( MAX_APPROX_VAL );
  isAllowed_.insert( NONLIN_DEPENDENCY );
}

ElectricConductionMaterial::~ElectricConductionMaterial() {

  materialDatabaseName_ = "ElectricConduction";

}

void ElectricConductionMaterial::Finalize() {
  ComputeFullMuTensor();
}


void ElectricConductionMaterial::SetScalar( Double param, MaterialType matType, 
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

void ElectricConductionMaterial::SetScalar( Complex param, MaterialType matType, 
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


void ElectricConductionMaterial::GetScalar( Double& param, MaterialType matType, 
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


void ElectricConductionMaterial::GetScalar( Complex& param, MaterialType matType, 
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

void ElectricConductionMaterial::SetTensor(const Matrix<Double>& param, 
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

void ElectricConductionMaterial::GetTensor( Matrix<Double>& param, 
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

void ElectricConductionMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
                                    const MaterialType& matType, 
                                    const SubTensorType& subTensor) const {

  tensorMap::const_iterator pos;
  pos = tensorParams_.find( matType );

  //2D tensor axi or plane is the same
  matMatrix.Resize(2,2);
  matMatrix.Init();
  pos->second.GetSubMatrix(matMatrix, 0, 0);

}


void ElectricConductionMaterial::ComputeFullMuTensor() {

  Matrix<Double> elecCond(3,3);
  Double scalarVal;

  // depending on symmetry, calculate full 3x3 tensor
  // check, if tensor was already set
  if( tensorParams_.find( ELEC_CONDUCTIVITY_TENSOR ) == tensorParams_.end() ) {

	  GetScalar( scalarVal, ELEC_CONDUCTIVITY, Global::REAL );
      elecCond[0][0] = scalarVal;
      elecCond[1][1] = scalarVal;
      elecCond[2][2] = scalarVal;
      SetTensor( elecCond, ELEC_CONDUCTIVITY_TENSOR, Global::REAL );
      //SetTensor( elecCond, ELEC_CONDUCTIVITY_TENSOR, Global::REAL );
  }
 }

PtrCoefFct ElectricConductionMaterial::GetScalCoefFncNonLin(MaterialType matType,
                                                           Global::ComplexPart matDataType,
                                                           PtrCoefFct fluxCoef ) {


  return BaseMaterial::GetScalCoefFncNonLin(matType, matDataType, fluxCoef);
}

PtrCoefFct ElectricConductionMaterial::GetScalCoefFncMultivariateNonLin(MaterialType matType,
                                                   NonLinType nlType,
                                                   Global::ComplexPart matDataType,
                                                   StdVector<PtrCoefFct> dependencies,
 					           StdVector<RegionIdType> & regs) {

  UInt ndeps = dependencies.GetSize(); // 1 no temp dep, 2: temp dep
  UInt nregs = regs.GetSize(); // 2 regs dipole, 3 regs tripole

  // Ensure that only real-valued parameters are used
  if( matDataType != Global::REAL ) {
    EXCEPTION( "Only real-valued nonlinear parameters are supported");
  }
  
  // check if isotropic material type is defined
  if( nonlinIsoParams_.find(matType) == nonlinIsoParams_.end() ) {
    EXCEPTION( "No nonlinear definition found for material type '"
        << MaterialTypeEnum.ToString(matType) << "'");
  }
  
  // check, if nonlinear curve was already calculated ??
  // using matType...
  MatDescriptorNl & matNl = nonlinIsoParams_[matType];

  // Check, if smooth spline approximation was already created 
  // and initialized
  if( !matNl.approxData ) {
    if ( matNl.approxType == SMOOTH_SPLINES ) {
      SmoothSpline * sp = new SmoothSpline( matNl.fileName, matType );
      sp->SetAccuracy( matNl.measAccuracy );
      sp->SetMaxY( matNl.maxVal );
      sp->CalcBestParameter();
      sp->CalcApproximation();
      sp->Print();
      matNl.approxData = sp;
    }
    else if ( matNl.approxType == LIN_INTERPOLATE ) {
      LinInterpolate * sp = new LinInterpolate( matNl.fileName, matType );
      //sp->SetAccuracy( matNl.measAccuracy );
      //sp->SetMaxY( matNl.maxVal );
      sp->Print();
      matNl.approxData = sp;
     }
     else if ( matNl.approxType == BILIN_INTERPOLATE ) {

      if ( ndeps != 2 )
        EXCEPTION("bilinear interpoltation requires two ptrCoefFcts");
      std::string val = stringParams_[NONLIN_DEPENDENCY];
      if  (  val == "voltage-temperature" ) {
        BiLinInterpolate * sp = new BiLinInterpolate( matNl.fileName, matType );
        //sp->SetAccuracy( matNl.measAccuracy );
        //sp->SetMaxY( matNl.maxVal );
        sp->Print();
        matNl.approxData = sp;
      } else {
        EXCEPTION("'dependency' string in material xml must be specified with value 'voltage-temperature' to force user check column ordering: | voltage | temperature | sigma |");
      }
    }
    else {
      EXCEPTION( "No nonlinear approx/interpolate type not available '"
          << ApproxCurveTypeEnum.ToString(matNl.approxType) << "'");
    }
  }

  ApproxData * sp = matNl.approxData;
  // get linear starting value
  Double startVal = 0.0;
  this->GetScalar( startVal, matType, Global::REAL );

  shared_ptr<CoefFunctionComposite> coef;
  if (nlType == NLELEC_BIPOLE) {
    assert ( ndeps == 1) ;
    assert ( nregs == 2) ;
    coef.reset(new CoefFunctionBipole());
    coef->SetRegion(ANODE,regs[0]);
    coef->SetRegion(CATHODE,regs[1]);
    coef->SetDependCoef(NLELEC_BIPOLE, dependencies[0] );
    coef->Init( startVal, sp);

  }
  else if (nlType == NLELEC_BIPOLE_TEMP_DEP) {
    assert ( ndeps == 2) ;
    assert ( nregs == 2) ;
    coef.reset(new CoefFunctionHeatBipole());
    coef->SetRegion(ANODE,regs[0]);
    coef->SetRegion(CATHODE,regs[1]);
    coef->SetDependCoef(NLELEC_BIPOLE, dependencies[0] );
    coef->SetDependCoef(NLELEC_CONDUCTIVITY, dependencies[1] ); // this means temperature dep
    coef->Init( startVal, sp);

  }
  else if (nlType == NLELEC_TRIPOLE) {
    assert ( ndeps == 1) ;
    assert ( nregs == 3) ;
    EXCEPTION("not implemented yet");
  }
  else if (nlType == NLELEC_TRIPOLE_TEMP_DEP) {
    assert ( ndeps == 2) ;
    assert ( nregs == 3) ;
    EXCEPTION("not implemented yet");
  }
  else {
    EXCEPTION("non linear type not recongized");
  }

    
  return coef;
}


}
