#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "ElectroMagneticMaterial.hh"

#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Preisach.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefXpr.hh"


namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElectroMagneticMaterial::ElectroMagneticMaterial(MathParser* mp,
                                                   CoordSystem * defaultCoosy)
  : BaseMaterial(mp, defaultCoosy) {

    materialDatabaseName_ = "Electromagnetics";

    //set the allowed material parameters
    isAllowed_.insert( MAG_PERMEABILITY );
    isAllowed_.insert( MAG_PERMEABILITY_1 );
    isAllowed_.insert( MAG_PERMEABILITY_2 );
    isAllowed_.insert( MAG_PERMEABILITY_3 );
    isAllowed_.insert( MAG_RELUCTIVITY );
    isAllowed_.insert( MAG_CONDUCTIVITY );
    isAllowed_.insert( ELEC_PERMITTIVITY );
    isAllowed_.insert( PREISACH_WEIGHTS );
    isAllowed_.insert( X_SATURATION );
    isAllowed_.insert( Y_SATURATION );
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
  }
  
  void ElectroMagneticMaterial::Finalize() {
    ComputeFullMuTensor();
  }

  void ElectroMagneticMaterial::SetScalar(const std::string& param, MaterialType matType) {


    if (( matType == HYST_MODEL ) || ( matType == MAG_RELUCTIVITY ) ||
        ( matType == MAG_RELUCTIVITY_DERIV )) {
      stringParams_[matType] = param;
      isSet_.insert( matType );
    }
    else {
      std::string dim = "string";
      matTypeNotAllowed( matType, dim );
    }
  }


  void ElectroMagneticMaterial::SetScalar( Double param, MaterialType matType, 
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


  void ElectroMagneticMaterial::SetScalar( Complex param, MaterialType matType, 
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


  void ElectroMagneticMaterial::SetTensor(const Matrix<Double>& param, MaterialType matType, 
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

  void ElectroMagneticMaterial::SetTensor(const Matrix<Complex>& param, MaterialType matType, 
                                          Global::ComplexPart dataType ) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
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

    if ( matType == ELEC_PERMITTIVITY ) {
      // to be consistent to old structure
      scalarParams_[matType] = param[2][2];
    }
  }


  void ElectroMagneticMaterial::GetScalar( Double& param, MaterialType matType, 
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

  void ElectroMagneticMaterial::GetScalar( Complex& param, MaterialType matType, 
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

  void ElectroMagneticMaterial::GetTensor( Matrix<Double>& param, 
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

  void ElectroMagneticMaterial::GetTensor( Matrix<Complex>& param, 
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
    if(matType==MAG_RELUCTIVITY){
       for(UInt i = 0; i<param.GetNumRows();i++){
         for(UInt j = 0; j<param.GetNumCols();j++){
           Complex tmp = param[i][j];
           param[i][j] = Complex(1.0,0.0) / tmp;
         }
       }
     }
  }


  void ElectroMagneticMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
                                                 MaterialType matType, 
                                                 SubTensorType subTensor) const {


    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor axi or plane is the same
    matMatrix.Resize(2,2);
    matMatrix.Init();
    pos->second.GetSubMatrix(matMatrix, 1, 1);
  }


//  void ElectroMagneticMaterial::InitApproxCurves() {
//
//    // check, if we need to approx BH curve
//    if (  needApproxMatCurves_.find( MAG_PERMEABILITY ) != needApproxMatCurves_.end() ) {
//      std::string nlfnc = GetNonlinFileName(MAG_PERMEABILITY);
//      nlinFncBH_ = new SmoothSpline(nlfnc, MAG_PERMEABILITY);
//
//      //get accuracy of approximation
//      Double dataAccuracy;
//      GetScalar( dataAccuracy, DATA_ACCURACY, Global::REAL );
//
//      //get maximal approximation value
//      Double maxApproxVal;
//      GetScalar( maxApproxVal, MAX_APPROX_VAL, Global::REAL );
//
//      nlinFncBH_->SetAccuracy( dataAccuracy );
//      nlinFncBH_->SetMaxY( maxApproxVal );   //maximal value of magnetic induction B
//      nlinFncBH_->CalcBestParameter();
//      nlinFncBH_->CalcApproximation();
//      nlinFncBH_->Print(); 
//    }
//  }

  void ElectroMagneticMaterial::InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                                          bool isInverse, bool computeHystInverse ) {

    isHystInverse_      = isInverse;
    computeHystInverse_ = computeHystInverse;

    std::string val = stringParams_[HYST_MODEL];
    if ( val != "preisach" ) {
      EXCEPTION( "Currently we just support Preisach Hysteresis Model" );
    }
    else {
      isHysteresis_ = true;
      dim_ = 2;

      // we use vector hysteresisi model
      hyst_ = NULL;

      std::cout << "computeHystInverse: " << computeHystInverse_  
          << " isHystInverse: " << isHystInverse_ << std::endl;  

      GetScalar(Xsat_, X_SATURATION, Global::REAL);
      GetScalar(Ysat_, Y_SATURATION, Global::REAL);
      Matrix<Double> weights;
      GetTensor(weights,  PREISACH_WEIGHTS, Global::REAL);
      bool isVirgin = true;

      hyst_ = new Preisach(numElemSD, Xsat_, Ysat_, weights, isVirgin);

      vecHyst_ = new Hysteresis* [dim_];
      vecHyst_[0] = new Preisach(numElemSD, Xsat_, Ysat_, weights, isVirgin);
      vecHyst_[1] = new Preisach(numElemSD, Xsat_, Ysat_, weights, isVirgin);

      // set map: global to local element number
      EntityIterator it = actSDList->GetIterator();
      UInt iel = 0;
      UInt globalElNr;
      for ( it.Begin(); !it.IsEnd(); it++, iel++) {
        globalElNr = it.GetElem()->elemNum;
        globalElem2Local_[globalElNr] = iel;
      }
    }

    //allocate memory for previous results, needed for the
    //effective material parameter formulation
    vecXprevious_.Resize(dim_,numElemSD);
    vecYprevious_.Resize(dim_,numElemSD);
    vecXprevious_.Init();
    vecYprevious_.Init();

    //
    Double startMatdiff;
    GetScalar( startMatdiff, MAG_RELUCTIVITY, Global::REAL );
    matDiffprevious_.Resize( numElemSD );
    matDiffprevious_.Init( startMatdiff );

    if (  computeHystInverse_ ) {
      //for inverse hysteresis
      vecXact_.Resize(dim_,numElemSD);
      vecYact_.Resize(dim_,numElemSD);
      vecXact_.Init();
      vecYact_.Init();
    }
  }


  void ElectroMagneticMaterial::SetPreviousHystVal( UInt nrElem, Vector<Double>& valVec) {

    UInt idx = globalElem2Local_[nrElem];

    if ( isHystInverse_ ) {
      vecYprevious_[0][idx] = valVec[0];
      vecYprevious_[1][idx] = valVec[1];
      vecXprevious_[0][idx] = vecHyst_[0]->computeValueAndUpdate( valVec[0], idx );
      vecXprevious_[1][idx] = vecHyst_[1]->computeValueAndUpdate( valVec[1], idx );
    }
    else if ( computeHystInverse_ ) {
      Vector<Double> newX(2), newY(2);
      ComputeInverseScalar( idx, 0, valVec[0], newX[0] );
      ComputeInverseScalar( idx, 1, valVec[1], newX[1]  );

      //! now perform also an updating
      newY[0] = vecHyst_[0]->computeValueAndUpdate( newX[0], idx);
      newY[1] = vecHyst_[1]->computeValueAndUpdate( newX[1], idx);

      Vector<Double> dX(2), dY(2);
      dX[0] = newX[0] - vecXprevious_[0][idx];
      dX[1] = newX[1] - vecXprevious_[1][idx];

      dY[0] = newY[0] - vecYprevious_[0][idx];
      dY[1] = newY[1] - vecYprevious_[1][idx];

      //Double dB = dY[0]*dY[0] + dY[1]*dY[1];
      Double dB = dY[1]*dY[1];
      if ( dB > 1e-5) {
        Double newMatDiff = ComputeMatDiff( dX, dY, idx );
        matDiffprevious_[idx] = newMatDiff;
      }

      vecXprevious_[0][idx] = newX[0];
      vecXprevious_[1][idx] = newX[1];
      vecYprevious_[0][idx] = newY[0];
      vecYprevious_[1][idx] = newY[1];

      vecXact_[0][idx] = vecXprevious_[0][idx];
      vecXact_[1][idx] = vecXprevious_[1][idx];
      vecYact_[0][idx] = vecYprevious_[0][idx];
      vecYact_[1][idx] = vecYprevious_[1][idx];
    }
    else {
      vecXprevious_[0][idx] = valVec[0];
      vecXprevious_[1][idx] = valVec[1];
      vecYprevious_[0][idx] = vecHyst_[0]->computeValueAndUpdate( valVec[0], idx );
      vecYprevious_[1][idx] = vecHyst_[1]->computeValueAndUpdate( valVec[1], idx );
    }
  }


  void ElectroMagneticMaterial::ComputeScalarDiffValues( UInt nrElem, 
                                                         Vector<Double>& valVec,
                                                         Vector<Double>& scalarValues ) {

    Vector<Double> Ycurrent(dim_);
    Vector<Double> Xcurrent(dim_);

    UInt idx = globalElem2Local_[nrElem];

    //   std::cout << "elNr=" << nrElem << "  idx=" << idx << std::endl;

    if ( isHystInverse_ ) {
      Ycurrent[0] = valVec[0];
      Ycurrent[1] = valVec[1];
      Xcurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
      Xcurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
      std::cout << "Ycurrent:\n" << Ycurrent << "\n Xcurrent: \n" << Ycurrent << std::endl;
    }
    else if ( computeHystInverse_ ) {
      ComputeInverseScalar( idx, 0, valVec[0], Xcurrent[0] );
      ComputeInverseScalar( idx, 1, valVec[1], Xcurrent[1] );
      Ycurrent[0] = vecHyst_[0]->computeValueAndUpdate(Xcurrent[0], idx);
      Ycurrent[1] = vecHyst_[1]->computeValueAndUpdate(Xcurrent[1], idx);
      //       Ycurrent[0] = valVec[0];
      //       Ycurrent[1] = valVec[1];
    }
    else {
      Xcurrent[0] = valVec[0];
      Xcurrent[1] = valVec[1];
      Ycurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
      Ycurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
    }

    //compute differential material parameter
    Vector<Double> dX(dim_), dY(dim_);
    dX[0] = Xcurrent[0] - vecXprevious_[0][idx];
    dX[1] = Xcurrent[1] - vecXprevious_[1][idx];

    dY[0] = Ycurrent[0] - vecYprevious_[0][idx];
    dY[1] = Ycurrent[1] - vecYprevious_[1][idx];

    std::cout << "dX \n" << dX << "  DdY:\n" << dY << std::endl;
    Double dB = dY[0]*dY[0] + dY[1]*dY[1];

    scalarValues.Init();
    if ( dB > 1e-10 ) {
      scalarValues[0] = ( dX[0] * dY[0] + dX[1] * dY[1] ) / dB ;
      scalarValues[1] = ( dX[0] * dY[1] + dX[1] * dY[0] ) / dB ;
    }
  }


  Double ElectroMagneticMaterial::ComputeScalarDiffVal( UInt nrElem, Vector<Double>& valVec) {

    // COMPWARNING: unused variable Double matDiff, eps;
    Vector<Double> Ycurrent(dim_);
    Vector<Double> Xcurrent(dim_);

    UInt idx = globalElem2Local_[nrElem];

    if ( isHystInverse_ ) {
      Ycurrent[0] = valVec[0];
      Ycurrent[1] = valVec[1];
      Xcurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
      Xcurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
    }
    else if ( computeHystInverse_ ) {
      ComputeInverseScalar( idx, 0, valVec[0], Xcurrent[0] );
      ComputeInverseScalar( idx, 1, valVec[1], Xcurrent[1] );
      Ycurrent[0] = valVec[0];
      Ycurrent[1] = valVec[1];
    }
    else {
      Xcurrent[0] = valVec[0];
      Xcurrent[1] = valVec[1];
      Ycurrent[0] = vecHyst_[0]->computeValueAndUpdate(valVec[0], idx);
      Ycurrent[1] = vecHyst_[1]->computeValueAndUpdate(valVec[1], idx);
    }

    //    std::cout << "Hx=" << Xcurrent[0] << "  Hy= " << Xcurrent[1] << std::endl;
    //    std::cout << "Bx=" << Ycurrent[0] << "  By= " << Ycurrent[1] << std::endl;
    //compute differential material parameter
    Vector<Double> dX(dim_), dY(dim_);
    dX[0] = Xcurrent[0] - vecXprevious_[0][idx];
    dX[1] = Xcurrent[1] - vecXprevious_[1][idx];

    dY[0] = Ycurrent[0] - vecYprevious_[0][idx];
    dY[1] = Ycurrent[1] - vecYprevious_[1][idx];
    Double newMatDiff = ComputeMatDiff( dX, dY, idx );

    return newMatDiff;
  }

  Double ElectroMagneticMaterial::ComputeMatDiff( Vector<Double>& dX, Vector<Double>& dY, UInt idx ) {

    Double matDiff, eps;
    Double dB = dY[0]*dY[0] + dY[1]*dY[1];
    //Double dB = dY[1]*dY[1];

    if ( dB < 1e-12 ) {
      matDiff = matDiffprevious_[idx];
      //matDiff = eps;
      //std::cout << "Startnu: " << matDiff << std::endl;
    }
    else {
      matDiff = ( dX[0] * dY[0] + dX[1] * dY[1] ) / dB;
      //matDiff = ( dX[1] * dY[1] ) / dB;
    }

    //    std::cout << "dB=" << dB  <<  "  dnu=" << matDiff << std::endl << std::endl;
    //    std::cout << "  dnu=" << matDiff << std::endl << std::endl;

    if ( matDiff < 50.0 ) {
      GetScalar(eps,MAG_RELUCTIVITY,Global::REAL);
      matDiff = eps;
    }

    return matDiff;
  }

  void ElectroMagneticMaterial::ComputeVectorHystVal( UInt nrElem, Vector<Double>& in, 
                                                      Vector<Double>& out ) {

    UInt idx = globalElem2Local_[nrElem];

    //    std::cout << "elNr=" << nrElem << "  idx=" << idx << std::endl;

    if ( isHystInverse_ ) {
      out[0] = vecHyst_[0]->computeValueAndUpdate( in[0], idx );
      out[1] = vecHyst_[1]->computeValueAndUpdate( in[1], idx );
    }
    else if ( computeHystInverse_ ) {
      ComputeInverseScalar( idx, 0, in[0], out[0] );
      ComputeInverseScalar( idx, 1, in[1], out[1] );
    }
    else {
      out[0] = vecHyst_[0]->computeValueAndUpdate( in[0], idx );
      out[1] = vecHyst_[1]->computeValueAndUpdate( in[1], idx );
    }
  }

  void ElectroMagneticMaterial::GetVectorHystVal( UInt nrElem, Vector<Double>& Val ) {

    UInt idx = globalElem2Local_[nrElem];

    if ( isHystInverse_ ) {
      Val[0] = vecHyst_[0]->getValue( idx );
      Val[1] = vecHyst_[1]->getValue( idx );
    }
    else if ( computeHystInverse_ ) {
      Val[0] = vecXprevious_[0][idx];
      Val[1] = vecXprevious_[1][idx];
    }
    else {
      Val[0] = vecHyst_[0]->getValue( idx );
      Val[1] = vecHyst_[1]->getValue( idx );
    }
  }

  Double ElectroMagneticMaterial::GetScalarHystVal( UInt nrElem ) {

    EXCEPTION("ElectroMagneticMaterial::GetScalarHystVal makes no sense");

    // COMPWARNING: unused variable UInt idx = globalElem2Local_[nrElem];
    Double Yval = 0.0; // = hyst_->getValue( idx );

    return Yval;
  }


  //====================================== INVERSE HYST =======================================


  void ElectroMagneticMaterial::ComputeInverseScalar( UInt idxEl, UInt comp, Double Yin, 
                                                      Double& Xout ) {


    Double eps = 1e-3;
    Double  dH = vecHyst_[0]->GetIncX();

    std::cout << "Yin= " << Yin << std::endl;

    if ( ( abs(Yin) + 0.05*Ysat_ ) > Ysat_ ) {
      Xout = Xsat_;
    }
    else {
      Double Hs, Ho, Hu, Hact, Bs, Bact, dB;
      bool found = false;

      //compute starting values
      Hs = vecXact_[comp][idxEl];
      Bs = vecHyst_[comp]->computeValueAndUpdate( Hs, idxEl, false); 

      std::cout << "Start Bs: " << Bs << "  Hs=" << Hs <<  std::endl;
      if  ( abs(Bs - Yin) < eps ) {
        found = true;
        Xout  = Hs;
        std::cout << "Direct found " << std::endl;
      }
      else if ( (Bs - Yin) > eps ) {
        Ho = Hs;
        Hu = Hs;
        std::cout << "Fix Ho= " << Ho << std::endl;
        do {
          Hu  -= dH;
          Bact = vecHyst_[comp]->computeValueAndUpdate( Hu, idxEl, false); 
          dB   = Bact - Yin; 
        } while ( dB > 0 ); 
      }
      else {
        Hu = Hs;
        Ho = Hs;
        std::cout << "Fix Hu=  " << Hu << std::endl;
        do {
          Ho  += dH;
          Bact = vecHyst_[comp]->computeValueAndUpdate( Ho, idxEl, false); 
          //  std::cout << "Compute Ho: " << Ho << "  Bact=" << Bact << std::endl;
          dB   = Bact - Yin;
        } while ( dB < 0 ); 
      }

      if ( found == false ) {
        std::cout << "Do iter: Bin=" << Yin << "  Bs=" << std::endl;
        do {
          Hact = ( Ho + Hu ) * 0.5;
          Bact = vecHyst_[comp]->computeValueAndUpdate( Hact, idxEl, false); 
          dB   = Bact - Yin;

          if ( dB < 0 ) 
            Hu = Hact;
          else 
            Ho = Hact;

          std::cout << "newB =" << Bact << "  Hact=" << Hact << "  Ho=" << Ho << "   Hu=" << Hu << std::endl; 

        } while ( abs(dB) > eps && abs(Ho-Hu) > abs(Ho)*1e-4 );

        Xout = Hact;
      }
    }

    vecXact_[comp][idxEl] = Xout;
    vecYact_[comp][idxEl] = Yin;

    // update
    //vecHyst_[comp]->updateMinMaxList( Xout, idxEl );

    //     if ( found ) 
    //       std::cout << " Hval = " << Xout << "  Bval=" << Yin << "   Bs=" << Bs << std::endl;
    //     else
    //       std::cout << " Hval = " << Xout << "  Bval=" << Yin << "  Bact=" << Bact << std::endl;
    //   }

  }
  
  void ElectroMagneticMaterial::ComputeFullMuTensor() {

    Matrix<Complex> muTensor(3,3);
    Complex mu1, mu2, mu3, isoMu;
    
    // depending on symmetry, calculate full 3x3 permeability tensor
    switch(GetSymmetryType(MAG_PERMEABILITY)) {
      
      case GENERAL:
        // in this case we have already the full permeability tensor
        
        //std::cout << "WTF?" << std::endl;
        GetTensor( muTensor, MAG_PERMEABILITY, Global::COMPLEX );
        break;
        
      case ISOTROPIC:
        GetScalar( isoMu, MAG_PERMEABILITY, Global::COMPLEX );
        muTensor[0][0] = isoMu;
        muTensor[1][1] = isoMu;
        muTensor[2][2] = isoMu;
        SetTensor( muTensor, MAG_PERMEABILITY, Global::COMPLEX );
        break;
        
      case ORTHOTROPIC:
        
        GetScalar( mu1, MAG_PERMEABILITY_1, Global::COMPLEX );
        GetScalar( mu2, MAG_PERMEABILITY_1, Global::COMPLEX );
        GetScalar( mu3, MAG_PERMEABILITY_1, Global::COMPLEX );
        muTensor[0][0] = mu1;
        muTensor[1][1] = mu2;
        muTensor[2][2] = mu3;
        SetTensor( muTensor, MAG_PERMEABILITY, Global::COMPLEX );
        break;
      default:
        EXCEPTION( "Calculation of full permeability matrix for symmetryType '"
            << GetSymmetryType(MAG_PERMEABILITY) << "' not implemented!" );
    }
    
    // Now we have the full mu-tensor, so we can invert the matrix
    // and store the reluctivity tensor
    Matrix<Double> nuTensor(3,3), temp;
    temp = muTensor.GetPart(Global::REAL);
    temp.Invert(nuTensor);
  
  /*  
    std::cout << "-----------" << std::endl;
    std::cout << "NU" << std::endl;
    std::cout << temp << std::endl;
    std::cout << "MU" << std::endl;   
    std::cout << nuTensor << std::endl;
    std::cout << "-----------" << std::endl;
   */
        
    SetTensor( nuTensor, MAG_RELUCTIVITY, Global::REAL );
    
    GetTensor( temp, MAG_RELUCTIVITY, Global::REAL );
  }

  PtrCoefFct ElectroMagneticMaterial::GetScalCoefFncNonLin(MaterialType matType,
                                                           Global::ComplexPart matDataType,
                                                           PtrCoefFct fluxCoef ) {
     //This method allocates the objects handling the nonlinear BH curve; thereby, we allow
     //approximation with smooth splines and analytically defined functions
     //
     //Please note: in the nonlinear bilinear form, we need the reluctivity (=1/permeability)
     //             therefore, we switch between permability and reluctivity quite often
     //             The analytic defined functions in the material file are
     //             reluctivity(magFluxDensity) = nu(B)
     //

     // Ensure that only MAG_RELUCTIVITY or MAG_RELUCTIVITY_DERIV are queried
     if( matType != MAG_RELUCTIVITY  ) {
       EXCEPTION("Scalar Nonlinearity for magnetic materials only allowed for MAG_RELUCTIVITY!");
     }
     
     // Ensure that only real-valued parameters are used
     if( matDataType != Global::REAL ) {
       EXCEPTION( "Only real-valued nonlinear parameters are supported");
     }
     PtrCoefFct ret;
     
     // check if material is isotropic or anisotropic
     if( nonlinIsoParams_.find(MAG_PERMEABILITY) != nonlinIsoParams_.end() ) {
       
       // ---------------------------
       // ISOTROPIC VERSION
       // ---------------------------
       // check, if nonlinear curve was already calculated
       MatDescriptorNl & matNl = nonlinIsoParams_[MAG_PERMEABILITY];

       //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
       if( matNl.approxType == SMOOTH_SPLINES ) {
         // Check, if smooth spline approximation was already created
         // and initialized
         if( !matNl.approxData ) {
           SmoothSpline * sp = new SmoothSpline( matNl.fileName, MAG_PERMEABILITY );
           sp->SetAccuracy( matNl.measAccuracy );
           sp->SetMaxY( matNl.maxVal );
           sp->CalcBestParameter();
           sp->CalcApproximation();
           sp->Print();
           matNl.approxData = sp;
         }

         ApproxData * sp = matNl.approxData;
         // get linear starting value
         Double startVal = 0.0;
         this->GetScalar( startVal, matType, Global::REAL );
         shared_ptr<CoefFunctionApprox> coef( new CoefFunctionApprox());
         coef->Init( startVal, sp, fluxCoef);
         ret = coef;

       }
       else if( matNl.approxType == ANALYTIC ) {
         // this is for describing the reluctivity directly in the xml as analytic formula
         // idea: the string from the xml describes a function with the same notation as
         // described in CoefFunctionCompound.hh
         // basically, all occurences of B_R are replaced with the CoefFunction fluxDensAbs
         // note: a good starting value for B->0 works miracles!

         // get Euclidean norm of B
         CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, fluxCoef, CoefXpr::OP_NORM );
         PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

         // get function of B
         std::string nuStr = stringParams_[matType];
         shared_ptr<CoefFunctionCompound<Double> > nuFnc(new CoefFunctionCompound<Double>(mp_));
         std::map<std::string,PtrCoefFct> symbolsNu;
         symbolsNu["B"] = fluxDensAbs;

         nuStr.insert(0,"( ");
         nuStr.append(" )");
         nuFnc->SetScalar(nuStr,symbolsNu);
         return(nuFnc);
       }

     } else if( nonlinAnisoParams_.find(MAG_PERMEABILITY) != nonlinAnisoParams_.end() ) {
       
       // ---------------------------
       // ANISOTROPIC VERSION: here we allow for different BH-curves as a function of the angle!
       // ---------------------------
       StdVector<MatDescriptorNl> & matNl = nonlinAnisoParams_[MAG_PERMEABILITY];
       UInt numCurves = matNl.GetSize();
       StdVector<Double> angles(numCurves);
       StdVector<Double> zScalings(numCurves);
       StdVector<shared_ptr<CoefFunction> > approx(numCurves);
       Double startValAveraged = 0.0;

       // Loop over all entries
       for( UInt i = 0; i < matNl.GetSize(); ++i ) {
         MatDescriptorNl & actNl = matNl[i];
         angles[i] = actNl.angle;
         zScalings[i] = actNl.zScaling;

         //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
         if( actNl.approxType == SMOOTH_SPLINES ) {
           // Check, if smooth spline approximation was already created
           // and initialized
           if( !actNl.approxData ) {
             SmoothSpline * sp = new SmoothSpline( actNl.fileName, MAG_PERMEABILITY );
             sp->SetAccuracy( actNl.measAccuracy );
             sp->SetMaxY( actNl.maxVal );
             sp->CalcBestParameter();
             sp->CalcApproximation();
             sp->Print();
             actNl.approxData = sp;
           }

           ApproxData * sp = actNl.approxData;
           // get linear starting value

           Double startVal;
           this->GetScalar( startVal, matType, Global::REAL );
           shared_ptr<CoefFunctionApprox> coef( new CoefFunctionApprox());
           coef->Init( startVal, sp, fluxCoef);

           //compute an averaged starting value
           startValAveraged += startVal / (Double)numCurves;

           //store in array
           approx[i] = coef;
         }
         else if( actNl.approxType == ANALYTIC ) {
           // this is for describing the reluctivity directly in the xml as analytic formula
           // idea: the string from the xml describes a function with the same notation as
           // described in CoefFunctionCompound.hh
           // basically, all occurences of B_R are replaced with the CoefFunction fluxDensAbs
           // note: a good starting value for B->0 works miracles!

           // get Euclidean norm of B
           CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, fluxCoef, CoefXpr::OP_NORM );
           PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

           // get function of B
           std::string nuStr = actNl.analyticExpr;
           shared_ptr<CoefFunctionCompound<Double> > nuFnc(new CoefFunctionCompound<Double>(mp_));
           std::map<std::string,PtrCoefFct> symbolsNu;
           symbolsNu["B"] = fluxDensAbs;

           nuStr.insert(0,"( ");
           nuStr.append(" )");
           nuFnc->SetScalar(nuStr,symbolsNu);

           //compute an averaged starting value directly from the string
           Double B_init = 0.0;
           MathParser::HandleType handle = mp_->GetNewHandle();
           mp_->RegisterExternalVar(handle,"B_R",&B_init);
           mp_->SetExpr(handle,nuStr);
           Double nuInit = mp_->Eval(handle);
           startValAveraged += nuInit / (Double)numCurves;

           //store in array
           approx[i] = nuFnc;
         }
       }
       
       // -------------------------
       // Insertion sort algorithm: we sort the BH-curves starting at smallest
       //                           specified angle
       // ------------------------
       Double compAngle;
       Double compZScaling;
       shared_ptr<CoefFunction> compApprox;
       UInt j;
       for( UInt i = 1; i < numCurves; i++ ) {
         compAngle = angles[i];
         compZScaling = zScalings[i];
         compApprox = approx[i];
         j = i;
         while( ( j > 0 ) && ( angles[j - 1] > compAngle ) ) {
           angles[j] = angles[j - 1];
           zScalings[j] = zScalings[j - 1];
           approx[j] = approx[j - 1];
           j = j - 1;
         }
         angles[j] = compAngle;
         zScalings[j] = compZScaling;
         approx[j] = compApprox;
       }

       // allocate the coef-Function for handling the ansiotropy
       shared_ptr<CoefFunctionApproxAniso> coef( new CoefFunctionApproxAniso());
       coef->Init( startValAveraged, approx, angles, zScalings, fluxCoef );
       ret = coef;
     }

     else {
       EXCEPTION( "No nonlinear definition found for material type '"
           << MaterialTypeEnum.ToString(matType) << "'");
     }

     
     return ret;
   }
  
  
  PtrCoefFct ElectroMagneticMaterial::GetTensorCoefFncNonLin( MaterialType matType,
                                                              SubTensorType type,
                                                              Global::ComplexPart matDataType,
                                                              PtrCoefFct dependency ) {
    //
    //This method allocates the objects handling the derivative of the reluctivity w.r.t.
    //the magnetic flux density ( nu'(B) ); therefore it is called to bulid up the nonlinear
    //bilinear form for the tangential stiffness matrix
    //
    //Please note: in the nonlinear bilinear form, we need the derivative of the
    //             reluctivity (=1/permeability); therefore, we switch between
    //             permability and reluctivity quite often
    //             The analytic defined functions in the material file diretcly define
    //             the derivative of the reluctivity as a function of mag. flux density!
    //

       // Ensure that only MAG_RELUCTIVITY or MAG_RELUCTIVITY_DERIV are queried
       if( matType != MAG_RELUCTIVITY && matType != MAG_RELUCTIVITY_DERIV ) {
         EXCEPTION("Nonlinearity for magnetic materials only allowed for MAG_RELUCTIVITY "
             << "or MAG_RELUCTIVITY_DERIV" );
       }
       
       // Ensure that only real-valued parameters are used
       if( matDataType != Global::REAL ) {
         EXCEPTION( "Only real-valued nonlinear parameters are supported");
       }
       PtrCoefFct ret;
       
       UInt dimDMat = (type == FULL) ? 3 : 2;
       
       // check if material is isotropic or anisotropic
       if( nonlinIsoParams_.find(MAG_PERMEABILITY) != nonlinIsoParams_.end() ) {
         
         // Check, if MAG_RELUCTIVITY is queried
         if( matType == MAG_RELUCTIVITY ) {
           EXCEPTION("An isotropic nonlinear MAG_RELUCTIVITY must be queried using "\
                     "GetScalCoefFncNonLin");
         }
         
         // ---------------------------
         // ISOTROPIC VERSION
         // ---------------------------
         // check, if nonlinear curve was already calculated
         MatDescriptorNl & matNl = nonlinIsoParams_[MAG_PERMEABILITY];

         if( matNl.approxType == SMOOTH_SPLINES ) {

           //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
           if( !matNl.approxData ) {
             SmoothSpline * sp = new SmoothSpline( matNl.fileName, MAG_PERMEABILITY );
             sp->SetAccuracy( matNl.measAccuracy );
             sp->SetMaxY( matNl.maxVal );
             sp->CalcBestParameter();
             sp->CalcApproximation();
             sp->Print();
             matNl.approxData = sp;
           }

           //now we need the object "CoefFunctionApproxDeriv", which returns by
           //calling "B^T [ e_B^T * nu' * |B| * e_B] B", so it is a tensor!!
           ApproxData * sp = matNl.approxData;
           shared_ptr<CoefFunctionApproxDeriv> coef( new CoefFunctionApproxDeriv());
           coef->Init( sp, dimDMat, dependency );
           ret = coef;

         } else if( matNl.approxType == ANALYTIC ) {
           //Here, we obtain " nu' " be evaluating the analytical defined function in
           //the material file, and then we have to build the tensor due to
           // "B^T [ e_B^T * nu' * |B| * e_B] B"

           // get Euclidean norm of B
           CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, dependency, CoefXpr::OP_NORM );
           PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

           // get function of B
           std::string dnudBStr = stringParams_[matType];
           shared_ptr<CoefFunctionCompound<Double> > scalFnc(new CoefFunctionCompound<Double>(mp_));
           std::map<std::string,PtrCoefFct> symbolsNu;
           symbolsNu["B"] = fluxDensAbs;

           // avoid divisions by zero, for very small B a fixed-step iteration is performed
           dnudBStr.insert(0,"( ( B_R lt 1e-6 ) ? ( 0.0 ) : ( ( ");
           dnudBStr.append(" ) / B_R ) )");
           scalFnc->SetScalar(dnudBStr,symbolsNu);

           shared_ptr<CoefFunctionCompound<Double> > dnudBTens(new CoefFunctionCompound<Double>(mp_));
           std::map<std::string,PtrCoefFct> symbolsTens;
           symbolsTens["a"] = dependency;
           StdVector<std::string> bbStr;
           if ( dimDMat == 3 ) {
             bbStr = "( a_0_R * a_0_R )" , "( a_0_R * a_1_R )" , "( a_0_R * a_2_R )" ,
                     "( a_0_R * a_1_R )" , "( a_1_R * a_1_R )" , "( a_1_R * a_2_R )" ,
                     "( a_0_R * a_2_R )" , "( a_1_R * a_2_R )" , "( a_2_R * a_2_R )" ;
           } else {
             bbStr = "( a_0_R * a_0_R )" , "( a_0_R * a_1_R )" ,
                     "( a_0_R * a_1_R )" , "( a_1_R * a_1_R )" ;
           }
           dnudBTens->SetTensor(bbStr,dimDMat,dimDMat,symbolsTens);

           CoefXprTensScalOp dnudBOp = CoefXprTensScalOp( mp_, dnudBTens, scalFnc, CoefXpr::OP_MULT );
           PtrCoefFct dnudBFnc = CoefFunction::Generate( mp_, Global::REAL, dnudBOp );
           return(dnudBFnc);
         }

       } else if( nonlinAnisoParams_.find(MAG_PERMEABILITY) != nonlinAnisoParams_.end() ) {
         // ---------------------------
         // ANISOTROPIC VERSION
         // ---------------------------
         StdVector<MatDescriptorNl> & matNl = nonlinAnisoParams_[MAG_PERMEABILITY];

         UInt numCurves = matNl.GetSize();
         StdVector<Double> angles(numCurves);
         StdVector<Double> zScalings(numCurves);
         StdVector<shared_ptr<CoefFunction> > approx(numCurves);
         // Loop over all entries
         for( UInt i = 0; i < matNl.GetSize(); ++i ) {
           MatDescriptorNl & actNl = matNl[i];
           angles[i] = actNl.angle;
           zScalings[i] = actNl.zScaling;

           //Here we really approximate H(B); see book Kaltenbacher, 2nd, 125ff
           if( actNl.approxType == SMOOTH_SPLINES ) {
             // Check, if smooth spline approximation was already created
             // and initialized
             if( actNl.approxData == NULL) {
               SmoothSpline * sp = new SmoothSpline( actNl.fileName, MAG_PERMEABILITY );
               sp->SetAccuracy( actNl.measAccuracy );
               sp->SetMaxY( actNl.maxVal );
               sp->CalcBestParameter();
               sp->CalcApproximation();
               sp->Print();
               actNl.approxData = sp;
             }
             //now we need the object "CoefFunctionApproxDeriv", which returns by
             //calling coef->getScalar( nuPrime, lmp) the derivative of the reluctivity;
             //Please note: In this case, we do not need the tensor (as in the isotropic case),
             //since the method "CoefFunctionApproxDerivAniso" (see below) will do the job;
             //That's why the object "CoefFunctionApproxDeriv" has the method "GetScalar"!
             //
             ApproxData * sp = actNl.approxData;
             shared_ptr<CoefFunctionApproxDeriv> coef( new CoefFunctionApproxDeriv());
             coef->Init( sp, dimDMat, dependency );

             approx[i] = coef;
           }
           else if( actNl.approxType == ANALYTIC ) {
             //Get the analytic expression for nu'(B)
             CoefXprUnaryOp fluxDensAbsOp = CoefXprUnaryOp( mp_, dependency, CoefXpr::OP_NORM );
             PtrCoefFct fluxDensAbs = CoefFunction::Generate( mp_, Global::REAL, fluxDensAbsOp );

             // get function of B
             std::string nuStr = actNl.analyticExprDeriv; //stringParams_[matType];
             shared_ptr<CoefFunctionCompound<Double> > nuFncDeriv(new CoefFunctionCompound<Double>(mp_));
             std::map<std::string,PtrCoefFct> symbolsNu;
             symbolsNu["B"] = fluxDensAbs;

             nuStr.insert(0,"( ");
             nuStr.append(" )");
             nuFncDeriv->SetScalar(nuStr,symbolsNu);

             //store in array
             approx[i] = nuFncDeriv;
           }
         }

         // -------------------------
         // Insertion sort algorithm: we sort the BH-curves starting at smallest
         //                           specified angle
         // ------------------------
         Double compAngle;
         Double compZScaling;
         shared_ptr<CoefFunction> compApprox;
         UInt j;
         for( UInt i = 1; i < numCurves; i++ ) {
           compAngle = angles[i];
           compZScaling = zScalings[i];
           compApprox = approx[i];
           j = i;
           while( ( j > 0 ) && ( angles[j - 1] > compAngle ) ) {
             angles[j] = angles[j - 1];
             zScalings[j] = zScalings[j - 1];
             approx[j] = approx[j - 1];
             j = j - 1;
           }
           angles[j] = compAngle;
           zScalings[j] = compZScaling;
           approx[j] = compApprox;
         }
         
         if( matType == MAG_RELUCTIVITY ) {
           // get linear starting value
           Double startVal = 0.0;
           this->GetScalar( startVal, matType, Global::REAL );
           shared_ptr<CoefFunctionApproxAniso> coef( new CoefFunctionApproxAniso());
           coef->Init( startVal, approx, angles, zScalings, dependency );
           ret = coef;
         }
         else if (matType == MAG_RELUCTIVITY_DERIV ) {
           //used for the bilinear form of the Newton method
           shared_ptr<CoefFunctionApproxDerivAniso> coef( new CoefFunctionApproxDerivAniso());
           coef->Init( approx, angles, zScalings, dimDMat, dependency );
           ret = coef;
         }

       } else {
         EXCEPTION( "No nonlinear definition found for material type '"
             << MaterialTypeEnum.ToString(matType) << "'");
       }

       return ret;
     }

}
