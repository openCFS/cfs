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
#include "Utils/StdVector.hh"
#include "MatVec/matrix.hh"
#include "Utils/coordSystem.hh"
#include "Utils/preisach.hh"
#include "Utils/simplePreisachInv.hh"
#include "Domain/entityList.hh"
#include "Utils/piezoMicroModel.hh"
#include "Utils/piezoMicroModelBK.hh"
#include "Domain/domain.hh"
#include "baseMaterial.hh"

using std::string;
using std::map;
using std::set;


namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  BaseMaterial::BaseMaterial() {


    materialDatabaseName_ = "";
    matFileName_ = "";
    nonlinFileName_ = "";

    coosy_ = NULL;
    hyst_  = NULL;

    symmetryType_  = GENERAL;
    isHysteresis_  = false;
    isHystInverse_ = false;

    piezoMicroModel_   = NULL;
    isPiezoMicroModel_ = false;
    
    mp_ = domain->GetMathParser();
  }

   BaseMaterial::~BaseMaterial() {


  }


  void BaseMaterial::NeedApproxMatCurve( ApproxMaterialCurves type ) {


    needApproxMatCurves_.insert( type );
  }


  void BaseMaterial::SetScalar(const std::string& param, MaterialType matType) {
        
      //check, if allowed
      if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
        std::string dim = "scalar";
        matTypeNotAllowed( matType, dim );
      }
      else {
        isSet_.insert( matType );
      }

      stringParams_[matType] = param;
      
    }
  
  void BaseMaterial::SetScalar(int param, MaterialType matType) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      integerParams_[matType] = param;
    }
  }
  void BaseMaterial::GetScalar( std::string& param, MaterialType matType)  const {


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

  
                 
  void BaseMaterial::GetScalar(Integer& param, MaterialType matType, Global::ComplexPart dataType ) const 
  {


    integerMap::const_iterator pos;
    pos = integerParams_.find( matType );

    if ( pos == integerParams_.end() ) {
      string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      if ( dataType == Global::INTEGER ) {	
        string msg = "GetScalar-Integer";
        dataTypeNotAllowed4SetGet( dataType, msg );
      }
      Integer val = pos->second;
      param = val;
    }
  }


  bool BaseMaterial::IsSet( MaterialType matType ) const {

    bool found = false;
    stringMap::const_iterator stringIt = stringParams_.find( matType );
    scalarMap::const_iterator scalarIt = scalarParams_.find( matType );
    tensorMap::const_iterator tensorIt = tensorParams_.find( matType );

    if( stringIt != stringParams_.end() ||
        scalarIt != scalarParams_.end() || 
        tensorIt != tensorParams_.end() ) {
      found = true;
    }

    return found;
  }

  void  BaseMaterial::matTypeNotAllowed(MaterialType matType, const string& dim ) const {

    string help = MaterialTypeEnum.ToString( matType );
    EXCEPTION( "Material type (" <<  dim <<  ") " << help << 
               " is not available for " << materialDatabaseName_ 
               << " Database" );
  }
  

  void  BaseMaterial::dataTypeNotAllowed4SetGet(Global::ComplexPart dataType, 
						 const string& msg ) const {

    string help;
    help = Global::complexPart.ToString( dataType );
    EXCEPTION( "Datatype " << help << " is not allowed in function " 
               << msg );
  }


  void BaseMaterial::dataTypeNotAllowed(Global::ComplexPart dataType, MaterialType matType ) const {

    string help1, help2;
    help1 = Global::complexPart.ToString( dataType );
    help2 = MaterialTypeEnum.ToString( matType );
    EXCEPTION( "Datatype " << help1 << " is not allowed for material type " 
               << help2 << " in material data base " << materialDatabaseName_ );
  }

  void BaseMaterial::matTypeNotInDataBase(MaterialType matType, const string& dim ) const {

    string help = MaterialTypeEnum.ToString( matType );
    EXCEPTION( "Material type (" << dim << ") " << help 
               << " was not read form/defined in material file" );
  }


  void  BaseMaterial::setMakesNoSense(Global::ComplexPart dataType, const string& msg ) const {
    EXCEPTION( "Set of " << msg << " makes no sense with datatype "
               << Global::complexPart.ToString( dataType ) );
  }


  void BaseMaterial::subTensorNotAvailable(MaterialType matType, SubTensorType subTensor )
  {
    string help1, help2;
    help1 = MaterialTypeEnum.ToString( matType );
    Enum2String(subTensor, help2);
    EXCEPTION("Subtensor " << help2 <<" not available for material type " << help1);
  }


  void BaseMaterial::ToInfo(PtrParamNode in)
  {
    set<MaterialType>::iterator iter;

    // in isSet the actually read material parameters are stored
    for(iter = isSet_.begin(); iter != isSet_.end(); iter++)
    {
      MaterialType mt = *iter;
      //check, if parameter is complex
      bool isComplex = isComplex_.find(mt) != isComplex_.end();

      //check for material data
      map<MaterialType, Matrix<Complex> >::const_iterator posTens = tensorParams_.find(mt);
      map<MaterialType, Complex >::const_iterator         posScal = scalarParams_.find(mt);
      map<MaterialType, string >::const_iterator          posStr  = stringParams_.find(mt);
      map<MaterialType, Integer >::const_iterator         posInt  = integerParams_.find(mt);

      PtrParamNode in_ = in->Get("property", ParamNode::APPEND);
      in_->Get("name")->SetValue(MaterialTypeEnum.ToString(mt));

      if(posTens != tensorParams_.end())
      {
        // get the tensor by default as complex
        if(isComplex)
          in_->Get("tensor")->SetValue(posTens->second);
        else
        {
          in_->Get("tensor")->SetValue(posTens->second.GetPart(Global::REAL));
        }
      }

      if(posScal != scalarParams_.end())
      {
        Complex val = posScal->second; // see tensor
        if(isComplex)
          in_->Get("value")->SetValue(val);
        else
          in_->Get("value")->SetValue(val.real());
      }

      if(posStr != stringParams_.end())
        in_->Get("value")->SetValue(posStr->second);

      if(posInt != integerParams_.end())
        in_->Get("value")->SetValue(posInt->second);
    }
  }


  void BaseMaterial::RotateTensorByRotationAngles( const Vector<Double> &rotAngle,
                                                   MaterialType matType,
                                                   bool persistent) {
    // Calculate rotation matrix( based on Kardan-Angles)
    // Ref.: C. Woernle, "Skript: Dynamik von Mehrkörpersystemen,
    // Kapitel 2 "Grundlagen der Kinematik", S. 12, Univ. Rostock
    // http://iamserver.fms.uni-rostock.de/studium/mehrkoerpersysteme/unterlagen.htm

    Double alpha = rotAngle[0] * PI / 180.0;
    Double beta  = rotAngle[1] * PI / 180.0;;
    Double gamma = rotAngle[2] * PI / 180.0;;
    Matrix<Double> R(3,3);
    R.Resize(3,3);

    R[0][0] =  cos(beta) * cos(gamma);
    R[0][1] = -cos(beta) * sin(gamma);
    R[0][2] =  sin(beta);
    R[1][0] =  cos(alpha)*sin(gamma) + sin(alpha)*sin(beta)*cos(gamma);
    R[1][1] =  cos(alpha)*cos(gamma) - sin(alpha)*sin(beta)*sin(gamma);
    R[1][2] = -sin(alpha)*cos(beta);
    R[2][0] =  sin(alpha)*sin(gamma) - cos(alpha)*sin(beta)*cos(gamma);
    R[2][1] =  sin(alpha)*cos(gamma) + cos(alpha)*sin(beta)*sin(gamma);
    R[2][2] =  cos(alpha)*cos(beta);
   
    RotateTensorByRotationMatrix( R, matType, persistent );
  }
  
  void BaseMaterial::RotateAllTensorsByRotationAngles( const Vector<Double>& rotAngle, 
                                                         bool persistent ) {
    BaseMaterial::tensorMap::iterator it = tensorParams_.begin();
    for( ; it != tensorParams_.end(); it++ ) {
      RotateTensorByRotationAngles( rotAngle, it->first, persistent );
    }
  }
  
  void BaseMaterial::RotateTensorByRotationMatrix( const Matrix<Double>& rotMatrix, 
                                                   MaterialType matType,
                                                   bool persistent ) {


    using namespace std;

    tensorMap::const_iterator pos;
    pos = this->tensorParams_.find( matType );

    tensorMap::const_iterator posOrig;
    posOrig = this->tensorParamsOrig_.find( matType );

    if ( pos == tensorParams_.end() ) {
      string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      if ( posOrig == tensorParamsOrig_.end() ) {
        string dim = "tensorOriginal";
        matTypeNotInDataBase( matType, dim );
      }

      Matrix<Complex> matTensor;
      matTensor = pos->second;
      Matrix<Complex> const matTensorOrig = posOrig->second;

      //compute complex rotation matrix 
      Matrix<Complex> RComplex (3,3);
      Matrix<Complex> helpTensor = matTensorOrig;
      RComplex.Init();
      RComplex.SetPart( Global::REAL, rotMatrix );
      PerformRotation(RComplex, matTensor, helpTensor); 
      tensorParams_[matType] = matTensor;
      
      // In the case of a persistent rotation, we override the
      // original tensor
      if ( persistent ) {
        tensorParamsOrig_[matType] = matTensor;
      }
    } // if parameter found
    
  }

    

  void BaseMaterial::RotateTensorByPointCoord( const Vector<Double>&  coord,
                                               MaterialType matType ) {


    // Determine rotation angles from attached coordinate system
    Matrix<Double> rotMatrix;
    coosy_->GetFullGlobRotationMatrix( rotMatrix, coord );

    // Calculate rotation. In this case this is always
    // non-persistent, as the orientation may change for each point in the
    // related coordinate system
    RotateTensorByRotationMatrix( rotMatrix, matType, false );

  }


  void BaseMaterial::PerformRotation( Matrix<Complex>& R,  
                                      Matrix<Complex>& matTensor,
                                      const Matrix<Complex>& matTensorOrig) {

    // get memory for transposed rotation matrix
    Matrix<Complex> RT;
    RT.Resize(3,3);
    R.Transpose(RT);

    //get dimension of matrix
    UInt rowSize = matTensorOrig.GetNumRows();
    UInt colSize = matTensorOrig.GetNumCols();

    Matrix<Complex> helpMat;

    if ( rowSize == 3 && colSize == 3) {
      // tensor is a 3x3 matrix: sol = R * matrixOrig * RT
      helpMat   = matTensorOrig * RT;
      matTensor = R * helpMat;
    }
    else {
      // we also need Q;
      Matrix<Complex> Q;

      // Composed Rotation Matrix
      // Ref.: M.Richter, "Entwicklung mechanischer Modelle zur analytischen
      // Beschreibung der Materialeigenschaften von textilbewehrtem Feinbeton",
      // Diss., Dresden, 2005, p. 27

      Q.Resize(6,6);  

      Q[0][0] = R[0][0]*R[0][0];
      Q[0][1] = R[0][1]*R[0][1];
      Q[0][2] = R[0][2]*R[0][2];
      Q[0][3] = 2.0*R[0][1]*R[0][2];
      Q[0][4] = 2.0*R[0][0]*R[0][2];
      Q[0][5] = 2.0*R[0][0]*R[0][1];

      Q[1][0] = R[1][0]*R[1][0];
      Q[1][1] = R[1][1]*R[1][1];
      Q[1][2] = R[1][2]*R[1][2];
      Q[1][3] = 2.0*R[1][1]*R[1][2];
      Q[1][4] = 2.0*R[1][0]*R[1][2];
      Q[1][5] = 2.0*R[1][0]*R[1][1];

      Q[2][0] = R[2][0]*R[2][0];
      Q[2][1] = R[2][1]*R[2][1];
      Q[2][2] = R[2][2]*R[2][2];
      Q[2][3] = 2.0*R[2][1]*R[2][2];
      Q[2][4] = 2.0*R[2][0]*R[2][2];
      Q[2][5] = 2.0*R[2][0]*R[2][1];

      Q[3][0] = R[1][0]*R[2][0];
      Q[3][1] = R[1][1]*R[2][1];
      Q[3][2] = R[1][2]*R[2][2];
      Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
      Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
      Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];

      Q[4][0] = R[0][0]*R[2][0];
      Q[4][1] = R[0][1]*R[2][1];
      Q[4][2] = R[0][2]*R[2][2];
      Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
      Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
      Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];

      Q[5][0] = R[0][0]*R[1][0];
      Q[5][1] = R[0][1]*R[1][1];
      Q[5][2] = R[0][2]*R[1][2];
      Q[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
      Q[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0];
      Q[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0];


      // 	std::cout << "R:\n" << R << std::endl;
      // 	std::cout << "Q:\n" << Q << std::endl;
      // 	std::cout << "Tensor orig:\n" << matTensor << std::endl;

      Matrix<Complex> QT;
      QT.Resize(6,6);
      Q.Transpose(QT);

      if ( rowSize == 3 && colSize == 6 ) {
        helpMat   = matTensorOrig * QT;
        matTensor = R * helpMat;
      }
      else if (rowSize == 6 && colSize == 6 ) {
        helpMat   = matTensorOrig * QT;
        matTensor = Q * helpMat;
      }
      // 	else {
      // 	  EXCEPTION("Cannot rotate tensor due to dimensions!");
      // 	}
    }
  }


  void BaseMaterial::InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                               bool isInverse, bool computeHystInverse ) {

    isHystInverse_      = isInverse;
    computeHystInverse_ = computeHystInverse;

    string val = stringParams_[HYST_MODEL];
    if ( val != "preisach" ) {
      EXCEPTION( "Currently we just support Preisach Hysteresis Model" );
    }
    else {
      isHysteresis_ = true;

      Double Xsat, Ysat;
      GetScalar(Xsat, X_SATURATION, Global::REAL);
      GetScalar(Ysat, Y_SATURATION, Global::REAL);
      Matrix<Double> weights;
      GetTensor(weights,  PREISACH_WEIGHTS, Global::REAL);
      bool isVirgin = true;   
      hyst_ = new Preisach(numElemSD, Xsat, Ysat, weights, isVirgin);

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
    Xprevious_.Resize(numElemSD);
    Yprevious_.Resize(numElemSD);
    Xprevious_.Init();
    Yprevious_.Init();
  }


  void BaseMaterial::InitVecHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                                  UInt dim ) {


    string val = stringParams_[HYST_MODEL];
    if ( val != "preisach" ) {
      EXCEPTION( "Currently we just support Preisach Hysteresis Model" );
    }
    else {
      isHysteresis_ = true;

      Double Xsat, Ysat;
      GetScalar(Xsat, X_SATURATION, Global::REAL);
      GetScalar(Ysat, Y_SATURATION, Global::REAL);
      Matrix<Double> weights;
      GetTensor(weights,  PREISACH_WEIGHTS, Global::REAL);
      bool isVirgin = true; 

      dimVecHyst_ = dim;
      hyst_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
      if ( dim == 2 ) 
        hystY_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
      else if ( dim == 3) {
        hystY_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
        hystZ_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
      }

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
    vecXprevious_.Resize(dim,numElemSD);
    vecYprevious_.Resize(dim,numElemSD);
    vecXprevious_.Init();
    vecYprevious_.Init();

    actDiffVal4VecHyst_.Resize(dim,numElemSD);
    previousDiffVal4VecHyst_.Resize(dim,numElemSD);
    Double nu;
    GetScalar(nu,MAG_RELUCTIVITY,Global::REAL);
    for (UInt i=0; i<dim; i++)
      for (UInt j=0; j<numElemSD; j++)
        previousDiffVal4VecHyst_[i][j] = nu;

  }


  Double BaseMaterial::GetScalarHystVal( UInt nrElem ) {

    UInt idx = globalElem2Local_[nrElem];
    return hyst_->getValue( idx );
  }


  Double BaseMaterial::GetScalarHystPrevVal( UInt nrElem ) {
    UInt idx = globalElem2Local_[nrElem];
    return Yprevious_[idx];
  }

  void BaseMaterial::ComputeRayleighDamping(string& alpha, string& beta,
                                            Double dampFreq, Double ratioDeltaF,
                                            bool adjustDamping, bool isHarmonic ) {

    // First, calculate ALPHA and BETA for the current dampFreq
    Double measuredFreq;
    std::string alphaOrig, betaOrig;
    GetScalar( measuredFreq, RAYLEIGH_FREQUENCY, Global::REAL ); 
    
    if ( IsSet( RAYLEIGH_ALPHA ) 
         && IsSet( RAYLEIGH_BETA ) 
         && IsSet(RAYLEIGH_FREQUENCY) ) {
      GetScalar( alphaOrig, RAYLEIGH_ALPHA ); 
      GetScalar( betaOrig, RAYLEIGH_BETA); 
      GetScalar( measuredFreq, RAYLEIGH_FREQUENCY, Global::REAL ); 

      if( abs(measuredFreq-dampFreq) > 0.001*measuredFreq ){
        alphaOrig = "(" + alphaOrig + "*" + lexical_cast<std::string>(dampFreq/measuredFreq) + ")";
        betaOrig= "(" + betaOrig + "*" + lexical_cast<std::string>(measuredFreq/dampFreq)+ ")";
//        SetScalar( alpha, RAYLEIGH_ALPHA, Global::REAL ); 
//        SetScalar( beta, RAYLEIGH_BETA, Global::REAL ); 
      }
    }
    else if ( IsSet(LOSS_TANGENS_DELTA) && IsSet(RAYLEIGH_FREQUENCY) ){

      std::string tanDelta, deltaFreq, omega1, omega2;

      GetScalar( tanDelta, LOSS_TANGENS_DELTA ); 
      // make sure to enclose the expression by brackets!
      tanDelta = "("+ tanDelta + ")";
      deltaFreq= lexical_cast<std::string>(ratioDeltaF)+"*"+ lexical_cast<std::string>(measuredFreq);

      omega1= "(("+lexical_cast<std::string>(measuredFreq)+"-"+deltaFreq+")*2.0*pi)";
      omega2= "(("+lexical_cast<std::string>(measuredFreq)+"+"+deltaFreq+")*2.0*pi)";

      //Computation of alpha and beta according to Habil.Kaltenbacher p. 50 ff
      // alpha + beta*omega_i*omega_i = omega_i*tanDelta_i
      //betaOrig=tanDelta*((omega2-omega1)/(omega2*omega2-omega1*omega1));
      betaOrig=tanDelta+"*((" + omega2 +"-" + omega1 +")/("+omega2+"*"+omega2+"-"+omega1+"*"+omega1+"))";
      //alphaOrig=(omega1*tanDelta)-(betaOrig*omega1*omega1);
      alphaOrig="("+omega1+"*"+tanDelta+")-("+betaOrig+"*"+omega1+"*"+omega1+")";
    }
    else
      EXCEPTION("Error in specification of Rayleigh damping!!!" );

    // If adjustDamping is true and we are in the frequency domain,
    // we can adjust alpha and beta to keep the loss factor (tangens Delta)
    // constant frequency.
    if( adjustDamping && isHarmonic ) {
      
      // calculate modified alpha' = alpha / measuredFreq * f
      alpha = "(" + alphaOrig + ")/" + lexical_cast<std::string>(measuredFreq)  + "* f";
      // calculate modified beta' = beta * measuredFreq / f
      beta =  "(" + betaOrig + ")*" + lexical_cast<std::string>(measuredFreq) + " / f";
    } else {
      alpha = alphaOrig;
      beta = betaOrig;
    }
  }


  void BaseMaterial::InitPiezoMicro( UInt numElemSD, shared_ptr<ElemList> actSDList, 
                                     BaseMaterial* mechMat, BaseMaterial* elecMat,
                                     SubTensorType tensorType, Double dt) {

      isPiezoMicroModel_ = true;

      piezoMicroModel_ = new PiezoMicroModelBK( numElemSD, this, 
                                                mechMat, elecMat, 
                                                tensorType, dt );

      // set map: global to local element number
      EntityIterator it = actSDList->GetIterator();
      UInt iel = 0;
      UInt globalElNr;
      for ( it.Begin(); !it.IsEnd(); it++, iel++) {
	globalElNr = it.GetElem()->elemNum;
	globalElem2Local_[globalElNr] = iel;
      }
  }
  

  void BaseMaterial::GetEffectiveTensors( Matrix<Double>& matMechC,
                                          Matrix<Double>& matMechS,
                                          Matrix<Double>& matElec,
                                          Matrix<Double>& matPiezo,
                                          Vector<Double>& stress, 
                                          Vector<Double>& elecField,
                                          UInt elemIdx, 
                                          bool recompute,
                                          bool previous ) {

    UInt idx = globalElem2Local_[elemIdx];

    piezoMicroModel_->GetEffectiveTensors( matMechC, matMechS, 
                                           matElec,  matPiezo, 
                                           stress, elecField, 
                                           idx, recompute, previous );
  }

  void BaseMaterial::GetEffectiveIrreversibleValues( Vector<Double>& Pirr,
                                                     Vector<Double>& Sirr,
                                                     UInt elemIdx,
                                                     bool recompute,
                                                     bool previous ) {

    UInt idx = globalElem2Local_[elemIdx];

    piezoMicroModel_->GetEffectiveIrreversibleValues( Pirr, Sirr, idx, recompute, previous );
  }

   void BaseMaterial::ComputeEffectiveCouplingTensor(Matrix<Double>& dMatEff, 
                                                     Vector<Double>& elecFieldAct,
                                                     Vector<Double>& elecFieldPrev,
                                                     UInt elemIdx) {

    UInt idx = globalElem2Local_[elemIdx];

    piezoMicroModel_->ComputeEffectiveCouplingTensor(dMatEff, elecFieldAct,
                                                     elecFieldPrev, idx);
   }

}
