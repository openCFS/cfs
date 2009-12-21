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
#include "Domain/entityList.hh"
#include "Utils/piezoMicroModel.hh"
#include "Utils/piezoMicroModelBK.hh"

#include "baseMaterial.hh"



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
  }

   BaseMaterial::~BaseMaterial() {


  }


  void BaseMaterial::NeedApproxMatCurve( ApproxMaterialCurves type ) {


    needApproxMatCurves_.insert( type );
  }


  void BaseMaterial::SetScalar(int param, MaterialType matType) {


    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      integerParams_[matType] = param;
    }
  }
                 
  void BaseMaterial::GetScalar(Integer& param, MaterialType matType, Global::ComplexPart dataType ) const 
  {


    integerMap::const_iterator pos;
    pos = integerParams_.find( matType );

    if ( pos == integerParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      if ( dataType == Global::INTEGER ) {	
        std::string msg = "GetScalar-Integer";
        dataTypeNotAllowed4SetGet( dataType, msg );
      }
      Integer val = pos->second;
      param = val;
    }
  }


  bool BaseMaterial::IsSet( MaterialType matType ) const {

    bool found = false;
    scalarMap::const_iterator scalarIt = scalarParams_.find( matType );
    tensorMap::const_iterator tensorIt = tensorParams_.find( matType );

    if( scalarIt != scalarParams_.end() || 
        tensorIt != tensorParams_.end() ) {
      found = true;
    }

    return found;
  }

  void  BaseMaterial::matTypeNotAllowed(MaterialType matType, const std::string& dim ) const {

    std::string help;
    Enum2String(matType, help);
    EXCEPTION( "Material type (" <<  dim <<  ") " << help << 
               " is not available for " << materialDatabaseName_ 
               << " Database" );
  }
  

  void  BaseMaterial::dataTypeNotAllowed4SetGet(Global::ComplexPart dataType, 
						 const std::string& msg ) const {

    std::string help;
    help = Global::complexPart.ToString( dataType );
    EXCEPTION( "Datatype " << help << " is not allowed in function " 
               << msg );
  }


  void BaseMaterial::dataTypeNotAllowed(Global::ComplexPart dataType, MaterialType matType ) const {

    std::string help1, help2;
    help1 = Global::complexPart.ToString( dataType );
    Enum2String( matType, help2 );
    EXCEPTION( "Datatype " << help1 << " is not allowed for material type " 
               << help2 << " in material data base " << materialDatabaseName_ );
  }

  void BaseMaterial::matTypeNotInDataBase(MaterialType matType, const std::string& dim ) const {

    std::string help;
    Enum2String(matType, help);
    EXCEPTION( "Material type (" << dim << ") " << help 
               << " was not read form/defined in material file" );
  }


  void  BaseMaterial::setMakesNoSense(Global::ComplexPart dataType, const std::string& msg ) const {

    std::string msgAll, help;
    help = Global::complexPart.ToString( dataType );
    EXCEPTION( "Set of " << msg << " makes no sense with datatype "
               << help );
  }


  void  BaseMaterial::subTensorNotAvailable(MaterialType matType, SubTensorType subTensor ) const {

    std::string msg, help1, help2;
    Enum2String(matType, help1);
    Enum2String(subTensor, help2);
    EXCEPTION( "Subtensor " << help2 <<" not available for material type " 
               << help1 );
  }

  std::ostream & operator << ( std::ostream & out, const BaseMaterial& matData)
  {

    std::set<MaterialType>::iterator iter;
    std::set<MaterialType> isSet = matData.GetIsSetInfo();
    std::set<MaterialType> isComplexData = matData.GetIsComplexInfo();

    std::map<MaterialType, Matrix<Complex> > tensorData  = matData.GetTensorParams();
    std::map<MaterialType, Complex >         scalarData  = matData.GetScalarParams();
    std::map<MaterialType, std::string >     stringData  = matData.GetStringParams();
    std::map<MaterialType, Integer >         integerData = matData.GetIntegerParams();

    std::string matTypeName;

    for ( iter = isSet.begin(); iter != isSet.end(); iter++ ) {
      //check, if parameter is complex
      bool isComplex = false;
      if (  isComplexData.find( *iter ) != isComplexData.end() ) {
	isComplex = true;
      }

      //get name of material type
      Enum2String(*iter, matTypeName);

      //check for material data
      std::map<MaterialType, Matrix<Complex> >::const_iterator posTens;
      posTens = tensorData.find( *iter );
      std::map<MaterialType, Complex >::const_iterator posScal;
      posScal = scalarData.find( *iter );
      std::map<MaterialType, std::string >::const_iterator posStr;
      posStr = stringData.find( *iter );
      std::map<MaterialType, Integer >::const_iterator posInt;
      posInt = integerData.find( *iter );



      if ( posTens != tensorData.end() ) {
	// tensor data
	Matrix<Complex> matTensor = posTens->second;
	Matrix<Double>  tensor = matTensor.GetPart( Global::REAL );

	out  << matTypeName << " (real part)" << ":\n\n" 
	     << tensor << std::endl;
      
	if ( isComplex ) {
	  tensor = matTensor.GetPart( Global::IMAG );
	  out  << matTypeName << " (imag. part)" << ":\n\n" 
	       << tensor << std::endl;
	}
      }
      else if ( posScal != scalarData.end() ) {
	// scalar data
	Complex val = posScal->second;
	out << matTypeName << " (real part) = " << val.real() << std::endl;
	
	if ( isComplex ) {
	  out << matTypeName << " (imag. part)" << val.imag() << std::endl;
	}
      }
      else if ( posStr != stringData.end() ) {
	// string data
	std::string val = posStr->second;
	out << matTypeName << ": " << val << std::endl;
      }
      else if ( posInt != integerData.end() ) {
	// integer data
	Integer val = posInt->second;
	out << matTypeName << " = " << val << std::endl;
      }
    }

    return out;
  }  


  void BaseMaterial::RotateTensorByRotationAngles( const Vector<Double>& rotAngle, 
                                                   MaterialType matType,
                                                   bool persistent ) {


    using namespace std;

    tensorMap::const_iterator pos;
    pos = this->tensorParams_.find( matType );

    tensorMap::const_iterator posOrig;
    posOrig = this->tensorParamsOrig_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      if ( posOrig == tensorParamsOrig_.end() ) {
	std::string dim = "tensorOriginal";
	matTypeNotInDataBase( matType, dim );
      }

      Matrix<Complex> matTensor;
      matTensor = pos->second;
      Matrix<Complex> const matTensorOrig = posOrig->second;

      // transfer to radiant
      Vector<Double> rotRadiant = rotAngle;
      for ( UInt i=0; i<rotAngle.GetSize(); i++ ) {
        rotRadiant[i] *= PI / 180.0;
      }
      // limit for angles used in special cases
      // COMPWARNING: unused variable Double eps = 1e-6;

      //compute rotation matrix; check also for special cases 
      Matrix<Complex> RComplex;
      Matrix<Double> R;
      R.Resize(3,3);  
      RComplex.Resize(3,3);

      Matrix<Complex> helpTensor = matTensorOrig;

      // Calculate rotation matrix( based on Kardan-Angles)
      // Ref.: C. Woernle, "Skript: Dynamik von Mehrkörpersystemen,
      // Kapitel 2 "Grundlagen der Kinematik", S. 12, Univ. Rostock
      // http://iamserver.fms.uni-rostock.de/studium/mehrkoerpersysteme/unterlagen.htm

      Double alpha = rotRadiant[0];
      Double beta  = rotRadiant[1];
      Double gamma = rotRadiant[2];
      
      R[0][0] =  cos(beta) * cos(gamma);
      R[0][1] = -cos(beta) * sin(gamma);
      R[0][2] =  sin(beta); 
      R[1][0] =  cos(alpha)*sin(gamma) + sin(alpha)*sin(beta)*cos(gamma);
      R[1][1] =  cos(alpha)*cos(gamma) - sin(alpha)*sin(beta)*sin(gamma);
      R[1][2] = -sin(alpha)*cos(beta);
      R[2][0] =  sin(alpha)*sin(gamma) - cos(alpha)*sin(beta)*cos(gamma);
      R[2][1] =  sin(alpha)*cos(gamma) + cos(alpha)*sin(beta)*sin(gamma);
      R[2][2] =  cos(alpha)*cos(beta);

      RComplex.Resize(3,3);
      RComplex.Init();
      RComplex.SetPart( Global::REAL, R );
      PerformRotation(RComplex, matTensor, helpTensor); 
      tensorParams_[matType] = matTensor;
      
      // In the case of a persistent rotation, we override the
      // original tensor
      if ( persistent ) {
        tensorParamsOrig_[matType] = matTensor;
      }
        

//       if ( abs(rotAngle[0]) > eps ) {
// 	// rotate around x-axis
// 	R.Init(Complex(0.0,0.0));
// 	R[0][0] =  Complex( 1.0, 0.0);
// 	R[1][1] =  Complex( std::cos(rotAngle[0]), 0.0 );
// 	R[1][2] =  -Complex( std::sin(rotAngle[0]), 0.0 );
// 	R[2][1] = -R[1][2];
// 	R[2][2] =  R[1][1];

//         //std::cerr << "matTensor before first rotation \n" << helpTensor << std::endl;
// 	PerformRotation(R, matTensor, helpTensor);  
// 	helpTensor = matTensor;
//         //std::cerr << "matTensor after first rotation \n" << matTensor << std::endl;
//       }


//       if ( abs(rotAngle[1]) > eps ) {
// 	// rotate around y-axis
// 	R.Init(Complex(0.0,0.0));
// 	R[0][0] = Complex( std::cos(rotAngle[1]), 0.0 );
// 	R[0][2] = Complex( std::sin(rotAngle[1]), 0.0 );
// 	R[1][1] = Complex( 1.0, 0.0);
// 	R[2][0] = -R[0][2];
// 	R[2][2] =  R[0][0];

// 	PerformRotation(R, matTensor, helpTensor);
// 	helpTensor = matTensor;  
//       }

//       if ( abs(rotAngle[2]) > eps ) {
// 	// rotate around z-axis
// 	R.Init(Complex(0.0,0.0));
// 	R[0][0] =  Complex( std::cos(rotAngle[2]), 0.0 );
// 	R[0][1] =  -Complex( std::sin(rotAngle[2]), 0.0 );
// 	R[1][0] = -R[0][1];
// 	R[1][1] =  R[0][0];
// 	R[2][2] =  Complex( 1.0, 0.0);
//         //std::cerr << "matTensor before 2nd rotation \n" << helpTensor << std::endl;
// 	PerformRotation(R, matTensor, helpTensor);  
//         //std::cerr << "matTensor after 2nd rotation \n" << matTensor << std::endl;
//       }
    
//       // save rotated matrix back
//       tensorParams_[matType] = matTensor;
    }
    
  }

  void BaseMaterial::RotateAllTensorsByRotationAngles( const Vector<Double>& rotAngle, 
                                                       bool persistent ) {

    BaseMaterial::tensorMap::iterator it = tensorParams_.begin();
    for( ; it != tensorParams_.end(); it++ ) {
      RotateTensorByRotationAngles( rotAngle, it->first, persistent );
    }
  }
    

  void BaseMaterial::RotateTensorByPointCoord( const Vector<Double>&  coord,
                                               MaterialType matType ) {


    // Determine rotation angles from attached coordinate system
    Vector<Double> angles;
    coosy_->GetGlobRotationAngles( angles, coord );

    // Calculate rotation. In this case this is always
    // non-persistent, as the orientation may change for each point in the
    // related coordinate system
    RotateTensorByRotationAngles( angles, matType, false );

  }


  void BaseMaterial::PerformRotation( Matrix<Complex>& R,  Matrix<Complex>& matTensor,
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
    computeHystInverse_ = computeHystInverse_;

    std::string val = stringParams_[HYST_MODEL];
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


  Double BaseMaterial::ConvertVec2ScalarWithSign( Vector<Double> & vecVal ) {

    Double absVal, maxVal;
    UInt idx;

    absVal = vecVal.NormL2();
    idx = 0;
    maxVal = abs( vecVal[0] );

    for ( UInt i=1; i<vecVal.GetSize(); i++ ) {
      if ( abs(vecVal[i]) > maxVal ) {
        maxVal = vecVal[i];
        idx++;
      }
    }
//     if ( vecVal[idx] < 0  )
//       absVal *= -1.0;

    return absVal;
  }


  void BaseMaterial::SetPreviousHystVal( UInt nrElem, Vector<Double>& Xvec) {

    UInt idx = globalElem2Local_[nrElem];

    Double Xval = ConvertVec2ScalarWithSign(Xvec);
    // first compute the differential material value
    Xprevious_[idx] = Xval;
    Yprevious_[idx] = hyst_->computeValue( Xval, idx );
  }

  Double BaseMaterial::ComputeScalarDiffVal( UInt nrElem, Vector<Double>& Xvec) {

    Double matDiff, eps;

    UInt idx = globalElem2Local_[nrElem];

    Double Xval     = ConvertVec2ScalarWithSign(Xvec);
    Double Ycurrent = hyst_->computeValue(Xval, idx);

    std::cout << "B=" << Xval << "  H=" << Ycurrent << std::endl;
    //compute differential material parameter
    Double dX = Xval - Xprevious_[idx];
    Double dY = Ycurrent -Yprevious_[idx];

    if ( (abs(dY) < 1e-12) || (abs(dX) < 1e-10) ) {
      //std::cout << "Use start eps" << std::endl;
      if ( materialDatabaseName_ == "Electrostatic" )
        GetScalar(eps,ELEC_PERMITTIVITY,Global::REAL);
      else if ( materialDatabaseName_ == "Electromagnetics" )
        GetScalar(eps,MAG_RELUCTIVITY,Global::REAL);

      matDiff = eps;
    }
    else {
      matDiff = dY / dX;
    }

    
    std::cout << "dB=" << dX << "  dH=" << dY <<  "  dnu=" << matDiff << std::endl << std::endl;

//     GetScalar(eps,MAG_RELUCTIVITY,Global::REAL);
//     matDiff = eps;

    return matDiff;
  }

  Double BaseMaterial::ComputeScalarHystVal( UInt nrElem, Vector<Double>& Xvec) {

    UInt idx = globalElem2Local_[nrElem];

    Double Xval = ConvertVec2ScalarWithSign(Xvec);
    Double Yval = hyst_->computeValue( Xval, idx );
    
    return Yval;
  }


  Double BaseMaterial::GetScalarHystVal( UInt nrElem ) {

    UInt idx = globalElem2Local_[nrElem];
    Double Yval = hyst_->getValue( idx );

    return Yval;
  }


  Double BaseMaterial::GetScalarHystPrevVal( UInt nrElem ) {
    UInt idx = globalElem2Local_[nrElem];
    return Yprevious_[idx];
  }

  void BaseMaterial::ComputeRayleighDamping(Double dampFreq, Double RatioDeltaF) {

    if ( IsSet( RAYLEIGH_ALPHA ) 
         && IsSet( RAYLEIGH_BETA ) 
         && IsSet(RAYLEIGH_FREQUENCY) ) {
      Double alpha, beta, freq;

      GetScalar( alpha, RAYLEIGH_ALPHA, Global::REAL ); 
      GetScalar( beta, RAYLEIGH_BETA, Global::REAL ); 
      GetScalar( freq, RAYLEIGH_FREQUENCY, Global::REAL ); 

      if( abs(freq-dampFreq) > 0.001*freq ){
        alpha*=(dampFreq/freq);
        beta*=(freq/dampFreq);
        SetScalar( alpha, RAYLEIGH_ALPHA, Global::REAL ); 
        SetScalar( beta, RAYLEIGH_BETA, Global::REAL ); 
      }
    }
    else if ( IsSet(LOSS_TANGENS_DELTA) && IsSet(RAYLEIGH_FREQUENCY) ){

      Double alpha, beta, tanDelta, deltaFreq, omega1, omega2;

      GetScalar( tanDelta, LOSS_TANGENS_DELTA, Global::REAL ); 

      deltaFreq=RatioDeltaF*dampFreq;

      omega1= (dampFreq-deltaFreq)*2.0*PI;
      omega2= (dampFreq+deltaFreq)*2.0*PI;

      //Computation of alpha and beta according to Habil.Kaltenbacher p. 50 ff
      // alpha + beta*omega_i*omega_i = omega_i*tanDelta_i
      beta=2.0*tanDelta*((omega2-omega1)/(omega2*omega2-omega1*omega1));
      alpha=(2.0*omega1*tanDelta)-(beta*omega1*omega1);

      SetScalar( alpha, RAYLEIGH_ALPHA, Global::REAL ); 
      SetScalar( beta, RAYLEIGH_BETA, Global::REAL ); 
    }
    else
      EXCEPTION("Error in specification of Rayleigh damping!!!" );

    //Info->PrintF()
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
