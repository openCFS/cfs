#include "BaseMaterial.hh"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "Utils/StdVector.hh"
#include "MatVec/Matrix.hh"

#include "Materials/Models/Preisach.hh"
#include "Materials/Models/SimplePreisachInv.hh"
#include "Materials/Models/PiezoMicroModelHF.hh"
#include "Materials/Models/PiezoMicroModelBK.hh"

#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
using std::string;
using std::map;
using std::set;


namespace CoupledField
{

  BaseMaterial::MatDescriptorNl::MatDescriptorNl() {
    measAccuracy = 0.0;
    maxVal = 0.0;
    angle = 0.0;
    factor = 0.0;
    approxData = NULL;
  }
  
  BaseMaterial::MatDescriptorNl::~MatDescriptorNl() {
    if (approxData) {
      delete approxData;
    }
  }



  // ***********************
  //   Default Constructor
  // ***********************
  BaseMaterial::BaseMaterial(MathParser* mp,
                             CoordSystem * defaultCoosy ) {


    materialDatabaseName_ = "";
    matFileName_ = "";
    //    nonlinFileName_ = "";

    coosy_ = NULL;
    hyst_  = NULL;

    isHysteresis_  = false;
    isHystInverse_ = false;

    piezoMicroModel_   = NULL;
    isPiezoMicroModel_ = false;
    
    mp_ = mp;
    
    // obtain standard coordinate system
    coosy_ = defaultCoosy;
  }

  BaseMaterial::~BaseMaterial() {
    handleMap::iterator it = scalarStringHandlesReal_.begin(),
                        itEnd = scalarStringHandlesReal_.end();
    for ( ; it != itEnd; ++it ) {
      mp_->ReleaseHandle(it->second);
    }

    it = scalarStringHandlesImag_.begin();
    itEnd = scalarStringHandlesImag_.end();
    for ( ; it != itEnd; ++it ) {
      mp_->ReleaseHandle(it->second);
    }
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

  void BaseMaterial::SetNonLinMatIso( MaterialType matType, MatDescriptorNl& data ) {
    
    if( nonlinIsoParams_.find(matType) != nonlinIsoParams_.end() ) {
      EXCEPTION( "Nonlinear material parameter '" << MaterialTypeEnum.ToString(matType)
                 << "' was already set");
    }
    nonlinIsoParams_[matType] = data;
  }

  void BaseMaterial::SetNonLinMatAniso( MaterialType matType, StdVector<MatDescriptorNl>& data ) {
    if( nonlinAnisoParams_.find(matType) != nonlinAnisoParams_.end() ) {
      EXCEPTION( "Nonlinear material parameter '" << MaterialTypeEnum.ToString(matType)
                 << "' was already set");
    }
    nonlinAnisoParams_[matType] = data;
  }
  
  void BaseMaterial::SetCoefFct( MaterialType matType, PtrCoefFct coef ) {
    
    // check, if material type is allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
       std::string dim = CoefFunction::CoefDimType_.ToString(coef->GetDimType());
       matTypeNotAllowed( matType, dim );
     }
    
    
    // switch depending on dimType of coefficient function
    switch(coef->GetDimType()) {
      case CoefFunction::SCALAR:
        scalarCoef_[matType] = coef;
        break;
      case CoefFunction::VECTOR:
        vectorCoef_[matType] = coef;
        break;
      case CoefFunction::TENSOR:
        tensorCoef_[matType] = coef;
        tensorOrigCoef_[matType] = coef;
        break;
      default:
        EXCEPTION("Unknown entry type of coefficient function");
        break;
    }
  }

  bool BaseMaterial::IsSet( MaterialType matType ) const {

    bool found = false;
    stringMap::const_iterator stringIt = stringParams_.find( matType );
    scalarMap::const_iterator scalarIt = scalarParams_.find( matType );
    tensorMap::const_iterator tensorIt = tensorParams_.find( matType );
    CoefMap::const_iterator coefScalIt = scalarCoef_.find( matType );
    CoefMap::const_iterator coefVecIt = vectorCoef_.find( matType );
    CoefMap::const_iterator coefTensIt = tensorCoef_.find( matType );
        
    if( stringIt != stringParams_.end() ||
        scalarIt != scalarParams_.end() || 
        tensorIt != tensorParams_.end() || 
        coefScalIt != scalarCoef_.end() ||
        coefVecIt != vectorCoef_.end()  ||
        coefTensIt != tensorCoef_.end() ) {
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
    EXCEPTION("Sub-tensor " << help2 <<" not available for material type " << help1);
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
    
    // 1) Constant-Valued tensors
    BaseMaterial::tensorMap::iterator it = tensorParams_.begin();
    for( ; it != tensorParams_.end(); it++ ) {
      RotateTensorByRotationAngles( rotAngle, it->first, persistent );
    }
    
    // 2) Coefficient tensors
    BaseMaterial::CoefMap::iterator cIt = tensorCoef_.begin();
    for( ; cIt != tensorCoef_.end(); cIt++ ) {
      RotateTensorByRotationAngles( rotAngle, cIt->first, persistent );
    }
  }
  
  void BaseMaterial::RotateTensorByRotationMatrix( const Matrix<Double>& rotMatrix, 
                                                   MaterialType matType,
                                                   bool persistent ) {


    using namespace std;

    
    
    tensorMap::iterator pos;
    pos = this->tensorParams_.find( matType );
    if( pos != tensorParams_.end() ) {
      // -----------------------------------------
      // 1) Tensor is defined as constant matrix
      // -----------------------------------------
    tensorMap::const_iterator posOrig;
    posOrig = this->tensorParamsOrig_.find( matType );

      if ( posOrig == tensorParamsOrig_.end() ) {
        string dim = "tensorOriginal";
        matTypeNotInDataBase( matType, dim );
      }

      Matrix<Complex> & matTensor = pos->second;
      const Matrix<Complex> & matTensorOrig = posOrig->second;

      //perform rotation 
      matTensorOrig.PerformRotation(rotMatrix, matTensor); 
      
      // In the case of a persistent rotation, we override the
      // original tensor as well
      if ( persistent ) {
        tensorParamsOrig_[matType] = matTensor;
      }
    }  else if( tensorCoef_.find( matType) != tensorCoef_.end() ) {
      // ----------------------------------------------
      // 2) Tensor is defined as coefficient function
      // ----------------------------------------------
      PtrCoefFct coef;
      PtrCoefFct coefOrig = tensorOrigCoef_[matType];

      //perform rotation
      PerformRotation( rotMatrix, coef, coefOrig ); 
      // store back rotated material
      tensorCoef_[matType] = coef;
      
      // In the case of a persistent rotation, we override the
      // original tensor as well
      if ( persistent ) {
        tensorOrigCoef_[matType] = coef;
      }
    } else {
      string dim = "tensor";
                  matTypeNotInDataBase( matType, dim );
    }
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


//  void BaseMaterial::PerformRotation( Matrix<Complex>& R,  
//                                      Matrix<Complex>& matTensor,
//                                      const Matrix<Complex>& matTensorOrig) {
//
//    // get memory for transposed rotation matrix
//    Matrix<Complex> RT;
//    RT.Resize(3,3);
//    R.Transpose(RT);
//
//    //get dimension of matrix
//    UInt rowSize = matTensorOrig.GetNumRows();
//    UInt colSize = matTensorOrig.GetNumCols();
//
//    Matrix<Complex> helpMat;
//
//    if ( rowSize == 3 && colSize == 3) {
//      // tensor is a 3x3 matrix: sol = R * matrixOrig * RT
//      helpMat   = matTensorOrig * RT;
//      matTensor = R * helpMat;
//    }
//    else {
//      // we also need Q;
//      Matrix<Complex> Q;
//
//      // Composed Rotation Matrix
//      // Ref.: M.Richter, "Entwicklung mechanischer Modelle zur analytischen
//      // Beschreibung der Materialeigenschaften von textilbewehrtem Feinbeton",
//      // Diss., Dresden, 2005, p. 27
//
//      Q.Resize(6,6);  
//
//      Q[0][0] = R[0][0]*R[0][0];
//      Q[0][1] = R[0][1]*R[0][1];
//      Q[0][2] = R[0][2]*R[0][2];
//      Q[0][3] = 2.0*R[0][1]*R[0][2];
//      Q[0][4] = 2.0*R[0][0]*R[0][2];
//      Q[0][5] = 2.0*R[0][0]*R[0][1];
//
//      Q[1][0] = R[1][0]*R[1][0];
//      Q[1][1] = R[1][1]*R[1][1];
//      Q[1][2] = R[1][2]*R[1][2];
//      Q[1][3] = 2.0*R[1][1]*R[1][2];
//      Q[1][4] = 2.0*R[1][0]*R[1][2];
//      Q[1][5] = 2.0*R[1][0]*R[1][1];
//
//      Q[2][0] = R[2][0]*R[2][0];
//      Q[2][1] = R[2][1]*R[2][1];
//      Q[2][2] = R[2][2]*R[2][2];
//      Q[2][3] = 2.0*R[2][1]*R[2][2];
//      Q[2][4] = 2.0*R[2][0]*R[2][2];
//      Q[2][5] = 2.0*R[2][0]*R[2][1];
//
//      Q[3][0] = R[1][0]*R[2][0];
//      Q[3][1] = R[1][1]*R[2][1];
//      Q[3][2] = R[1][2]*R[2][2];
//      Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
//      Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
//      Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];
//
//      Q[4][0] = R[0][0]*R[2][0];
//      Q[4][1] = R[0][1]*R[2][1];
//      Q[4][2] = R[0][2]*R[2][2];
//      Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
//      Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
//      Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];
//
//      Q[5][0] = R[0][0]*R[1][0];
//      Q[5][1] = R[0][1]*R[1][1];
//      Q[5][2] = R[0][2]*R[1][2];
//      Q[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
//      Q[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0];
//      Q[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0];
//
//
//      // 	std::cout << "R:\n" << R << std::endl;
//      // 	std::cout << "Q:\n" << Q << std::endl;
//      // 	std::cout << "Tensor orig:\n" << matTensor << std::endl;
//
//      Matrix<Complex> QT;
//      QT.Resize(6,6);
//      Q.Transpose(QT);
//
//      if ( rowSize == 3 && colSize == 6 ) {
//        helpMat   = matTensorOrig * QT;
//        matTensor = R * helpMat;
//      }
//      else if (rowSize == 6 && colSize == 6 ) {
//        helpMat   = matTensorOrig * QT;
//        matTensor = Q * helpMat;
//      }
//      // 	else {
//      // 	  EXCEPTION("Cannot rotate tensor due to dimensions!");
//      // 	}
//    }
//  }

  void BaseMaterial::PerformRotation( const Matrix<Double>& R, PtrCoefFct& rotatedCoef,
                                       PtrCoefFct origCoef) {
    
    // determine entry type (double / complex) of coefficient
    Global::ComplexPart part = 
        origCoef->IsComplex() ? Global::COMPLEX : Global::REAL;
    
    // get memory for transposed rotation matrix
    Matrix<Double> RT;
    RT.Resize(3,3);
    R.Transpose(RT);
    
    // convert both to coefficient functions
    shared_ptr<CoefFunctionConst<Double> > cR(new CoefFunctionConst<Double>());
    shared_ptr<CoefFunctionConst<Double> > cRT(new CoefFunctionConst<Double>());
    cR->SetTensor( R );
    cRT->SetTensor( RT );

    //get dimension of matrix
    UInt rowSize, colSize;
    origCoef->GetTensorSize(rowSize, colSize );

    Matrix<Complex> helpMat;

    // perform rotation of 3 x 3 matrix
    if ( rowSize == 3 && colSize == 3) {
      // tensor is a 3x3 matrix: sol = R * matrixOrig * RT
      CoefXprBinOp tmp( mp_, origCoef, cRT, CoefXpr::OP_MULT );
      CoefXprBinOp final( mp_, cR, tmp, CoefXpr::OP_MULT );
      rotatedCoef = CoefFunction::Generate( mp_, part, final );
    }
    else {
      // we also need Q;
      Matrix<Double> Q;

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


      //  std::cout << "R:\n" << R << std::endl;
      //  std::cout << "Q:\n" << Q << std::endl;
      //  std::cout << "Tensor orig:\n" << matTensor << std::endl;

      Matrix<Double> QT;
      QT.Resize(6,6);
      Q.Transpose(QT);
      
      // convert both to coefficient functions
      shared_ptr<CoefFunctionConst<Double> > cQ(new CoefFunctionConst<Double>());
      shared_ptr<CoefFunctionConst<Double> > cQT(new CoefFunctionConst<Double>());
      cQ->SetTensor( Q );
      cQT->SetTensor( QT );

      if ( rowSize == 3 && colSize == 6 ) {
        CoefXprBinOp tmp( mp_, origCoef, cQT, CoefXpr::OP_MULT );
        UInt nRows, nCols;
        StdVector<std::string>t1, t2;
        tmp.GetTensorXpr(nRows, nCols, t1, t2 );
        CoefXprBinOp final( mp_, cR, tmp, CoefXpr::OP_MULT );
        final.GetTensorXpr(nRows, nCols, t1, t2 );
        rotatedCoef = CoefFunction::Generate(mp_,  part, final );
      }
      else if (rowSize == 6 && colSize == 6 ) {
        CoefXprBinOp tmp( mp_, origCoef, cQT, CoefXpr::OP_MULT );
        PtrCoefFct tmp2 = CoefFunction::Generate(mp_,  part, tmp );
        CoefXprBinOp final( mp_, cQ, tmp, CoefXpr::OP_MULT );
        rotatedCoef = CoefFunction::Generate(mp_, part, final );
      
      }
      //  else {
      //    EXCEPTION("Cannot rotate tensor due to dimensions!");
      //  }
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


   PtrCoefFct  BaseMaterial::GetTensorCoefFnc(MaterialType matType, SubTensorType type,
                                              Global::ComplexPart matDataType, 
                                              bool transpose ) {
     PtrCoefFct mFunct;
     
     if( tensorCoef_.find(matType) !=  tensorCoef_.end() ) {
       // --------------------------------------
       //  Coefficient Function already defined
       // --------------------------------------
       PtrCoefFct tmp = GetSubTensorCoefFnc( matType, type, transpose  );
       PtrCoefFct ret = tmp->GetComplexPart( matDataType );
       return ret;
       
       
     } else {
       // -------------------------------------------
       //  Create CoefFunction from constant entries
       // -------------------------------------------
       if(matDataType == Global::REAL){
         CoefFunctionConst<Double>* tmpFnc = new CoefFunctionConst<Double>();
         Matrix<Double> coefMat;
         GetTensor(coefMat,matType,matDataType,type);
         // transpose if flag is true
         if( transpose ) {
           Matrix<Double> temp;
           coefMat.Transpose(temp);
           coefMat = temp;
         }
         tmpFnc->SetTensor(coefMat);
         mFunct.reset(tmpFnc);
       }else if(matDataType == Global::COMPLEX){
         CoefFunctionConst<Complex>* tmpFnc = new CoefFunctionConst<Complex>();
         Matrix<Complex> coefMat;
         GetTensor(coefMat,matType,matDataType,type);
         // transpose if flag is true
         if( transpose ) {
           Matrix<Complex> temp;
           coefMat.Transpose(temp);
           coefMat = temp;
         }
         tmpFnc->SetTensor(coefMat);
         mFunct.reset(tmpFnc);
       }else{
         EXCEPTION("Material Data Type not supported");
       }
     }
     mFunct->SetCoordinateSystem(this->coosy_);
     return mFunct;
   }
   
   PtrCoefFct  BaseMaterial::GetVectorCoefFnc(MaterialType matType,
                                              Global::ComplexPart matDataType) {
     //EXCEPTION("Not implemented")
     PtrCoefFct mFunct;

     if( vectorCoef_.find(matType) !=  vectorCoef_.end() ) {
       // --------------------------------------
       //  Coefficient Function already defined
       // --------------------------------------
       mFunct = vectorCoef_[matType]->GetComplexPart( matDataType );

     } else {
       // -------------------------------------------
       //  Create CoefFunction from constant entries
       // -------------------------------------------
       if(matDataType == Global::REAL){
         CoefFunctionConst<Double>* tmpFnc = new CoefFunctionConst<Double>();
         Vector<Double> real;
         GetVector(real,matType,matDataType);
         tmpFnc->SetVector(real);
         mFunct.reset(tmpFnc);
       }else if(matDataType == Global::COMPLEX){
         CoefFunctionConst<Complex>* tmpFnc = new CoefFunctionConst<Complex>();
         Vector<Complex> val;
         GetVector(val,matType,matDataType);
         tmpFnc->SetVector(val);
         mFunct.reset(tmpFnc);
       }else{
         EXCEPTION("Material Data Type not supported");
       }
     }
     mFunct->SetCoordinateSystem(this->coosy_);
     return mFunct;
   }

   PtrCoefFct  BaseMaterial::GetScalCoefFnc(MaterialType matType,
                                            Global::ComplexPart matDataType) {
     PtrCoefFct mFunct;

     if( scalarCoef_.find(matType) !=  scalarCoef_.end() ) {
       // --------------------------------------
       //  Coefficient Function already defined
       // --------------------------------------
       mFunct = scalarCoef_[matType]->GetComplexPart( matDataType );

     } else {
       // -------------------------------------------
       //  Create CoefFunction from constant entries
       // -------------------------------------------
       if(matDataType == Global::REAL){
         CoefFunctionConst<Double>* tmpFnc = new CoefFunctionConst<Double>();
         Double real;
         GetScalar(real,matType,matDataType);
         tmpFnc->SetScalar(real);
         mFunct.reset(tmpFnc);
       }else if(matDataType == Global::COMPLEX){
         CoefFunctionConst<Complex>* tmpFnc = new CoefFunctionConst<Complex>();
         Complex val;
         // transpose if flag is true
         GetScalar(val,matType,matDataType);
         tmpFnc->SetScalar(val);
         mFunct.reset(tmpFnc);
       }else{
         EXCEPTION("Material Data Type not supported");
       }
     }
     mFunct->SetCoordinateSystem(this->coosy_);
     return mFunct;
   }

   PtrCoefFct BaseMaterial::GetSubTensorCoefFnc( MaterialType matType, 
                                                 SubTensorType tensorType,
                                                 bool transposed  ) {

     PtrCoefFct mFunct;
     if( tensorCoef_.find(matType) !=  tensorCoef_.end() ) {
       CoefXprSubTensor subTensorXpr(mp_,  tensorCoef_[matType] );

       subTensorXpr.SetSubTensorType( tensorType, transposed );
       mFunct = CoefFunction::Generate( mp_, Global::COMPLEX, subTensorXpr );
     } else {
       EXCEPTION( "Material tensor not found" );
     }
     return mFunct;
   }

   PtrCoefFct BaseMaterial::GetTensorCoefFncNonLin( MaterialType matType,
                                                    SubTensorType type,
                                                    Global::ComplexPart matDataType,
                                                    PtrCoefFct dependency ) {
     EXCEPTION("Currently only implemented for ElectroMagnetic material")
   }
   
   

   PtrCoefFct BaseMaterial::GetScalCoefFncNonLin(MaterialType matType,
                                                 Global::ComplexPart matDataType,
                                                 PtrCoefFct dependency) {
     shared_ptr<CoefFunctionApprox> coef;
     
     // Ensure that only real-valued parameters are used
     if( matDataType != Global::REAL ) {
       EXCEPTION( "Only real-valued nonlinear parameters are supported");
     }
     
     // check if isotropic material type is defined
     if( nonlinIsoParams_.find(matType) == nonlinIsoParams_.end() ) {
       EXCEPTION( "No nonlinear definition found for material type '"
           << MaterialTypeEnum.ToString(matType) << "'");
     }
     
     // check, if nonlinear curve was already calculated
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
       else {
         EXCEPTION( "No nonlinear approx/interpolate type not available '"
             << ApproxCurveTypeEnum.ToString(matNl.approxType) << "'");
       }
     }

     ApproxData * sp = matNl.approxData;
     
     // get linear starting value
     Double startVal = 0.0;
     this->GetScalar( startVal, matType, Global::REAL );
     coef.reset(new CoefFunctionApprox());
     coef->Init( startVal, sp, dependency );

     return coef;
   }
} // end of namespace
