#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>
#include "Utils/StdVector.hh"
#include "Matrix/matrix.hh"

#include "baseMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  BaseMaterial::BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");

    materialDatabaseName_ = "";
    matFileName_ = "";
    nonlinFileName_ = "";
    isScalar    = FALSE;
    isIsotrop   = FALSE;
    isOrthotrop = FALSE;
    isTensor    = FALSE;

  }


   BaseMaterial::~BaseMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }


  void BaseMaterial::SetScalar( Integer& param, const MaterialType& matType) {

    ENTER_FCN( "Basematerial::SetScalar" );

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

  void BaseMaterial::GetScalar( Integer& param, const MaterialType& matType, 
				const DataType& dataType ) const {

    ENTER_FCN( "BaseMaterial::GetScalar" );

    integerMap::const_iterator pos;
    pos = integerParams_.find( matType );

    if ( pos == integerParams_.end() ) {
      std::string dim = "scalar";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      if ( dataType == INTEGER ) {	
	std::string msg = "GetScalar-Integer";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
      Integer val = pos->second;
      param = val;
    }
  }


  bool BaseMaterial::IsSet( MaterialType matType ) const {
    ENTER_FCN(" BaseMaterial::IsSet" );

    bool found = false;
    scalarMap::const_iterator scalarIt = scalarParams_.find( matType );
    tensorMap::const_iterator tensorIt = tensorParams_.find( matType );

    if( scalarIt != scalarParams_.end() || 
        tensorIt != tensorParams_.end() ) {
      found = true;
    }

    return found;
  }

  void  BaseMaterial::matTypeNotAllowed( const MaterialType& matType, 
					 std::string dim ) const {
    ENTER_FCN("BaseMaterial::matTypeNotAllowed");

    std::string msg, help;
    Enum2String(matType, help);
    msg = "Material type (" + dim + ") " + help + 
          " is not available for " + materialDatabaseName_ + " Database";
    Error(msg.c_str(), __FILE__, __LINE__);
  }
  

  void  BaseMaterial::dataTypeNotAllowed4SetGet( const DataType& dataType, 
						 const std::string msg ) const {
    ENTER_FCN("BaseMaterial::dataTypeNotAllowed4Set");

    std::string msgAll, help;
    Enum2String( dataType, help );
    msgAll = "Datatype " + help + " is not allowed in function " 
           + msg;
    Error(msgAll.c_str(), __FILE__, __LINE__);
  }


  void BaseMaterial::dataTypeNotAllowed( const DataType& dataType, 
					  const MaterialType& matType ) const {
    ENTER_FCN("BaseMaterial::dataTypeNotAllowed");

    std::string msg, help1, help2;
    Enum2String( dataType, help1 );
    Enum2String( matType, help2 );
    msg = "Datatype " + help1 + " is not allowed for material type " 
           + help2 + " in material data base " + materialDatabaseName_;
    Error(msg.c_str(), __FILE__, __LINE__);
  }

  void BaseMaterial::matTypeNotInDataBase( const MaterialType& matType,
					   std::string dim ) const {
    ENTER_FCN("BaseMaterial::matTypeNotInFile");

    std::string msg, help;
    Enum2String(matType, help);
    msg = "Material type (" + dim + ") " + help + " was not read form/defined in material file";
    Error(msg.c_str(), __FILE__, __LINE__);
  }


  void  BaseMaterial::setMakesNoSense( const DataType& dataType, 
				       std::string msg ) const {
    ENTER_FCN("BaseMaterial::setMakesNoSense");

    std::string msgAll, help;
    Enum2String( dataType, help );
    msgAll = "Set of " + msg + " makes no sense with datatype " + help;
    Error(msgAll.c_str(), __FILE__, __LINE__);
  }


  void  BaseMaterial::subTensorNotAvailable( const MaterialType& matType, 
					     const SubTensorType subTensor ) const {
    ENTER_FCN("BaseMaterial::subTensorNotAvailable") 

    std::string msg, help1, help2;
    Enum2String(matType, help1);
    Enum2String(subTensor, help2);
    msg = "Subtensor " + help2 + " not available for material type " + help1;
    Error(msg.c_str(), __FILE__, __LINE__);
  }

  std::ostream & operator << ( std::ostream & out, const BaseMaterial& matData)
  {
    ENTER_IFCN( "operator <<(BaseMaterial)" );

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
      Boolean isComplex = FALSE;
      if (  isComplexData.find( *iter ) != isComplexData.end() ) {
	isComplex = TRUE;
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
	Matrix<Double>  tensor = matTensor.GetPart( REAL );

	out  << matTypeName << " (real part)" << ":\n\n" 
	     << tensor << std::endl;
      
	if ( isComplex ) {
	  tensor = matTensor.GetPart( IMAG );
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


  void BaseMaterial::RotateTensorByRotationAngles( StdVector<Double>& rotAngle, 
						   const MaterialType& matType) {

    ENTER_FCN( "BaseMaterial::RotateTensorByRotationAngles" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    tensorMap::const_iterator posOrig;
    posOrig = tensorParamsOrig_.find( matType );

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

      Matrix<Complex> matTensorOrig;
      matTensorOrig = posOrig->second;

      // transfer to radiant
      for ( UInt i=0; i<rotAngle.GetSize(); i++ ) {
	rotAngle[i] *= PI / 180.0;
      }

      // limit for angles used in special cases
      Double eps = 1e-6;

      //compute rotation matrix; check also for special cases 
      Matrix<Complex> R;
      R.Resize(3,3);  

      Matrix<Complex> helpTensor = matTensorOrig;

      if ( abs(rotAngle[0]) > eps ) {
	// rotate around x-axis
	R.Init(Complex(0.0,0.0));
	R[0][0] =  Complex( 1.0, 0.0);
	R[1][1] =  Complex( std::cos(rotAngle[0]), 0.0 );
	R[1][2] =  Complex( std::sin(rotAngle[0]), 0.0 );
	R[2][1] = -R[1][2];
	R[2][2] =  R[1][1];

	PerformRotation(R, matTensor, helpTensor);  
	helpTensor = matTensor;
      }


      if ( abs(rotAngle[1]) > eps ) {
	// rotate around y-axis
	R.Init(Complex(0.0,0.0));
	R[0][0] =  Complex( std::cos(rotAngle[1]), 0.0 );
	R[0][2] = -Complex( std::sin(rotAngle[1]), 0.0 );
	R[1][1] =  Complex( 1.0, 0.0);
	R[2][0] = -R[0][2];
	R[2][2] =  R[0][0];

	PerformRotation(R, matTensor, helpTensor);
	helpTensor = matTensor;  
      }

      if ( abs(rotAngle[2]) > eps ) {
	// rotate around z-axis
	R.Init(Complex(0.0,0.0));
	R[0][0] =  Complex( std::cos(rotAngle[2]), 0.0 );
	R[0][1] =  Complex( std::sin(rotAngle[2]), 0.0 );
	R[1][0] = -R[0][1];
	R[1][1] =  R[0][0];
	R[2][2] =  Complex( 1.0, 0.0);

	PerformRotation(R, matTensor, helpTensor);  
      }

    }
  }


  void BaseMaterial::PerformRotation( Matrix<Complex>& R,  Matrix<Complex>& matTensor,
				      Matrix<Complex>& matTensorOrig) {

    ENTER_FCN( "BaseMaterial::PerformRotation" );

      // get memory for transposed rotation matrix
      Matrix<Complex> RT;
      RT.Resize(3,3);
      R.Transpose(RT);

      //get dimension of matrix
      UInt rowSize = matTensorOrig.GetSizeRow();
      UInt colSize = matTensorOrig.GetSizeCol();

      Matrix<Complex> helpMat;

      if ( rowSize == 3 && colSize == 3) {
	// tensor is a 3x3 matrix: sol = R * matrixOrig * RT
	helpMat   = matTensorOrig * RT;
	matTensor = R * helpMat;
      }
      else {
	// we also need Q;
	Matrix<Complex> Q;
	Q.Resize(6,6);  // Composed Rotation Matrix

	Q[0][0] = R[0][0]*R[0][0];
	Q[0][1] = R[0][1]*R[0][1];
	Q[0][2] = R[0][2]*R[0][2];
	Q[0][3] = R[0][1]*R[0][2];
	Q[0][4] = R[0][0]*R[0][2];
	Q[0][5] = R[0][0]*R[0][1];

	Q[1][0] = R[1][0]*R[1][0];
	Q[1][1] = R[1][1]*R[1][1];
	Q[1][2] = R[1][2]*R[1][2];
	Q[1][3] = R[1][1]*R[1][2];
	Q[1][4] = R[1][0]*R[1][2];
	Q[1][5] = R[1][0]*R[1][1];

	Q[2][0] = R[2][0]*R[2][0];
	Q[2][1] = R[2][1]*R[2][1];
	Q[2][2] = R[2][2]*R[2][2];
	Q[2][3] = R[2][1]*R[2][2];
	Q[2][4] = R[2][0]*R[2][2];
	Q[2][5] = R[2][0]*R[2][1];

	Q[3][0] = 2.0*R[1][0]*R[2][0];
	Q[3][1] = 2.0*R[1][1]*R[2][1];
	Q[3][2] = 2.0*R[1][2]*R[2][2];
	Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
	Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
	Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];

	Q[4][0] = 2.0*R[0][0]*R[2][0];
	Q[4][1] = 2.0*R[0][1]*R[2][1];
	Q[4][2] = 2.0*R[0][2]*R[2][2];
	Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
	Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
	Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];

	Q[5][0] = 2.0*R[0][0]*R[1][0];
	Q[5][1] = 2.0*R[0][1]*R[1][1];
	Q[5][2] = 2.0*R[0][2]*R[1][2];
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
	else {
	  Error("Cannot rotate tensor due to dimensions!",__FILE__,__LINE__);
	}
      }
  }

}
