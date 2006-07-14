 #include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "mechanicMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  MechanicMaterial::MechanicMaterial() : BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");
    materialDatabaseName_ = "Mechanics";

    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( MECH_STIFFNESS_TENSOR );
    isAllowed_.insert( MECH_EMODULUS );
    isAllowed_.insert( MECH_EMODULUS_X );
    isAllowed_.insert( MECH_EMODULUS_Y );
    isAllowed_.insert( MECH_EMODULUS_Z );
    isAllowed_.insert( MECH_POISSON );
    isAllowed_.insert( MECH_POISSON_XY );
    isAllowed_.insert( MECH_POISSON_YZ );
    isAllowed_.insert( MECH_POISSON_XZ );
    isAllowed_.insert( MECH_GMODULUS_YZ );
    isAllowed_.insert( MECH_GMODULUS_ZX );
    isAllowed_.insert( MECH_GMODULUS_XY );
    isAllowed_.insert( RAYLEIGH_ALPHA );
    isAllowed_.insert( RAYLEIGH_BETA );
    isAllowed_.insert( RAYLEIGH_FREQUENCY);
    isAllowed_.insert( LOSS_TANGENS_DELTA);
    isAllowed_.insert( ACOU_ALPHA );
    isAllowed_.insert( FRACTIONAL_EXPONENT );

  }

  MechanicMaterial::~MechanicMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }

  void MechanicMaterial::Finalize() {
    ENTER_FCN( "MechanicMaterial::Finalize" );


    // Trigger calculation of stiffness tensor
    ComputeFullStiffTensor();

  }

  void MechanicMaterial::SetScalar( std::string& param, const MaterialType& matType) {

    ENTER_FCN( "AcousticMaterial::SetScalar" );
    Error("SetScalar for 'String' not implemented",__FILE__,__LINE__);
  }


  void MechanicMaterial::SetScalar( Double& param, const MaterialType& matType, 
				    const DataType& dataType ) {

    ENTER_FCN( "MechanicMaterial::SetScalar" );

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
  }


  void MechanicMaterial::SetScalar( Complex& param, const MaterialType& matType, 
				    const DataType& dataType ) {

    ENTER_FCN( "MechanicMaterial::SetScalar" );

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
  }

  void MechanicMaterial::SetTensor( Matrix<Double>& param, const MaterialType& matType, 
				    const DataType& dataType ) {
    
    ENTER_FCN( "MechanicMaterial::SetTensor" );

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

	if ( dataType == IMAG ) {
	  isComplex_.insert( matType );
	}
      }
      else {
	std::string msg = "SetTensor-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
  }

  void MechanicMaterial::SetTensor( Matrix<Complex>& param, const MaterialType& matType, 
				    const DataType& dataType ) {
    
    ENTER_FCN( "MechanicMaterial::SetTensor" );

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
  }


  void MechanicMaterial::GetScalar( Double& param, const MaterialType& matType, 
				    const DataType& dataType )  const {

    ENTER_FCN( "MechanicMaterial::GetScalar" );

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

  void MechanicMaterial::GetScalar( Complex& param, const MaterialType& matType, 
				    const DataType& dataType )  const {

    ENTER_FCN( "MechanicMaterial::GetScalar" );
    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "tensor";
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

  void MechanicMaterial::GetTensor( Matrix<Double>& param, 
				    const MaterialType& matType, 
				    const DataType& dataType,
				    const SubTensorType subTensor ) const {
    
    ENTER_FCN( "MechanicMaterial::GetTensor" );

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

  void MechanicMaterial::GetTensor( Matrix<Complex>& param, 
				    const MaterialType& matType, 
				    const DataType& dataType,
				    const SubTensorType subTensor ) const {	
    
    ENTER_FCN( "MechanicMaterial::GetTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "scalar";
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
        param.Init();
	param.SetPart( dataType, help );
      }
      else if ( dataType == COMPLEX ) {
	param = matTensor;
      }
    }
  }
  

  void MechanicMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
					  const MaterialType& matType, 
					  const SubTensorType& subTensor) const {

    ENTER_FCN( "MechanicMaterial::ComputeSubTensor" );

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    Matrix<Complex> const &mat = pos->second;

    if ( subTensor == AXI ) {
      UInt nrElemsAxi = 4;
      matMatrix.Resize( nrElemsAxi, nrElemsAxi );
      matMatrix.Init();

      UInt rowPtr[] = {2,3,4,1};  // indices of rows and lines for yz-plane
      for ( UInt i=0; i<nrElemsAxi; i++ )
	for ( UInt j=0; j<nrElemsAxi; j++ )
	  matMatrix[i][j] = mat[rowPtr[i]-1][rowPtr[j]-1];
    }
    else if ( subTensor == PLANE_STRAIN ) {
      UInt nrElems = 3;
      matMatrix.Resize( nrElems, nrElems    );
      matMatrix.Init();

      UInt rowPtr[] = {2,3,4};  // indices of rows and lines for yz-plane
      for ( UInt i=0; i<nrElems; i++ )
	for ( UInt j=0; j<nrElems; j++ )
	  matMatrix[i][j] = mat[rowPtr[i]-1][rowPtr[j]-1];
    }
    else if ( subTensor == PLANE_STRESS ) {
      UInt nrElems = 3;
      matMatrix.Resize( nrElems, nrElems );
      matMatrix.Init();

      if ( abs(mat[0][0]) < 1.09E-15 ) {
	Error("Singular material tensor when computing plane stress case",
	      __FILE__,__LINE__);
      }

      //explicite computation
      matMatrix[0][0] = mat[1][1] - mat[0][1]*mat[1][0]/mat[0][0];
      matMatrix[0][1] = mat[1][2] - mat[0][2]*mat[1][0]/mat[0][0];
      matMatrix[0][2] = mat[1][3] - mat[0][3]*mat[1][0]/mat[0][0];
      matMatrix[1][0] = mat[2][1] - mat[0][1]*mat[2][0]/mat[0][0];
      matMatrix[1][1] = mat[2][2] - mat[0][2]*mat[2][0]/mat[0][0];
      matMatrix[1][2] = mat[2][3] - mat[0][3]*mat[2][0]/mat[0][0];
      matMatrix[2][0] = mat[3][1] - mat[0][1]*mat[3][0]/mat[0][0];
      matMatrix[2][1] = mat[3][2] - mat[0][2]*mat[3][0]/mat[0][0];
      matMatrix[2][2] = mat[3][3] - mat[0][3]*mat[3][0]/mat[0][0];
      //std::cout << "MatMatrix:\n" << matMatrix << std::endl;
    }
    else {
      subTensorNotAvailable( matType, subTensor );
    }
  }
 

  void MechanicMaterial::ComputeFullStiffTensor() {
    ENTER_FCN( "MechanicMaterial::ComputeFullStiffnessTensor" ) 

      Matrix<Complex> elasticityTensor;

    
    if( symmetryType_ == GENERAL ) {

      // check that general stiffness tensor is present
 
    } else if( symmetryType_ == ISOTROPIC ) {

      // get complex valued values
      Complex EModul, poisson;
        
      GetScalar( EModul, MECH_EMODULUS, COMPLEX ); 
      GetScalar( poisson, MECH_POISSON, COMPLEX ); 
  
      // calculate isothropic case        
      Complex LameLambda, LameMu;
      elasticityTensor.Resize(6,6);
      elasticityTensor.Init();
      LameLambda = (poisson*EModul)/
        ((Complex(1.0,0) + poisson)*
         (Complex(1.0,0)  - Complex(2.0,0)*poisson));
      LameMu = (EModul)/(Complex(2.0,0)*(Complex(1.0)+poisson));
      
      elasticityTensor[0][0]=LameLambda+Complex(2.0,0)*LameMu;
      elasticityTensor[1][1]=LameLambda+Complex(2.0,0)*LameMu;
      elasticityTensor[2][2]=LameLambda+Complex(2.0,0)*LameMu;
      
      elasticityTensor[0][1]=LameLambda;
      elasticityTensor[0][2]=LameLambda;
      elasticityTensor[1][0]=LameLambda;
      elasticityTensor[1][2]=LameLambda;
      elasticityTensor[2][0]=LameLambda;
      elasticityTensor[2][1]=LameLambda;
      
      elasticityTensor[3][3]=LameMu;
      elasticityTensor[4][4]=LameMu;
      elasticityTensor[5][5]=LameMu;

      SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, COMPLEX ); 



    } else if( symmetryType_ == ORTHOTROPIC ) {
        
      Complex EX, EY, EZ, nuXY, nuYZ, nuXZ, GYZ, GZX, GXY;
      GetScalar( EX, MECH_EMODULUS_X, COMPLEX ); 
      GetScalar( EY, MECH_EMODULUS_Y, COMPLEX ); 
      GetScalar( EZ, MECH_EMODULUS_Z, COMPLEX ); 
      GetScalar( nuXY, MECH_POISSON_XY, COMPLEX ); 
      GetScalar( nuYZ, MECH_POISSON_YZ, COMPLEX ); 
      GetScalar( nuXZ, MECH_POISSON_XZ, COMPLEX ); 
      GetScalar( GYZ, MECH_GMODULUS_YZ, COMPLEX ); 
      GetScalar( GZX, MECH_GMODULUS_ZX, COMPLEX ); 
      GetScalar( GXY, MECH_GMODULUS_XY, COMPLEX ); 
        
      Complex nuYX, nuZY, nuZX, aux;
      nuYX=(EY/EX)*nuXY;
      nuZY=(EZ/EY)*nuYZ;
      nuZX=(EZ/EX)*nuXZ;
        
      aux= (Complex(1,0)-nuXY*nuYX-nuYZ*nuZY-nuXZ*nuZX-
            Complex(2.0,0)*nuYX*nuZY*nuXZ)/(EX*EY*EZ);
        
      elasticityTensor.Resize(6,6);
      elasticityTensor.Init();
        
      elasticityTensor[0][0]=(Complex(1,0)-nuYZ*nuZY)/(EY*EZ*aux);
      elasticityTensor[1][1]=(Complex(1,0)-nuXZ*nuZX)/(EX*EZ*aux);
      elasticityTensor[2][2]=(Complex(1,0)-nuXY*nuYX)/(EX*EY*aux);
        
      elasticityTensor[0][1]=(nuYX+nuZX*nuYZ)/(EY*EZ*aux);
      elasticityTensor[0][2]=(nuZX+nuYX*nuZY)/(EY*EZ*aux);
      elasticityTensor[1][0]=(nuYX+nuZX*nuYZ)/(EY*EZ*aux);
      elasticityTensor[1][2]=(nuZY+nuXY*nuZX)/(EX*EZ*aux);
      elasticityTensor[2][0]=(nuZX+nuYX*nuZY)/(EY*EZ*aux);
      elasticityTensor[2][1]=(nuZY+nuXY*nuZX)/(EX*EZ*aux);
        
      elasticityTensor[3][3]=GYZ;
      elasticityTensor[4][4]=GZX;
      elasticityTensor[5][5]=GXY;
      SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, COMPLEX );
    } else {
      *error << "Calculation of full stiffness matrix for symmetryType '"
             << symmetryType_ << "' not implemented!";
      Error( __FILE__, __LINE__ );
    }
    
  }
}
