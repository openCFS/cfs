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

#include "mechanicMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  MechanicMaterial::MechanicMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Mechanics";

    //set the allowed material parameters
    isAllowed_.insert( DENSITY );
    isAllowed_.insert( MECH_STIFFNESS_TENSOR );
    isAllowed_.insert( COEFF_STRAIN_IRREVERSIBLE );
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
    isAllowed_.insert( MECH_LAME_LAMBDA );
    isAllowed_.insert( MECH_LAME_MU );
    isAllowed_.insert( RAYLEIGH_ALPHA );
    isAllowed_.insert( RAYLEIGH_BETA );
    isAllowed_.insert( RAYLEIGH_FREQUENCY);
    isAllowed_.insert( LOSS_TANGENS_DELTA);
    isAllowed_.insert( ACOU_ALPHA );
    isAllowed_.insert( FRACTIONAL_EXPONENT );
    isAllowed_.insert( NONLIN_COEFFICIENT );
    isAllowed_.insert( NONLIN_DEPENDENCY );
    isAllowed_.insert( NONLIN_APPROXIMATION_TYPE );
    isAllowed_.insert( NONLIN_DATA_NAME );

  }

  MechanicMaterial::~MechanicMaterial() {


  }

  void MechanicMaterial::Finalize() {


    // Trigger calculation of stiffness tensor
    ComputeFullStiffTensor();

  }

  void MechanicMaterial::SetScalar(const std::string& param, MaterialType matType) {


      
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
  

  void MechanicMaterial::SetScalar(Double param, MaterialType matType, 
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


  void MechanicMaterial::SetScalar( Complex param, MaterialType matType, 
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


  void MechanicMaterial::SetVector(const Vector<Double>& param, MaterialType matType, 
				    Global::ComplexPart dataType ) {
    

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "vector";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      if ( dataType == Global::REAL || dataType == Global::IMAG ) {
	if ( vectorParams_[matType].GetSize() == 0 ) {
	  vectorParams_[matType].Resize( param.GetSize() );
          vectorParams_[matType].Init();
	}

	vectorParams_[matType].SetPart( dataType, param );

	if ( dataType == Global::IMAG ) {
	  isComplex_.insert( matType );
	}
      }
      else {
	std::string msg = "SeVector-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
  }


  void MechanicMaterial::SetTensor(const Matrix<Double>& param, MaterialType matType, 
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

        if ( dataType == Global::IMAG ) {
          isComplex_.insert( matType );
        }
      }
      else {
        std::string msg = "SetTensor-Double";
        dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
  }

  void MechanicMaterial::SetTensor(const Matrix<Complex>& param, MaterialType matType, 
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
  }

  void MechanicMaterial::GetScalar( std::string& param, MaterialType matType)  const {


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
 
   void MechanicMaterial::GetScalar( Integer& param, MaterialType matType)  const {
    
    
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
  


  void MechanicMaterial::GetScalar( Double& param, MaterialType matType, 
				    Global::ComplexPart dataType )  const {


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

  void MechanicMaterial::GetScalar( Complex& param, MaterialType matType, 
				    Global::ComplexPart dataType )  const {

    scalarMap::const_iterator pos;
    pos = scalarParams_.find( matType );

    if ( pos == scalarParams_.end() ) {
      std::string dim = "tensor";
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


  void MechanicMaterial::GetVector( Vector<Double>& param, 
				    MaterialType matType, 
				    Global::ComplexPart dataType ) const {
    

    vectorMap::const_iterator pos;
    pos = vectorParams_.find( matType );

    if ( pos == vectorParams_.end() ) {
      std::string dim = "vector";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Vector<Complex> matVector;
      matVector = pos->second;

      if ( dataType == Global::REAL || dataType == Global::IMAG) {
	param = matVector.GetPart( dataType );
      }
      else {
	std::string msg = "GetVector-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }


  void MechanicMaterial::GetTensor( Matrix<Double>& param, 
				    MaterialType matType, 
				    Global::ComplexPart dataType,
				    SubTensorType subTensor ) const {
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

  void MechanicMaterial::GetTensor( Matrix<Complex>& param, 
				    MaterialType matType, 
				    Global::ComplexPart dataType,
				    SubTensorType subTensor ) const {	
    

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

      if ( dataType == Global::REAL || dataType == Global::IMAG) {
	Matrix<Double> help; 
	help = matTensor.GetPart( dataType );
	param.Resize( matTensor.GetNumRows(), matTensor.GetNumCols() );
        param.Init();
	param.SetPart( dataType, help );
      }
      else if ( dataType == Global::COMPLEX ) {
	param = matTensor;
      }
    }
  }
  

  void MechanicMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
					  MaterialType matType, 
					  SubTensorType subTensor) const {


    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    Matrix<Complex> const &mat = pos->second;

    if ( subTensor == AXI ) {
      UInt nrElemsAxi = 4;
      matMatrix.Resize( nrElemsAxi, nrElemsAxi );
      matMatrix.Init();

      // indices of rows and lines for xy-plane (rr, zz, rz, phiphi)
      UInt rowPtr[] = {1,2,6,3};  
      for ( UInt i=0; i<nrElemsAxi; i++ )
	for ( UInt j=0; j<nrElemsAxi; j++ )
	  matMatrix[i][j] = mat[rowPtr[i]-1][rowPtr[j]-1];
    }
    else if ( subTensor == PLANE_STRAIN ) {
      UInt nrElems = 3;
      matMatrix.Resize( nrElems, nrElems    );
      matMatrix.Init();

      // indices of rows and lines for xy-plane (xx, yy, xy)
      UInt rowPtr[] = {1,2,6}; 
      for ( UInt i=0; i<nrElems; i++ )
	for ( UInt j=0; j<nrElems; j++ )
	  matMatrix[i][j] = mat[rowPtr[i]-1][rowPtr[j]-1];
    }
    else if ( subTensor == PLANE_STRESS ) {
      UInt nrElems = 3;
      matMatrix.Resize( nrElems, nrElems );
      matMatrix.Init();

      if ( abs(mat[0][0]) < 1.09E-15 ) {
	EXCEPTION("Singular material tensor when computing plane stress case" );
      }

      // calculate plane stress matrix for xy-plane
      matMatrix[0][0] = mat[0][0] - mat[2][0]*mat[0][2]/mat[2][2];
      matMatrix[0][1] = mat[0][1] - mat[2][1]*mat[0][2]/mat[2][2];
      matMatrix[0][2] = mat[0][5];
      matMatrix[1][0] = mat[1][0] - mat[2][0]*mat[1][2]/mat[2][2];
      matMatrix[1][1] = mat[1][1] - mat[2][1]*mat[1][2]/mat[2][2];
      matMatrix[1][2] = mat[1][5];
      matMatrix[2][0] = mat[5][0];
      matMatrix[2][1] = mat[5][1];
      matMatrix[2][2] = mat[5][5];


      //explicite computation (old case for yz-plane
//       matMatrix[0][0] = mat[1][1] - mat[0][1]*mat[1][0]/mat[0][0];
//       matMatrix[0][1] = mat[1][2] - mat[0][2]*mat[1][0]/mat[0][0];
//       matMatrix[0][2] = mat[1][3] - mat[0][3]*mat[1][0]/mat[0][0];
//       matMatrix[1][0] = mat[2][1] - mat[0][1]*mat[2][0]/mat[0][0];
//       matMatrix[1][1] = mat[2][2] - mat[0][2]*mat[2][0]/mat[0][0];
//       matMatrix[1][2] = mat[2][3] - mat[0][3]*mat[2][0]/mat[0][0];
//       matMatrix[2][0] = mat[3][1] - mat[0][1]*mat[3][0]/mat[0][0];
//       matMatrix[2][1] = mat[3][2] - mat[0][2]*mat[3][0]/mat[0][0];
//       matMatrix[2][2] = mat[3][3] - mat[0][3]*mat[3][0]/mat[0][0]
        ;
      
    }
    else {
      subTensorNotAvailable( matType, subTensor );
    }
  }
 

  void MechanicMaterial::ComputeFullStiffTensor() {

      Matrix<Complex> elasticityTensor;

    
    if( symmetryType_ == GENERAL ) {

      // check that general stiffness tensor is present
 
    } else if( symmetryType_ == ISOTROPIC ) {

      // get complex valued values
      Complex EModul, poisson;
        
      GetScalar( EModul, MECH_EMODULUS, Global::COMPLEX ); 
      GetScalar( poisson, MECH_POISSON, Global::COMPLEX ); 
  
      // calculate isotropic case        
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

      SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, Global::COMPLEX );
      SetScalar(LameLambda, MECH_LAME_LAMBDA, Global::COMPLEX);
      SetScalar(LameMu, MECH_LAME_MU, Global::COMPLEX);



    } else if( symmetryType_ == ORTHOTROPIC ) {
        
      Complex EX, EY, EZ, nuXY, nuYZ, nuXZ, GYZ, GZX, GXY;
      GetScalar( EX, MECH_EMODULUS_X, Global::COMPLEX ); 
      GetScalar( EY, MECH_EMODULUS_Y, Global::COMPLEX ); 
      GetScalar( EZ, MECH_EMODULUS_Z, Global::COMPLEX ); 
      GetScalar( nuXY, MECH_POISSON_XY, Global::COMPLEX ); 
      GetScalar( nuYZ, MECH_POISSON_YZ, Global::COMPLEX ); 
      GetScalar( nuXZ, MECH_POISSON_XZ, Global::COMPLEX ); 
      GetScalar( GYZ, MECH_GMODULUS_YZ, Global::COMPLEX ); 
      GetScalar( GZX, MECH_GMODULUS_ZX, Global::COMPLEX ); 
      GetScalar( GXY, MECH_GMODULUS_XY, Global::COMPLEX ); 
        
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
      SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, Global::COMPLEX );
    } else {
      EXCEPTION( "Calculation of full stiffness matrix for symmetryType '"
                 << symmetryType_ << "' not implemented!" );
    }
    
  }
}
