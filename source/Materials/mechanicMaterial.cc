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
#include <boost/bind.hpp> // TODO what do we need bind here?? - Fabian
#include <boost/lexical_cast.hpp>

#include "mechanicMaterial.hh"
#include "Domain/domain.hh"
#include "DataInOut/Logging/cfslog.hh"

DECLARE_LOG(mat)
DEFINE_LOG(mat, "mat")


namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  MechanicMaterial::MechanicMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Mechanics";
    mHandle_ = mp_->GetNewHandle(true);

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

    mp_->ReleaseHandle(mHandle_);
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
  

  void MechanicMaterial::SetScalar(const std::string& param, MaterialType matType, 
                                   Global::ComplexPart dataType ) {

    MathParser::HandleType actHandle;
    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      if ( dataType == Global::REAL ) {
        actHandle = mp_->GetNewHandle(true);
        scalarStringParamsReal_[matType] = param;
        scalarStringHandlesReal_[matType] = actHandle; 
        mp_->SetExpr(actHandle, param );        
        // Register call-back function:
        // Every time some variable within the expression of the elasticity modulus
        // changes, we have to recalculate the tensor
        mp_->AddExpChangeCallBack(boost::bind(&MechanicMaterial::ComputeFullStiffTensor,
                                              this), actHandle );
      }
      else if (dataType == Global::IMAG ) {
        actHandle = mp_->GetNewHandle(true);
        scalarStringParamsImag_[matType] = param;
        scalarStringHandlesImag_[matType] = actHandle;
        mp_->SetExpr(actHandle, param );
        // Register call-back function:
        // Every time some variable within the expression of the elasticity modulus
        // changes, we have to recalculate the tensor
        mp_->AddExpChangeCallBack(boost::bind(&MechanicMaterial::ComputeFullStiffTensor,
                                              this), actHandle );
      }
      else {
        std::string msg = "SetScalar-Double";
        dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
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

      Complex val, oldVal;
      oldVal = scalarParams_[matType];
      if ( dataType == Global::REAL ) {
        val = Complex ( param, oldVal.imag() );
      }
      else if (dataType == Global::IMAG ) {
        val = Complex ( oldVal.real(), param );
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

  void MechanicMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(Matrix<Double>& out, Double emod, Double poi)
  {
    Complex EModul(emod);
    Complex poisson(poi); 
    Complex LameLambda = (poisson*EModul)/
          ((Complex(1.0,0) + poisson)*(Complex(1.0,0)  - Complex(2.0,0)*poisson));
    Complex LameMu = (EModul)/(Complex(2.0,0)*(Complex(1.0)+poisson));
    
    Matrix<Complex> elasticityTensor;
    CalcComplexIsotropicStiffnessTensor(elasticityTensor, LameLambda, LameMu);
    out = elasticityTensor.GetPart(Global::REAL);
  }

  void MechanicMaterial::CalcIsotropicStiffnessTensorFromLame(Matrix<Double>& out, Double lambda, Double mu)
  { 
    Matrix<Complex> elasticityTensor;
    CalcComplexIsotropicStiffnessTensor(elasticityTensor, static_cast<Complex>(lambda), static_cast<Complex>(mu));
    
    out = elasticityTensor.GetPart(Global::REAL);
  }

  double MechanicMaterial::CalcIsotropyError(const Matrix<double>& tensor, SubTensorType stt)
  {
    double v = CalcIsotropicPoissonsRatio(tensor, stt);
    double E = CalcIsotropicYoungsModulus(tensor, stt);

    // this is FULL
    Matrix<double> full_hom;
    CalcIsotropicStiffnessTensorFromEAndPoisson(full_hom, E, v);
    // eventually reduce
    Matrix<double> hom;
    ComputeSubTensor(hom, stt, full_hom);

    LOG_DBG(mat) << "MM::CIE E=" << E << " v=" << v << " err=" << hom.DiffNormL1(tensor);
    LOG_DBG2(mat) << "MM::CIE tensor=" << tensor.ToString(0,1);
    LOG_DBG2(mat) << "MM::CIE full_hom=" << full_hom.ToString(0,1);
    LOG_DBG2(mat) << "MM::CIE hom=" << hom.ToString(0,1);

    return hom.DiffNormL1(tensor);
  }

  double MechanicMaterial::CalcIsotropicPoissonsRatio(const Matrix<double>& tensor, SubTensorType subTensor)
  {
    assert(tensor.GetNumCols() == 3 || tensor.GetNumCols() == 6);

    double E11 = tensor[0][0];
    double E12 = tensor[0][1];

    switch(subTensor)
    {
    case FULL:
    case PLANE_STRAIN:
      return E12 / (E11 + E12);

    case PLANE_STRESS:
      return E12 / E11;

    default:
      EXCEPTION("fail");
    }
  }


  double MechanicMaterial::CalcIsotropicYoungsModulus(const Matrix<double>& tensor, SubTensorType subTensor)
  {
    double E11 = tensor[0][0];
    double v = CalcIsotropicPoissonsRatio(tensor, subTensor);

    switch(subTensor)
    {
    case FULL:
    case PLANE_STRAIN:
      E11 *= (1.0 + v) * (1.0 - 2.0 * v) / (1.0 - v);
      break;

    case PLANE_STRESS:
      E11 *= (1.0 - v*v);
      break;

    default:
      assert(false);
    }

    return E11;
  }

  

  void MechanicMaterial::CalcComplexIsotropicStiffnessTensor(Matrix<Complex>& out, Complex LameLambda, Complex LameMu)
  {
    out.Resize(6);
    out.Init();
        
    out[0][0] = LameLambda + Complex(2.0,0) * LameMu;
    out[1][1] = LameLambda + Complex(2.0,0) * LameMu;
    
    out[0][1] = LameLambda;
    out[1][0] = LameLambda;
    
    out[2][2] = LameMu;

    out[0][2] = LameLambda;
    out[1][2] = LameLambda;
    out[2][0] = LameLambda;
    out[2][1] = LameLambda;

    out[2][2] = LameLambda + Complex(2.0,0) * LameMu;
    out[3][3] = LameMu;
    out[4][4] = LameMu;
    out[5][5] = LameMu;
  }
  

  void MechanicMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix, MaterialType matType, SubTensorType subTensor) const
  {
    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    Matrix<Complex> const &mat = pos->second;

    ComputeSubTensor<Complex>(matMatrix, subTensor, mat);
  }

  template<class T>
  void MechanicMaterial::ComputeSubTensor(Matrix<T>& matMatrix, SubTensorType subTensor, const Matrix<T>& mat)
  {
    switch(subTensor)
    {
    case AXI:
    {
      UInt nrElemsAxi = 4;
      matMatrix.Resize( nrElemsAxi, nrElemsAxi );
      matMatrix.Init();

      // indices of rows and lines for xy-plane (rr, zz, rz, phiphi)
      UInt rowPtr[] = {1,2,6,3};  
      for ( UInt i=0; i<nrElemsAxi; i++ )
        for ( UInt j=0; j<nrElemsAxi; j++ )
          matMatrix[i][j] = mat[rowPtr[i]-1][rowPtr[j]-1];
      break;
    }
    case PLANE_STRAIN:
    {
      UInt nrElems = 3;
      matMatrix.Resize(nrElems, nrElems);
      matMatrix.Init();

      // indices of rows and lines for xy-plane (xx, yy, xy)
      UInt rowPtr[] = {1,2,6}; 
      for ( UInt i=0; i<nrElems; i++ )
        for ( UInt j=0; j<nrElems; j++ )
          matMatrix[i][j] = mat[rowPtr[i]-1][rowPtr[j]-1];
      break;
    }
    case PLANE_STRESS:
    {
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
      break;
    }
    case FULL:
      matMatrix = mat; // copy as nothing changes
      break;
    default:
      // PLAIN is unspecific
      subTensorNotAvailable(NO_MATERIAL, subTensor); // shall be clear
    }
  }
 

  void MechanicMaterial::ComputeFullStiffTensor() {

      Matrix<Complex> elasticityTensor;

    switch(symmetryType_)
    {
    case GENERAL:
      // check that general stiffness tensor is present
      break;
    case ISOTROPIC:
    {

      // ====================================================================
      // Hard coded section for isotropic material data with variable 
      // coefficients
      // ====================================================================
      
      // Loop over real / imag stiffnes parameters, evaluate string
      // representation and store it in double map
      handleMap::iterator it = scalarStringHandlesReal_.begin();
      for ( ; it != scalarStringHandlesReal_.end(); ++it ) {
        SetScalar( mp_->Eval(it->second),it->first, Global::REAL );
      }
      it = scalarStringHandlesImag_.begin();
      for ( ; it != scalarStringHandlesImag_.end(); ++it ) {
        SetScalar( mp_->Eval(it->second),it->first, Global::IMAG );
      }
      // ====================================================================
      
      // get complex valued values
      Complex EModul, poisson;

      GetScalar( EModul, MECH_EMODULUS, Global::COMPLEX ); 
      GetScalar( poisson, MECH_POISSON, Global::COMPLEX );
      Complex LameLambda = (poisson*EModul)/
          ((Complex(1.0,0) + poisson)*
              (Complex(1.0,0)  - Complex(2.0,0)*poisson));
      Complex LameMu = (EModul)/(Complex(2.0,0)*(Complex(1.0)+poisson));

      CalcComplexIsotropicStiffnessTensor(elasticityTensor, LameLambda, LameMu);
      SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, Global::COMPLEX );
      SetScalar(LameLambda, MECH_LAME_LAMBDA, Global::COMPLEX);
      SetScalar(LameMu, MECH_LAME_MU, Global::COMPLEX);
      break;
    }
    case ORTHOTROPIC:
    {
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

      break;
    }
    default:
      EXCEPTION( "Calculation of full stiffness matrix for symmetryType '"
                       << symmetryType_ << "' not implemented!" );
    }
    
  }

  StdVector<double> MechanicMaterial::CalcOrthotropeYoungsModulus(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol)
  {
    Matrix<double> D;
    tensor.Invert(D);

    assert((tensor.GetNumRows() == 6 && stt == FULL) || (tensor.GetNumRows() == 3 && stt == PLANE_STRAIN)); // no plain-stress yet!

    StdVector<double> res(3);

    if(stt == FULL)
    {
    for(UInt i = 0; i < 3; i++)
      res[i] = 1.0/D[i][i];
    }
    else
    {
      assert(mat != NULL && vol > 0 && vol <= 1.0);
      //E1 = 1/(B(1,1) + nusteg^2/rho/Esteg )
      //E2 = 1/(B(2,2) + nusteg^2/rho/Esteg )
      //E3 = rho*Esteg
      // assert(CalcIsotropyError(tens, stt) < 1e-3);
      // core properties
      // core properties
      double E_core, v_core;
      mat->GetScalar(E_core, MECH_EMODULUS, Global::REAL);
      mat->GetScalar(v_core, MECH_POISSON, Global::REAL);

      double E3 = E_core * vol;

      res[0] = 1.0 / (D[0][0] + v_core*v_core / E3);
      res[1] = 1.0 / (D[1][1] + v_core*v_core / E3);
      res[2] = E3;
    }

    return res;
  }

  StdVector<double> MechanicMaterial::CalcOrthotropePoissonsRatio(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol)
  {
    Matrix<double> D;
    tensor.Invert(D);

    StdVector<double> res;
    StdVector<double> E = CalcOrthotropeYoungsModulus(tensor, mat, stt, vol); // does the asserts!
    // v_21, v_12, v_31, v_13, v_32, v_23
    switch(stt)
    {
    case FULL:
      res.Push_back(-1.0 * D[1-1][2-1] * E[2-1]); // v_21
      res.Push_back(-1.0 * D[2-1][1-1] * E[1-1]); // v_12
      res.Push_back(-1.0 * D[1-1][3-1] * E[3-1]); // v_31
      res.Push_back(-1.0 * D[3-1][1-1] * E[1-1]); // v_13
      res.Push_back(-1.0 * D[2-1][3-1] * E[3-1]); // v_32
      res.Push_back(-1.0 * D[3-1][2-1] * E[2-1]); // v_23
      break;

    case PLANE_STRAIN:
    {
      // core properties
      double E_core, v_core;
      mat->GetScalar(E_core, MECH_EMODULUS, Global::REAL);
      mat->GetScalar(v_core, MECH_POISSON, Global::REAL);

      // % the two 2D poisson ratios and Es have been derived by Bastian and Fabian S.
      // v12 = (-B(1,2)-nusteg^2/rho/Esteg)*E1
      // v21 = E2/E1*v12
      double v_12 = (-D[1-1][2-1] - v_core * v_core / E[3-1]) * E[1-1];
      double v_21 = (E[2-1] / E[1-1]) * v_12;

      // % the remaining poisson ratios
      // v13 = E1*nusteg/rho/Esteg
      // v31 = E3/E1*v13
      double v_13 = E[1-1] * v_core / E[3-1];
      double v_31 = (E[3-1] / E[1-1]) * v_13;
      // v23 = E2*nusteg/rho/Esteg
      // v32 = E3/E2*v23
      double v_23 = E[2-1] * v_core / E[3-1];
      double v_32 = (E[3-1] / E[2-1]) * v_23;

      res.Push_back(v_21);
      res.Push_back(v_12);
      res.Push_back(v_31);
      res.Push_back(v_13);
      res.Push_back(v_32);
      res.Push_back(v_23);
      break;
    }

    default:
      assert(false);
    }

    return res;
  }

  double MechanicMaterial::CalcOrthotropeError(const Matrix<double>& tensor)
  {
    // measure zeros in 1-norm
    double err = 0.0;
    if(tensor.GetNumRows() == 3)
    {
       err += abs(tensor[1-1][3-1]);
       err += abs(tensor[2-1][3-1]);
    }
    else
    {
      for(UInt r = 1; r <= 5; r++)
      {
        for(UInt c = 4; c <= 6; c++)
        {
          // only above the diagonal
          if(r >= c) continue;
          err += abs(tensor[r-1][c-1]);
        }
      }
    }
    return err;
  }

  StdVector<std::pair<std::string, double> > MechanicMaterial::CalcOrthotropeProperties(const Matrix<double>& tensor, BaseMaterial* mat, SubTensorType stt, double vol)
  {
    Matrix<double> D;
    tensor.Invert(D);

    UInt dim = tensor.GetNumRows() == 3 ? 2 : 3; // rank 3 or 6

    StdVector<std::pair<std::string, double> > res;

    StdVector<double> E = CalcOrthotropeYoungsModulus(tensor, mat, stt, vol);

    res.Push_back(std::make_pair("E_1", E[0]));
    res.Push_back(std::make_pair("E_2", E[1]));
    res.Push_back(std::make_pair("E_3", E[2]));

    StdVector<double> v = CalcOrthotropePoissonsRatio(tensor, mat, stt, vol);
    // v_21=0, v_12=1, v_31=2, v_13=3, v_32=4, v_23=5
    res.Push_back(std::make_pair("v_21", v[0]));
    res.Push_back(std::make_pair("v_12", v[1]));
    res.Push_back(std::make_pair("v_31", v[2]));
    res.Push_back(std::make_pair("v_13", v[3]));
    res.Push_back(std::make_pair("v_32", v[4]));
    res.Push_back(std::make_pair("v_23", v[5]));

    if(dim == 2)
    {
      res.Push_back(std::make_pair("G_12", 1.0 / D[3-1][3-1]));
    }
    else
    {
      res.Push_back(std::make_pair("G_23", 1.0 / D[4-1][4-1]));
      res.Push_back(std::make_pair("G_13", 1.0 / D[5-1][5-1]));
      res.Push_back(std::make_pair("G_12", 1.0 / D[6-1][6-1]));
    }

    double err = CalcOrthotropeError(tensor);
    res.Push_back(std::make_pair("err", err));

    return res;
  }

  // required in ErsatzMaterial. The complex version is explictely called here
  template void MechanicMaterial::ComputeSubTensor<Double>(Matrix<Double>& matMatrix, SubTensorType subTensor, const Matrix<Double>& mat);

  
} // end of namespace
