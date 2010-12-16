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
  
  void MechanicMaterial::CalcIsotropicStiffnessTensorFromEAndPoisson(Matrix<Double>& out, Double emod, Double poi, UInt dim)
  {
    Complex EModul(emod);
    Complex poisson(poi); 
    Complex LameLambda = (poisson*EModul)/
          ((Complex(1.0,0) + poisson)*(Complex(1.0,0)  - Complex(2.0,0)*poisson));
    Complex LameMu = (EModul)/(Complex(2.0,0)*(Complex(1.0)+poisson));
    
    Matrix<Complex> elasticityTensor;
    CalcComplexIsotropicStiffnessTensor(elasticityTensor, LameLambda, LameMu, dim);
    out = elasticityTensor.GetPart(Global::REAL);
  }

  void MechanicMaterial::CalcIsotropicStiffnessTensorFromLame(Matrix<Double>& out, Double lambda, Double mu, UInt dim)
  { 
    Matrix<Complex> elasticityTensor;
    CalcComplexIsotropicStiffnessTensor(elasticityTensor, static_cast<Complex>(lambda), static_cast<Complex>(mu), dim);
    
    out = elasticityTensor.GetPart(Global::REAL);
  }

  double MechanicMaterial::CalcIsotropyError(const Matrix<double>& E, bool normed)
  {
     assert((E.GetNumRows() == 3 && E.GetNumCols() == 3)
         || (E.GetNumRows() == 6 && E.GetNumCols() == 6));

     bool D3 = E.GetNumRows() == 6;

     double E11 = E[0][0];
     double E22 = E[1][1];
     double E33 = E[2][2];
     double E44 = D3 ? E[3][3] : 0.0;
     double E55 = D3 ? E[4][4] : 0.0;
     double E66 = D3 ? E[5][5] : 0.0;

     double E12 = E[0][1];
     double E13 = E[0][2];
     double E14 = D3 ? E[0][3] : 0.0;
     double E15 = D3 ? E[0][4] : 0.0;
     double E16 = D3 ? E[0][5] : 0.0;

     double E23 = E[1][2];
     double E24 = D3 ? E[1][3] : 0.0;
     double E25 = D3 ? E[1][4] : 0.0;
     double E26 = D3 ? E[1][5] : 0.0;

     double E34 = D3 ? E[2][3] : 0.0;
     double E35 = D3 ? E[2][4] : 0.0;
     double E36 = D3 ? E[2][5] : 0.0;

     double E45 = D3 ? E[3][4] : 0.0;
     double E46 = D3 ? E[3][5] : 0.0;

     double E56 = D3 ? E[4][5] : 0.0;

     double result = 0;

     if(!D3)
     {
       result =  abs(E11 - E22);          // non-shear diagonal const
       result += abs(E11 - E12 - 2.0*E33);// combination
       result += abs(E13);                // rest is zero
       result += abs(E23);
     }
     else
     {
       // non-shear diagonal is constant E11 = E22 = E33
       result =  abs(E11 - E22);
       result += abs(E22 - E33);

       // upper non-shear triangle is constant
       // E12 = E13 = E23 -> E12 - E13 = 0, E13 - E23 = 0
       result += abs(E12 - E13);
       result += abs(E13 - E23);

       // shear diagonal is constant
       result += abs(E44 - E55);
       result += abs(E55 - E66);

       // relationship of the three unique values
       // E11 - E12 - 2E66 = E11 - E12 - E66 - E66 = 0
       result += abs(E11 - E12 - E66 - E66);

       // the rest is zero
       // E14 = E15 = E16 = E24 = E25 = E26 = E34 = E35 = E36 = E45 = E46 = E56 = 0
       result += abs(E14);
       result += abs(E15);
       result += abs(E16);
       result += abs(E24);
       result += abs(E25);
       result += abs(E26);
       result += abs(E34);
       result += abs(E35);
       result += abs(E36);
       result += abs(E45);
       result += abs(E46);
       result += abs(E56);
     }

     return normed ? result / E11 : result;
  }

  double MechanicMaterial::CalcIsotropicPoissonsRatio(const Matrix<double>& tensor)
  {
    // isotropic values are calculated differently than orthotropic values
    // as a result, the values can differ significantly from the orthotropic ones
    // if the isotropy-error is positive (this is also true for small isotropy-errors!)
    double E11 = tensor[0][0];
    double E12 = tensor[0][1];

    return E12 / (E11 + E12);
  }


  double MechanicMaterial::CalcIsotropicYoungsModulus(const Matrix<double>& tensor)
  {
    // see comment in MechanicMaterial::CalcIsotropicPoissonsRatio
    double E11 = tensor[0][0];
    double v = CalcIsotropicPoissonsRatio(tensor);

    return (E11 * (1.0 + v) * (1.0 - 2.0 * v)) / (1.0 - v);
  }

  

  void MechanicMaterial::CalcComplexIsotropicStiffnessTensor(Matrix<Complex>& out, 
      const Complex &LameLambda, const Complex &LameMu, UInt dim)
  {
    out.Resize(dim == 2 ? 3 : 6);
    out.Init();
        
    out[0][0] = LameLambda + Complex(2.0,0) * LameMu;
    out[1][1] = LameLambda + Complex(2.0,0) * LameMu;
    
    out[0][1] = LameLambda;
    out[1][0] = LameLambda;
    
    out[2][2] = LameMu;
    
    if(dim == 3)
    {
      out[0][2] = LameLambda;
      out[1][2] = LameLambda;
      out[2][0] = LameLambda;
      out[2][1] = LameLambda;

      out[2][2] = LameLambda + Complex(2.0,0) * LameMu;
      out[3][3] = LameMu;
      out[4][4] = LameMu;
      out[5][5] = LameMu;
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

      CalcComplexIsotropicStiffnessTensor(elasticityTensor, LameLambda, LameMu, 3);
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

  StdVector<double> MechanicMaterial::CalcOrthotropeYoungsModulus(const Matrix<double>& tensor)
  {
    Matrix<double> D;
    tensor.Invert(D);

    UInt dim = tensor.GetNumRows() == 3 ? 2 : 3; // rank is 3 or 6

    StdVector<double> res(dim);

    for(UInt i = 0; i < dim; i++)
      res[i] = 1.0/D[i][i];

    return res;
  }

  StdVector<double> MechanicMaterial::CalcOrthotropePoissonsRatio(const Matrix<double>& tensor)
  {
    Matrix<double> D;
    tensor.Invert(D);

    UInt dim = tensor.GetNumRows() == 3 ? 2 : 3; // rank is 3 or 6

    StdVector<double> res;
    StdVector<double> E = CalcOrthotropeYoungsModulus(tensor);

    if(dim == 2)
    {
      res.Push_back(-1.0 * D[1-1][2-1] * E[2-1]); // v_21
      res.Push_back(-1.0 * D[2-1][1-1] * E[1-1]); // v_12
    }
    if(dim == 3)
    {
      res.Push_back(-1.0 * D[2-1][3-1] * E[3-1]); // v_32
      res.Push_back(-1.0 * D[1-1][3-1] * E[3-1]); // v_31
      res.Push_back(-1.0 * D[1-1][2-1] * E[2-1]); // v_21
      res.Push_back(-1.0 * D[3-1][2-1] * E[2-1]); // v_23
      res.Push_back(-1.0 * D[3-1][1-1] * E[1-1]); // v_13
      res.Push_back(-1.0 * D[2-1][1-1] * E[1-1]); // v_12
    }

    return res;
  }

  double MechanicMaterial::CalcOrthotropeError(const Matrix<double>& tensor)
  {


    Matrix<double> D;
    tensor.Invert(D);

    StdVector<double> E = CalcOrthotropeYoungsModulus(tensor);
    StdVector<double> v = CalcOrthotropePoissonsRatio(tensor);

    double err = 0.0;
    if(tensor.GetNumRows() == 3)
    {
       err += abs(D[3-1][2-1]);
       err += abs(D[3-1][1-1]);
       err += abs(abs(E[0]*v[1]) - abs(E[1]*v[0]));
    }
    else
    {
      for(UInt r = 1; r <= 5; r++)
        for(UInt c = 4; c <= 6; c++)
        {
          // only above the diagonal
          if(r >= c) continue;
          err += abs(D[r-1][c-1]);
        }
      err += abs(abs(E[0]*v[2]) - abs(E[1]*v[5]));  // E_1 v_21 = E_2 v_12
      err += abs(abs(E[1]*v[0]) - abs(E[2]*v[3]));  // E_2 v_32 = E_3 v_23
      err += abs(abs(E[2]*v[4]) - abs(E[0]*v[1]));  // E_3 v_13 = E_1 v_31
    }
    return err;
  }

  StdVector<std::pair<std::string, double> > MechanicMaterial::CalcOrthotropeProperties(
      const Matrix<double>& tensor)
  {
    Matrix<double> D;
    tensor.Invert(D);

    UInt dim = tensor.GetNumRows() == 3 ? 2 : 3; // rank 3 or 6

    StdVector<std::pair<std::string, double> > res;

    StdVector<double> E = CalcOrthotropeYoungsModulus(tensor);

    res.Push_back(std::make_pair("E_1", E[0]));
    res.Push_back(std::make_pair("E_2", E[1]));
    if(dim == 3)
      res.Push_back(std::make_pair("E_3", E[2]));

    StdVector<double> v = CalcOrthotropePoissonsRatio(tensor);

    if(dim == 2)
    {
      res.Push_back(std::make_pair("v_21", v[0]));
      res.Push_back(std::make_pair("v_12", v[1]));
    }
    else
    {
      res.Push_back(std::make_pair("v_32", v[0]));
      res.Push_back(std::make_pair("v_31", v[1]));
      res.Push_back(std::make_pair("v_21", v[2]));
      res.Push_back(std::make_pair("v_23", v[3]));
      res.Push_back(std::make_pair("v_13", v[4]));
      res.Push_back(std::make_pair("v_12", v[5]));
    }

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



  
} // end of namespace
