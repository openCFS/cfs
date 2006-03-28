#include "XMLMaterialHandler.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

namespace CoupledField {

  XMLMaterialHandler::XMLMaterialHandler( const std::string & fileName )
    : MaterialHandler( fileName) {
    ENTER_FCN( "XMLMaterialHandler::XMLMaterialHandler" );

    // Create new XMLParamHandler and parse material file
    std::string schema = commandLine->GetSchemaPath();
    schema += "/CFS-Material/CFS_Material.xsd";
    parser_ = new XMLParamHandler( fileName.c_str(), schema.c_str() );
    
    
  }
  
  XMLMaterialHandler::~XMLMaterialHandler() {
    ENTER_FCN( "XMLMaterialHandler::~XMLMaterialHandler" ) 
      
      delete parser_;
  }
  
  void XMLMaterialHandler::GetMaterial( BaseMaterial * material, 
                                        const std::string matName, 
                                        const MaterialClass matClass ) {
    ENTER_FCN( "XMLMaterialHandler::GetMaterial");

    if ( matClass == PIEZO ) {
      ReadPiezo( material, matName);
    }
    else if ( matClass == MECHANIC ) {
      ReadMechanic( material, matName );
    }    
    else if ( matClass == FLUID ) {
      ReadAcoustic( material, matName );
    }
    else if ( matClass == ELECTROMAGNETIC ) {
      ReadMagnetic( material, matName );
    }
    else if ( matClass == ELECTROSTATIC ) {
      ReadElectrostatic( material, matName );
    }
    else if ( matClass == THERMIC ) {
      ReadThermic( material, matName );
    }
    else if ( matClass == FLOW ) {
      ReadFlow( material, matName );
    }
    else {
      (*error) << "Warning: material type:" << matClass << " not defined ";
      Error( __FILE__, __LINE__ );
    }
  }
//**********************************************************************
//*************  READ PIEZO ********************************************
//**********************************************************************
  void XMLMaterialHandler::ReadPiezo(BaseMaterial *material,
                                     const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadPiezo" );

    Integer     inteValue;
    std::string striValue;

    Matrix<Double> couplingTensor(3,6);

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    //read real piezo coupling tensor
    const unsigned int dim1=3, dim2=6;
    keyVec = "material" ,"piezo" , "piezoCouplingTensor", "real";
    attrVec= "name"     , ""     , "dim1";
    valVec =  matName   , ""     , "3";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, couplingTensor );
      material->SetTensor( couplingTensor, PIEZO_TENSOR, REAL ); 
      std::cerr << "real couplingTensor=" << std::endl << couplingTensor << std::endl;
    }

    //read imaginary piezo coupling tensor
    keyVec = "material" ,"piezo" , "piezoCouplingTensor", "imag";
    attrVec= "name"     , ""     , "dim1";
    valVec =  matName   , ""     , "3";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, couplingTensor );
      material->SetTensor( couplingTensor, PIEZO_TENSOR, IMAG ); 
      std::cerr << "imaginary couplingTensor=" << std::endl << couplingTensor << std::endl;
    }

    //read nonlinearity of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "entry";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, NONLIN_COEFFICIENT ); 
      std::cerr << "entry=" << inteValue << std::endl;
    }

    //read non linear dependency of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "dependency";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      std::cerr << "dependency=" << striValue << std::endl;
    }

    //read non linear approxType of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "approxType";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      std::cerr << "approxType=" << striValue << std::endl;
    }

    //read non linear data name of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "dataName";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      std::cerr << "dataName=" << striValue << std::endl;
    }
  }

//**********************************************************************
//*************  READ MECHANICS ****************************************
//**********************************************************************
  void XMLMaterialHandler::ReadMechanic(BaseMaterial *material,
                                        const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadMechanic" );
    
    Double      doubValue;
    Integer     inteValue;
    std::string striValue;

    Matrix<Double> elasticityTensor(6,6);

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    //read material density
    keyVec = "material","mechanical","density";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, DENSITY, REAL ); 
      std::cerr << "density=" << doubValue << std::endl;
    }

    //read real elasticity tensor
    const unsigned int dim1=6, dim2=6;
    keyVec = "material","mechanical","elasticity","tensor","real";
    attrVec= "name"    ,""          ,""          ,"dim1";
    valVec =  matName  ,""          ,""          ,"6";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, elasticityTensor );
      material->SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, REAL ); 
      std::cerr << "real elasticityTensor=" << std::endl << elasticityTensor << std::endl;
    }

    //read imaginary elasticity tensor
    keyVec = "material","mechanical","elasticity","tensor","imag";
    attrVec= "name"    ,""          ,""          ,"dim1";
    valVec =  matName  ,""          ,""          ,"6";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, elasticityTensor );
      material->SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, IMAG ); 
      std::cerr << "imaginary elasticityTensor=" << std::endl << elasticityTensor << std::endl;
    }

    //read elasticity modulus
    keyVec = "material","mechanical","elasticity","isotropic","real","elasticityModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_EMODULUS, REAL ); 
      std::cerr << "elasticityModulus=" << doubValue << std::endl;
    }

    //read imaginary elasticity modulus
    keyVec = "material","mechanical","elasticity","isotropic","imag","elasticityModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_EMODULUS, IMAG ); 
      std::cerr << "imaginary elasticityModulus=" << doubValue << std::endl;
    }

    //read shear modulus
    keyVec = "material","mechanical","elasticity","isotropic","real","shearModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_GMODULUS, REAL ); 
      std::cerr << "shearModulus=" << doubValue << std::endl;
    }

    //read imaginary shear modulus
    keyVec = "material","mechanical","elasticity","isotropic","imag","shearModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_GMODULUS, IMAG ); 
      std::cerr << "imaginary shearModulus=" << doubValue << std::endl;
    }

    //read compression modulus
    keyVec = "material","mechanical","elasticity","isotropic","real","compressionModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_KMODULUS, REAL ); 
      std::cerr << "compressionModulus=" << doubValue << std::endl;
    }

    //read imaginary compression modulus
    keyVec = "material","mechanical","elasticity","isotropic","imag","compressionModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_KMODULUS, IMAG ); 
      std::cerr << "imaginary compressionModulus=" << doubValue << std::endl;
    }

    //read Lame parameter mu
    keyVec = "material","mechanical","elasticity","isotropic","real","lameParameterMu";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_LAME_MU, REAL ); 
      std::cerr << "Lame parameter mu=" << doubValue << std::endl;
    }

    //read imaginary Lame parameter mu
    keyVec = "material","mechanical","elasticity","isotropic","imag","lameParameterMu";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_LAME_MU, IMAG ); 
      std::cerr << "imaginary lameParameterMu=" << doubValue << std::endl;
    }

    //read Lame parameter lambda
    keyVec = "material","mechanical","elasticity","isotropic","real","lameParameterLamda";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_LAME_LAMBDA, REAL ); 
      std::cerr << "Lame parameter lamda=" << doubValue << std::endl;
    }

    //read imaginary Lame parameter lamda
    keyVec = "material","mechanical","elasticity","isotropic","imag","lameParameterLambda";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_LAME_LAMBDA, IMAG ); 
      std::cerr << "imaginary lameParameterLambda=" << doubValue << std::endl;
    }

    //read Poisson number
    keyVec = "material","mechanical","elasticity","isotropic","real","poissonNumber";
    attrVec= "name"     ,""         ,""          ,""         ,"";
    valVec =  matName   ,""         ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_POISSON, REAL ); 
      std::cerr << "poissonNumber=" <<  doubValue << std::endl;
    }

    //read imaginary Poisson number
    keyVec = "material","mechanical","elasticity","isotropic","imag","poissonNumber";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_POISSON, IMAG ); 
      std::cerr << "imaginary poissonNumber=" << doubValue << std::endl;
    }

    //read nonlinearity of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "entry";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, NONLIN_COEFFICIENT ); 
      std::cerr << "entry=" << inteValue << std::endl;
    }

    //read non linear dependency of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "dependency";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      std::cerr << "dependency=" << striValue << std::endl;
    }

    //read non linear approxType of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "approxType";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      std::cerr << "approxType=" << striValue << std::endl;
    }

    //read non linear data name of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "dataName";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      std::cerr << "dataName=" << striValue << std::endl;
    }

    //read alpha of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","alpha";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_ALPHA, REAL ); 
      std::cerr << "Rayleigh damping alpha=" << doubValue << std::endl;
    }

    //read beta of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","beta";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_BETA, REAL ); 
      std::cerr << "Rayleigh damping beta=" << doubValue << std::endl;
    }

    //read loss tangens delta of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","lossTangensDelta";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, LOSS_TANGENS_DELTA, REAL ); 
      std::cerr << "Rayleigh lossTangensDelta=" << doubValue << std::endl;
    }

    //read frequency of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","freq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_FREQUENCY, REAL ); 
      std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read delta frequency of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","deltaFreq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, RAYLEIGH_DELTA_FREQ, REAL ); 
      std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read algorithm of fractional damping
    keyVec = "material","mechanical","mechanicalDamping","fractional","alg";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, FRACTIONAL_ALG ); 
      std::cerr << "Fractional damping algorithm=" << striValue << std::endl;
    }

    //read memory of fractional damping
    keyVec = "material","mechanical","mechanicalDamping","fractional","memory";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, FRACTIONAL_MEMORY ); 
      std::cerr << "Fractional damping memory=" << inteValue << std::endl;
    }

    //read interpolation of fractional damping
    keyVec = "material","mechanical","mechanicalDamping","fractional","interpolation";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, FRACTIONAL_INTERPOL ); 
      std::cerr << "Fractional damping interpolation=" << striValue << std::endl;
    }
  }

//**********************************************************************
//*************  READ ACOUSTICS ****************************************
//**********************************************************************
  void XMLMaterialHandler::ReadAcoustic(BaseMaterial *material,
                                        const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadAcoustic" );

    Double      doubValue;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    //read density
    keyVec = "material","acoustic","density";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, DENSITY, REAL ); 
      std::cerr << "density=" << doubValue << std::endl;
    }

    //read compression modulus
    keyVec = "material","acoustic","compressionModulus";
    attrVec= "name"    , "";
    valVec =  matName  , "";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, ACOU_BULK_MODULUS, REAL );  ???
      std::cerr << "compressionModulus=" << doubValue << std::endl;
    }

    //read alpha of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","alpha";
    attrVec= "name"    ,""        ,""               ,"";
    valVec =  matName  ,""        ,""               ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_ALPHA, REAL ); 
      std::cerr << "Rayleigh damping alpha=" << doubValue << std::endl;
    }

    //read bata of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","beta";
    attrVec= "name"    ,""        ,""               ,"";
    valVec =  matName  ,""        ,""               ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_BETA, REAL ); 
      std::cerr << "Rayleigh damping beta=" << doubValue << std::endl;
    }

    //read loss tangens delta of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","lossTangensDelta";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, LOSS_TANGENS_DELTA, REAL ); 
      std::cerr << "Rayleigh lossTangensDelta=" << doubValue << std::endl;
    }

    //read frequency of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","freq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_FREQUENCY, REAL ); 
      std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read delta frequency of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","deltaFreq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, RAYLEIGH_DELTA_FREQ, REAL ); 
      std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read alpha0 of thermo viscous damping
    keyVec = "material","acoustic","acousticDamping","thermoViscous","alpha0";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, ACOU_ALPHA, REAL ); 
      std::cerr << "thermo viscous alpha0=" << doubValue << std::endl;
    }

    //read alpha0 of fractional damping
    keyVec = "material","acoustic","acousticDamping","fractional","alpha0";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, ACOU_ALPHA, REAL ); 
      std::cerr << "fractional alpha0=" << doubValue << std::endl;
    }

    //read exponent of fractional damping
    keyVec = "material","acoustic","acousticDamping","fractional","y";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, FRACTIONAL_EXPONENT, REAL ); 
      std::cerr << "fractional exponent=" << doubValue << std::endl;
    }

    //read acoustic non linearity
    keyVec = "material","acoustic","acousticNonlinear","bOverA";
    attrVec= "name"    ,""          ,"";
    valVec =  matName  ,""          ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, BOVERA, REAL ); 
      std::cerr << "bOverA=" << doubValue << std::endl;
    }
  }

//**********************************************************************
//*************  READ ELECTROSTATICS ************************************
//**********************************************************************
  void XMLMaterialHandler::ReadElectrostatic(BaseMaterial *material,
                                             const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadElectrostatic" );

    Double      doubValue;
    Integer     inteValue;
    std::string striValue;
    Matrix<Double> permittivityTensor(3,3);

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    //read real permittivity tensor
    const unsigned int dim1=3, dim2=3;
    keyVec = "material","electric","permittivity","tensor","real";
    attrVec= "name"    ,""        ,""            ,"dim1";
    valVec =  matName  ,""        ,""            ,"3";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, permittivityTensor );
      material->SetTensor( permittivityTensor, ELEC_PERMITTIVITY, REAL ); 
      std::cerr << "real permittivityTensor=" << std::endl << permittivityTensor << std::endl;
    }

    //read imaginary permittivity tensor
    keyVec = "material","electric","permittivity","tensor","imag";
    attrVec= "name"    ,""        ,""            ,"dim1";
    valVec =  matName  ,""        ,""            ,"3";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, permittivityTensor );
      material->SetTensor( permittivityTensor, ELEC_PERMITTIVITY, IMAG ); 
      std::cerr << "imaginary permittivityTensor=" << std::endl << permittivityTensor << std::endl;
    }

    //read isotropic permittivity
    keyVec = "material","electric","permittivity","isotropic","real";
    attrVec= "name"    ,""        ,""            ,"";
    valVec =  matName  ,""        ,""            ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, ELEC_PERMITTIVITY, REAL ); 
      std::cerr << "isotropic electric permittivity=" << doubValue << std::endl;
    }

    //read imaginary isotropic permittivity
    keyVec = "material","electric","permittivity","isotropic","imag";
    attrVec= "name"    ,""        ,""            ,"";
    valVec =  matName  ,""        ,""            ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, ELEC_PERMITTIVITY, IMAG ); 
      std::cerr << "imaginary isotropic electric permittivity=" << doubValue << std::endl;
    }

    //read nonlinearity of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient", "entry";
    attrVec= "name"     , ""      ,"nonlinear";
    valVec =  matName   , ""      ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, NONLIN_COEFFICIENT ); 
      std::cerr << "entry=" << inteValue << std::endl;
    }

    //read non linear dependency of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","dependency";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      std::cerr << "dependency=" << striValue << std::endl;
    }

    //read non linear approxType of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","approxType";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      std::cerr << "approxType=" << striValue << std::endl;
    }

    //read non linear data name of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","dataName";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      std::cerr << "dataName=" << striValue << std::endl;
    }

    //read non linear hysterese model of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","hystModel";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"hysteresis";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, HYST_MODEL ); 
      std::cerr << "hysterese model=" << striValue << std::endl;
    }

    //read e saturation of non linear hysterese model for a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","eSat";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"hysteresis";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, E_SATURATION, REAL ); 
      std::cerr << "eSat=" << doubValue << std::endl;
    }

    //read p saturation of non linear hysterese model for a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","pSat";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"hysteresis";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, P_SATURATION, REAL ); 
      std::cerr << "eSat=" << doubValue << std::endl;
    }

    //read p function of non linear hysterese model for a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","pFunction";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"hysteresis";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, P_FUNCTION ); 
      std::cerr << "pFunction=" << striValue << std::endl;
    }


  }

//**********************************************************************
//*************  READ MAGNETIC *****************************************
//**********************************************************************
  void XMLMaterialHandler::ReadMagnetic(BaseMaterial *material,
                                        const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadMagnetic" );

    Double      doubValue;
    std::string striValue;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    //read electric conductivity
    keyVec = "material","magnetic","electricConductivity";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, MAG_CONDUCTIVITY, REAL ); ???
      std::cerr << "electricConductivity=" << doubValue << std::endl;
    }

    //read magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","linear","isotropic";
    attrVec= "name"    ,""        ,""                    ,"";
    valVec =  matName  ,""        ,""                    ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MAG_PERMEABILITY, REAL ); 
      std::cerr << "magneticPermeability=" << doubValue << std::endl;
    }

    //read nonlinear dependency of magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","nonlinear","isotropic","dependency";
    attrVec= "name"    ,""        ,""                    ,""         ,"";
    valVec =  matName  ,""        ,""                    ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      std::cerr << "nonlinear dependency of magneticPermeability=" << striValue << std::endl;
    }

    //read nonlinear approxType of magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","nonlinear","isotropic","approxType";
    attrVec= "name"    ,""        ,""                    ,""         ,"";
    valVec =  matName  ,""        ,""                    ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      std::cerr << "nonlinear approxType of magneticPermeability=" << striValue << std::endl;
    }

    //read nonlinear dataName of magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","nonlinear","isotropic","dataName";
    attrVec= "name"    ,""        ,""                    ,""         ,"";
    valVec =  matName  ,""        ,""                    ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      std::cerr << "nonlinear dataName of magneticPermeability=" << striValue << std::endl;
    }
  }

//**********************************************************************
//*************  READ THERMIC ******************************************
//**********************************************************************
  void XMLMaterialHandler::ReadThermic(BaseMaterial *material,
                                       const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadThermic" );

    Double      doubValue;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    //read density
    keyVec = "material","heatConduction","density";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, DENSITY, REAL );
      std::cerr << "density=" << doubValue << std::endl;
    }

    //read heat capacity
    keyVec = "material","heatConduction","heatCapacity";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, HEAT_CAPACITY, REAL );
      std::cerr << "heatCapacity=" << doubValue << std::endl;
    }

    //read thermal conductivity
    keyVec = "material","heatConduction","thermalConductivity";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, HEAT_CONDUCTIVITY, REAL );
      std::cerr << "thermalConductivity=" << doubValue << std::endl;
    }
  }

//**********************************************************************
//*************  READ FLOW *********************************************
//**********************************************************************
  void XMLMaterialHandler::ReadFlow(BaseMaterial *material,
                                    const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadFlow" );

    Double      doubValue;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    //read density
    keyVec = "material","flow","density";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, DENSITY, REAL );
      std::cerr << "density=" << doubValue << std::endl;
    }

    //read density
    keyVec = "material","flow","dynamicViscosity";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, DYNAMIC_VISCOSITY, REAL );
      std::cerr << "dynamicViscosity=" << doubValue << std::endl;
    }

  }

}
