#include "XMLMaterialHandler.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

// header for materials
#include "Materials/electroMagneticMaterial.hh"
#include "Materials/electrostaticMaterial.hh"
#include "Materials/heatMaterial.hh"
#include "Materials/acousticMaterial.hh"
#include "Materials/mechanicMaterial.hh"
#include "Materials/piezoMaterial.hh"
#include "Materials/flowMaterial.hh"

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
  
  BaseMaterial * XMLMaterialHandler::
  LoadMaterial( const std::string matName, 
               const MaterialClass matClass ) {
    ENTER_FCN( "XMLMaterialHandler::LoadMaterial");
    
    BaseMaterial * material = NULL;
    
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    
    std::string strMatClass;

    Enum2String(matClass,strMatClass);
    
    //keyVec = "material", "electrostatic";
    keyVec = "material", strMatClass;
    attrVec= "name","";
    valVec =  matName,"";

    if (!parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      (*error) << "Error: " << strMatClass << " material " << matName 
               <<" does not exist.";
      Error( __FILE__, __LINE__ );
    }

    if ( matClass == PIEZO ) {
      material = new PiezoMaterial();
      ReadPiezo( material, matName);
    }
    else if ( matClass == MECHANIC ) {
      material = new MechanicMaterial();
      ReadMechanic( material, matName );
    }    
    else if ( matClass == FLUID ) {
      material = new AcousticMaterial();
      ReadAcoustic( material, matName );
    }
    else if ( matClass == ELECTROMAGNETIC ) {
      material = new ElectroMagneticMaterial();
      ReadMagnetic( material, matName );
    }
    else if ( matClass == ELECTROSTATIC ) {
      material = new ElectroStaticMaterial();
      ReadElectrostatic( material, matName );
    }
    else if ( matClass == THERMIC ) {
      material = new HeatMaterial();
      ReadThermic( material, matName );
    }
    else if ( matClass == FLOW ) {
      material = new FlowMaterial();
      ReadFlow( material, matName );
    }
    else {
      (*error) << "Error: material type:" << matClass << " not defined ";
      Error( __FILE__, __LINE__ );
    }
    // Finalize setup of material
    material->Finalize();

    return material;
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
      //std::cerr << "real couplingTensor=" << std::endl << couplingTensor << std::endl;
    }

    //read imaginary piezo coupling tensor
    keyVec = "material" ,"piezo" , "piezoCouplingTensor", "imag";
    attrVec= "name"     , ""     , "dim1";
    valVec =  matName   , ""     , "3";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, couplingTensor );
      material->SetTensor( couplingTensor, PIEZO_TENSOR, IMAG ); 
      //std::cerr << "imaginary couplingTensor=" << std::endl << couplingTensor << std::endl;
    }

    //read nonlinearity of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "entry";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, NONLIN_COEFFICIENT ); 
      //std::cerr << "entry=" << inteValue << std::endl;
    }

    //read non linear dependency of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "dependency";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      //std::cerr << "dependency=" << striValue << std::endl;
    }

    //read non linear approxType of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "approxType";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      //std::cerr << "approxType=" << striValue << std::endl;
    }

    //read non linear data name of a coupling coefficient
    keyVec = "material" ,"piezo" , "piezoCouplingCoefficient", "dataName";
    attrVec= "name"     , ""     , "nonlinear";
    valVec =  matName   , ""     , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      //std::cerr << "dataName=" << striValue << std::endl;
    }

    // Print material information to .info-file
    Info->PrintMaterial(material );
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

    bool     flagEModulReal=false;
    bool     flagPoissonReal=false;
    bool     flagEModulImag=false;
    bool     flagPoissonImag=false;

    bool     flagEModulXReal=false;
    bool     flagEModulYReal=false;
    bool     flagEModulZReal=false;
    bool     flagPoissonXYReal=false;
    bool     flagPoissonYZReal=false;
    bool     flagPoissonXZReal=false;
    bool     flagShearModulYZReal=false;
    bool     flagShearModulZXReal=false;
    bool     flagShearModulXYReal=false;

    bool     flagElastTensorReal=false;
    bool     flagElastTensorImag=false;

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
      //std::cerr << "density=" << doubValue << std::endl;
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
      flagElastTensorReal=true;
      //std::cerr << "real elasticityTensor=" << std::endl << elasticityTensor << std::endl;
    }

    //read imaginary elasticity tensor
    keyVec = "material","mechanical","elasticity","tensor","imag";
    attrVec= "name"    ,""          ,""          ,"dim1";
    valVec =  matName  ,""          ,""          ,"6";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, elasticityTensor );
      material->SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, IMAG ); 
      flagElastTensorImag=true;
      // std::cerr << "imaginary elasticityTensor=" << std::endl << elasticityTensor << std::endl;
    }

    //read elasticity modulus
    keyVec = "material","mechanical","elasticity","isotropic","real","elasticityModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_EMODULUS, REAL ); 
      flagEModulReal=true;
      // std::cerr << "elasticityModulus=" << doubValue << std::endl;
    }

    //read imaginary elasticity modulus
    keyVec = "material","mechanical","elasticity","isotropic","imag","elasticityModulus";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_EMODULUS, IMAG ); 
      flagEModulImag=true;
      // std::cerr << "imaginary elasticityModulus=" << doubValue << std::endl;
    }

    //read Poisson number
    keyVec = "material","mechanical","elasticity","isotropic","real","poissonNumber";
    attrVec= "name"     ,""         ,""          ,""         ,"";
    valVec =  matName   ,""         ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_POISSON, REAL ); 
      flagPoissonReal=true;
      // std::cerr << "poissonNumber=" <<  doubValue << std::endl;
    }

    //read imaginary Poisson number
    keyVec = "material","mechanical","elasticity","isotropic","imag","poissonNumber";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_POISSON, IMAG ); 
      flagPoissonImag=true;
      // std::cerr << "imaginary poissonNumber=" << doubValue << std::endl;
    }

    //read orthotropic elasticity modulus
    keyVec = "material","mechanical","elasticity","orthotropic","real","elasticityModulus_X";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_EMODULUS_X, REAL ); 
      flagEModulXReal=true;
      std::cerr << "elasticityModulus_X=" << doubValue << std::endl;
    }

    keyVec = "material","mechanical","elasticity","orthotropic","real","elasticityModulus_Y";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_EMODULUS_Y, REAL ); 
      flagEModulYReal=true;
      std::cerr << "elasticityModulus_Y=" << doubValue << std::endl;
    }

    keyVec = "material","mechanical","elasticity","orthotropic","real","elasticityModulus_Z";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_EMODULUS_Z, REAL ); 
      flagEModulZReal=true;
      std::cerr << "elasticityModulus_Z=" << doubValue << std::endl;
    }

    //read orthotropic Poisson numbers
    keyVec = "material","mechanical","elasticity","orthotropic","real","poissonNumber_XY";
    attrVec= "name"     ,""         ,""          ,""         ,"";
    valVec =  matName   ,""         ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_POISSON_XY, REAL ); 
      flagPoissonXYReal=true;
      //std::cerr << "poissonNumber_XY=" <<  doubValue << std::endl;
    }

    keyVec = "material","mechanical","elasticity","orthotropic","real","poissonNumber_YZ";
    attrVec= "name"     ,""         ,""          ,""         ,"";
    valVec =  matName   ,""         ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_POISSON_YZ, REAL ); 
      flagPoissonYZReal=true;
      //std::cerr << "poissonNumber_YZ=" <<  doubValue << std::endl;
    }

    keyVec = "material","mechanical","elasticity","orthotropic","real","poissonNumber_XZ";
    attrVec= "name"     ,""         ,""          ,""         ,"";
    valVec =  matName   ,""         ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_POISSON_XZ, REAL ); 
      flagPoissonXZReal=true;
      //std::cerr << "poissonNumber_XZ=" <<  doubValue << std::endl;
    }

    //read orthotropic shear modulus
    keyVec = "material","mechanical","elasticity","orthotropic","real","shearModulus_YZ";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_GMODULUS_YZ, REAL ); 
      flagShearModulYZReal=true;
      //std::cerr << "shearModulus_YZ=" << doubValue << std::endl;
    }

    keyVec = "material","mechanical","elasticity","orthotropic","real","shearModulus_ZX";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_GMODULUS_ZX, REAL ); 
      flagShearModulZXReal=true;
      //std::cerr << "shearModulus_ZX=" << doubValue << std::endl;
    }

    keyVec = "material","mechanical","elasticity","orthotropic","real","shearModulus_XY";
    attrVec= "name"    ,""          ,""          ,""         ,"";
    valVec =  matName  ,""          ,""          ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MECH_GMODULUS_XY, REAL ); 
      flagShearModulXYReal=true;
      //std::cerr << "shearModulus_XY=" << doubValue << std::endl;
    }


    if (flagEModulReal==true && 
        flagPoissonReal==true && 
        flagElastTensorReal==false) {
      Double EModul, PoissonNumber;
      material->SetSymmetryType(BaseMaterial::ISOTROPIC);
    }
    else if (flagEModulReal==false && 
        flagPoissonReal==false && 
        flagElastTensorReal==true) {
      //stiffness tensor is already set
    }
    else if (flagEModulXReal==true && 
             flagEModulYReal==true && 
             flagEModulZReal==true && 
             flagPoissonXYReal==true && 
             flagPoissonYZReal==true && 
             flagPoissonXZReal==true && 
             flagShearModulYZReal==true &&
             flagShearModulZXReal==true &&
             flagShearModulXYReal==true) {
      material->SetSymmetryType(BaseMaterial::ORTHOTROPIC);
    }
    else if (flagEModulReal==true && 
        flagPoissonReal==true && 
        flagElastTensorReal==true) {
      (*error) << "Error: mechanical stiffness tensor is over determined."
               << " You specified the tensor as well as E-Modul and Poisson number";
      Error( __FILE__, __LINE__ );
    }
    else if (flagEModulReal==false && 
             flagEModulXReal==false && 
             flagElastTensorReal==false) {
      (*error) << "Error: mechanical stiffness must be specified somehow.";
      Error( __FILE__, __LINE__ );
    }
    else {
      (*error) << "Error: mechanical stiffness tensor can not be computed.";
      Error( __FILE__, __LINE__ );
    }


    //*********************************************************************
    //******* end of stiffness tensor definition **************************
    //*********************************************************************

    //read nonlinearity of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "entry";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, NONLIN_COEFFICIENT ); 
      // std::cerr << "entry=" << inteValue << std::endl;
    }

    //read non linear dependency of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "dependency";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      // std::cerr << "dependency=" << striValue << std::endl;
    }

    //read non linear approxType of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "approxType";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      // std::cerr << "approxType=" << striValue << std::endl;
    }

    //read non linear data name of a elasticity coefficient
    keyVec = "material" ,"mechanical" , "elasticityCoefficient", "dataName";
    attrVec= "name"     , ""          , "nonlinear";
    valVec =  matName   , ""          , "function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      // std::cerr << "dataName=" << striValue << std::endl;
    }

    //read alpha of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","alpha";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_ALPHA, REAL ); 
      // std::cerr << "Rayleigh damping alpha=" << doubValue << std::endl;
    }

    //read beta of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","beta";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_BETA, REAL ); 
      // std::cerr << "Rayleigh damping beta=" << doubValue << std::endl;
    }

    //read loss tangens delta of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","lossTangensDelta";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, LOSS_TANGENS_DELTA, REAL ); 
      // std::cerr << "Rayleigh lossTangensDelta=" << doubValue << std::endl;
    }

    //read frequency of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","freq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_FREQUENCY, REAL ); 
      // std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read delta frequency of Rayleigh damping
    keyVec = "material","mechanical","mechanicalDamping","rayleigh","deltaFreq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, RAYLEIGH_DELTA_FREQ, REAL ); 
      // std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read algorithm of fractional damping
    keyVec = "material","mechanical","mechanicalDamping","fractional","alg";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, FRACTIONAL_ALG ); 
      // std::cerr << "Fractional damping algorithm=" << striValue << std::endl;
    }

    //read memory of fractional damping
    keyVec = "material","mechanical","mechanicalDamping","fractional","memory";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, FRACTIONAL_MEMORY ); 
      // std::cerr << "Fractional damping memory=" << inteValue << std::endl;
    }

    //read interpolation of fractional damping
    keyVec = "material","mechanical","mechanicalDamping","fractional","interpolation";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, FRACTIONAL_INTERPOL ); 
      // std::cerr << "Fractional damping interpolation=" << striValue << std::endl;
    }

    // Print information to info file
    Info->PrintMaterial( material);
  }

//   void XMLMaterialHandler::ComputeIsoMechStiffnesTensor(Double EModul, 
//                                                         Double PoissonNumber,
//                                                         Matrix<Double>& elasticityTensor){
//     Double LameLambda, LameMu;
//     elasticityTensor.Resize(6,6);
//     elasticityTensor.Init();
//     LameLambda = (PoissonNumber*EModul)/((1.0 + PoissonNumber)*(1.0 - 2.0*PoissonNumber));
//     LameMu    = (EModul)/(2.0*(1.0+PoissonNumber));
//     // std::cerr << "LameLambda=" << LameLambda << "\nLameMu=" << LameMu << std::endl;

//     elasticityTensor[0][0]=LameLambda+2.0*LameMu;
//     elasticityTensor[1][1]=LameLambda+2.0*LameMu;
//     elasticityTensor[2][2]=LameLambda+2.0*LameMu;

//     elasticityTensor[0][1]=LameLambda;
//     elasticityTensor[0][2]=LameLambda;
//     elasticityTensor[1][0]=LameLambda;
//     elasticityTensor[1][2]=LameLambda;
//     elasticityTensor[2][0]=LameLambda;
//     elasticityTensor[2][1]=LameLambda;

//     elasticityTensor[3][3]=LameMu;
//     elasticityTensor[4][4]=LameMu;
//     elasticityTensor[5][5]=LameMu;
//   }

//   void XMLMaterialHandler::ComputeOrthoMechStiffnesTensor(Double EX, Double EY, Double EZ, 
//                                                           Double nuXY, Double nuYZ, Double nuXZ,
//                                                           Double GYZ, Double GZX, Double GXY,
//                                                           Matrix<Double>& elasticityTensor){
//     Double nuYX, nuZY, nuZX, aux;
//     nuYX=(EY/EX)*nuXY;
//     nuZY=(EZ/EY)*nuYZ;
//     nuZX=(EZ/EX)*nuXZ;

//     aux=(1-nuXY*nuYX-nuYZ*nuZY-nuXZ*nuZX-2.0*nuYX*nuZY*nuXZ)/(EX*EY*EZ);

//     elasticityTensor.Resize(6,6);
//     elasticityTensor.Init();

//     elasticityTensor[0][0]=(1-nuYZ*nuZY)/(EY*EZ*aux);
//     elasticityTensor[1][1]=(1-nuXZ*nuZX)/(EX*EZ*aux);
//     elasticityTensor[2][2]=(1-nuXY*nuYX)/(EX*EY*aux);

//     elasticityTensor[0][1]=(nuYX+nuZX*nuYZ)/(EY*EZ*aux);
//     elasticityTensor[0][2]=(nuZX+nuYX*nuZY)/(EY*EZ*aux);
//     elasticityTensor[1][0]=(nuYX+nuZX*nuYZ)/(EY*EZ*aux);
//     elasticityTensor[1][2]=(nuZY+nuXY*nuZX)/(EX*EZ*aux);
//     elasticityTensor[2][0]=(nuZX+nuYX*nuZY)/(EY*EZ*aux);
//     elasticityTensor[2][1]=(nuZY+nuXY*nuZX)/(EX*EZ*aux);

//     elasticityTensor[3][3]=GYZ;
//     elasticityTensor[4][4]=GZX;
//     elasticityTensor[5][5]=GXY;
//   }

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
      // std::cerr << "density=" << doubValue << std::endl;
    }

    //read compression modulus
    keyVec = "material","acoustic","compressionModulus";
    attrVec= "name"    , "";
    valVec =  matName  , "";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, ACOU_BULK_MODULUS, REAL );
      // std::cerr << "compressionModulus=" << doubValue << std::endl;
    }

    //read alpha of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","alpha";
    attrVec= "name"    ,""        ,""               ,"";
    valVec =  matName  ,""        ,""               ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_ALPHA, REAL ); 
      // std::cerr << "Rayleigh damping alpha=" << doubValue << std::endl;
    }

    //read bata of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","beta";
    attrVec= "name"    ,""        ,""               ,"";
    valVec =  matName  ,""        ,""               ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_BETA, REAL ); 
      // std::cerr << "Rayleigh damping beta=" << doubValue << std::endl;
    }

    //read loss tangens delta of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","lossTangensDelta";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, LOSS_TANGENS_DELTA, REAL ); 
      // std::cerr << "Rayleigh lossTangensDelta=" << doubValue << std::endl;
    }

    //read frequency of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","freq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, RAYLEIGH_FREQUENCY, REAL ); 
      // std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read delta frequency of Rayleigh damping
    keyVec = "material","acoustic","acousticDamping","rayleigh","deltaFreq";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, RAYLEIGH_DELTA_FREQ, REAL ); 
      // std::cerr << "Rayleigh damping freq=" << doubValue << std::endl;
    }

    //read alpha0 of thermo viscous damping
    keyVec = "material","acoustic","acousticDamping","thermoViscous","alpha0";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, ACOU_ALPHA, REAL ); 
      // std::cerr << "thermo viscous alpha0=" << doubValue << std::endl;
    }

    //read alpha0 of fractional damping
    keyVec = "material","acoustic","acousticDamping","fractional","alpha0";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, ACOU_ALPHA, REAL ); 
      // std::cerr << "fractional alpha0=" << doubValue << std::endl;
    }

    //read exponent of fractional damping
    keyVec = "material","acoustic","acousticDamping","fractional","y";
    attrVec= "name"    ,""          ,""                 ,"";
    valVec =  matName  ,""          ,""                 ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, FRACTIONAL_EXPONENT, REAL ); 
      // std::cerr << "fractional exponent=" << doubValue << std::endl;
    }

    //read acoustic non linearity
    keyVec = "material","acoustic","acousticNonlinear","bOverA";
    attrVec= "name"    ,""          ,"";
    valVec =  matName  ,""          ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, BOVERA, REAL ); 
      // std::cerr << "bOverA=" << doubValue << std::endl;
    }

    // Print material information to info-file
    Info->PrintMaterial( material );
  }

//**********************************************************************
//*************  READ ELECTROSTATICS ************************************
//**********************************************************************
  void XMLMaterialHandler::ReadElectrostatic(BaseMaterial *material,
                                             const std::string matName) {
    ENTER_FCN( "XMLMaterialHandler::ReadElectrostatic" );

    Double      doubValue;
    Integer     inteValue, dim;
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
      // std::cerr << "real permittivityTensor=" << std::endl << permittivityTensor << std::endl;
    }

    //read imaginary permittivity tensor
    keyVec = "material","electric","permittivity","tensor","imag";
    attrVec= "name"    ,""        ,""            ,"dim1";
    valVec =  matName  ,""        ,""            ,"3";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim1, dim2, permittivityTensor );
      material->SetTensor( permittivityTensor, ELEC_PERMITTIVITY, IMAG ); 
      // std::cerr << "imaginary permittivityTensor=" << std::endl << permittivityTensor << std::endl;
    }

    //read isotropic permittivity
    keyVec = "material","electric","permittivity","isotropic","real";
    attrVec= "name"    ,""        ,""            ,"";
    valVec =  matName  ,""        ,""            ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, ELEC_PERMITTIVITY, REAL ); 
      // std::cerr << "isotropic electric permittivity=" << doubValue << std::endl;
    }

    //read imaginary isotropic permittivity
    keyVec = "material","electric","permittivity","isotropic","imag";
    attrVec= "name"    ,""        ,""            ,"";
    valVec =  matName  ,""        ,""            ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, ELEC_PERMITTIVITY, IMAG ); 
      // std::cerr << "imaginary isotropic electric permittivity=" << doubValue << std::endl;
    }

    //read nonlinearity of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient", "entry";
    attrVec= "name"     , ""      ,"nonlinear";
    valVec =  matName   , ""      ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      material->SetScalar( inteValue, NONLIN_COEFFICIENT ); 
      // std::cerr << "entry=" << inteValue << std::endl;
    }

    //read non linear dependency of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","dependency";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      // std::cerr << "dependency=" << striValue << std::endl;
    }

    //read non linear approxType of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","approxType";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      // std::cerr << "approxType=" << striValue << std::endl;
    }

    //read non linear data name of a permittivity coefficient
    keyVec = "material","electric","permittivityCoefficient","dataName";
    attrVec= "name"    ,""        ,"nonlinear";
    valVec =  matName  ,""        ,"function";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      // std::cerr << "dataName=" << striValue << std::endl;
    }

    //read Preisach hysterese model
    keyVec = "material","electric","hystModel","preisach";
    attrVec= "name"    ,""        ,"";
    valVec =  matName  ,""        ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      striValue="preisach";
      material->SetScalar(striValue, HYST_MODEL ); 
    }

    //read E saturation of Preisach hysterese model
    keyVec = "material","electric","hystModel","preisach","eSat";
    attrVec= "name"    ,""        ,""         ,"";
    valVec =  matName  ,""        ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, E_SATURATION, REAL ); 
    }

    //read P saturation of Preisach hysterese model
    keyVec = "material","electric","hystModel","preisach","pSat";
    attrVec= "name"    ,""        ,""         ,"";
    valVec =  matName  ,""        ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, P_SATURATION, REAL ); 
    }

    //read weight dimension of Preisach hysterese model
    keyVec = "material","electric","hystModel","preisach","dim";
    attrVec= "name"    ,""        ,""         ,"";
    valVec =  matName  ,""        ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, inteValue );
      dim=inteValue;
    }

    Matrix<Double> preisachWeightTensor(dim,dim);

    //read real permittivity tensor
    keyVec = "material","electric","hystModel","preisach","weights";
    attrVec= "name"    ,""        ,""            ,"";
    valVec =  matName  ,""        ,""            ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
                                   dim, dim, preisachWeightTensor );
      material->SetTensor( preisachWeightTensor, PREISACH_WEIGHTS, REAL); 
    }

    // Print information to info file
    Info->PrintMaterial( material );

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
      material->SetScalar( doubValue, MAG_CONDUCTIVITY, REAL );
      std::cerr << matName << "electricConductivity=" << doubValue << std::endl;
    }

    //read magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","linear","isotropic";
    attrVec= "name"    ,""        ,""                    ,"";
    valVec =  matName  ,""        ,""                    ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, MAG_PERMEABILITY, REAL ); 
      std::cerr << matName << "magneticPermeability=" << doubValue << std::endl;
    }

    //read nonlinear dependency of magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","nonlinear","isotropic","dependency";
    attrVec= "name"    ,""        ,""                    ,""         ,"";
    valVec =  matName  ,""        ,""                    ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      //material->SetScalar( striValue, NONLIN_DEPENDENCY ); 
      // std::cerr << "nonlinear dependency of magneticPermeability=" << striValue << std::endl;
    }

    //read nonlinear approxType of magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","nonlinear","isotropic","approxType";
    attrVec= "name"    ,""        ,""                    ,""         ,"";
    valVec =  matName  ,""        ,""                    ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      //material->SetScalar( striValue, NONLIN_APPROXIMATION_TYPE ); 
      // std::cerr << "nonlinear approxType of magneticPermeability=" << striValue << std::endl;
    }

    //read nonlinear dataName of magnetic permeability
    keyVec = "material","magnetic","magneticPermeability","nonlinear","isotropic","dataName";
    attrVec= "name"    ,""        ,""                    ,""         ,"";
    valVec =  matName  ,""        ,""                    ,""         ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, striValue );
      //material->SetScalar( striValue, NONLIN_DATA_NAME ); 
      material->SetNonlinFileName( striValue.c_str() );
      // std::cerr << "nonlinear dataName of magneticPermeability=" << striValue << std::endl;
    }

    // Print information to info file
    Info->PrintMaterial( material ); 
    
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
      // std::cerr << "density=" << doubValue << std::endl;
    }

    //read heat capacity
    keyVec = "material","heatConduction","heatCapacity";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, HEAT_CAPACITY, REAL );
      // std::cerr << "heatCapacity=" << doubValue << std::endl;
    }

    //read thermal conductivity
    keyVec = "material","heatConduction","thermalConductivity";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      material->SetScalar( doubValue, HEAT_CONDUCTIVITY, REAL );
      // std::cerr << "thermalConductivity=" << doubValue << std::endl;
    }

    // Print information to info file
    Info->PrintMaterial( material );
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
      // std::cerr << "density=" << doubValue << std::endl;
    }

    //read density
    keyVec = "material","flow","dynamicViscosity";
    attrVec= "name"    ,"";
    valVec =  matName  ,"";
    if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
      parser_->Get( keyVec, attrVec, valVec, doubValue );
      //material->SetScalar( doubValue, DYNAMIC_VISCOSITY, REAL );
      // std::cerr << "dynamicViscosity=" << doubValue << std::endl;
    }

    // Print information to info file
    Info->PrintMaterial( material );
  }

}
