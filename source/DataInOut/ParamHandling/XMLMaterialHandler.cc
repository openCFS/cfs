#include "XMLMaterialHandler.hh"

#include "Domain/CoefFunction/CoefFunction.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ProgramOptions.hh"

// header for materials
#include "Materials/ElectroMagneticMaterial.hh"
#include "Materials/ElectroStaticMaterial.hh"
#include "Materials/HeatMaterial.hh"
#include "Materials/AcousticMaterial.hh"
#include "Materials/MechanicMaterial.hh"
#include "Materials/PiezoMaterial.hh"
#include "Materials/FlowMaterial.hh"
#include "Materials/TestMaterial.hh"
#include "Materials/ElectricConductionMaterial.hh"
//#include "Materials/thermoelasticMaterial.hh"
//#include "Materials/pyroelectricMaterial.hh"
#include "Materials/magStrictMaterial.hh"
#include "Utils/tools.hh"

// Note, that the methods ComputeIso/OrthoMechStiffnesTensor were commented out
// in revision 7562 and are not in the code -> check the repository!

namespace CoupledField {
  
  XMLMaterialHandler::XMLMaterialHandler()
  : MaterialHandler() {
  }
  
  XMLMaterialHandler::~XMLMaterialHandler()
  {
  }
  
  void XMLMaterialHandler::LoadFromFile( const std::string& fileName ) {
    
    this->fileName_ = fileName;
    
    // Create a ParamNode and parse the material file
    std::string schema = progOpts->GetSchemaPathStr();
    schema += "/CFS-Material/CFS_Material.xsd";
    
    parser_ = XmlReader::ParseFile(fileName, schema);
  }
  
  void XMLMaterialHandler::LoadFromString( const std::string& str ) {
    // Create a ParamNode and parse the material file
    std::string schema = progOpts->GetSchemaPathStr();
    schema += "/CFS-Material/CFS_Material.xsd";
    
    parser_ = XmlReader::ParseString(str, schema);
  }
  
  BaseMaterial * XMLMaterialHandler::LoadMaterial( const std::string matName,
          const MaterialClass matClass ) {
    
    BaseMaterial * material = NULL;
    
    std::string strMatClass;
    
    Enum2String(matClass,strMatClass);
    if(!parser_->HasByVal("material", "name", matName))
      EXCEPTION("Cannot find material '" << matName << "'");
    
    // first get the material element:  <material name="iron">
    PtrParamNode pn;
    pn = parser_->GetByVal("material", "name", matName);
    
    if( !pn ) {
      EXCEPTION( "Material with name '" << matName 
              << "' could not be found in material file!" );
    }
    // the the requested material class: <mechanical>
    pn = pn->Get(strMatClass);  
    MathParser * mp = domain_->GetMathParser();
    CoordSystem * cs = domain_->GetCoordSystem();
    try {
      if ( matClass == PIEZO ) {
        material = new PiezoMaterial(mp, cs);
        ReadPiezo( material, pn);
      }
      else if ( matClass == MECHANIC ) {
        material = new MechanicMaterial(mp, cs);
        ReadMechanic( dynamic_cast<MechanicMaterial *>(material), pn );
      }    
      else if ( matClass == FLUID ) {\
        material = new AcousticMaterial(mp, cs);
      ReadAcoustic( material, pn );
      }
      else if ( matClass == ELECTROMAGNETIC ) {
        material = new ElectroMagneticMaterial(mp, cs);
        ReadMagnetic( material, pn );
      }
      else if ( matClass == ELECTROSTATIC ) {
        material = new ElectroStaticMaterial(mp, cs);
        ReadElectrostatic( material, pn );
      }
      else if ( matClass == THERMIC ) {
        material = new HeatMaterial(mp, cs);
        ReadThermic( material, pn );
      }
      else if ( matClass == FLOW ) {
        material = new FlowMaterial(mp, cs);
        ReadFlow( material, pn );
      }
      else if ( matClass == TESTMAT ) {
        material = new TestMaterial(mp, cs);
        ReadTest( material, pn );
      }
      else if ( matClass == PYROELECTRIC ) {
        REFACTOR;
        //material = new PyroelectricMaterial();
        //ReadPyroelectric( material,pn );
      }
      else if ( matClass == THERMOELASTIC ) {
        REFACTOR;//
        //material = new ThermoelasticMaterial();
        //ReadThermoelastic( material, pn );
      }
      else if ( matClass == MAGNETOSTRICTIVE ) {
        //REFACTOR;
        material = new MagStrictMaterial(mp,cs);
        ReadMagStrict( material, pn );
      }
      else if ( matClass == ELECTRICCONDUCTION ) {
        material = new ElectricConductionMaterial(mp, cs);
        ReadElectricConduction( material, pn );
      }
      else {
        EXCEPTION( "material type:" << matClass << " not defined" );
      }
      
      // Finalize setup of material
      material->Finalize();
    } catch (Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not load material '" << matName  
              << "' of class '" << matClass << "'" );
    }
    return material;
  }
  
  //**********************************************************************
  //*************  READ PIEZO ********************************************
  //**********************************************************************
  
  void XMLMaterialHandler::ReadPiezo(BaseMaterial *material, PtrParamNode pn) 
  {
    //    std::cout << "ReadPiezo" << std::endl;
    //read real piezo coupling tensor
    if(pn->Has("piezoCouplingTensor"))
    {
      StdVector<std::string> realVals(18), imagVals(18);
      realVals.Init("0.0");
      imagVals.Init("0.0");
      PtrCoefFct pctCoef;
      PtrParamNode pct = pn->Get("piezoCouplingTensor");
      if(pct->Has("real")){
        ParamTools::AsStringTensor( pct->Get("real"), 18, realVals );
      }
      if(pct->Has("imag"))
      {
        ParamTools::AsStringTensor( pct->Get("imag"), 18, imagVals );
      }
      pctCoef = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 6,
              realVals, imagVals );
      material->SetCoefFct( PIEZO_TENSOR, pctCoef);
    } 
    
    //read nonlinearity of a coupling coefficient
    if(pn->HasByVal("piezoCouplingCoefficient", "nonlinear", "function"))
    {
      PtrParamNode pcc = pn->GetByVal("piezoCouplingCoefficient", 
              std::string("nonlinear"), 
              "function");
      if(pcc->Has("entry"))
        material->SetScalar(pcc->Get("entry")->As<Integer>(), NONLIN_COEFFICIENT);
      
      if(pcc->Has("dependency"))
        material->SetScalar(pcc->Get("dependency")->As<std::string>(), NONLIN_DEPENDENCY);
      
      if(pcc->Has("approxType"))
        material->SetScalar(pcc->Get("approxType")->As<std::string>(), NONLIN_APPROXIMATION_TYPE);
      
      if(pcc->Has("dataName"))
        material->SetScalar(pcc->Get("dataName")->As<std::string>(), NONLIN_DATA_NAME);
    }
    
    if ( pn->Has("piezoMicroData"))
    {
      if( pn->Get("piezoMicroData")->Has("HuberFleck"))
      {
        PtrParamNode pcc = pn->Get("piezoMicroData")->Get("HuberFleck");
        
        // force name
        //        material->SetScalar("BelovKreher", PIEZO_MICRO_MODEL);
        
        // read remanent polarisation
        if(pcc->Has("sponPolarization"))
          material->SetScalar(pcc->Get("sponPolarization")->As<Double>(), SPON_POLARIZATION, Global::REAL ); 
        
        // read remanent strain
        if(pcc->Has("sponStrain")) 
          material->SetScalar(pcc->Get("sponStrain")->As<Double>(), SPON_STRAIN, Global::REAL ); 
        
        // 
        if(pcc->Has("Efield0")) 
          material->SetScalar(pcc->Get("Efield0")->As<Double>(), EFIELD0, Global::REAL ); 
        
        // 
        if(pcc->Has("Stress0")) 
          material->SetScalar(pcc->Get("Stress0")->As<Double>(), STRESS0, Global::REAL ); 
        
        // 
        if(pcc->Has("dCouple0")) 
          material->SetScalar(pcc->Get("dCouple0")->As<Double>(), DCOUPLE0, Global::REAL ); 
        
        // read rate constant
        if(pcc->Has("rateConstant")) 
          material->SetScalar(pcc->Get("rateConstant")->As<Double>(), RATE_CONSTANT, Global::REAL ); 
        
        // read visco-plasti index
        if(pcc->Has("viscoPlasticIndex")) 
          material->SetScalar(pcc->Get("viscoPlasticIndex")->As<Double>(), VISCO_PLASTIC_INDEX, Global::REAL ); 
        
        // read saturation index
        if(pcc->Has("saturationIndex")) 
          material->SetScalar(pcc->Get("saturationIndex")->As<Double>(), SATURATION_INDEX, Global::REAL ); 
        
        // read init value for volume fraction
        if(pcc->Has("volumeFracInit")) 
          material->SetScalar(pcc->Get("volumeFracInit")->As<Double>(), VOLUME_FRAC_INIT, Global::REAL ); 
        
        // 
        if(pcc->Has("scaleForceElec")) 
          material->SetScalar(pcc->Get("scaleForceElec")->As<Double>(), SCALE_FORCE_ELEC, Global::REAL ); 
        
        // 
        if(pcc->Has("scaleForceMech")) 
          material->SetScalar(pcc->Get("scaleForceMech")->As<Double>(), SCALE_FORCE_MECH, Global::REAL ); 
        
        // 
        if(pcc->Has("scaleForceCouple")) 
          material->SetScalar(pcc->Get("scaleForceCouple")->As<Double>(), SCALE_FORCE_COUPLE, Global::REAL ); 
        
        // read mean temperatute
        if(pcc->Has("Tmean"))
          material->SetScalar(pcc->Get("Tmean")->As<Double>(), MEAN_TEMPERATURE, Global::REAL ); 
        
      }
    }
    
  }
  
  //**********************************************************************
  //*************  READ MECHANICS ****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadMechanic(BaseMaterial *material, PtrParamNode mech)
  {
    //    std::cout << "ReadMechanics" << std::endl;
    bool     flagElasticitySet=false;
    
    // bool     flagElastTensorImag=false;
    
    
    //! [Read PtrParamNode]
    //read material density
    if(mech->Has("density")) {
      PtrCoefFct densFct =
              CoefFunction::Generate(mp_, Global::REAL, 
              mech->Get("density")->As<std::string>() );
      material->SetCoefFct( DENSITY, densFct );
    }
    //! [Read PtrParamNode]
    
    // quite a lot is elasitcity
    if(mech->Has("elasticity"))
    {
      PtrParamNode elast = mech->Get("elasticity");
      
      if(elast->HasByVal("tensor", "dim1", "6"))
      {
        PtrParamNode tens = elast->GetByVal("tensor", std::string("dim1"), "6");
        StdVector<std::string> realVals(36), imagVals(36);
        realVals.Init("0.0");
        imagVals.Init("0.0");
        PtrCoefFct elastCoef;
        //read real elasticity tensor   
        if(tens->Has("real")) {
          ParamTools::AsStringTensor( tens->Get("real"), 36, realVals );
        }
        if(tens->Has("imag")) {
          ParamTools::AsStringTensor( tens->Get("imag"), 36, imagVals );
        }
        elastCoef = CoefFunction::Generate(mp_,  Global::COMPLEX, 6, 6,
                realVals, imagVals );
        material->SetCoefFct( MECH_STIFFNESS_TENSOR, elastCoef);
        material->SetSymmetryType(MECH_STIFFNESS_TENSOR,BaseMaterial::GENERAL);
        flagElasticitySet = true;
      } // end tensor  
      
      // check values for isotropic
      if(elast->Has("isotropic"))
      {
        PtrParamNode iso = elast->Get("isotropic");
        StdVector<std::string> vRe(36), vIm(36);
        vRe.Init("0");
        vIm.Init("0");
        if ( iso->Has("real") ) // read the real part
        {
          vRe = ReadMechanicIsotropic(iso->Get("real"));
        }
        if ( iso->Has("imag") ) // read the imaginary part
        {
          vIm = ReadMechanicIsotropic(iso->Get("imag"));
        }
        PtrCoefFct C = CoefFunction::Generate(mp_, Global::IMAG,6,6,vRe,vIm);
        material->SetCoefFct(MECH_STIFFNESS_TENSOR,C);
        material->SetSymmetryType(MECH_STIFFNESS_TENSOR,BaseMaterial::ISOTROPIC);
        flagElasticitySet = true;
      } // end of isotropic
      
      // Note, in revision 7565 there were some details in the code now deleted 
      if(elast->Has("transversalIsotropic"))
        throw Exception("transversalIsotropic for elasticity in mechanical is not implemented!");
      
      // check orthotropic stuff
      if(elast->Has("orthotropic"))
      {
        // onyl real is implemented
        if(elast->Get("orthotropic")->Has("imag"))
          throw Exception("imaginary orthotropic elasitcity parameters for mechanical not implemented");
        
        PtrParamNode real = elast->Get("orthotropic")->Get("real");
        
        //read orthotropic elasticity modulus
        if(real->Has("elasticityModulus_1"))
        {
          material->SetScalar(real->Get("elasticityModulus_1")->As<Double>(), MECH_EMODULUS_X, Global::REAL ); 
        }
        
        if(real->Has("elasticityModulus_2"))
        {
          material->SetScalar(real->Get("elasticityModulus_2")->As<Double>(), MECH_EMODULUS_Y, Global::REAL ); 
        }
        
        if(real->Has("elasticityModulus_3"))
        {
          material->SetScalar(real->Get("elasticityModulus_3")->As<Double>(), MECH_EMODULUS_Z, Global::REAL ); 
        }
        
        // read orthotropic Poisson numbers
        if(real->Has("poissonNumber_12"))
        {
          material->SetScalar(real->Get("poissonNumber_12")->As<Double>(), MECH_POISSON_XY, Global::REAL ); 
        }
        
        if(real->Has("poissonNumber_23"))
        {
          material->SetScalar(real->Get("poissonNumber_23")->As<Double>(), MECH_POISSON_YZ, Global::REAL ); 
        }
        
        if(real->Has("poissonNumber_13"))
        {
          material->SetScalar(real->Get("poissonNumber_13")->As<Double>(), MECH_POISSON_XZ, Global::REAL ); 
        }
        
        // read orthotropic shear modulus
        if(real->Has("shearModulus_23"))
        {
          material->SetScalar(real->Get("shearModulus_23")->As<Double>(), MECH_GMODULUS_YZ, Global::REAL ); 
        }
        
        if(real->Has("shearModulus_31"))
        {
          material->SetScalar(real->Get("shearModulus_31")->As<Double>(), MECH_GMODULUS_ZX, Global::REAL ); 
        }
        
        if(real->Has("shearModulus_12"))
        {
          material->SetScalar(real->Get("shearModulus_12")->As<Double>(), MECH_GMODULUS_XY, Global::REAL ); 
        }
      }  // orthotropic      
    } // end of elasticity
    
    
    if (flagElasticitySet==false) {
      EXCEPTION ("mechanical stiffness must be specified somehow.");
    }
    
    //read viscoelastic parameters
    if(mech->Has("viscoelasticity")) {
    	PtrParamNode visco = mech->Get("viscoelasticity");
    	if(visco->Has("isotropic")) {
    		PtrParamNode tens = visco->Get("isotropic");
    		UInt dim2;
    		if ( tens->Has("dim2") )
    			dim2 = tens->Get("dim2")->As<Integer>();
    		else
    			EXCEPTION("Viscoelasticity parameters: dim2 has to be specified");
        
    		Matrix<Double> realMat(3,dim2);
        realMat.Init();
        PtrCoefFct viscoelastCoef;
        //read real viscoelastic coefficients
        if ( tens->Has("real") ) {
          ParamTools::AsTensor( tens->Get("real"), 3, dim2, realMat );
        }
        Vector<Double> dummy(dim2);
        for (UInt i=0; i<dim2; i++)
          dummy[i] = realMat[0][i];
        material->SetVector( dummy, MECH_VISCOALPHA_VECTOR, Global::REAL);
        for (UInt i=0; i<dim2; i++)
          dummy[i] = realMat[1][i];
        material->SetVector( dummy, MECH_VISCOK_VECTOR, Global::REAL);
        for (UInt i=0; i<dim2; i++)
          dummy[i] = realMat[2][i];
        material->SetVector( dummy, MECH_VISCOG_VECTOR, Global::REAL);
    	}
    }
    
    // elasticityCoefficient of type <elasticityCoefficient nonlinear="function">
    if(mech->HasByVal("elasticityCoefficient", "nonlinear", "function"))
    {
      PtrParamNode ec = mech->GetByVal("elasticityCoefficient", 
              std::string("nonlinear"),
              "function");
      if(ec->Has("entry")) 
        material->SetScalar(ec->Get("entry")->As<Integer>(), NONLIN_COEFFICIENT);
      
      if(ec->Has("dependency"))
        material->SetScalar(ec->Get("dependency")->As<std::string>(), NONLIN_DEPENDENCY );
      
      if(ec->Has("approxType"))
        material->SetScalar(ec->Get("approxType")->As<std::string>(), NONLIN_APPROXIMATION_TYPE );
      
      if(ec->Has("dataName"))
        material->SetScalar(ec->Get("dataName")->As<std::string>(), NONLIN_DATA_NAME );
    }; // end of <elasticityCoefficient nonlinear="function">
    
    
    //read coefficients for irreversible mechanical strain
    if(mech->Has("irreversibleStrainCoefficient"))
    {
      PtrParamNode isc = mech->Get("irreversibleStrainCoefficient");
      // the dimension is only printed in the old param handler version 7562
      //if(isc->Has("dim")) std::cout << "dim=" << isc->Get("dim")->As<Integer>() << std::endl;
      
      if(isc->Has("coeffs"))
      {
        // read matrix
        Matrix<Double> matrixCoeffs(5,1);
        ParamTools::AsTensor<double>(isc->Get("coeffs"),1,5,matrixCoeffs); 
        
        // transform to vector
        Vector<Double> coeffs;
        coeffs.Resize( matrixCoeffs.GetNumCols());
        for( UInt i=0; i<matrixCoeffs.GetNumCols(); i++)
          coeffs[i] = matrixCoeffs[0][i];
        
        material->SetVector( coeffs, COEFF_STRAIN_IRREVERSIBLE, Global::REAL ); 
      }
    } // end of irreversibleStrainCoefficient
    
    
    // check and read thermal expansion coefficients (TECs)
    if (mech->Has("thermalExpanison")) {
      PtrParamNode tec = mech->Get("thermalExpanison");
      
      // read reference temperature
      PtrCoefFct tRef;// = CoefFunction::Generate( mp_, Global::IMAG,tecR,tecI);
      if (tec->Has("refTemperature")) {
        PtrParamNode refT = tec->Get("refTemperature");
        tRef = ParamTools::AsScalarCoefFct(mp_,refT);
      }
      else { // set ref temperature to zero
        tRef = CoefFunction::Generate( mp_, Global::REAL, "0.0");
      }
      material->SetCoefFct(MECH_TE_REFTEMPERATURE,tRef);
      
      // read thermal expansion coefficient to create a vector valued coef function (in Voigt notation)
      StdVector<std::string> tecR(6), tecI(6);
      tecR.Init("0.0");
      tecI.Init("0.0");
      if (tec->Has("isotropic")) {
        if (tec->Get("isotropic")->Has("real")) {
          std::string coef = tec->Get("isotropic")->Get("real")->As<std::string>();
          tecR[0] = coef;
          tecR[1] = coef;
          tecR[2] = coef;
        }
        if (tec->Get("isotropic")->Has("imag")) {
          std::string coef = tec->Get("isotropic")->Get("imag")->As<std::string>();
          tecI[0] = coef;
          tecI[1] = coef;
          tecI[2] = coef;
        }
      }
      else if (tec->Has("orthotropic")) {
        PtrParamNode node = tec->Get("orthotropic");
        StdVector<std::string> vals(3);
        vals.Init("0.0");
        if (node->Has("real")) {
          ParamTools::AsStringTensor( node->Get("real"), 3, vals );
          for (UInt i = 0; i < 3; ++i) {
            tecR[i] = vals[i];
          }
        }
        if (node->Has("imag")) {
          ParamTools::AsStringTensor( node->Get("imag"), 3, vals );
          for (UInt i = 0; i < 3; ++i) {
            tecI[i] = vals[i];
          }
        }
      }
      else if(tec->Has("anisotropic")) {
        PtrParamNode node = tec->Get("anisotropic");
        StdVector<std::string> vals(6);
        vals.Init("0.0");
        if (node->Has("real")) {
          ParamTools::AsStringTensor( node->Get("real"), 6, vals );
          for (UInt i = 0; i < 6; ++i) {
            tecR[i] = vals[i];
          }
        }
        if (node->Has("imag")) {
          ParamTools::AsStringTensor( node->Get("imag"), 6, vals );
          for (UInt i = 0; i < 6; ++i) {
            tecI[i] = vals[i];
          }
        }
      }
      PtrCoefFct tecVect = CoefFunction::Generate( mp_, Global::COMPLEX,tecR,tecI);
      material->SetCoefFct(MECH_TE_TENSOR,tecVect);
    }
    
    // read mechanical damping
    if(mech->Has("mechanicalDamping"))
    {
      // first rayleigh damping
      if(mech->Get("mechanicalDamping")->Has("rayleigh"))
      {
        PtrParamNode r = mech->Get("mechanicalDamping")->Get("rayleigh");
        
        if(r->Has("alpha"))
          material->SetScalar(r->Get("alpha")->As<std::string>(), RAYLEIGH_ALPHA);
        
        if(r->Has("beta"))
          material->SetScalar(r->Get("beta")->As<std::string>(), RAYLEIGH_BETA);
        
        if(r->Has("lossTangensDelta"))
          material->SetScalar(r->Get("lossTangensDelta")->As<std::string>(), LOSS_TANGENS_DELTA);
        
        if(r->Has("measuredFreq"))
          material->SetScalar(r->Get("measuredFreq")->As<Double>(), RAYLEIGH_FREQUENCY, Global::REAL);
      }
      if(mech->Get("mechanicalDamping")->Has("fractional"))
      {
        PtrParamNode f = mech->Get("mechanicalDamping")->Get("fractional");
        
        if(f->Has("alg"))        
          material->SetScalar(f->Get("alg")->As<std::string>(), FRACTIONAL_ALG );
        
        if(f->Has("memory"))        
          material->SetScalar(f->Get("memory")->As<Integer>(), FRACTIONAL_MEMORY );
        
        if(f->Has("interpolation"))        
          material->SetScalar(f->Get("interpolation")->As<std::string>(), FRACTIONAL_INTERPOL );
      }
    }
    
    //read real magmech coupling tensor
    if(mech->Has("magnetoStrictionTensor_h_mech"))
    {
      StdVector<std::string> realVals(18), imagVals(18);
      realVals.Init("0.0");
      imagVals.Init("0.0");
      PtrCoefFct pctCoef;
      PtrParamNode pct = mech->Get("magnetoStrictionTensor_h_mech");
      if(pct->Has("real")){
        ParamTools::AsStringTensor( pct->Get("real"), 18, realVals );
      }
      if(pct->Has("imag"))
      {
        ParamTools::AsStringTensor( pct->Get("imag"), 18, imagVals );
      }
      pctCoef = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 6,
              realVals, imagVals );
      material->SetCoefFct( MAGNETOSTRICTION_TENSOR_h_mech, pctCoef);
    }
    /* old style
     if(mech->Has("magnetoStrictionTensor_h_mech")) {
     Matrix<Double> couplingTensor(3,6);
     
     PtrParamNode mst = mech->Get("magnetoStrictionTensor_h_mech");
     if(mst->Has("real")) {
     ParamTools::AsTensor<double>(mst->Get("real"), 3, 6, couplingTensor);
     material->SetTensor( couplingTensor, MAGNETOSTRICTION_TENSOR_h_mech, Global::REAL );
     }
     if(mst->Has("imag")) {
     ParamTools::AsTensor<double>(mst->Get("imag"), 3, 6, couplingTensor);
     material->SetTensor( couplingTensor, MAGNETOSTRICTION_TENSOR_h_mech, Global::IMAG );
     }
     }
     */
  }
  
  StdVector<std::string> XMLMaterialHandler::ReadMechanicIsotropic(PtrParamNode node)
  {
    StdVector<std::string> vals(36);
    vals.Init("0");
    std::string mu, lam;
    // read as strings and add formulas if lame parameters are not given directly
    if(node->Has("elasticityModulus") && node->Has("poissonNumber"))
    {
      std::string E = Bracket(node->Get("elasticityModulus")->As<std::string>());
      std::string nu = Bracket(node->Get("poissonNumber")->As<std::string>());
      std::stringstream sLameMu,sLameLambda;
      sLameMu << Bracket(E) << "/(2*(1+"<<Bracket(nu)<<"))";
      sLameLambda << nu<<"*"<<E<<"/((1+"<<nu<<")*(1-2*"<<nu<<"))";
      mu = E+"/(2*(1+"+nu+"))";// E/(2*(1+nu))
      lam = nu+"*"+E+"/((1+"+nu+")*(1-2*"+nu+"))";// nu*E/((1+nu)*(1-2*nu))
      //material->SetScalar(sLameMu.str(),MECH_LAME_MU,dataType);
      //material->SetScalar(sLameLambda.str(),MECH_LAME_LAMBDA,dataType);
    }
    else if (node->Has("compressionModulus") && node->Has("shearModulus"))
    {
      std::string K = Bracket(node->Get("compressionModulus")->As<std::string>());
      std::string G = Bracket(node->Get("shearModulus")->As<std::string>());
      std::stringstream sLameLambda;
      sLameLambda <<"("<< K << "-2*"<< G << "/3)";
      mu = G;
      lam = Bracket( K+"-2/3*"+G );//sLameLambda.str();
      //material->SetScalar(G,MECH_LAME_MU,dataType);
      //material->SetScalar(sLameLambda.str(),MECH_LAME_LAMBDA,dataType);
    }
    else if (node->Has("lameParameterMu") && node->Has("lameParameterLambda"))
    {
      mu = Bracket(node->Get("lameParameterMu")->As<std::string>());
      lam = Bracket(node->Get("lameParameterLambda")->As<std::string>());
      //material->SetScalar( node->Get("lameParameterMu")->As<std::string>(), MECH_LAME_MU, dataType);
      //material->SetScalar( node->Get("lameParameterLambda")->As<std::string>(), MECH_LAME_LAMBDA, dataType);
    }
    else {
      EXCEPTION ("Specify E+nu, K+B, or lambda+mu!" );
    }
    // build the tensor elements
    std::string diag = Bracket("2*"+mu+"+"+lam); // 2*mu+lam
    vals[0] = diag;
    vals[1] = lam;
    vals[2] = lam;
    vals[6] = lam;
    vals[7] = diag;
    vals[8] = lam;
    vals[12] = lam;
    vals[13] = lam;
    vals[14] = diag;
    vals[21] = mu;
    vals[28] = mu;
    vals[35] = mu;
    return vals;
  }
  
  
  //**********************************************************************
  //*************  READ ACOUSTICS ****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadAcoustic(BaseMaterial *material, PtrParamNode acou)
  {
    
    //read density
    if(acou->Has("density")) {
      PtrCoefFct densFct =
              CoefFunction::Generate(mp_, Global::REAL,
              acou->Get("density")->As<std::string>() );
      material->SetCoefFct( DENSITY, densFct );
    }
    
    //check for complex valued density
    if ( acou->Has("densityComplex") ) {
    	PtrParamNode densNode = acou->Get("densityComplex");
      
    	// read the real part
    	std::string realStr = densNode->Get("real")->As<std::string>();
      
    	// read the imaginary part
    	std::string imagStr = densNode->Get("imag")->As<std::string>();
      
    	PtrCoefFct densFct = CoefFunction::Generate( mp_, Global::COMPLEX,
              realStr, imagStr );
      
      material->SetCoefFct( ACOU_DENSITY_COMPLEX, densFct );
    }
    
    //read dauabatic exponent
    if(acou->Has("adiabaticExponent")) {
      PtrCoefFct fct =
              CoefFunction::Generate(mp_, Global::REAL,
              acou->Get("adiabaticExponent")->As<std::string>() );
      material->SetCoefFct( ADIABATIC_EXPONENT, fct );
    }
    
    //read compression modulus
    if(acou->Has("compressionModulus")) { 
    	PtrCoefFct blkFct =
              CoefFunction::Generate(mp_, Global::REAL, 
              acou->Get("compressionModulus")->As<std::string>() );
    	material->SetCoefFct( ACOU_BULK_MODULUS, blkFct );
    }
    
    //check for complex valued density
    if ( acou->Has("compressionModulusComplex") ) {
    	PtrParamNode compNode = acou->Get("compressionModulusComplex");
      
    	// read the real part
    	std::string realStr = compNode->Get("real")->As<std::string>();
      
    	// read the imaginary part
    	std::string imagStr = compNode->Get("imag")->As<std::string>();
      
    	PtrCoefFct compFct = CoefFunction::Generate( mp_, Global::COMPLEX,
              realStr, imagStr );
      
      material->SetCoefFct( ACOU_BULK_MODULUS_COMPLEX, compFct );
    }
    
    //read kinematic viscosity
    if(acou->Has("kinematicViscosity")) {
      PtrCoefFct kinVisc =
              CoefFunction::Generate(mp_, Global::REAL,
              acou->Get("kinematicViscosity")->As<std::string>() );
      material->SetCoefFct( KINEMATIC_VISCOSITY, kinVisc );
    }
    // check for acousticDamping
    if(acou->Has("acousticDamping"))
    {
      PtrParamNode ad = acou->Get("acousticDamping");
      
      // check rayleigh
      if(ad->Has("rayleigh"))
      {
        PtrParamNode r = ad->Get("rayleigh");
        
        if(r->Has("alpha"))
          material->SetScalar(r->Get("alpha")->As<std::string>(), RAYLEIGH_ALPHA);
        
        if(r->Has("beta"))
          material->SetScalar(r->Get("beta")->As<std::string>(), RAYLEIGH_BETA);
        
        if(r->Has("lossTangensDelta"))
          material->SetScalar(r->Get("lossTangensDelta")->As<std::string>(), LOSS_TANGENS_DELTA);
        
        if(r->Has("measuredFreq"))
          material->SetScalar(r->Get("measuredFreq")->As<Double>(), RAYLEIGH_FREQUENCY, Global::REAL );
      } // end of acousticDamping:rayleigh
      
      // read alpha0 of thermo viscous damping
      if(ad->Has("thermoViscous"))
      {
        if(ad->Get("thermoViscous")->Has("alpha0"))
          material->SetScalar(ad->Get("thermoViscous")->Get("alpha0")->As<Double>(), ACOU_ALPHA, Global::REAL );
      }
      // read fractional damping
      if(ad->Has("fractional"))
      {
        PtrParamNode f = ad->Get("fractional");
        
        if(f->Has("alpha0")) 
          material->SetScalar(f->Get("alpha0")->As<Double>(), ACOU_ALPHA, Global::REAL );
        
        // read exponent of fractional damping      
        if(f->Has("y")) 
          material->SetScalar(f->Get("y")->As<Double>(), FRACTIONAL_EXPONENT, Global::REAL );
      }
    } // end of acousticDamping
    
    //read acoustic non linearity
    if(acou->Has("acousticNonlinear"))
    {
      if(acou->Get("acousticNonlinear")->Has("bOverA"))
        material->SetScalar(acou->Get("acousticNonlinear")->Get("bOverA")->As<Double>(), BOVERA, Global::REAL );
    }  
  }
  
  
  //**********************************************************************
  //*************  READ ELECTROSTATICS ************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadElectrostatic(BaseMaterial *material, PtrParamNode elec)
  {
    //    std::cout << "ReadElectrostatic" << std::endl;
    
    // check for permittivity
    if(elec->Has("permittivity"))
    {
      PtrParamNode p = elec->Get("permittivity");
      
      // check for tensor with dim1 = 3 <tensor dim1="3">
      if(p->HasByVal("tensor", "dim1", "3")) {
        StdVector<std::string> realVals(9), imagVals(9);
        realVals.Init("0.0");
        imagVals.Init("0.0");
        PtrCoefFct epsCoef;
        
        // read real permittivity tensor 
        if(p->GetByVal("tensor", std::string("dim1"), "3")->Has("real"))
        {
          PtrParamNode tensor =  p->GetByVal("tensor","dim1","3")->Get("real");        
          ParamTools::AsStringTensor( tensor, 9, realVals );
        }
        
        // read imaginary permittivity tensor
        if(p->GetByVal("tensor", "dim1", "3")->Has("imag"))
        {
          PtrParamNode tensor =  p->GetByVal("tensor","dim1","3")->Get("imag");
          ParamTools::AsStringTensor( tensor, 9, imagVals );
        }
        epsCoef = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 3,
                realVals, imagVals );
        material->SetCoefFct( ELEC_PERMITTIVITY, epsCoef);
      } // end of <tensor dim1="3">
      
      
      
    } // end of permittivity
    
    // check for <permittivityCoefficient nonlinear="function">
    if(elec->HasByVal("permittivityCoefficient", "nonlinear", "function"))
    {
      PtrParamNode pc = 
              elec->GetByVal("permittivityCoefficient","nonlinear","function");
      
      // read nonlinearity of a permittivity coefficient
      if(pc->Has("entry"))
        material->SetScalar(pc->Get("entry")->As<Integer>(), NONLIN_COEFFICIENT);
      
      // read non linear dependency of a permittivity coefficient
      if(pc->Has("dependency"))
        material->SetScalar(pc->Get("dependency")->As<std::string>(), NONLIN_DEPENDENCY);
      
      // read non linear approxType of a permittivity coefficient
      if(pc->Has("approxType"))
        material->SetScalar(pc->Get("approxType")->As<std::string>(), NONLIN_APPROXIMATION_TYPE);
      
      // read non linear data name of a permittivity coefficient        
      if(pc->Has("dataName"))
        material->SetScalar(pc->Get("dataName")->As<std::string>(), NONLIN_DATA_NAME);
    } // end of permittivityCoefficient
    
    //read hysterese model
    if(elec->Has("hystModel"))
    {
      PtrParamNode hystNode = elec->Get("hystModel");
      ReadHysteresis(material, hystNode);
    }
  }
  
  //**********************************************************************
  //*************  READ MAGNETIC *****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadMagnetic(BaseMaterial *material, PtrParamNode mag)
  {
    // flags for symmetry type
    bool isIsotropic = false;
    bool isOrthotropic = false;
    bool isTensor = false;
    
    // read electric conductivity
    if(mag->Has("electricConductivity")){
      
      if (mag->Get("electricConductivity")->Has("linear")) {
        PtrParamNode lin = mag->Get("electricConductivity")->Get("linear");
        
        // === ISOTROPIC ===
        if (lin->Has("isotropic")) {
          material->SetScalar(lin->Get("isotropic")->As<Double>(), MAG_CONDUCTIVITY, Global::REAL);
          material->SetSymmetryType(MAG_CONDUCTIVITY,BaseMaterial::ISOTROPIC);
        }
        else if (lin->Has("tensor")){
          EXCEPTION("For magnetic simulation, no tensor-valued el. conductivity allowed");          }
      }
      
      // we know only nonlinear isotropic material
      if(mag->Get("electricConductivity")->Has("nonlinear") &&
              mag->Get("electricConductivity")->Get("nonlinear")->Has("isotropic"))
      {
        PtrParamNode iso = mag->Get("electricConductivity")->Get("nonlinear")->Get("isotropic");
        BaseMaterial::MatDescriptorNl info;
        info. approxType = NO_APPROX_TYPE;
        info.measAccuracy = 0.01;
        info.maxVal = 1000;
        info.fileName = "";
        info.factor = 1.;
        
        // read dependency: can be temperature
        if(iso->Has("dependency"))
          material->SetScalar(iso->Get("dependency")->As<std::string>(), NONLIN_DEPENDENCY);
        
        // read approximation type
        if(iso->Has("approxType")) {
          std::string type =  iso->Get("approxType")->As<std::string>();
          info.approxType = ApproxCurveTypeEnum.Parse(type );
        }
        
        // read measurement accuracy
        if(iso->Has("measAccuracy"))
          info.measAccuracy = iso->Get("measAccuracy")->As<Double>();
        
        // read maximum value for approximation
        if(iso->Has("maxApproxVal"))
          info.maxVal = iso->Get("maxApproxVal")->As<Double>();
        
        // read name of function file
        if(iso->Has("dataName"))
          info.fileName = iso->Get("dataName")->As<std::string>().c_str();
        
        // read factor
        if(iso->Has("factor"))
          info.factor = iso->Get("factor")->As<Double>();
        
        //set info to material class
        material->SetNonLinMatIso(MAG_CONDUCTIVITY, info);
      } // nonlinear isotropic material
      
    }
    
    
    // read nonlinear reluctivity for magnetostrictive strains
    if(mag->Has("magneticReluctivity_MagStrict"))
    {
      
      // we know only nonlinear isotropic material
      if(mag->Get("magneticReluctivity_MagStrict")->Has("nonlinear") ) {
        if (mag->Get("magneticReluctivity_MagStrict")->Get("nonlinear")->Has("isotropic")) {
          PtrParamNode iso = mag->Get("magneticReluctivity_MagStrict")->Get("nonlinear")->Get("isotropic");
          BaseMaterial::MatDescriptorNl info;
          info.approxType = NO_APPROX_TYPE;
          info.measAccuracy = 0.01;
          info.maxVal = 2.5;
          info.fileName = "";
          
          // read approximation type  
          if(iso->Has("approxType")) {
            std::string type =  iso->Get("approxType")->As<std::string>();
            info.approxType = ApproxCurveTypeEnum.Parse(type );
          }
          
          // read nonlinear approxType of magnetic permeability
          if(iso->Has("measAccuracy"))
            info.measAccuracy = iso->Get("measAccuracy")->As<Double>();
          
          // read nonlinear approxType of magnetic permeability
          if(iso->Has("maxApproxVal"))
            info.maxVal = iso->Get("maxApproxVal")->As<Double>();
          
          // read nonlinear dataName of magnetic permeability
          if(iso->Has("dataName"))
            info.fileName = iso->Get("dataName")->As<std::string>();
          
          // pass info to material class
          material->SetNonLinMatIso( MAGSTRICT_RELUCTIVITY, info);
        } // nonlinear isotropic material   
        
      } // end of nonlinear section
    } // end of magneticReluctivity_MagStrict  
    
    
    // read magnetic permeability
    if(mag->Has("magneticPermeability"))
    {
      if(mag->Get("magneticPermeability")->Has("linear"))
      {
        PtrParamNode lin = mag->Get("magneticPermeability")->Get("linear");
        double eps = 1e-10;
        
        // === ISOTROPIC ===
        if(lin->Has("isotropic"))
        {
          if(lin->Get("isotropic")->As<Double>() < eps)
            EXCEPTION("Magnetic permeability is near zero. Check material database");
          material->SetScalar(lin->Get("isotropic")->As<Double>(), MAG_PERMEABILITY, Global::REAL );
          isIsotropic = true;
        }
        
        // === ORTHOTROPIC ===
        if(lin->Has("orthotropic"))
        {
          PtrParamNode ortho = lin->Get("orthotropic");
          bool permOrtho_1=false, permOrtho_2=false, permOrtho_3=false;
          
          if(ortho->Has("permeability_1"))
          {
            if(ortho->Get("permeability_1")->As<Double>() < eps)
              EXCEPTION("Magnetic permeability is near zero; Check material database");
            material->SetScalar(ortho->Get("permeability_1")->As<Double>(), MAG_PERMEABILITY_1, Global::REAL); 
            permOrtho_1 = true;  
          }
          
          if(ortho->Has("permeability_2"))
          {
            if(ortho->Get("permeability_2")->As<Double>() < eps)
              EXCEPTION("Magnetic permeability is near zero; Check material database");
            material->SetScalar(ortho->Get("permeability_2")->As<Double>(), MAG_PERMEABILITY_2, Global::REAL); 
            permOrtho_2 = true;  
          }
          
          if(ortho->Has("permeability_3"))
          {
            if(ortho->Get("permeability_3")->As<Double>() < eps)
              EXCEPTION("Magnetic permeability is near zero; Check material database");
            material->SetScalar(ortho->Get("permeability_3")->As<Double>(), MAG_PERMEABILITY_3, Global::REAL); 
            permOrtho_3 = true;  
          }
          
          // check, if there is an orthotropic permeability!!
          if (permOrtho_1 == true && permOrtho_2 == true && permOrtho_3 == true)
            isOrthotropic = true;
        } // end of linear orthotropic
        
        // === TENSOR / GENERAL ===
        if(lin->Has("tensor")) {
          
          Matrix<Double> muTensor(3,3);
          
          // read permeability tensor (real part)
          if(lin->GetByVal("tensor", std::string("dim1"), "3")->Has("real")) {
            PtrParamNode tens = 
                    lin->GetByVal("tensor", std::string("dim1"),"3" )->Get("real");
            ParamTools::AsTensor<double>(tens, 3, 3, muTensor); 
            material->SetTensor(muTensor, MAG_PERMEABILITY, Global::REAL);
            isTensor = true;
          }
          
          // read permeability tensor (imaginary part)
          if(lin->GetByVal("tensor", std::string("dim1"), "3")->Has("imag")) {
            PtrParamNode tens = 
                    lin->GetByVal("tensor", std::string("dim1"),"3" )->Get("imag");
            ParamTools::AsTensor<double>(tens, 3, 3, muTensor);
            material->SetTensor(muTensor, MAG_PERMEABILITY, Global::IMAG);
          }
        } // tensor
      } // end of linear
      
      // Try to determine, if a unique symmetry type can be obtained
      if( isIsotropic && !isOrthotropic && !isTensor ) {
        material->SetSymmetryType(MAG_PERMEABILITY,BaseMaterial::ISOTROPIC);
        
      } else if( !isIsotropic && isOrthotropic && !isTensor ) {
        material->SetSymmetryType(MAG_PERMEABILITY,BaseMaterial::ORTHOTROPIC );
        
      } else if( !isIsotropic && !isOrthotropic && isTensor ) {
        material->SetSymmetryType(MAG_PERMEABILITY,BaseMaterial::GENERAL );
        
      } else {
        EXCEPTION("Could not determine unique material symmetry. "
                "Please ensure, that only ONE of isotropic, orthotropic or"
                "the material tensor of the permeability are given!")
      }
      
      // we know only nonlinear isotropic material
      if(mag->Get("magneticPermeability")->Has("nonlinear") ) {
        if (mag->Get("magneticPermeability")->Get("nonlinear")->Has("isotropic")) {
          PtrParamNode iso = mag->Get("magneticPermeability")->Get("nonlinear")->Get("isotropic");
          BaseMaterial::MatDescriptorNl info;
          info.approxType = NO_APPROX_TYPE;
          info.measAccuracy = 0.01;
          info.maxVal = 2.5;
          info.fileName = "";
          
          // read approximation type  
          if(iso->Has("approxType")) {
            std::string type =  iso->Get("approxType")->As<std::string>();
            info.approxType = ApproxCurveTypeEnum.Parse(type );
          }
          
          // read nonlinear approxType of magnetic permeability
          if(iso->Has("measAccuracy"))
            info.measAccuracy = iso->Get("measAccuracy")->As<Double>();
          
          // read nonlinear approxType of magnetic permeability
          if(iso->Has("maxApproxVal"))
            info.maxVal = iso->Get("maxApproxVal")->As<Double>();
          
          // read nonlinear dataName of magnetic permeability
          if(iso->Has("dataName"))
            info.fileName = iso->Get("dataName")->As<std::string>();
          
          // read analytic function of material parameter
          if(iso->Has("nuExpr"))
            info.analyticExpr = iso->Get("nuExpr")->As<std::string>();
          
          // read analytic derivative of material parameter
          if(iso->Has("nuDerivExpr"))
            info.analyticExprDeriv = iso->Get("nuDerivExpr")->As<std::string>();
          
          // pass info to material class
          material->SetNonLinMatIso( MAG_PERMEABILITY, info);
        } // nonlinear isotropic material   
        
        else if (mag->Get("magneticPermeability")->Get("nonlinear")->Has("anisotropic")) {
          
          //anisotropic case: bundle of nonlinear curves
          PtrParamNode nonLin = mag->Get("magneticPermeability")->Get("nonlinear")->Get("anisotropic");
          
          // fetch paramnodes for hdbc
          ParamNodeList anIsoNodes = nonLin->GetList("data");
          
          if ( anIsoNodes.GetSize() > 0 ) {
            StdVector<BaseMaterial::MatDescriptorNl> nlData;
            nlData.Resize(anIsoNodes.GetSize());
            
            // iterate over all parameter nodes
            for( UInt i = 0; i < anIsoNodes.GetSize(); i++ ) {
              // read parameters
              BaseMaterial::MatDescriptorNl info;
              info.angle = 0.0;
              info.approxType = NO_APPROX_TYPE;
              info.measAccuracy = 0.01;
              info.maxVal = 2.5;
              info.zScaling= 1.0;
              info.fileName = "";              
              info.analyticExpr = "";
              info.analyticExprDeriv = "";
              
              // read angle  
              if(anIsoNodes[i]->Has("angle")) {
                info.angle =  anIsoNodes[i]->Get("angle")->As<Double>();
              }
              
              // read approximation type  
              if(anIsoNodes[i]->Has("approxType")) {
                std::string type =  anIsoNodes[i]->Get("approxType")->As<std::string>();
                info.approxType = ApproxCurveTypeEnum.Parse(type );
              }
              
              // read measurement accuracy
              if(anIsoNodes[i]->Has("measAccuracy")) 
                info.measAccuracy = anIsoNodes[i]->Get("measAccuracy")->As<Double>();
              
              // read name of function file
              if(anIsoNodes[i]->Has("dataName"))
                info.fileName = anIsoNodes[i]->Get("dataName")->As<std::string>().c_str();
              
              // read maximum value for approximation
              if(anIsoNodes[i]->Has("maxApproxVal")) 
                info.maxVal = anIsoNodes[i]->Get("maxApproxVal")->As<Double>();
              
              // read z-scaling factor  
              if(anIsoNodes[i]->Has("angle"))
                info.zScaling =  anIsoNodes[i]->Get("zScaling")->As<Double>();
              
              // read analytic function of material parameter
              if(anIsoNodes[i]->Has("nuExpr")) {
                info.analyticExpr = anIsoNodes[i]->Get("nuExpr")->As<std::string>().c_str();
              }
              // read analytic derivative of material parameter
              if(anIsoNodes[i]->Has("nuDerivExpr"))
                info.analyticExprDeriv = anIsoNodes[i]->Get("nuDerivExpr")->As<std::string>().c_str();
              
              nlData[i].angle        = info.angle;
              nlData[i].fileName     = info.fileName;
              nlData[i].approxType   = info.approxType;
              nlData[i].measAccuracy = info.measAccuracy;
              nlData[i].maxVal       = info.maxVal;
              nlData[i].analyticExpr = info.analyticExpr;
              nlData[i].analyticExprDeriv = info.analyticExprDeriv;
              nlData[i].zScaling          = info.zScaling;
            }
            material->SetNonLinMatAniso( MAG_PERMEABILITY, nlData );
          }        
          
        } // end of anisotropic nonlinear material
      } // end of nonlinear section
    } // end of magneticPermeability  
    
    //read hysterese model
    if(mag->Has("hystModel"))
    {
      PtrParamNode hystNode = mag->Get("hystModel");
      ReadHysteresis(material, hystNode);
    }
    
    //read real magmech coupling tensor
    if(mag->Has("magnetoStrictionTensor_h_mag"))
    {
      StdVector<std::string> realVals(18), imagVals(18);
      realVals.Init("0.0");
      imagVals.Init("0.0");
      PtrCoefFct pctCoef;
      PtrParamNode pct = mag->Get("magnetoStrictionTensor_h_mag");
      if(pct->Has("real")){
        ParamTools::AsStringTensor( pct->Get("real"), 18, realVals );
      }
      if(pct->Has("imag"))
      {
        ParamTools::AsStringTensor( pct->Get("imag"), 18, imagVals );
      }
      pctCoef = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 6,
              realVals, imagVals );
      material->SetCoefFct( MAGNETOSTRICTION_TENSOR_h_mag, pctCoef);
    }
    /* old style
     if(mag->Has("magnetoStrictionTensor_h_mag")) {
     Matrix<Double> couplingTensor(3,6);
     
     PtrParamNode mst = mag->Get("magnetoStrictionTensor_h_mag");
     if(mst->Has("real")) {
     ParamTools::AsTensor<double>(mst->Get("real"), 3, 6, couplingTensor);
     material->SetTensor( couplingTensor, MAGNETOSTRICTION_TENSOR_h_mag, Global::REAL );
     }
     if(mst->Has("imag")) {
     ParamTools::AsTensor<double>(mst->Get("imag"), 3, 6, couplingTensor);
     material->SetTensor( couplingTensor, MAGNETOSTRICTION_TENSOR_h_mag, Global::IMAG );
     }
     }
     */
    // read core loss
    if(mag->Has("coreLoss")){
      PtrParamNode clParam = mag->Get("coreLoss");
      BaseMaterial::MatDescriptorNl info;
      info.approxType = LIN_INTERPOLATE;
      info.fileName = "";
      
      if(clParam->Has("approxType")) {
        std::string type =  clParam->Get("approxType")->As<std::string>();
        info.approxType = ApproxCurveTypeEnum.Parse( type );
      }
      
      if(clParam->Has("measAccuracy"))
        info.measAccuracy = clParam->Get("measAccuracy")->As<Double>();
      
      if(clParam->Has("maxApproxVal"))
        info.maxVal = clParam->Get("maxApproxVal")->As<Double>();
      
      if(clParam->Has("dataName"))
        info.fileName = clParam->Get("dataName")->As<std::string>();
      
      material->SetNonLinMatIso( CORE_LOSS, info );
    }
    
    // read density needed for core loss
    PtrCoefFct densFct;
    if(mag->Has("density")){
      densFct = CoefFunction::Generate( mp_, Global::REAL,
              mag->Get("density")->As<std::string>() );
      material->SetCoefFct( DENSITY, densFct );
    } else {
      densFct = CoefFunction::Generate( mp_, Global::REAL, "0.0");
      material->SetCoefFct( DENSITY, densFct );
    }
    
  }
  
  //**********************************************************************
  //*************  READ THERMIC ******************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadThermic(BaseMaterial *material, PtrParamNode therm)
  {
    // read density
    if(therm->Has("density")) {
      PtrCoefFct densFct =
              CoefFunction::Generate(mp_, Global::REAL, 
              therm->Get("density")->As<std::string>() );
      material->SetCoefFct( DENSITY, densFct );
    }
    
    // read heat capacity
    if(therm->Has("heatCapacity")) {
      if(therm->Get("heatCapacity")->Has("linear")) {
        PtrParamNode lin = therm->Get("heatCapacity")->Get("linear");
        
        // === ISOTROPIC ===
        if(lin->Has("isotropic")) {
          material->SetScalar(lin->Get("isotropic")->As<Double>(),  HEAT_CAPACITY, Global::REAL );
          //          isIsotropic = true;
        }
      }
      // we know only nonlinear isotropic material
      if(therm->Get("heatCapacity")->Has("nonlinear") && 
              therm->Get("heatCapacity")->Get("nonlinear")->Has("isotropic")) {
        PtrParamNode iso = therm->Get("heatCapacity")->Get("nonlinear")->Get("isotropic");
        
        BaseMaterial::MatDescriptorNl info;
        info. approxType = NO_APPROX_TYPE;
        info.measAccuracy = 0.01;
        info.maxVal = 1000;
        info.fileName = "";
        // read approximation type  
        if(iso->Has("approxType")) {
          std::string type =  iso->Get("approxType")->As<std::string>();
          info.approxType = ApproxCurveTypeEnum.Parse(type);
        }
        // read measurement accuracy
        if(iso->Has("measAccuracy")) 
          info.measAccuracy = iso->Get("measAccuracy")->As<Double>();
        
        // read maximum value for approximation
        if(iso->Has("maxApproxVal")) 
          info.maxVal = iso->Get("maxApproxVal")->As<Double>();
        
        // read name of function file 
        if(iso->Has("dataName")) 
          info.fileName = iso->Get("dataName")->As<std::string>().c_str();
        
        //set info to material class
        material->SetNonLinMatIso(HEAT_CAPACITY, info);
        
      } // nonlinear isotropic material  
    }
    
    // read thermal conductivity
    if(therm->Has("heatConductivity")) {
      if (therm->Get("heatConductivity")->Has("linear")) {
        PtrParamNode lin = therm->Get("heatConductivity")->Get("linear");
        
        // === ISOTROPIC ===
        if (lin->Has("isotropic")) {
          material->SetScalar(lin->Get("isotropic")->As<Double>(), HEAT_CONDUCTIVITY, Global::REAL);
        }
        else if (lin->Has("tensor")){
          
          // can only be a real 3x3 tensor
          Matrix<double> tensor(3,3);
          PtrParamNode tens_pn = 
                  lin->GetByVal("tensor","dim1","3")->Get("real");
          
          ParamTools::AsTensor<double>(tens_pn, 3, 3, tensor);
          material->SetTensor(tensor, HEAT_CONDUCTIVITY_TENSOR, Global::REAL);
          material->SetSymmetryType(HEAT_CONDUCTIVITY_TENSOR,BaseMaterial::GENERAL);
          
          //set mean value of diagonal entries
          Double meanCond = 0.0;
          for ( UInt i=0; i<tensor.GetNumRows(); i++)
            meanCond += tensor[i][i];
          meanCond /= (Double)tensor.GetNumRows();
          material->SetScalar(meanCond, HEAT_CONDUCTIVITY, Global::REAL);
        }
      }
      
      // we know only nonlinear isotropic material
      if(therm->Get("heatConductivity")->Has("nonlinear") && 
              therm->Get("heatConductivity")->Get("nonlinear")->Has("isotropic"))
      {
        PtrParamNode iso = therm->Get("heatConductivity")->Get("nonlinear")->Get("isotropic");
        BaseMaterial::MatDescriptorNl info;
        info. approxType = NO_APPROX_TYPE;
        info.measAccuracy = 0.01;
        info.maxVal = 1000;
        info.fileName = "";
        
        // read approximation type  
        if(iso->Has("approxType")) {
          std::string type =  iso->Get("approxType")->As<std::string>();
          info.approxType = ApproxCurveTypeEnum.Parse(type );
        }
        
        // read measurement accuracy
        if(iso->Has("measAccuracy")) 
          info.measAccuracy = iso->Get("measAccuracy")->As<Double>();
        
        // read maximum value for approximation
        if(iso->Has("maxApproxVal")) 
          info.maxVal = iso->Get("maxApproxVal")->As<Double>();
        
        // read name of function file 
        if(iso->Has("dataName")) 
          info.fileName = iso->Get("dataName")->As<std::string>().c_str();
        
        //set info to material class
        material->SetNonLinMatIso(HEAT_CONDUCTIVITY, info);
        //material->SetSymmetryType(HEAT_CONDUCTIVITY,BaseMaterial::ISOTROPIC);
      } // nonlinear isotropic material  
      
    }
  }
  
  //**********************************************************************
  //*************  READ FLOW *********************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadFlow(BaseMaterial *material, PtrParamNode flow)
  {    
    // read density
    if(flow->Has("density")) {
      PtrCoefFct densFct =
              CoefFunction::Generate(mp_, Global::REAL, 
              flow->Get("density")->As<std::string>() );
      material->SetCoefFct( DENSITY, densFct );
    }
    
    // read dynamicViscosity 
    if( flow->Has("dynamicViscosity") &&
            flow->Has("kinematicViscosity") ) {
      EXCEPTION("Please specify either dynamic or kinematic viscosity but not both!");
    }
    
    if(flow->Has("dynamicViscosity")) {
      PtrCoefFct dynVisc =
              CoefFunction::Generate(mp_, Global::REAL, 
              flow->Get("dynamicViscosity")->As<std::string>() );
      material->SetCoefFct( DYNAMIC_VISCOSITY, dynVisc );
    }
    
    // read kinematicViscosity 
    if(flow->Has("kinematicViscosity")) {
      PtrCoefFct kinVisc =
              CoefFunction::Generate(mp_, Global::REAL, 
              flow->Get("kinematicViscosity")->As<std::string>() );
      material->SetCoefFct( KINEMATIC_VISCOSITY, kinVisc );
    }
  }
  
  //**********************************************************************
  //*************  READ TEST *********************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadTest(BaseMaterial *material, PtrParamNode test)
  {
    // read alpha
    if(test->Has("alpha")) {
      PtrCoefFct alphaFct =
              CoefFunction::Generate(mp_, Global::REAL,
              test->Get("alpha")->As<std::string>() );
      material->SetCoefFct( TEST_ALPHA, alphaFct );
    }
    
    // read beta
    if(test->Has("beta")) {
      PtrCoefFct betaFnc =
              CoefFunction::Generate(mp_, Global::REAL,
              test->Get("beta")->As<std::string>() );
      material->SetCoefFct( TEST_BETA, betaFnc );
    }
  }
  
  
  //**********************************************************************
  //*************  READ PYROELECTRIC *************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadPyroelectric(BaseMaterial *material, 
          PtrParamNode pyro){
    
    if (pyro->Has("pyrocoefficient")){
      PtrParamNode py = pyro->Get("pyrocoefficient");
      if(py->Has("tensor"))
      {
        // can only be a real 3x3 tensor
        Matrix<double> tensor(3,3);
        PtrParamNode tens_pn = 
                py->GetByVal("tensor","dim1","3")->Get("real");
        ParamTools::AsTensor<double>(tens_pn, 3, 3, tensor);
        material->SetTensor(tensor,PYROCOEFFICIENT_TENSOR,Global::REAL);
      }
    }
  }
  
  //**********************************************************************
  //*************  READ THERMOELASTIC ************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadThermoelastic(BaseMaterial *material,
          PtrParamNode thermExp) {
    
    if(thermExp->Has("thermalExpansion")){
      PtrParamNode te = thermExp->Get("thermalExpansion");
      if(te->Has("tensor"))
      {
        // can only be a real 3x3 tensor
        Matrix<double> tensor(3,3);
        PtrParamNode tens_pn = 
                te->GetByVal("tensor","dim1","3")->Get("real");
        ParamTools::AsTensor<double>(tens_pn, 3, 3, tensor);
        material->SetTensor(tensor,THERMAL_EXPANSION_TENSOR,Global::REAL);
      }
    }
  }
  
  //**********************************************************************
  //*************  READ MAGNETOSTRICTIVE  ********************************
  //**********************************************************************
  void XMLMaterialHandler::ReadMagStrict(BaseMaterial *material,
          PtrParamNode pn) {
    //read real magmech coupling tensor
    if(pn->Has("magnetoStrictionTensor_h"))
    {
      StdVector<std::string> realVals(18), imagVals(18);
      realVals.Init("0.0");
      imagVals.Init("0.0");
      PtrCoefFct pctCoef;
      PtrParamNode pct = pn->Get("magnetoStrictionTensor_h");
      if(pct->Has("real")){
        ParamTools::AsStringTensor( pct->Get("real"), 18, realVals );
      }
      if(pct->Has("imag"))
      {
        ParamTools::AsStringTensor( pct->Get("imag"), 18, imagVals );
      }
      pctCoef = CoefFunction::Generate( mp_, Global::COMPLEX, 3, 6,
              realVals, imagVals );
      material->SetCoefFct( MAGNETOSTRICTION_TENSOR_h, pctCoef);
    }
    
    /* old way
     if(pn->Has("magnetoStrictionTensor_h")) {
     Matrix<Double> couplingTensor(3,6);
     
     PtrParamNode mst = pn->Get("magnetoStrictionTensor_h");
     if(mst->Has("real")) {
     ParamTools::AsTensor<double>(mst->Get("real"), 3, 6, couplingTensor);
     material->SetTensor( couplingTensor, MAGNETOSTRICTION_TENSOR_h, Global::REAL );
     }
     if(mst->Has("imag")) {
     ParamTools::AsTensor<double>(mst->Get("imag"), 3, 6, couplingTensor);
     material->SetTensor( couplingTensor, MAGNETOSTRICTION_TENSOR_h, Global::IMAG );
     }
     }
     */
  }
  
  //**********************************************************************
  //*************  READ ELECTRIC CONDUCTIVITY ******************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadElectricConduction(BaseMaterial *material, PtrParamNode elec)
  {
    // read electric conductivity
    if(elec->Has("elecConductivity")) {
      if (elec->Get("elecConductivity")->Has("linear")) {
        PtrParamNode lin = elec->Get("elecConductivity")->Get("linear");
        
        // === ISOTROPIC ===
        if (lin->Has("isotropic")) {
          material->SetScalar(lin->Get("isotropic")->As<Double>(), ELEC_CONDUCTIVITY, Global::REAL);
          material->SetSymmetryType(ELEC_CONDUCTIVITY,BaseMaterial::ISOTROPIC);
        }
        else if (lin->Has("tensor")){
          // can only be a real 3x3 tensor
          Matrix<double> tensor(3,3);
          PtrParamNode tens_pn =
                  lin->GetByVal("tensor","dim1","3")->Get("real");
          
          ParamTools::AsTensor<double>(tens_pn, 3, 3, tensor);
          material->SetTensor(tensor, ELEC_CONDUCTIVITY_TENSOR, Global::REAL);
          material->SetSymmetryType(ELEC_CONDUCTIVITY_TENSOR,BaseMaterial::GENERAL);
        }
      }
      
      // we know only nonlinear isotropic material
      if(elec->Get("elecConductivity")->Has("nonlinear") &&
              elec->Get("elecConductivity")->Get("nonlinear")->Has("isotropic"))
      {
        PtrParamNode iso = elec->Get("elecConductivity")->Get("nonlinear")->Get("isotropic");
        BaseMaterial::MatDescriptorNl info;
        info. approxType = NO_APPROX_TYPE;
        info.measAccuracy = 0.01;
        info.maxVal = 1000;
        info.fileName = "";
        info.factor = 1.;
        
        // read dependency: can be "temperature" or voltage
        if(iso->Has("dependency"))
          material->SetScalar(iso->Get("dependency")->As<std::string>(), NONLIN_DEPENDENCY);
        
        // read approximation type
        if(iso->Has("approxType")) {
          std::string type =  iso->Get("approxType")->As<std::string>();
          info.approxType = ApproxCurveTypeEnum.Parse(type );
        }
        
        // read measurement accuracy
        if(iso->Has("measAccuracy"))
          info.measAccuracy = iso->Get("measAccuracy")->As<Double>();
        
        // read maximum value for approximation
        if(iso->Has("maxApproxVal"))
          info.maxVal = iso->Get("maxApproxVal")->As<Double>();
        
        // read name of function file
        if(iso->Has("dataName"))
          info.fileName = iso->Get("dataName")->As<std::string>().c_str();
        
        // read factor
        if(iso->Has("factor"))
          info.factor = iso->Get("factor")->As<Double>();
        
        //set info to material class
        material->SetNonLinMatIso(ELEC_CONDUCTIVITY, info);
      } // nonlinear isotropic material
      
    }
  }
  
  //**********************************************************************
  //*************  READ HYSTERESIS****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadHysteresis(BaseMaterial *material, PtrParamNode hystNode){
    PtrParamNode operatorNode = NULL;
    
    if(hystNode->Has("elecPolarization")){
      operatorNode = hystNode->Get("elecPolarization");
      ReadHystOperator(material, operatorNode, false);
    } else if(hystNode->Has("magPolarization")){
      operatorNode = hystNode->Get("magPolarization");
      ReadHystOperator(material, operatorNode, false);
    } else {
      EXCEPTION("Either elecPolarization or magPolarization has to be defined for hysteresis node.")
    }
    
    PtrParamNode couplingNode = NULL;
    int usePolarization = 0;
    if(hystNode->Has("piezoCouplingAndStrains")){
      couplingNode = hystNode->Get("piezoCouplingAndStrains");
    } else if(hystNode->Has("magstrictCouplingAndStrains")){
      couplingNode = hystNode->Get("magstrictCouplingAndStrains");
    }
    
    int couplingDefined = 0;
    if(couplingNode != NULL){
      couplingDefined = 1;
      PtrParamNode strainModelingNode = couplingNode->Get("strainModeling");
      int irrStrainImplementation = -1;
      Double strainSat = 0.0;
      Double irrStrains_c1 = 0.0;
      Double irrStrains_c2 = 0.0;
      Double irrStrains_c3 = 0.0;
      int coefdim = 1;
      Matrix<Double> betaCoefs = Matrix<Double>(1,1);
      
      if(strainModelingNode != NULL){
        PtrParamNode innerNode = NULL;
        if(strainModelingNode->Has("muDatRelated")){
          irrStrainImplementation = 0;
          innerNode = strainModelingNode->Get("muDatRelated");
          strainSat = innerNode->Get("strainSat")->As<Double>();
          irrStrains_c1 = innerNode->Get("c1")->As<Double>();
          irrStrains_c2 = innerNode->Get("c2")->As<Double>();
          irrStrains_c3 = innerNode->Get("c3")->As<Double>();
        } else if(strainModelingNode->Has("higherOrderPolynomial")){
          irrStrainImplementation = 1;
          strainSat = innerNode->Get("strainSat")->As<Double>();
          coefdim = innerNode->Get("dim_betaCoefs")->As<Integer>();
          betaCoefs.Resize(1,coefdim);
          ParamTools::AsTensor<double>(innerNode->Get("betaCoefs"),1, coefdim, betaCoefs);
        }
      }
      material->SetScalar(strainSat, S_SATURATION, Global::REAL ); 
      material->SetScalar(irrStrainImplementation, HYST_IRRSTRAINS);
      material->SetScalar(irrStrains_c1, HYST_IRRSTRAIN_C1, Global::REAL);
      material->SetScalar(irrStrains_c2, HYST_IRRSTRAIN_C2, Global::REAL);
      material->SetScalar(irrStrains_c3, HYST_IRRSTRAIN_C3, Global::REAL);
      material->SetScalar(coefdim, DIM_BETA_COEFS);
      material->SetTensor(betaCoefs, HYST_BETA_COEFS, Global::REAL);
            
      PtrParamNode strainOperatorNode = couplingNode->Get("hystOperatorForStrains");
      if(strainOperatorNode->Has("usePolarization")){
        usePolarization = 1;
      } else if(strainOperatorNode->Has("separateHystOperator")){
        operatorNode = strainOperatorNode->Get("separateHystOperator");
        ReadHystOperator(material, operatorNode, true);
      }
      material->SetScalar(usePolarization, IRRSTRAIN_REUSE_P); 
      
      PtrParamNode smallSignalNode = couplingNode->Get("smallSignalForm");
      /*
       * -1 : coupled but without small signal tensors
       *  0 : coupled e-form/h-form (piezo/magstrict)
       *  1 : coupled d-form (piezo)
       *  2 : coupled g-form (magstrict)
       */
      int strainForm = -1;
      if(smallSignalNode != NULL){
        
        if(smallSignalNode->Has("noSmallSignalCoupling")){
          strainForm = -1;
        }
        if(smallSignalNode->Has("piezo_eform")){
          strainForm = 0;
        }
        if(smallSignalNode->Has("piezo_dform")){
          strainForm = 1;
        }
        if(smallSignalNode->Has("magstrict_hform")){
          strainForm = 0;
        }
        if(smallSignalNode->Has("magstrict_gform")){
          strainForm = 2;
        }
      }
      material->SetScalar(strainForm, HYST_STRAIN_FORM);
    }
    material->SetScalar(couplingDefined, HYST_COUPLING_DEFINED ); 
  }
  
  void XMLMaterialHandler::ReadHystOperator(BaseMaterial *material, PtrParamNode operatorNode, bool setStrains){
//    std::cout << "ReadHystOperator" << std::endl;
    PtrParamNode model;
    PtrParamNode pWeight = NULL;
    PtrParamNode pAnhyst = NULL;
    PtrParamNode pInversion = NULL;
    PtrParamNode pStrainForm = NULL;
    PtrParamNode initialState = NULL;
    PtrParamNode irrStrainNode = NULL;
    
    /*
     * Instead of setting each parameter twice in the form
     *    if(setStrains){
     *      material->SetScalar(model->Get("inputSat")->As<Double>(), X_SATURATION_STRAIN, Global::REAL ); 
     *    } else {
     *      material->SetScalar(model->Get("inputSat")->As<Double>(), X_SATURATION, Global::REAL ); 
     *    }
     * 
     * we define the enums for strains like X_SATURATION_STRAIN with a fixed offset towards the equivalent 
     * parameter for polarization X_SATURATION
     * > see Environment.hh
     */
    int enumOffset = 0;
    if(setStrains){
      enumOffset = 100;
    }
    
    if(operatorNode->Has("scalarPreisach")){
      model = operatorNode->Get("scalarPreisach");
      if(setStrains){
        material->SetScalar("scalarPreisach", HYST_MODEL_STRAIN);
        material->SetScalar("SCALAR", PREISACH_DIM_STRAIN);
      } else {
        material->SetScalar("scalarPreisach", HYST_MODEL);
        material->SetScalar("SCALAR", PREISACH_DIM);
      }
      
      // read input/output saturation of Preisach hysterese model
      material->SetScalar(model->Get("inputSat")->As<Double>(), MaterialType(X_SATURATION+enumOffset), Global::REAL ); 
      material->SetScalar(model->Get("outputSat")->As<Double>(), MaterialType(Y_SATURATION+enumOffset), Global::REAL ); 
      
      Matrix<Double> directionVector;
      if(model->Has("dirPolarization"))
      {
        //std::cout << "InitialState found" << std::endl;
        ParamTools::AsTensor<double>(model->Get("dirPolarization"),1, 3, directionVector);
        //std::cout << "IntialState: " << initialStateTensor.ToString() << std::endl;
      } else {
        //std::cout << "NO InitialState found" << std::endl;
        directionVector.Resize(1,3);
        directionVector.Init();
      }
      
      material->SetScalar( directionVector[0][0], MaterialType(P_DIRECTION_X+enumOffset), Global::REAL);
      material->SetScalar( directionVector[0][1], MaterialType(P_DIRECTION_Y+enumOffset), Global::REAL);
      material->SetScalar( directionVector[0][2], MaterialType(P_DIRECTION_Z+enumOffset), Global::REAL);
      
      if(model->Has("weights")){
        pWeight = model->Get("weights");
        // Read in weights separately below as they require the same steps for VectorHysteresis, too
      }
      
      if(model->Has("initialState")){
        initialState = model->Get("initialState");
      }
    }
    else if(operatorNode->Has("vectorPreisach_Sutor")){
//      std::cout << "vectorPreisach_Sutor" << std::endl;
      model = operatorNode->Get("vectorPreisach_Sutor");
      
      material->SetScalar("vectorPreisach_Sutor", MaterialType(HYST_MODEL+enumOffset));
      material->SetScalar("VECTOR", MaterialType(PREISACH_DIM+enumOffset));
      
      // read input/output saturation of Preisach hysterese model
      material->SetScalar(model->Get("inputSat")->As<Double>(), MaterialType(X_SATURATION+enumOffset), Global::REAL ); 
      material->SetScalar(model->Get("outputSat")->As<Double>(), MaterialType(Y_SATURATION+enumOffset), Global::REAL ); 
      
      /*
       * new numbering: 1 -> classical vector model (sutor2012)
       *                2 -> revised model (sutor2015)
       *                10 -> classical vector model, matrix based
       *                20 -> revised model, matrix based
       */
      int evalVersion = 2;
      if(model->Has("evalVersion")){
//        std::cout << "evalVersion" << std::endl;
        if(model->Get("evalVersion")->Has("Classical_Version_2012__list_implementation")){
          evalVersion = 1;
        } 
        if(model->Get("evalVersion")->Has("Revised_Version_2015__list_implementation")){
          evalVersion = 2;
        }
        if(model->Get("evalVersion")->Has("Classical_Version_2012__matrix_implementation")){
          evalVersion = 10;
        }
        if(model->Get("evalVersion")->Has("Revised_Version_2015__matrix_implementation")){
          evalVersion = 20;
        }
      }
      
      material->SetScalar(evalVersion, MaterialType(EVAL_VERSION+enumOffset));
      
      Double rotResistance = 1.0;
      Double angularDistance = 0.0;
      bool enforceSatOutputAtSatInput = true;
      if(model->Has("rotResistance")){
//        std::cout << "rotResistance" << std::endl;
        rotResistance = model->Get("rotResistance")->As<Double>();
      } 
      
      material->SetScalar(rotResistance, MaterialType(ROT_RESISTANCE+enumOffset), Global::REAL);
      
      if(model->Has("angularDistance")){
//                std::cout << "angularDistance" << std::endl;
        angularDistance = model->Get("angularDistance")->As<Double>();
      } 
      
      material->SetScalar(angularDistance, MaterialType(ANG_DISTANCE+enumOffset), Global::REAL);
      
      if(model->Has("enforceSatOutputAtSatInput")){
//        std::cout << "enforceSatOutputAtSatInput" << std::endl;
        enforceSatOutputAtSatInput = model->Get("enforceSatOutputAtSatInput")->As<bool>();
      } 
      
      if(enforceSatOutputAtSatInput){
        material->SetScalar(1, MaterialType(SCALETOSAT+enumOffset));
      } else {
        material->SetScalar(0, MaterialType(SCALETOSAT+enumOffset));
      }

      if(model->Has("weights")){
        pWeight = model->Get("weights");
        // Read in weights separately below as they require the same steps for VectorHysteresis, too
      }
      
      if(model->Has("hystInversion")){
        pInversion = model->Get("hystInversion");
      }
      
      if(model->Has("initialState")){
        initialState = model->Get("initialState");
      }
      
      PtrParamNode debugNode = NULL;
      Double angClip = 0.0;
      Double angRes = 1e-16;
      Double ampRes = 1e-16;
      int printOut = 0;
      int bmpRes = 1000;
      
      if(model->Has("debuggingParameter")){
        debugNode = model->Get("debuggingParameter");
        if(debugNode->Has("angularClipping")){
          angClip = debugNode->Get("angularClipping")->As<Double>();
        } 
        
        if(debugNode->Has("angularResolution")){
          debugNode->Get("angularResolution")->As<Double>();
        } 
        
        if(debugNode->Has("amplitudeResolution")){
          debugNode->Get("amplitudeResolution")->As<Double>();
        } 
        
        /*
         * if printOut > 0, the overlaid rotation and switching state of each printOut timestep will be
         * written to a bmp file of resolution bmpResolution
         */
        if(debugNode->Has("printOut")){
          debugNode->Get("printOut")->As<Integer>();
        } 
        
        if(debugNode->Has("bmpResolution")){
          debugNode->Get("bmpResolution")->As<Integer>();
        } 
      }
      
      material->SetScalar(angClip, MaterialType(ANG_CLIPPING+enumOffset), Global::REAL);
      material->SetScalar(angRes, MaterialType(ANG_RESOLUTION+enumOffset), Global::REAL);
      material->SetScalar(ampRes, MaterialType(AMP_RESOLUTION+enumOffset), Global::REAL);
      // extra parameter that will not be set for strains
      if(!setStrains){
        material->SetScalar(printOut, PRINT_PREISACH);
        material->SetScalar(bmpRes, PRINT_PREISACH_RESOLUTION);
      }
    }
    
    else if(operatorNode->Has("vectorPreisach_Mayergoyz")){
      model = operatorNode->Get("vectorPreisach_Mayergoyz");
      
      material->SetScalar("vectorPreisach_Mayergoyz", MaterialType(HYST_MODEL+enumOffset));
      material->SetScalar("VECTOR", MaterialType(PREISACH_DIM+enumOffset));
      
      PtrParamNode innerModel = NULL;
      PtrParamNode singleModel = NULL;
      if(model->Has("isotropic")){
        // read isotropic Mayergoyz vector model
        innerModel = model->Get("isotropic");
        
        material->SetScalar(1, MaterialType(PREISACH_MAYERGOYZ_ISOTROPIC+enumOffset) );
        
        
        int numDir = 11;
        if(innerModel->Has("numDirections")){
          numDir = innerModel->Get("numDirections")->As<Integer>();
        } 
        
        material->SetScalar(numDir, MaterialType(PREISACH_MAYERGOYZ_NUM_DIR+enumOffset) );
        
        if(innerModel->Has("ScalarModel")){
          singleModel = innerModel->Get("ScalarModel");
        } else {
          EXCEPTION("Single scalar Preisach model required for isotropic vector model");
        }
        
        // read input/output saturation of Preisach hysterese model
        material->SetScalar(model->Get("inputSat")->As<Double>(), MaterialType(X_SATURATION+enumOffset), Global::REAL ); 
        material->SetScalar(model->Get("outputSat")->As<Double>(), MaterialType(Y_SATURATION+enumOffset), Global::REAL ); 
        
        if(singleModel->Has("weights")){
          pWeight = singleModel->Get("weights");
          // Read in weights separately below as they require the same steps for VectorHysteresis, too
        }
        
        int adaptedToVectorCase = 0;
        if(singleModel->Has("weightsAdaptedToMayergoyzVectorModel")){
          bool adapted = singleModel->Get("weightsAdaptedToMayergoyzVectorModel")->As<bool>();
          if(adapted){
            adaptedToVectorCase = 1;
          }
        }
        material->SetScalar(adaptedToVectorCase, MaterialType(PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR+enumOffset));
        
      } else if(model->Has("anIsotropic")){
        EXCEPTION("Anisotropic Mayergoyz vector hysteresis model not yet supported");
      }
      
      int clipOutput = 2;
      if(model->Has("clipOutputToSat")){
        if(model->Get("clipOutputToSat")->Has("noClipping")){
          clipOutput = 0;
        }
        if(model->Get("clipOutputToSat")->Has("clipAmplitude")){
          clipOutput = 1;
        }
        if(model->Get("clipOutputToSat")->Has("clipComponentParallelToInput")){
          // leads to most reasonable results
          clipOutput = 2;
        }
      }
      material->SetScalar(clipOutput, MaterialType(PREISACH_MAYERGOYZ_CLIPOUTPUT+enumOffset));
      
      if(model->Has("hystInversion")){
        pInversion = model->Get("hystInversion");
      }
      
      if(model->Has("initialState")){
        initialState = model->Get("initialState");
      }      
    }
    else {
      material->SetScalar("none", MaterialType(HYST_MODEL+enumOffset));  
      
      return;
    }
    
    if(pWeight != NULL){
      // Read in weights
      int dim = -1;
      if(pWeight->Has("dim_weights")) dim = pWeight->Get("dim_weights")->As<Integer>();
      
      material->SetScalar( dim, MaterialType(PREISACH_WEIGHTS_DIM+enumOffset));
      
      int weightType = 0; // 0 = const, 1 = muDat, 2 = muDatExtended, 3 = givenTensor
      PtrParamNode pWeightInner;
      if(pWeight->Has("weightType")){
        pWeightInner = pWeight->Get("weightType");
      }
      
      if(pWeightInner->Has("const")){
        weightType = 0;
        Double constValue = pWeightInner->Get("const")->As<Double>();
        
        material->SetScalar(constValue, MaterialType(PREISACH_WEIGHTS_CONSTVALUE+enumOffset), Global::REAL );
        material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));
        
      } else if(pWeightInner->Has("muDat")){
        weightType = 1;
        PtrParamNode muDat = pWeightInner->Get("muDat");
        Double A = muDat->Get("A")->As<Double>();
        Double h = muDat->Get("h")->As<Double>();
        Double sigma = muDat->Get("sigma")->As<Double>();
        Double eta = muDat->Get("eta")->As<Double>();
        
        material->SetScalar(A, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL );
        material->SetScalar(h, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL );
        material->SetScalar(sigma, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL );
        material->SetScalar(eta, MaterialType(PREISACH_WEIGHTS_MUDAT_ETA+enumOffset), Global::REAL );
        material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));
        
      } else if(pWeightInner->Has("muDatExtended")){
        weightType = 2;
        PtrParamNode muDatExt = pWeightInner->Get("muDatExtended");
        Double A = muDatExt->Get("A")->As<Double>();
        Double h1 = muDatExt->Get("h1")->As<Double>();
        Double h2 = muDatExt->Get("h2")->As<Double>();
        Double sigma1 = muDatExt->Get("sigma1")->As<Double>();
        Double sigma2 = muDatExt->Get("sigma2")->As<Double>();
        Double eta = muDatExt->Get("eta")->As<Double>();        
        
        material->SetScalar(A, MaterialType(PREISACH_WEIGHTS_MUDATEXT_A+enumOffset), Global::REAL );
        material->SetScalar(h1, MaterialType(PREISACH_WEIGHTS_MUDATEXT_H1+enumOffset), Global::REAL );
        material->SetScalar(h2, MaterialType(PREISACH_WEIGHTS_MUDATEXT_H2+enumOffset), Global::REAL );
        material->SetScalar(sigma1, MaterialType(PREISACH_WEIGHTS_MUDATEXT_SIGMA1+enumOffset), Global::REAL );
        material->SetScalar(sigma2, MaterialType(PREISACH_WEIGHTS_MUDATEXT_SIGMA2+enumOffset), Global::REAL );
        material->SetScalar(eta, MaterialType(PREISACH_WEIGHTS_MUDATEXT_ETA+enumOffset), Global::REAL );
        material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));
        
      } else if(pWeightInner->Has("weightTensor")){
        weightType = 3;
        Matrix<Double> preisachWeightTensor(dim,dim);
        ParamTools::AsTensor<double>(pWeightInner->Get("weightTensor"), dim, dim, preisachWeightTensor);
        
        material->SetTensor( preisachWeightTensor, MaterialType(PREISACH_WEIGHTS_TENSOR+enumOffset), Global::REAL);
        material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));
        
      } else {
        EXCEPTION("No valid Preisach weights found");
      }
      if(pWeightInner->Has("anhystereticParameter")){
        pAnhyst = pWeightInner->Get("anhystereticParameter");
      } 
    } else {
      EXCEPTION("No valid Preisach weights found");
    }
    
    Double a,b,c;
    bool onlyanhyst;
    if(pAnhyst != NULL){
      if(pAnhyst->Has("a")){
        a = pAnhyst->Get("a")->As<Double>();
      } else {
        a = 0.0;
      }
      if(pAnhyst->Has("b")){
        b = pAnhyst->Get("b")->As<Double>();
      } else {
        b = 0.0;
      }
      if(pAnhyst->Has("c")){
        c = pAnhyst->Get("c")->As<Double>();
      } else {
        c = 0.0;
      }
      if(pAnhyst->Has("onlyAnhyst")){
        onlyanhyst = pAnhyst->Get("onlyAnhyst")->As<bool>();
      } else {
        onlyanhyst = false;
      }
    } else {
      a = 0;
      b = 0;
      c = 0;
      onlyanhyst = false;
    }
    
    material->SetScalar(a, MaterialType(PREISACH_WEIGHTS_ANHYST_A+enumOffset), Global::REAL);
    material->SetScalar(b, MaterialType(PREISACH_WEIGHTS_ANHYST_B+enumOffset), Global::REAL);
    material->SetScalar(c, MaterialType(PREISACH_WEIGHTS_ANHYST_C+enumOffset), Global::REAL);
    if(onlyanhyst){
      material->SetScalar(1, MaterialType(PREISACH_WEIGHTS_ANHYST_ONLY+enumOffset));
    } else {
      material->SetScalar(0, MaterialType(PREISACH_WEIGHTS_ANHYST_ONLY+enumOffset));
    }
    
    Matrix<Double> initialStateTensor = Matrix<Double>(1,3);
    initialStateTensor.Init();
    bool prescribeOutput = false;
    bool scaleBySaturation = false;
    
    if(initialState != NULL){
      PtrParamNode InOutState = NULL;
      if(initialState->Has("initialOutput")){
        InOutState = initialState->Get("initialOutput");
        prescribeOutput = true;
      } else if(initialState->Has("initialInput")){
        InOutState = initialState->Get("initialInput");
        prescribeOutput = false;
      } else {
        EXCEPTION("Either input or output must be given at this point");
      }
      if(InOutState != NULL){
        if(InOutState->Has("Vector")){
          //std::cout << "InitialState found" << std::endl;
          ParamTools::AsTensor<double>(InOutState->Get("Vector"),1, 3, initialStateTensor);
          //std::cout << "IntialState: " << initialStateTensor.ToString() << std::endl;
        }
        if(InOutState->Has("scaleVectorBySaturation")){
          scaleBySaturation = InOutState->Get("scaleVectorBySaturation")->As<bool>();
        } 
      }
    }
    
    // strain hyst operator will automatically get the same initial input as the polarization
    // if output is prescribed for polarization, the resulting input will be passed to strain operator
    if(!setStrains){
      material->SetScalar( initialStateTensor[0][0], INITIAL_STATE_X, Global::REAL);
      material->SetScalar( initialStateTensor[0][1], INITIAL_STATE_Y, Global::REAL);
      material->SetScalar( initialStateTensor[0][2], INITIAL_STATE_Z, Global::REAL);
      if(scaleBySaturation){
        material->SetScalar(1, PREISACH_SCALEINITIALSTATE);
      } else {
        material->SetScalar(0, PREISACH_SCALEINITIALSTATE);
      }
      
      if(prescribeOutput){
        material->SetScalar(1, PREISACH_PRESCRIBEOUTPUT);
      } else {
        material->SetScalar(0, PREISACH_PRESCRIBEOUTPUT);
      }
    }
    int inversionMethod = 0;
    int maxNumIts = 35;
    int maxNumLSIts = 100;
    double tolH = 1e-12;
    double tolB = 1e-12;
    double jacRes = 1e-12;
    double alphaLSStart = 0.25;
    double alphaLSMin = 1.0/512.0;
    double alphaLSMax = 8192.0;
    bool setInversion = false;
    if(pInversion != NULL){
      setInversion = true;
      if(pInversion->Has("InversionMethod"))
      {
        if(pInversion->Get("InversionMethod")->Has("LevenbergMarquardt")){
          inversionMethod = 0;
        } else if(pInversion->Get("InversionMethod")->Has("Newton")){
          inversionMethod = 1;
        } else if(pInversion->Get("InversionMethod")->Has("JacobianFreeNewtonKrylov")){
          inversionMethod = 2;
        }
      }

      if(pInversion->Has("maxNumberOuterIterations"))
      {
        maxNumIts = pInversion->Get("maxNumberOuterIterations")->As<Integer>();
      }
      
      if(pInversion->Has("residualTolH"))
      {
        tolH = pInversion->Get("residualTolH")->As<double>();
      }
      
      if(pInversion->Has("residualTolB"))
      {
        tolB = pInversion->Get("residualTolB")->As<double>();
      }
      
      if(pInversion->Has("jacobiResolution"))
      {
        jacRes = pInversion->Get("jacobiResolution")->As<double>();
      }
      
      if(pInversion->Has("maxNumberLinesearchIterations"))
      {
        maxNumLSIts = pInversion->Get("maxNumberLinesearchIterations")->As<Integer>();
      }
            
      if(pInversion->Has("alphaRegStart"))
      {
        alphaLSStart = pInversion->Get("alphaRegStart")->As<double>();
      }
      
      if(pInversion->Has("alphaRegMin"))
      {
        alphaLSMin = pInversion->Get("alphaRegMin")->As<double>();
      }
      
      if(pInversion->Has("alphaRegMax"))
      {
        alphaLSMax = pInversion->Get("alphaRegMax")->As<double>();
      }
    }
    //Hyst operator for strains does not use inversion! Only forward mode used
    if((setInversion==true) && (setStrains==false)){
      material->SetScalar(maxNumIts, MAX_NUM_IT_HYST_INV);
      material->SetScalar(tolH, RES_TOL_H_HYST_INV, Global::REAL);
      material->SetScalar(tolB, RES_TOL_B_HYST_INV, Global::REAL);
      material->SetScalar(jacRes, JAC_RESOLUTION_HYST_INV, Global::REAL);
      material->SetScalar(inversionMethod, VEC_HYST_INV_METHOD);
      material->SetScalar(maxNumLSIts, MAX_NUM_LS_IT_HYST_INV);
      material->SetScalar(alphaLSStart, ALPHA_LS_HYST_INV, Global::REAL);
      material->SetScalar(alphaLSMin, ALPHA_LS_MIN_HYST_INV, Global::REAL);
      material->SetScalar(alphaLSMax, ALPHA_LS_MAX_HYST_INV, Global::REAL);
    }
//    std::cout << "ReadHystOperator - done" << std::endl;
  }
  
} // end of namespace
