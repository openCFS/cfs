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
    
    //read Preisach hysterese model
    if(elec->Has("hystModel"))
    {
      if(elec->Get("hystModel")->Has("preisach"))
      {
        PtrParamNode p = elec->Get("hystModel")->Get("preisach");
        
        // force name
        material->SetScalar("preisach", HYST_MODEL);
        
        // read E saturation of Preisach hysterese model
        if(p->Has("eSat"))
          material->SetScalar(p->Get("eSat")->As<Double>(), X_SATURATION, Global::REAL ); 
        
        // read P saturation of Preisach hysterese model
        if(p->Has("pSat"))
          material->SetScalar(p->Get("pSat")->As<Double>(), Y_SATURATION, Global::REAL ); 
        
        // never used
        //        // read P saturation of Preisach hysterese model
        //        if(p->Has("Pr"))
        //          material->SetScalar(p->Get("Pr")->As<Double>(), Y_REMANENCE, Global::REAL );
        
        // read direction of polarization
        if(p->Has("dirP"))
        {
          int dir = p->Get("dirP")->As<Integer>();
          
          if(dir == 1) material->SetScalar("X", P_DIRECTION );
          if(dir == 2) material->SetScalar("Y", P_DIRECTION );
          if(dir == 3) material->SetScalar("Z", P_DIRECTION );
          
          if(dir != 1 && dir != 2 && dir != 3)
            EXCEPTION(dir << " is valid coordinate direction for electric preisach "
                    << " hysteresis model polarization");
        }
        // not needed anymore -> preisach is always scalar; vector model has its own name now
        material->SetScalar("SCALAR", PREISACH_DIM);
        
        // read weight dimension of Preisach hysterese model for weights
        int dim = -1;
        if(p->Has("dim_weights")) dim = p->Get("dim_weights")->As<Integer>();
        
        // read real permittivity tensor    
        if(p->Has("weights"))
        {
          Matrix<Double> preisachWeightTensor(dim,dim);
          ParamTools::AsTensor<double>(p->Get("weights"), dim, dim, preisachWeightTensor);
          material->SetTensor( preisachWeightTensor, PREISACH_WEIGHTS, Global::REAL);
        }
        
        if(p->Has("strainForm")){
          material->SetScalar(p->Get("strainForm")->As<Integer>(), HYST_STRAIN_FORM);
        } else {
          material->SetScalar(0, HYST_STRAIN_FORM);
        }
        
        int coefdim = -1;
        if(p->Has("dim_betaCoefs")) coefdim = p->Get("dim_betaCoefs")->As<Integer>();
        material->SetScalar(coefdim, DIM_BETA_COEFS);
        
        Matrix<Double> betaCoef;
        if(p->Has("betaCoefs")&&(coefdim != -1))
        {
          //std::cout << "beta coefs found" << std::endl;
          ParamTools::AsTensor<double>(p->Get("betaCoefs"),1, coefdim, betaCoef);
          //std::cout << "beta coefs: " << betaCoef.ToString() << std::endl;
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        } else {
          //std::cout << "NO betaCoef found" << std::endl;
          betaCoef = Matrix<Double>(1,1);
          betaCoef.Init();
          //std::cout << "Beta coefs set to: " << betaCoef.ToString() << std::endl;
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        }
                
      }
      else if (elec->Get("hystModel")->Has("vectorPreisach"))
      {
        PtrParamNode p = elec->Get("hystModel")->Get("vectorPreisach");
        
        // force name
        material->SetScalar("vectorPreisach", HYST_MODEL);
        
        // read E saturation of Preisach hysterese model
        if(p->Has("eSat"))
          material->SetScalar(p->Get("eSat")->As<Double>(), X_SATURATION, Global::REAL );
        // read P saturation of Preisach hysterese model
        if(p->Has("pSat"))
          material->SetScalar(p->Get("pSat")->As<Double>(), Y_SATURATION, Global::REAL );
        
        /*
         * new numbering: 1 -> classical vector model (sutor2012)
         *                2 -> revised model (sutor2015)
         *                10 -> classical vector model, matrix based
         *                20 -> revised model, matrix based
         */
        if(p->Has("evalVersion")){
          material->SetScalar(p->Get("evalVersion")->As<Integer>(), EVAL_VERSION);
        } else {
          material->SetScalar(2, EVAL_VERSION);
        }
        
        //not needed anymore as vectorPreisach always is vectorial
        material->SetScalar("VECTOR", PREISACH_DIM);
        
        if(p->Has("rotResistance")){
          material->SetScalar(p->Get("rotResistance")->As<Double>(), ROT_RESISTANCE, Global::REAL);
        } else {
          material->SetScalar(1.0, ROT_RESISTANCE, Global::REAL);
        }
        
        if(p->Has("angularDistance")){
          material->SetScalar(p->Get("angularDistance")->As<Double>(), ANG_DISTANCE, Global::REAL);
        } else {
          material->SetScalar(0.0, ANG_DISTANCE, Global::REAL);
        }
        
        if(p->Has("angularResolution")){
          material->SetScalar(p->Get("angularResolution")->As<Double>(), ANG_RESOLUTION, Global::REAL);
        } else {
          material->SetScalar(0.0001, ANG_RESOLUTION, Global::REAL);
        }
        
        if(p->Has("angularClipping")){
          material->SetScalar(p->Get("angularClipping")->As<Double>(), ANG_CLIPPING, Global::REAL);
        } else {
          material->SetScalar(0.0001, ANG_CLIPPING, Global::REAL);
        }
        
        if(p->Has("amplitudeResolution")){
          material->SetScalar(p->Get("amplitudeResolution")->As<Double>(), AMP_RESOLUTION, Global::REAL);
        } else {
          material->SetScalar(0.0, AMP_RESOLUTION, Global::REAL);
        }
        
        // read weight dimension of Preisach hysterese model for weights
        int dim = -1;
        if(p->Has("dim_weights")) dim = p->Get("dim_weights")->As<Integer>();
        
        // read real permittivity tensor
        if(p->Has("weights"))
        {
          Matrix<Double> preisachWeightTensor(dim,dim);
          ParamTools::AsTensor<double>(p->Get("weights"), dim, dim, preisachWeightTensor);
          material->SetTensor( preisachWeightTensor, PREISACH_WEIGHTS, Global::REAL);
        }
        
        Matrix<Double> initialStateTensor;
        if(p->Has("initialState"))
        {
          //std::cout << "InitialState found" << std::endl;
          ParamTools::AsTensor<double>(p->Get("initialState"),1, 3, initialStateTensor);
          //std::cout << "IntialState: " << initialStateTensor.ToString() << std::endl;
        } else {
          //std::cout << "NO InitialState found" << std::endl;
          initialStateTensor.Resize(1,3);
          initialStateTensor.Init();
        }
        
        material->SetScalar( initialStateTensor[0][0], INITIAL_STATE_X, Global::REAL);
        material->SetScalar( initialStateTensor[0][1], INITIAL_STATE_Y, Global::REAL);
        material->SetScalar( initialStateTensor[0][2], INITIAL_STATE_Z, Global::REAL);
        
        /*
         * if printOut > 0, the overlaid rotation and switching state of each printOut timestep will be
         * written to a bmp file of resolution bmpResolution
         */
        if(p->Has("printOut")){
          material->SetScalar(p->Get("printOut")->As<Integer>(), PRINT_PREISACH);
        } else {
          material->SetScalar(0, PRINT_PREISACH);
        }
        
        if(p->Has("bmpResolution")){
          material->SetScalar(p->Get("bmpResolution")->As<Integer>(), PRINT_PREISACH_RESOLUTION);
        } else {
          material->SetScalar(1000, PRINT_PREISACH_RESOLUTION);
        }
        
        if(p->Has("strainForm")){
          material->SetScalar(p->Get("strainForm")->As<Integer>(), HYST_STRAIN_FORM);
        } else {
          material->SetScalar(0, HYST_STRAIN_FORM);
        }
        
        int coefdim = -1;
        if(p->Has("dim_betaCoefs")) coefdim = p->Get("dim_betaCoefs")->As<Integer>();
        material->SetScalar(coefdim, DIM_BETA_COEFS);
        
        Matrix<Double> betaCoef;
        if(p->Has("betaCoefs")&&(coefdim != -1))
        {
          //std::cout << "beta coefs found" << std::endl;
          ParamTools::AsTensor<double>(p->Get("betaCoefs"),1, coefdim, betaCoef);
          //std::cout << "beta coefs: " << betaCoef.ToString() << std::endl;
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        } else {
          //std::cout << "NO betaCoef found" << std::endl;
          betaCoef = Matrix<Double>(1,1);
          betaCoef.Init();
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        }
        
      }
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
    
    //read Preisach hysterese model
    if(mag->Has("hystModel"))
    {
      if(mag->Get("hystModel")->Has("preisach"))
      {
        PtrParamNode p = mag->Get("hystModel")->Get("preisach");
        
        // force name
        material->SetScalar("preisach", HYST_MODEL);
        
        // read E saturation of Preisach hysterese model
        if(p->Has("hSat"))
          material->SetScalar(p->Get("hSat")->As<Double>(), X_SATURATION, Global::REAL );
        // read P saturation of Preisach hysterese model
        if(p->Has("jSat"))
          material->SetScalar(p->Get("jSat")->As<Double>(), Y_SATURATION, Global::REAL );
        
        // never used
        //        // read M saturation of Preisach hysterese model
        //        if(p->Has("Mr"))
        //          material->SetScalar(p->Get("Mr")->As<Double>(), Y_REMANENCE, Global::REAL );
        
        // read direction of polarization
        if(p->Has("dirJ"))
        {
          int dir = p->Get("dirJ")->As<Integer>();
          
          if(dir == 1) material->SetScalar("X", P_DIRECTION );
          if(dir == 2) material->SetScalar("Y", P_DIRECTION );
          if(dir == 3) material->SetScalar("Z", P_DIRECTION );
          
          if(dir != 1 && dir != 2 && dir != 3)
            EXCEPTION(dir << " is valid coordinate direction for electric preisach "
                    << " hysteresis model polarization");
        }
        // not needed anymore -> preisach is always scalar; vector model has its own name now
        material->SetScalar("SCALAR", PREISACH_DIM);
        
        // read weight dimension of Preisach hysterese model for weights
        int dim = -1;
        if(p->Has("dim_weights")) dim = p->Get("dim_weights")->As<Integer>();
        
        // read real permittivity tensor
        if(p->Has("weights"))
        {
          Matrix<Double> preisachWeightTensor(dim,dim);
          ParamTools::AsTensor<double>(p->Get("weights"), dim, dim, preisachWeightTensor);
          material->SetTensor( preisachWeightTensor, PREISACH_WEIGHTS, Global::REAL);
        }
        
        if(p->Has("strainForm")){
          material->SetScalar(p->Get("strainForm")->As<Integer>(), HYST_STRAIN_FORM);
        } else {
          material->SetScalar(0, HYST_STRAIN_FORM);
        }
        
        int coefdim = -1;
        if(p->Has("dim_betaCoefs")) coefdim = p->Get("dim_betaCoefs")->As<Integer>();
        material->SetScalar(coefdim, DIM_BETA_COEFS);
        
        Matrix<Double> betaCoef;
        if(p->Has("betaCoefs")&&(coefdim != -1))
        {
          //std::cout << "beta coefs found" << std::endl;
          ParamTools::AsTensor<double>(p->Get("betaCoefs"),1, coefdim, betaCoef);
          //std::cout << "beta coefs: " << betaCoef.ToString() << std::endl;
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        } else {
          //std::cout << "NO betaCoef found" << std::endl;
          betaCoef = Matrix<Double>(1,1);
          betaCoef.Init();
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        }
        
      }
      
      else if (mag->Get("hystModel")->Has("vectorPreisach"))
      {
        PtrParamNode p = mag->Get("hystModel")->Get("vectorPreisach");
        
        // force name
        material->SetScalar("vectorPreisach", HYST_MODEL);
        
        // read E saturation of Preisach hysterese model
        if(p->Has("hSat"))
          material->SetScalar(p->Get("hSat")->As<Double>(), X_SATURATION, Global::REAL );
        // read P saturation of Preisach hysterese model
        if(p->Has("jSat"))
          material->SetScalar(p->Get("jSat")->As<Double>(), Y_SATURATION, Global::REAL );
        
        /*
         * new numbering: 1 -> classical vector model (sutor2012)
         *                2 -> revised model (sutor2015)
         *                10 -> classical vector model, matrix based
         *                20 -> revised model, matrix based
         */
        if(p->Has("evalVersion")){
          material->SetScalar(p->Get("evalVersion")->As<Integer>(), EVAL_VERSION);
        } else {
          material->SetScalar(2, EVAL_VERSION);
        }
        
        //not needed anymore as vectorPreisach always is vectorial
        material->SetScalar("VECTOR", PREISACH_DIM);
        
        if(p->Has("rotResistance")){
          material->SetScalar(p->Get("rotResistance")->As<Double>(), ROT_RESISTANCE, Global::REAL);
        } else {
          material->SetScalar(1.0, ROT_RESISTANCE, Global::REAL);
        }
        
        if(p->Has("angularDistance")){
          material->SetScalar(p->Get("angularDistance")->As<Double>(), ANG_DISTANCE, Global::REAL);
        } else {
          material->SetScalar(0.0, ANG_DISTANCE, Global::REAL);
        }
        
        if(p->Has("angularResolution")){
          material->SetScalar(p->Get("angularResolution")->As<Double>(), ANG_RESOLUTION, Global::REAL);
        } else {
          material->SetScalar(0.0, ANG_RESOLUTION, Global::REAL);
        }
        
        if(p->Has("angularClipping")){
          material->SetScalar(p->Get("angularClipping")->As<Double>(), ANG_CLIPPING, Global::REAL);
        } else {
          material->SetScalar(0.0001, ANG_CLIPPING, Global::REAL);
        }
        
        if(p->Has("amplitudeResolution")){
          material->SetScalar(p->Get("amplitudeResolution")->As<Double>(), AMP_RESOLUTION, Global::REAL);
        } else {
          material->SetScalar(0.0, AMP_RESOLUTION, Global::REAL);
        }
        
        // read weight dimension of Preisach hysterese model for weights
        int dim = -1;
        if(p->Has("dim_weights")) dim = p->Get("dim_weights")->As<Integer>();
        
        // read real permittivity tensor
        if(p->Has("weights"))
        {
          Matrix<Double> preisachWeightTensor(dim,dim);
          ParamTools::AsTensor<double>(p->Get("weights"), dim, dim, preisachWeightTensor);
          material->SetTensor( preisachWeightTensor, PREISACH_WEIGHTS, Global::REAL);
        }
        
        Matrix<Double> initialStateTensor;
        if(p->Has("initialState"))
        {
          //std::cout << "InitialState found" << std::endl;
          ParamTools::AsTensor<double>(p->Get("initialState"),1, 3, initialStateTensor);
          //std::cout << "IntialState: " << initialStateTensor.ToString() << std::endl;
        } else {
          //std::cout << "NO InitialState found" << std::endl;
          initialStateTensor.Resize(1,3);
          initialStateTensor.Init();
        }
        
        material->SetScalar( initialStateTensor[0][0], INITIAL_STATE_X, Global::REAL);
        material->SetScalar( initialStateTensor[0][1], INITIAL_STATE_Y, Global::REAL);
        material->SetScalar( initialStateTensor[0][2], INITIAL_STATE_Z, Global::REAL);
        
        /*
         * if printOut > 0, the overlaid rotation and switching state of each printOut timestep will be
         * written to a bmp file of resolution bmpResolution
         */
        if(p->Has("printOut")){
          material->SetScalar(p->Get("printOut")->As<Integer>(), PRINT_PREISACH);
        } else {
          material->SetScalar(0, PRINT_PREISACH);
        }
        
        if(p->Has("bmpResolution")){
          material->SetScalar(p->Get("bmpResolution")->As<Integer>(), PRINT_PREISACH_RESOLUTION);
        } else {
          material->SetScalar(1000, PRINT_PREISACH_RESOLUTION);
        }
        
        if(p->Has("strainForm")){
          material->SetScalar(p->Get("strainForm")->As<Integer>(), HYST_STRAIN_FORM);
        } else {
          material->SetScalar(0, HYST_STRAIN_FORM);
        }
        
        int coefdim = -1;
        if(p->Has("dim_betaCoefs")) coefdim = p->Get("dim_betaCoefs")->As<Integer>();
        material->SetScalar(coefdim, DIM_BETA_COEFS);
        
        Matrix<Double> betaCoef;
        if(p->Has("betaCoefs")&&(coefdim != -1))
        {
          //std::cout << "beta coefs found" << std::endl;
          ParamTools::AsTensor<double>(p->Get("betaCoefs"),1, coefdim, betaCoef);
          //std::cout << "beta coefs: " << betaCoef.ToString() << std::endl;
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        } else {
          //std::cout << "NO betaCoef found" << std::endl;
          betaCoef = Matrix<Double>(1,1);
          betaCoef.Init();
          material->SetTensor(betaCoef, HYST_BETA_COEFS, Global::REAL);
        }
        
      }
      
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

    // read reference temperature
    if(therm->Has("refTemperature")) {
          PtrCoefFct refTemp =
                  CoefFunction::Generate(mp_, Global::REAL,
                		  therm->Get("refTemperature")->As<std::string>() );
          material->SetCoefFct( REF_TEMPERATURE, refTemp );
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

    // read bulk viscosity
    if(flow->Has("bulkViscosity")) {
          PtrCoefFct bulkVisc =
                  CoefFunction::Generate(mp_, Global::REAL,
                  flow->Get("bulkViscosity")->As<std::string>() );
          material->SetCoefFct( BULK_VISCOSITY, bulkVisc );
    }

    // read adiabatic exponent
    if(flow->Has("adiabaticExponent")) {
    	PtrCoefFct exp =
    		CoefFunction::Generate(mp_, Global::REAL,
    					flow->Get("adiabaticExponent")->As<std::string>() );
    	material->SetCoefFct( ADIABATIC_EXPONENT, exp );
    }

    // read reference pressure
    if(flow->Has("refPressure")) {
          PtrCoefFct refPres =
                  CoefFunction::Generate(mp_, Global::REAL,
                  flow->Get("refPressure")->As<std::string>() );
          material->SetCoefFct( REF_PRESSURE, refPres );
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
  
} // end of namespace
