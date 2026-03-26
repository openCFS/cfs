#include "XMLMaterialHandler.hh"

#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ProgramOptions.hh"

// header for logging
#include "DataInOut/Logging/LogConfigurator.hh"

// header for materials
#include "Materials/ElectroMagneticMaterial.hh"
#include "Materials/ElectroStaticMaterial.hh"
#include "Materials/ElecQuasistaticMaterial.hh"
#include "Materials/HeatMaterial.hh"
#include "Materials/AcousticMaterial.hh"
#include "Materials/MechanicMaterial.hh"
#include "Materials/SmoothMaterial.hh"
#include "Materials/PiezoMaterial.hh"
#include "Materials/FlowMaterial.hh"
#include "Materials/TestMaterial.hh"
#include "Materials/ElectricConductionMaterial.hh"
//#include "Materials/thermoelasticMaterial.hh"
//#include "Materials/pyroelectricMaterial.hh"
#include "Materials/MagStrictMaterial.hh"
#include "Utils/tools.hh"

// Note, that the methods ComputeIso/OrthoMechStiffnesTensor were commented out
// in revision 7562 and are not in the code -> check the repository!

// define shorthand notation
typedef BaseMaterial BM;

namespace CoupledField {
  
  DEFINE_LOG(xmlmathandler, "xmlmathandler")

  // Path to the material XML schema relative to share/xml
  const std::string XMLMaterialHandler::schemaFile_ = "/CFS-Material/CFS_Material.xsd";

  // XML Schema URL
  const std::string XMLMaterialHandler::schemaUrl_ = "http://www.cfs++.org/material";

  XMLMaterialHandler::XMLMaterialHandler()
  : MaterialHandler() {
  }
  
  XMLMaterialHandler::~XMLMaterialHandler()
  {
  }
  
  void XMLMaterialHandler::LoadFromFile( const std::string& fileName ) {
    
    this->fileName_ = fileName;
    
    // Create a ParamNode and parse the material file
    std::string schema = progOpts->GetSchemaPathStr() + schemaFile_;
    
    rootNode_ = XmlReader::ParseFile(fileName, schema, schemaUrl_);
  }
  
  void XMLMaterialHandler::LoadFromString( const std::string& str ) {
    // Create a ParamNode and parse the material file
    std::string schema = progOpts->GetSchemaPathStr() + schemaFile_;
    
    rootNode_ = XmlReader::ParseString(str, schema, schemaUrl_);
  }
  
  BaseMaterial * XMLMaterialHandler::LoadMaterial( const std::string &matName,
                                                   MaterialClass matClass )
  {
    BaseMaterial * material = NULL;
    
    std::string strMatClass;
    
    Enum2String(matClass,strMatClass);
    if (!rootNode_->HasByVal("material", "name", matName)) {
      EXCEPTION("Cannot find material '" << matName << "'");
    }
    
    // first get the material element:  <material name="iron">
    PtrParamNode pn;
    pn = rootNode_->GetByVal("material", "name", matName);
    
    if ( !pn ) {
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
      else if (matClass == SMOOTH) {
        material = new SmoothMaterial(mp, cs);
        ReadSmooth( material, pn );
      }
      else if ( matClass == ACOUSTIC ) {\
        material = new AcousticMaterial(mp, cs);
        ReadAcoustic( material, pn );
      }
      else if ( matClass == ELECTROMAGNETIC ) {
        material = new ElectroMagneticMaterial(mp, cs);
        ReadMagnetic( material, pn );
      }
      else if ( matClass == ELECTROMAGNETIC_DARWIN ) {
        material = new ElectroMagneticMaterial(mp, cs, true);
        ReadMagnetic( material, pn );
      }
      else if ( matClass == ELECTROSTATIC ) {
        material = new ElectroStaticMaterial(mp, cs);
        ReadElectrostatic( material, pn );
      }
      else if ( matClass == ELECQUASISTATIC ) {
        material = new ElecQuasistaticMaterial(mp, cs);
        ReadElecQuasistatic( material, pn );
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
        REFACTOR;
        //material = new ThermoelasticMaterial();
        //ReadThermoelastic( material, pn );
      }
      else if ( matClass == MAGNETOSTRICTIVE ) {
        material = new MagStrictMaterial(mp,cs);
        ReadMagStrict( material, pn );
      }
      else if ( matClass == ELECTRICCONDUCTION ) {
        material = new ElectricConductionMaterial(mp, cs);
        ReadElectricConduction( material, pn );
      }
      else {
        EXCEPTION( "material type" << matClass << " not defined" );
      }
      
      // Finalize setup of material
      material->Finalize();
      material->SetName(matName);

    }
    catch (Exception& ex ) {
      RETHROW_EXCEPTION(ex, "Could not load material '" << matName  
          << "' of class '" << MaterialClassEnum.ToString(matClass) << "'" );
    }

    return material;
  }
  
  //**********************************************************************
  //*************  READ PIEZO ********************************************
  //**********************************************************************
  
  void XMLMaterialHandler::ReadPiezo(BaseMaterial *material, PtrParamNode pn) 
  {
    PtrParamNode cpl = pn->Get("piezoCoupling", ParamNode::PASS);

    if (cpl) {
      // read piezo coupling tensor
      if (cpl->Has("linear")) {
        BM::SymmetryType symType = BM::NOSYMMETRY;

        if (cpl->Get("linear")->Has("tensor")) {
          PtrParamNode pcTensor = cpl->Get("linear")->Get("tensor");
          PtrCoefFct piezoCoef = ReadTensor(pcTensor, Global::COMPLEX );
          material->SetCoefFct(PIEZO_TENSOR, piezoCoef);
          symType = BM::GENERAL;
        }
        material->SetSymmetryType(PIEZO_TENSOR, symType);
      }

      // read nonlinearity of a coupling coefficient
      if (cpl->Has("nonlinear") && cpl->Get("nonlinear")->Has("isotropic")) {
        BaseMaterial::MatDescriptorNl nlInfo =
            ReadNonlinDescriptor(cpl->Get("nonlinear")->Get("isotropic"), material);
        material->SetNonLinMatIso(PIEZO_TENSOR, nlInfo);
      }
    }
    
    if (pn->Has("piezoMicroData")) {
      if (pn->Get("piezoMicroData")->Has("HuberFleck")) {
        PtrParamNode pcc = pn->Get("piezoMicroData")->Get("HuberFleck");
        
        // force name
        //        material->SetScalar("BelovKreher", PIEZO_MICRO_MODEL);
        
        // read remanent polarisation
        if (pcc->Has("sponPolarization")) {
          material->SetScalar(pcc->Get("sponPolarization")->As<Double>(), SPON_POLARIZATION, Global::REAL ); 
        }
        // read remanent strain
        if (pcc->Has("sponStrain")) {
          material->SetScalar(pcc->Get("sponStrain")->As<Double>(), SPON_STRAIN, Global::REAL ); 
        }
        // 
        if (pcc->Has("Efield0")) {
          material->SetScalar(pcc->Get("Efield0")->As<Double>(), EFIELD0, Global::REAL ); 
        }
        // 
        if (pcc->Has("Stress0")) {
          material->SetScalar(pcc->Get("Stress0")->As<Double>(), STRESS0, Global::REAL ); 
        }
        // 
        if (pcc->Has("dCouple0")) {
          material->SetScalar(pcc->Get("dCouple0")->As<Double>(), DCOUPLE0, Global::REAL ); 
        }
        // read rate constant
        if (pcc->Has("rateConstant")) {
          material->SetScalar(pcc->Get("rateConstant")->As<Double>(), RATE_CONSTANT, Global::REAL ); 
        }
        // read visco-plasti index
        if (pcc->Has("viscoPlasticIndex")) {
          material->SetScalar(pcc->Get("viscoPlasticIndex")->As<Double>(), VISCO_PLASTIC_INDEX, Global::REAL ); 
        }
        // read saturation index
        if (pcc->Has("saturationIndex")) {
          material->SetScalar(pcc->Get("saturationIndex")->As<Double>(), SATURATION_INDEX, Global::REAL ); 
        }
        // read init value for volume fraction
        if (pcc->Has("volumeFracInit")) {
          material->SetScalar(pcc->Get("volumeFracInit")->As<Double>(), VOLUME_FRAC_INIT, Global::REAL ); 
        }
        // 
        if (pcc->Has("scaleForceElec")) {
          material->SetScalar(pcc->Get("scaleForceElec")->As<Double>(), SCALE_FORCE_ELEC, Global::REAL ); 
        }
        // 
        if (pcc->Has("scaleForceMech")) {
          material->SetScalar(pcc->Get("scaleForceMech")->As<Double>(), SCALE_FORCE_MECH, Global::REAL ); 
        }
        //
        if (pcc->Has("scaleForceCouple")) {
          material->SetScalar(pcc->Get("scaleForceCouple")->As<Double>(), SCALE_FORCE_COUPLE, Global::REAL ); 
        }
        // read mean temperature
        if (pcc->Has("Tmean")) {
          material->SetScalar(pcc->Get("Tmean")->As<Double>(), MEAN_TEMPERATURE, Global::REAL ); 
        }
      }
    }
    
  }
  
  //**********************************************************************
  //*************  READ MECHANICS ****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadMechanic(BaseMaterial *material, PtrParamNode mech)
  {
    // ---------
    //  density
    // ---------
    if (mech->Has("density")) {
      PtrCoefFct densCoef = ReadScalarLin( mech, "density", Global::COMPLEX );
      material->SetCoefFct(DENSITY, densCoef);
    }
    
    PtrParamNode elast = mech->Get("elasticity");
    if (elast) {
      // -------------------
      //  linear elasticity
      // -------------------
      if (elast->Has("linear")) {
        PtrParamNode lin = elast->Get("linear");
        BM::SymmetryType symType = BM::NOSYMMETRY;
        BM::CoefMap coefMap;

        symType = ReadStiffnessTensor(lin, Global::COMPLEX, coefMap);

        BM::CoefMap::iterator coefIt = coefMap.begin(), coefEnd = coefMap.end();
        for ( ; coefIt != coefEnd; ++coefIt) {
          material->SetCoefFct(coefIt->first, coefIt->second);
        }

        material->SetSymmetryType( MECH_STIFFNESS_TENSOR, symType );
      }

      // -------------------
      //  nonlinear elasticity
      // -------------------
      if (elast->Has("nonlinear") && elast->Get("nonlinear")->Has("isotropic")) {
        BaseMaterial::MatDescriptorNl nlInfo =
            ReadNonlinDescriptor(elast->Get("nonlinear")->Get("isotropic"), material);
        material->SetNonLinMatIso(MECH_EMODULUS, nlInfo);
      }
    }
    
    // -------------------
    //  viscoelasticity
    // -------------------
    if (mech->Has("viscoElasticity")) {
      PtrParamNode visco = mech->Get("viscoElasticity");

      if (visco->Has("isotropic")) {
        PtrParamNode viscoIso = visco->Get("isotropic");

        UInt pronySize;
        Vector<Double> relaxTimes, relaxModuli;

        if (viscoIso->Has("bulkModulus")) {
          // read the initial bulk modulus
          PtrCoefFct bulkInit = ReadScalarLin(viscoIso, "bulkModulus", Global::REAL);
          material->SetCoefFct(MECH_VISCO_BULK_INITIAL, bulkInit);

          // read the prony terms
          ParamNodeList pronyList = viscoIso->Get("bulkModulus")->GetList("pronyTerm");
          pronySize = pronyList.size();
          if (pronySize > 0) {
            relaxTimes.Resize(pronySize);
            relaxModuli.Resize(pronySize);
            for (UInt i = 0; i < pronySize; ++i) {
              relaxTimes[i] = pronyList[i]->Get("relaxationTime")->As<Double>();
              relaxModuli[i] = pronyList[i]->Get("relaxationModulus")->As<Double>();
            }
            material->SetVector(relaxTimes, MECH_BULK_RELAX_TIMES, Global::REAL);
            material->SetVector(relaxModuli, MECH_BULK_RELAX_MODULI, Global::REAL);
          }
        }

        if (viscoIso->Has("shearModulus")) {
          // read the initial bulk modulus
          PtrCoefFct shearInit = ReadScalarLin(viscoIso, "shearModulus", Global::REAL);
          material->SetCoefFct(MECH_VISCO_SHEAR_INITIAL, shearInit);

          // read the prony terms
          ParamNodeList pronyList = viscoIso->Get("shearModulus")->GetList("pronyTerm");
          pronySize = pronyList.size();
          if (pronySize > 0) {
            relaxTimes.Resize(pronySize);
            relaxModuli.Resize(pronySize);
            for (UInt i = 0; i < pronySize; ++i) {
              relaxTimes[i] = pronyList[i]->Get("relaxationTime")->As<Double>();
              relaxModuli[i] = pronyList[i]->Get("relaxationModulus")->As<Double>();
            }
            material->SetVector(relaxTimes, MECH_SHEAR_RELAX_TIMES, Global::REAL);
            material->SetVector(relaxModuli, MECH_SHEAR_RELAX_MODULI, Global::REAL);
          }
        }
      }
      else {
        EXCEPTION("Only isotropic viscoelasticity is supported");
      }
    }
    
    
    //read coefficients for irreversible mechanical strain
    /*if (mech->Has("irreversibleStrainCoefficient")) {
      PtrParamNode isc = mech->Get("irreversibleStrainCoefficient");
      // the dimension is only printed in the old param handler version 7562
      //if (isc->Has("dim")) std::cout << "dim=" << isc->Get("dim")->As<Integer>() << std::endl;
      
      if (isc->Has("coeffs")) {
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
    }*/ // end of irreversibleStrainCoefficient
    
    
    // check and read thermal expansion coefficients (TECs)
    if (mech->Has("thermalExpansion")) {
      PtrParamNode tec = mech->Get("thermalExpansion");
      
      // read reference temperature
      PtrCoefFct tRef;
      if (tec->Has("refTemperature")) {
        tRef = ReadScalar(tec, "refTemperature", Global::REAL);
      }
      else { // set ref temperature to zero
        tRef = CoefFunction::Generate( mp_, Global::REAL, "0.0");
      }
      material->SetCoefFct(MECH_TE_REFTEMPERATURE,tRef);
      
      // read thermal expansion coefficient to create a vector valued coef function (in Voigt notation)
      std::string coef;
      StdVector<std::string> tecR(6), tecI(6);
      tecR.Init("0.0");
      tecI.Init("0.0");

      if (ReadScalar(tec, coef, "isotropic", "real")) {
        tecR[0] = coef;
        tecR[1] = coef;
        tecR[2] = coef;
        if (ReadScalar(tec, coef, "isotropic", "imag")) {
          tecI[0] = coef;
          tecI[1] = coef;
          tecI[2] = coef;
        }
        PtrCoefFct fct = CoefFunction::Generate(mp_, Global::COMPLEX, tecR[0], tecI[0]);
        material->SetCoefFct(MECH_THERMAL_EXPANSION_SCALAR, fct);
      }
      else if (tec->Has("transversalIsotropic")) {
        PtrParamNode node = tec->Get("transversalIsotropic");
        if (ReadScalar(node, coef, "value", "real")) {
          tecR[0] = coef;
          tecR[1] = coef;
        }
        if (ReadScalar(node, coef, "value", "imag")) {
          tecI[0] = coef;
          tecI[1] = coef;
        }
        if (ReadScalar(node, coef, "value_3", "real")) {
          tecR[2] = coef;
        }
        if (ReadScalar(node, coef, "value_3", "imag")) {
          tecI[2] = coef;
        }
        PtrCoefFct fct = CoefFunction::Generate(mp_, Global::COMPLEX, tecR[0], tecI[0]);
        material->SetCoefFct(MECH_THERMAL_EXPANSION_SCALAR, fct);
        fct = CoefFunction::Generate(mp_, Global::COMPLEX, tecR[2], tecI[2]);
        material->SetCoefFct(MECH_THERMAL_EXPANSION_3, fct);
      }
      else if (tec->Has("orthotropic")) {
        PtrParamNode node = tec->Get("orthotropic");
        if (ReadScalar(node, coef, "value_1", "real")) {
          tecR[0] = coef;
        }
        if (ReadScalar(node, coef, "value_1", "imag")) {
          tecI[0] = coef;
        }
        if (ReadScalar(node, coef, "value_2", "real")) {
          tecR[1] = coef;
        }
        if (ReadScalar(node, coef, "value_2", "imag")) {
          tecI[1] = coef;
        }
        if (ReadScalar(node, coef, "value_3", "real")) {
          tecR[2] = coef;
        }
        if (ReadScalar(node, coef, "value_3", "imag")) {
          tecI[2] = coef;
        }
        PtrCoefFct fct = CoefFunction::Generate(mp_, Global::COMPLEX, tecR[0], tecI[0]);
        material->SetCoefFct(MECH_THERMAL_EXPANSION_1, fct);
        fct = CoefFunction::Generate(mp_, Global::COMPLEX, tecR[1], tecI[1]);
        material->SetCoefFct(MECH_THERMAL_EXPANSION_2, fct);
        fct = CoefFunction::Generate(mp_, Global::COMPLEX, tecR[2], tecI[2]);
        material->SetCoefFct(MECH_THERMAL_EXPANSION_3, fct);
      }
      else if(tec->Has("tensor")) {
        PtrParamNode node = tec->Get("tensor");
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

      PtrCoefFct tecVect = CoefFunction::Generate( mp_, Global::COMPLEX, tecR, tecI);
      material->SetCoefFct(MECH_THERMAL_EXPANSION_TENSOR,tecVect);
    }

    // read mechanical damping
    if (mech->Has("damping")) {
      // first rayleigh damping
      if (mech->Get("damping")->Has("rayleigh")) {
        ReadRayleighDamping(mech->Get("damping")->Get("rayleigh"), material);
      }
      if (mech->Get("damping")->Has("adaptedLossTangensDelta")) {
        ReadLossTanDeltaDamping(mech->Get("damping")->Get("adaptedLossTangensDelta"), material);
      }
      // Read Kelvin-Voigt viscous tensor
      if (mech->Get("damping")->Has("kelvinVoigt")) {
        ReadKelvinVoigtDamping(mech->Get("damping")->Get("kelvinVoigt"), material);
      }
    }

    // read real magmech coupling tensor
    if (mech->Has("magnetoStrictionTensor_h_mech")) {
      PtrCoefFct pctCoef = ReadTensor(mech->Get("magnetoStrictionTensor_h_mech"),
                                      Global::COMPLEX);
      material->SetCoefFct( MAGNETOSTRICTION_TENSOR_h_mech, pctCoef);
    }

    if (mech->Has("piezoCoupling")){
      LOG_DBG3(xmlmathandler)<< "It has piezoCoupling" << std::endl;
      PtrParamNode cpl = mech->Get("piezoCoupling", ParamNode::PASS);
      if (cpl) {
        // read piezo coupling tensor
        if (cpl->Has("linear")) {
          BM::SymmetryType symType = BM::NOSYMMETRY;

          if (cpl->Get("linear")->Has("tensor")) {
            PtrParamNode pcTensor = cpl->Get("linear")->Get("tensor");
            PtrCoefFct piezoCoef = ReadTensor(pcTensor, Global::COMPLEX );
            material->SetCoefFct(PIEZO_TENSOR, piezoCoef);
            symType = BM::GENERAL;
          }
          material->SetSymmetryType(PIEZO_TENSOR, symType);
        }

        // read nonlinearity of a coupling coefficient
        if (cpl->Has("nonlinear") && cpl->Get("nonlinear")->Has("isotropic")) {
          BaseMaterial::MatDescriptorNl nlInfo =
              ReadNonlinDescriptor(cpl->Get("nonlinear")->Get("isotropic"), material);
          material->SetNonLinMatIso(PIEZO_TENSOR, nlInfo);
        }
      }
    }
  }


  //**********************************************************************
  //*************  READ SMOOTH ****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadSmooth(BaseMaterial *material, PtrParamNode smooth)
  {
    // ---------
    //  density
    // ---------
    if (smooth->Has("density")) {
      PtrCoefFct densCoef = ReadScalarLin( smooth, "density", Global::COMPLEX );
      material->SetCoefFct(DENSITY, densCoef);
    }
    
    PtrParamNode elast = smooth->Get("elasticity");
    if (elast) {
      // -------------------
      //  linear elasticity
      // -------------------
      if (elast->Has("linear")) {
        PtrParamNode lin = elast->Get("linear");
        BM::SymmetryType symType = BM::NOSYMMETRY;
        BM::CoefMap coefMap;

        symType = ReadStiffnessTensorSmooth(lin, Global::COMPLEX, coefMap);

        BM::CoefMap::iterator coefIt = coefMap.begin(), coefEnd = coefMap.end();
        for ( ; coefIt != coefEnd; ++coefIt) {
          material->SetCoefFct(coefIt->first, coefIt->second);
        }

        material->SetSymmetryType( SMOOTH_STIFFNESS_TENSOR, symType );
      }
    }  
  }
  
  //**********************************************************************
  //*************  READ ACOUSTICS ****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadAcoustic(BaseMaterial *material, PtrParamNode acou)
  {
    //! [Read PtrParamNode]
    // read density
    if (acou->Has("density")) { // check if PtrParamNode acou has <density> tag
      PtrParamNode dens = acou->Get("density");
      if (dens->Has("linear")) {
        PtrCoefFct densFct;
        // generate a complex valued coefficient function, if only real part is given, the imaginary part will be set to 0
        densFct = ReadScalarLin(acou, "density", Global::COMPLEX);
        material->SetCoefFct(DENSITY, densFct); // set it to the material
      }
      // create coef functions for the coefficients of the rational function approximation representing the INVERSE complex frequency-dependent density
      if (dens->Has("inverseFormApproxTDEF")) {
        std::cout << std::endl
                  << "Rational function approximation provided for the inverse density (specific volume)" << std::endl;

        if (dens->Get("inverseFormApproxTDEF")->Has("rationalFunctionApprox")) {

          // read constant part of rational function approximation (high-freq limit)
          PtrCoefFct densInf;
          PtrParamNode rationalInv = dens->Get("inverseFormApproxTDEF");
          densInf = ReadScalarLin(rationalInv, "rationalFunctionApprox", Global::COMPLEX);

          material->SetCoefFct(DENSITY, CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, "1.0", densInf, CoefXpr::OP_DIV))); // set it to the material (required for acoustic environment)
          material->SetCoefFct(ACOU_TDEF_INVDENS_CONST, densInf);

          //  read the variable number of REAL poles (denominator: pole=alpha, numerator: residue=A)
          if (dens->Get("inverseFormApproxTDEF")->Get("rationalFunctionApprox")->Has("poleListReal")) {
            PtrParamNode invDensPoleRealNode;
            ParamNodeList invDensPolesReal;

            invDensPoleRealNode = rationalInv->Get("rationalFunctionApprox")->Get("poleListReal", ParamNode::PASS);
            invDensPolesReal = invDensPoleRealNode->GetChildren();

            UInt numRealP = invDensPolesReal.GetSize();
            std::cout << "Number of real poles: " << numRealP << std::endl;
            Vector<Double> parNumerRe(numRealP);
            Vector<Double> parDenomRe(numRealP);
            for (UInt i = 0; i < numRealP; i++) {
              parDenomRe[i] = invDensPolesReal[i]->Get("pole")->Get("real")->As<Double>();
              parNumerRe[i] = invDensPolesReal[i]->Get("residue")->Get("real")->As<Double>();
            }
            material->SetVector(parNumerRe, ACOU_TDEF_INVDENS_A, Global::REAL);
            material->SetVector(parDenomRe, ACOU_TDEF_INVDENS_ALPHA, Global::REAL);
          }
          else {
            // set empty vector to avoid errors
            Vector<Double> emptyVec;
            material->SetVector(emptyVec, ACOU_TDEF_INVDENS_A, Global::REAL);
          }

          //  read the variable number of COMPLEX-conjugated poles (denominator: pole=beta+j gamma, numerator: residue=B+jC)
          if (dens->Get("inverseFormApproxTDEF")->Get("rationalFunctionApprox")->Has("poleListComplex")) {
            PtrParamNode invDensPoleComplNode;
            ParamNodeList invDensPolesCompl;

            invDensPoleComplNode = rationalInv->Get("rationalFunctionApprox")->Get("poleListComplex", ParamNode::PASS);
            invDensPolesCompl = invDensPoleComplNode->GetChildren();

            UInt numComplP = invDensPolesCompl.GetSize();
            std::cout << "Number of complex-conjugated pole pairs: " << numComplP << std::endl;

            Vector<Double> parNumerRe(numComplP);
            Vector<Double> parDenomRe(numComplP);
            Vector<Double> parNumerIm(numComplP);
            Vector<Double> parDenomIm(numComplP);
            for (UInt i = 0; i < numComplP; i++) {
              parDenomRe[i] = invDensPolesCompl[i]->Get("pole")->Get("real")->As<Double>();
              parNumerRe[i] = invDensPolesCompl[i]->Get("residue")->Get("real")->As<Double>();
              parDenomIm[i] = invDensPolesCompl[i]->Get("pole")->Get("imag")->As<Double>();
              parNumerIm[i] = invDensPolesCompl[i]->Get("residue")->Get("imag")->As<Double>();
            }

            material->SetVector(parNumerRe, ACOU_TDEF_INVDENS_B, Global::REAL);
            material->SetVector(parDenomRe, ACOU_TDEF_INVDENS_BETA, Global::REAL);
            material->SetVector(parNumerIm, ACOU_TDEF_INVDENS_C, Global::REAL);
            material->SetVector(parDenomIm, ACOU_TDEF_INVDENS_GAMMA, Global::REAL);
          }
          else {
            // set empty vector to avoid errors
            Vector<Double> emptyVec;
            material->SetVector(emptyVec, ACOU_TDEF_INVDENS_B, Global::REAL);
          }
        }
      }
    }

    // read adiabatic exponent
    if (acou->Has("adiabaticExponent")) {
      PtrCoefFct fct = ReadScalarLin(acou, "adiabaticExponent", Global::REAL);
      material->SetCoefFct( FLUID_ADIABATIC_EXPONENT, fct );
    }
    
    // read compression modulus
    if (acou->Has("compressionModulus")) {
      PtrParamNode blk = acou->Get("compressionModulus");
      if (blk->Has("linear")) {
        PtrCoefFct blkFct;
        // generate a complex valued coefficient function, if only real part is given, the imaginary part will be set to 0
        blkFct = ReadScalarLin(acou, "compressionModulus", Global::COMPLEX);
        material->SetCoefFct(ACOU_BULK_MODULUS, blkFct); // set it to the material
      }
      // create coef functions for the coefficients of the rational function approximation representing the INVERSE complex frequency-dependent bulk modulus
      if (blk->Has("inverseFormApproxTDEF")) {
        std::cout << std::endl
                  << "Rational function approximation provided for the inverse bulk modulus (compressibility)" << std::endl;

        // read constant part of rational function approximation (high-freq limit)
        PtrCoefFct BlkInf;
        PtrParamNode rationalInv = blk->Get("inverseFormApproxTDEF");
        BlkInf = ReadScalarLin(rationalInv, "rationalFunctionApprox", Global::COMPLEX);

        material->SetCoefFct(ACOU_BULK_MODULUS, CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, "1.0", BlkInf, CoefXpr::OP_DIV))); // set it to the material (required for acoustic environment)
        material->SetCoefFct(ACOU_TDEF_INVBLK_CONST, BlkInf);                                                                                  // set it to the material (used in coef function)

        //  read the variable number of REAL poles (denominator: pole=alpha, numerator: residue=A)
        if (blk->Get("inverseFormApproxTDEF")->Get("rationalFunctionApprox")->Has("poleListReal")) {
          PtrParamNode invBlkPoleRealNode;
          ParamNodeList invBlkPolesReal;

          invBlkPoleRealNode = rationalInv->Get("rationalFunctionApprox")->Get("poleListReal", ParamNode::PASS);
          invBlkPolesReal = invBlkPoleRealNode->GetChildren();

          UInt numRealP = invBlkPolesReal.GetSize();
          std::cout << "Number of real poles: " << numRealP << std::endl;
          Vector<Double> parNumerRe(numRealP);
          Vector<Double> parDenomRe(numRealP);
          for (UInt i = 0; i < numRealP; i++) {
            parDenomRe[i] = invBlkPolesReal[i]->Get("pole")->Get("real")->As<Double>();
            parNumerRe[i] = invBlkPolesReal[i]->Get("residue")->Get("real")->As<Double>();
          }

          material->SetVector(parNumerRe, ACOU_TDEF_INVBLK_A, Global::REAL);
          material->SetVector(parDenomRe, ACOU_TDEF_INVBLK_ALPHA, Global::REAL);
        }
        else {
          // set empty vector to avoid errors
          Vector<Double> emptyVec;
          material->SetVector(emptyVec, ACOU_TDEF_INVBLK_A, Global::REAL);
        }

        //  read the variable number of COMPLEX-conjugated poles (denominator: pole=beta+j gamma, numerator: residue=B+jC)
        if (blk->Get("inverseFormApproxTDEF")->Get("rationalFunctionApprox")->Has("poleListComplex")) {
          PtrParamNode invBLKPoleComplNode;
          ParamNodeList invBLKPolesCompl;

          invBLKPoleComplNode = rationalInv->Get("rationalFunctionApprox")->Get("poleListComplex", ParamNode::PASS);
          invBLKPolesCompl = invBLKPoleComplNode->GetChildren();

          UInt numComplP = invBLKPolesCompl.GetSize();
          std::cout << "Number of complex-conjugated pole pairs: " << numComplP << std::endl;
          Vector<Double> parNumerRe(numComplP);
          Vector<Double> parDenomRe(numComplP);
          Vector<Double> parNumerIm(numComplP);
          Vector<Double> parDenomIm(numComplP);
          for (UInt i = 0; i < numComplP; i++) {
            parDenomRe[i] = invBLKPolesCompl[i]->Get("pole")->Get("real")->As<Double>();
            parNumerRe[i] = invBLKPolesCompl[i]->Get("residue")->Get("real")->As<Double>();
            parDenomIm[i] = invBLKPolesCompl[i]->Get("pole")->Get("imag")->As<Double>();
            parNumerIm[i] = invBLKPolesCompl[i]->Get("residue")->Get("imag")->As<Double>();
          }

          // set material params: Vector (see e.g. "pronyList")
          material->SetVector(parNumerRe, ACOU_TDEF_INVBLK_B, Global::REAL);
          material->SetVector(parDenomRe, ACOU_TDEF_INVBLK_BETA, Global::REAL);
          material->SetVector(parNumerIm, ACOU_TDEF_INVBLK_C, Global::REAL);
          material->SetVector(parDenomIm, ACOU_TDEF_INVBLK_GAMMA, Global::REAL);
        }
        else {
          // set empty vector to avoid errors
          Vector<Double> emptyVec;
          material->SetVector(emptyVec, ACOU_TDEF_INVBLK_B, Global::REAL);
        }
      }
    }

    // read kinematic viscosity
    if (acou->Has("kinematicViscosity")) {
      PtrCoefFct kinVisc = ReadScalarLin(acou, "kinematicViscosity", Global::REAL);
      material->SetCoefFct( FLUID_KINEMATIC_VISCOSITY, kinVisc );
    }

    // check for acousticDamping
    if (acou->Has("damping")) {
      // first rayleigh damping
      if (acou->Get("damping")->Has("rayleigh")) {
        ReadRayleighDamping(acou->Get("damping")->Get("rayleigh"), material);
      }
      if (acou->Get("damping")->Has("adaptedLossTangensDelta")) {
        ReadLossTanDeltaDamping(acou->Get("damping")->Get("adaptedLossTangensDelta"), material);
      }
    }

    // read acoustic non linearity
    if (acou->Has("nonlinear"))
    {
      if (acou->Get("nonlinear")->Has("bOverA")) {
        PtrCoefFct bovera = ReadScalar(acou->Get("nonlinear"), "bOverA", Global::REAL);
        material->SetCoefFct(ACOU_BOVERA, bovera);
      }
    }
  }
  
  
  //**********************************************************************
  //*************  READ ELECTROSTATICS ************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadElectrostatic(BaseMaterial *material, PtrParamNode elec)
  {
    // check for permittivity
    PtrParamNode permit = elec->Get("permittivity", ParamNode::PASS);
    PtrParamNode model;

    if (permit) {
      if (permit->Has("linear")) {
        MaterialType orthoProp[3] = {
            ELEC_PERMITTIVITY_1, ELEC_PERMITTIVITY_2, ELEC_PERMITTIVITY_3
        };
        ReadSquare3x3Tensor(permit->Get("linear"), material, ELEC_PERMITTIVITY_SCALAR,
                            orthoProp, ELEC_PERMITTIVITY_TENSOR, Global::COMPLEX);
      }
      
      if (permit->Has("nonlinear") && permit->Get("nonlinear")->Has("isotropic")) {
        BaseMaterial::MatDescriptorNl nlInfo =
            ReadNonlinDescriptor(permit->Get("nonlinear")->Get("isotropic"), material);
        material->SetNonLinMatIso(ELEC_PERMITTIVITY_SCALAR, nlInfo);
      }

      if (permit->Has("model") && permit->Get("model")->Has("isotropic")){
        if (permit->Get("model")->Get("isotropic")->Has("JilesAthertonModel")){

          model = permit->Get("model")->Get("isotropic")->Get("JilesAthertonModel");

          material->SetScalar(model->Get("Ps")->As<Double>(), MaterialType(ELEC_PS_JILES), Global::REAL );
          material->SetScalar(model->Get("a")->As<Double>(), MaterialType(ELEC_A_JILES), Global::REAL );
          material->SetScalar(model->Get("alpha")->As<Double>(), MaterialType(ELEC_ALPHA_JILES), Global::REAL );
          material->SetScalar(model->Get("k")->As<Double>(), MaterialType(ELEC_K_JILES), Global::REAL );
          material->SetScalar(model->Get("c")->As<Double>(), MaterialType(ELEC_C_JILES), Global::REAL );

        }
      }
    } // end of permittivity
    
    // read hysteresis model
    if (elec->Has("hystModel")) {
      PtrParamNode hystNode = elec->Get("hystModel");
      bool isMagnetic = false;
      ReadHysteresis(material, hystNode, isMagnetic);
//      ReadHysteresis(material, hystNode);
    }
  }
  

  //**********************************************************************
  //*************  READ ELECTRO-QUASISTATICS ************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadElecQuasistatic(BaseMaterial *material, PtrParamNode elec)
  {
    // read electric conductivity
    PtrParamNode cond = elec->Get("electricConductivity", ParamNode::PASS);
    if ( cond ) {
      if (cond->Has("linear")) {        
          MaterialType orthoProps[3] = { ELEC_CONDUCTIVITY_1, ELEC_CONDUCTIVITY_2, ELEC_CONDUCTIVITY_3 };
          ReadSquare3x3Tensor(cond->Get("linear"), material, ELEC_CONDUCTIVITY_SCALAR,
                                  orthoProps, ELEC_CONDUCTIVITY_TENSOR, Global::COMPLEX);
      }
    }
    
    // check for permittivity
    PtrParamNode permit = elec->Get("permittivity", ParamNode::PASS);
    if ( permit ) {      
      if (permit->Has("linear")) {
        MaterialType orthoProps[3] = { ELEC_PERMITTIVITY_1, ELEC_PERMITTIVITY_2, ELEC_PERMITTIVITY_3 };
          ReadSquare3x3Tensor(permit->Get("linear"), material, ELEC_PERMITTIVITY_SCALAR,
                                  orthoProps, ELEC_PERMITTIVITY_TENSOR, Global::COMPLEX);
      }
    }    
  }


  //**********************************************************************
  //*************  READ MAGNETIC *****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadMagnetic(BaseMaterial *material, PtrParamNode mag)
  {
    bool hasFixedMagnetization = false;

    // read electric conductivity
    if (mag->Has("electricConductivity")) {
      PtrParamNode cond = mag->Get("electricConductivity");
      
      if (cond->Has("linear")) {
        PtrParamNode lin = cond->Get("linear");
        MaterialType orthoProps[3] = {
            MAG_CONDUCTIVITY_1, MAG_CONDUCTIVITY_2, MAG_CONDUCTIVITY_3
        };
        ReadSquare3x3Tensor(lin, material, MAG_CONDUCTIVITY_SCALAR, orthoProps,
                            MAG_CONDUCTIVITY_TENSOR, Global::COMPLEX);
      }
      
      // we know only nonlinear isotropic material
      if (cond->Has("nonlinear")) {
        PtrParamNode condNl = cond->Get("nonlinear");

        if (condNl->Has("isotropic")) {
          BaseMaterial::MatDescriptorNl info =
              ReadNonlinDescriptor(condNl->Get("isotropic"), material);
          material->SetNonLinMatIso(MAG_CONDUCTIVITY_SCALAR, info);
        }
      }
    }


    // read permittivity
    if (mag->Has("permittivity")) {
      PtrParamNode perm = mag->Get("permittivity");

      if (perm->Has("linear")) {
        MaterialType orthoProps[3] = {
            MAG_PERMITTIVITY_1, MAG_PERMITTIVITY_2, MAG_PERMITTIVITY_3
        };
        ReadSquare3x3Tensor(perm->Get("linear"), material, MAG_PERMITTIVITY_SCALAR,
                            orthoProps, MAG_PERMITTIVITY_TENSOR, Global::COMPLEX);
      }

      // we know only nonlinear isotropic material
      if (perm->Has("nonlinear") ) {
        EXCEPTION("Nonlinearity for permittivity in DarwinPDE not yet implemented!")
      } // end of nonlinear section
    } // end of permittivity

    // read magnetic permeability
    if (mag->Has("permeability"))
    {
      // check if there is more than one permeability definition. If so, then we have the case of temperature-dependent BH-curves
      ParamNodeList permNodes = mag->GetList("permeability");
      UInt numPermNodes = permNodes.GetSize();

      if (numPermNodes > 1)
      {
        // case of temperature-dependent permeabilities
        tempDependBH_ = true;
        StdVector<BaseMaterial::MatDescriptorNl> nlData;
        nlData.Resize(numPermNodes);
        // loop over every permeability node and extract relevant data
        for (UInt i = 0; i < numPermNodes; i++)
        {
          // Let's start with some checks
          if (!permNodes[i]->Has("temperature"))
          {
            EXCEPTION("In the case of temperature-dependent permeabilities, every permeability needs a temperature attribute");
          }
          if (!permNodes[i]->Has("nonlinear"))
          {
            EXCEPTION("In the case of temperature-dependent permeabilities, every permeability needs a nonlinear node");
          }
          if (!permNodes[i]->Get("nonlinear")->Has("isotropic"))
          {
            EXCEPTION("In the case of temperature-dependent permeabilities, every nonlinear permeability definition needs an isotropic definition");
          }
          if (permNodes[i]->Get("nonlinear")->Has("anisotropic"))
          {
            EXCEPTION("In the case of temperature-dependent permeabilities, only isotropic definition are implemented");
          }

          // TODO: Currently, every permeability needs a linear node. In the case of temperature-dependent permeabilities we only
          // need to define this once, however, every linear definition is read in
          PtrParamNode perm = permNodes[i];
          if (perm->Has("linear"))
          {
            MaterialType orthoProps[3] = {
                MAG_PERMEABILITY_1, MAG_PERMEABILITY_2, MAG_PERMEABILITY_3};
            ReadSquare3x3Tensor(perm->Get("linear"), material, MAG_PERMEABILITY_SCALAR,
                                orthoProps, MAG_PERMEABILITY_TENSOR, Global::COMPLEX);
          }
          
          // read parameters
          BaseMaterial::MatDescriptorNl &info = nlData[i];
          info = ReadNonlinDescriptor(permNodes[i]->Get("nonlinear/isotropic"), material);
          info.temperature = 0.0;
          info.analyticExpr = "";
          info.analyticExprDeriv = "";

          // read temperature for which the current BH curve is defined
          info.temperature = permNodes[i]->Get("temperature")->As<Double>();
          // read analytic function of material parameter
          if(permNodes[i]->Has("nonlinear/isotropic/nuExpr")){
            info.analyticExpr = permNodes[i]->Get("nonlinear/isotropic/nuExpr")->As<std::string>().c_str();
          }
          // read analytic derivative of material parameter
          if(permNodes[i]->Has("nonlinear/isotropic/nuDerivExpr")){
            info.analyticExprDeriv = permNodes[i]->Get("nonlinear/isotropic/nuDerivExpr")->As<std::string>().c_str();
          }
        }
        material->SetNonLinMatIsoTempDependBH(MAG_PERMEABILITY_SCALAR, nlData);
      }
      else
      {
        // Classical case with just one permeability
        PtrParamNode perm = mag->Get("permeability");
        PtrParamNode model;

        if (perm->Has("linear"))
        {
          MaterialType orthoProps[3] = {
              MAG_PERMEABILITY_1, MAG_PERMEABILITY_2, MAG_PERMEABILITY_3};
          ReadSquare3x3Tensor(perm->Get("linear"), material, MAG_PERMEABILITY_SCALAR,
                              orthoProps, MAG_PERMEABILITY_TENSOR, Global::COMPLEX);
        }

        // we know only nonlinear isotropic material
        if (perm->Has("nonlinear"))
        {
          PtrParamNode permNl = perm->Get("nonlinear");
          if (permNl->Has("isotropic"))
          {
            PtrParamNode iso = permNl->Get("isotropic");
            BaseMaterial::MatDescriptorNl info = ReadNonlinDescriptor(iso, material);
            // read analytic function of material parameter
            if (iso->Has("nuExpr"))
            {
              info.analyticExpr = iso->Get("nuExpr")->As<std::string>();
            }
            // read analytic derivative of material parameter
            if (iso->Has("nuDerivExpr"))
            {
              info.analyticExprDeriv = iso->Get("nuDerivExpr")->As<std::string>();
            }
            if (iso->Has("muExpr"))
            {
              info.analyticExpr = iso->Get("muExpr")->As<std::string>();
            }
            // read analytic derivative of material parameter
            if (iso->Has("muDerivExpr"))
            {
              info.analyticExprDeriv = iso->Get("muDerivExpr")->As<std::string>();
            }
            // pass info to material class
            material->SetNonLinMatIso(MAG_PERMEABILITY_SCALAR, info);
          } // nonlinear isotropic material
          else if (permNl->Has("anisotropic"))
          {
            // anisotropic case: bundle of nonlinear curves
            PtrParamNode aniso = permNl->Get("anisotropic");
            // fetch paramnodes for hdbc
            ParamNodeList anIsoNodes = aniso->GetList("data");
            if (anIsoNodes.GetSize() > 0)
            {
              StdVector<BaseMaterial::MatDescriptorNl> nlData;
              nlData.Resize(anIsoNodes.GetSize());
              // iterate over all parameter nodes
              for (UInt i = 0; i < anIsoNodes.GetSize(); i++)
              {
                // read parameters
                BaseMaterial::MatDescriptorNl &info = nlData[i];
                info = ReadNonlinDescriptor(anIsoNodes[i], material);
                info.angle = 0.0;
                info.zScaling = 1.0;
                info.analyticExpr = "";
                info.analyticExprDeriv = "";

                // read angle
                if (anIsoNodes[i]->Has("angle"))
                {
                  info.angle = anIsoNodes[i]->Get("angle")->As<Double>();
                }
                // read z-scaling factor
                if (anIsoNodes[i]->Has("zScaling"))
                {
                  info.zScaling = anIsoNodes[i]->Get("zScaling")->As<Double>();
                }
                // read analytic function of material parameter
                if (anIsoNodes[i]->Has("nuExpr"))
                {
                  info.analyticExpr = anIsoNodes[i]->Get("nuExpr")->As<std::string>().c_str();
                }
                // read analytic derivative of material parameter
                if (anIsoNodes[i]->Has("nuDerivExpr"))
                {
                  info.analyticExprDeriv = anIsoNodes[i]->Get("nuDerivExpr")->As<std::string>().c_str();
                }
              }
              material->SetNonLinMatAniso(MAG_PERMEABILITY_SCALAR, nlData);
            }
          } // end of anisotropic nonlinear material
        }   // end of nonlinear section
        if (perm->Has("model") && perm->Get("model")->Has("isotropic"))
        {
          if (perm->Get("model")->Get("isotropic")->Has("EBHysteresisModel"))
          {
            model = perm->Get("model")->Get("isotropic")->Get("EBHysteresisModel");

            material->SetScalar(model->Get("Ps")->As<Double>(), MaterialType(MAG_PS_EB), Global::REAL);
            material->SetScalar(model->Get("A")->As<Double>(), MaterialType(MAG_A_EB), Global::REAL);
            material->SetScalar(model->Get("mu0")->As<Double>(), MaterialType(MAG_MU0_EB), Global::REAL);
            material->SetScalar(model->Get("numS")->As<Double>(), MaterialType(MAG_NUMS_EB), Global::REAL);
            material->SetScalar(model->Get("chi_factor")->As<Double>(), MaterialType(MAG_CHI_FACTOR_EB), Global::REAL);
          }
        }
      }
    } // end of permeability

  // read in prescribed magnetization
    if(mag->Has("prescribedMagnetization"))
    {
      hasFixedMagnetization = true;
      
      Matrix<Double> prescribedMagnetization = Matrix<Double>(1,3);
      prescribedMagnetization.Init();
      
      if(mag->Get("prescribedMagnetization")->Has("MagnetizationVector")){
        ParamTools::AsTensor<double>(mag->Get("prescribedMagnetization")->Get("MagnetizationVector"),1, 3, prescribedMagnetization);
      }
      
      material->SetScalar( prescribedMagnetization[0][0], PRESCRIBED_MAGNETIZATION_X, Global::REAL);
      material->SetScalar( prescribedMagnetization[0][1], PRESCRIBED_MAGNETIZATION_Y, Global::REAL);
      material->SetScalar( prescribedMagnetization[0][2], PRESCRIBED_MAGNETIZATION_Z, Global::REAL);
    }
    
    if(hasFixedMagnetization){
      material->SetScalar( 1, PRESCRIBED_MAGNETIZATION);
    } else {
      material->SetScalar( 0, PRESCRIBED_MAGNETIZATION);
    }

    // read nonlinear reluctivity for magnetostrictive strains
    if (mag->Has("reluctivity_MagStrict")) {
      // we know only nonlinear isotropic material
      if (mag->Get("magneticReluctivity_MagStrict")->Has("nonlinear") ) {
        if (mag->Get("magneticReluctivity_MagStrict")->Get("nonlinear")->Has("isotropic")) {
          PtrParamNode iso = mag->Get("magneticReluctivity_MagStrict")->Get("nonlinear")->Get("isotropic");
          BaseMaterial::MatDescriptorNl info = ReadNonlinDescriptor(iso, material);
          material->SetNonLinMatIso(MAGSTRICT_RELUCTIVITY, info);
        }
      } // end of nonlinear section
    } // end of magneticReluctivity_MagStrict

    // read real magmech coupling tensor
    if (mag->Has("magnetoStrictionTensor_h_mag")) {
      PtrCoefFct pctCoef = ReadTensor(mag->Get("magnetoStrictionTensor_h_mag"),
                                      Global::COMPLEX);
      material->SetCoefFct( MAGNETOSTRICTION_TENSOR_h_mag, pctCoef);
    }

    // read hysteresis model
    if (mag->Has("hystModel")) {
      PtrParamNode hystNode = mag->Get("hystModel");

      if(hasFixedMagnetization == false){
        bool isMagnetic = true;
        ReadHysteresis(material, hystNode, isMagnetic);
      } else {
        WARN("Hysteresis model cannot be used as magnetization is prescribed explicitely. No hysteresis model will be added.");
      }      
    }

    // read core loss
    if (mag->Has("coreLoss")){
      PtrParamNode clParam = mag->Get("coreLoss");
      BaseMaterial::MatDescriptorNl info = ReadNonlinDescriptor(clParam, material);
      material->SetNonLinMatIso( MAG_CORE_LOSS_PER_MASS, info );
    }
    
    // read density needed for core loss
    PtrCoefFct densFct;
    if (mag->Has("density")){
      densFct = ReadScalarLin(mag, "density", Global::REAL);
    } else {
      densFct = CoefFunction::Generate( mp_, Global::REAL, "0.0");
    }
    material->SetCoefFct( DENSITY, densFct );
    
  }
  
  //**********************************************************************
  //*************  READ THERMIC ******************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadThermic(BaseMaterial *material, PtrParamNode therm)
  {
    // read density
    if (therm->Has("density")) {
      PtrParamNode dens = therm->Get("density");
      
      if (dens->Has("linear")) {
        PtrCoefFct densFct = ReadScalarLin(therm, "density", Global::REAL);
        material->SetCoefFct( DENSITY, densFct );
      }

      // we know only nonlinear isotropic material
      if (dens->Has("nonlinear") && dens->Get("nonlinear")) {
        PtrParamNode nl = dens->Get("nonlinear");
        BaseMaterial::MatDescriptorNl info = ReadNonlinDescriptor(nl, material);
        material->SetNonLinMatIso(DENSITY, info);
      }
    }
    
    // read heat capacity
    if (therm->Has("heatCapacity")) {
      PtrParamNode capa = therm->Get("heatCapacity");

      if (capa->Has("linear")) {
        PtrCoefFct capaFct = ReadScalarLin(therm, "heatCapacity", Global::REAL);
        material->SetCoefFct( HEAT_CAPACITY, capaFct );
      }

      // we know only nonlinear isotropic material
      if (capa->Has("nonlinear") && capa->Get("nonlinear")) {
        PtrParamNode nl = capa->Get("nonlinear");
        BaseMaterial::MatDescriptorNl info = ReadNonlinDescriptor(nl, material);
        material->SetNonLinMatIso(HEAT_CAPACITY, info);
      }
    }
    
    // read thermal conductivity
    if (therm->Has("heatConductivity")) {
      PtrParamNode cond = therm->Get("heatConductivity");

      if (cond->Has("linear")) {
        MaterialType orthoProps[3] = {
            HEAT_CONDUCTIVITY_1, HEAT_CONDUCTIVITY_2, HEAT_CONDUCTIVITY_3
        };
        ReadSquare3x3Tensor(cond->Get("linear"), material, HEAT_CONDUCTIVITY_SCALAR,
                            orthoProps, HEAT_CONDUCTIVITY_TENSOR, Global::COMPLEX);
      }
      
      // we know only nonlinear isotropic material
      if (cond->Has("nonlinear") && cond->Get("nonlinear")->Has("isotropic")) {
        PtrParamNode iso = cond->Get("nonlinear")->Get("isotropic");
        BaseMaterial::MatDescriptorNl info = ReadNonlinDescriptor(iso, material);
        material->SetNonLinMatIso(HEAT_CONDUCTIVITY_SCALAR, info);
      }
      
    }

    // read reference temperature
    if (therm->Has("refTemperature")) {
      PtrCoefFct refTemp = ReadScalarLin(therm, "refTemperature", Global::REAL);
      material->SetCoefFct( HEAT_REF_TEMPERATURE, refTemp );
    }
  }
  
  //**********************************************************************
  //*************  READ FLOW *********************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadFlow(BaseMaterial *material, PtrParamNode flow)
  {    
    // read density
    if (flow->Has("density")) {
      PtrCoefFct densFct = ReadScalarLin(flow, "density", Global::REAL);
      material->SetCoefFct( DENSITY, densFct );
    }
    
    // read dynamicViscosity 
    if (flow->Has("dynamicViscosity") && flow->Has("kinematicViscosity")) {
      EXCEPTION("Please specify either dynamic or kinematic viscosity but not both!");
    }
    
    if (flow->Has("dynamicViscosity")) {
      PtrCoefFct dynVisc = ReadScalarLin(flow, "dynamicViscosity", Global::REAL);
      material->SetCoefFct( FLUID_DYNAMIC_VISCOSITY, dynVisc );
    }
    
    // read kinematicViscosity 
    if (flow->Has("kinematicViscosity")) {
      PtrCoefFct kinVisc = ReadScalarLin(flow, "kinematicViscosity", Global::REAL);
      material->SetCoefFct( FLUID_KINEMATIC_VISCOSITY, kinVisc );
    }

    // read bulk viscosity
    if (flow->Has("bulkViscosity")) {
      PtrCoefFct bulkVisc = ReadScalarLin(flow, "bulkViscosity", Global::REAL);
      material->SetCoefFct( FLUID_BULK_VISCOSITY, bulkVisc );
    }

    // read adiabatic exponent
    if (flow->Has("adiabaticExponent")) {
    	PtrCoefFct exp = ReadScalarLin(flow, "adiabaticExponent", Global::REAL);
    	material->SetCoefFct( FLUID_ADIABATIC_EXPONENT, exp );
    }

    // read compression modulus
    if (flow->Has("compressionModulus")) {
        PtrCoefFct blkFct = ReadScalarLin(flow, "compressionModulus", Global::REAL);
        material->SetCoefFct( FLUID_BULK_MODULUS, blkFct );
    }
  }
  
  //**********************************************************************
  //*************  READ TEST *********************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadTest(BaseMaterial *material, PtrParamNode test)
  {
    // read alpha
    if (test->Has("alpha")) {
      PtrCoefFct alphaFct =
              CoefFunction::Generate(mp_, Global::REAL,
              test->Get("alpha")->As<std::string>() );
      material->SetCoefFct( TEST_ALPHA, alphaFct );
    }
    
    // read beta
    if (test->Has("beta")) {
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
    REFACTOR
    /*if (pyro->Has("pyroCoefficient")) {
      if (pyro->Get("pyroCoefficient")->Has("linear")) {
        if (pyro->Get("pyroCoefficient")->Get("linear")->Has("tensor")) {
          PtrCoefFct py = ReadTensor(pyro->Get("pyroCoefficient")->Get("linear")->Get("tensor"),
                                     Global::REAL);
          material->SetCoefFct(PYROCOEFFICIENT_TENSOR, py)
        }
      }
    }*/
  }
  
  //**********************************************************************
  //*************  READ THERMOELASTIC ************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadThermoelastic(BaseMaterial *material,
          PtrParamNode thermExp) {
    REFACTOR
    /*if (thermExp->Has("thermalExpansion")) {
      if (thermExp->Get("thermalExpansion")->Has("tensor")) {
        MaterialType orthoProps[3] = {
            MECH_THERMAL_EXPANSION_1,
            MECH_THERMAL_EXPANSION_2,
            MECH_THERMAL_EXPANSION_3
        };
        PtrCoefFct te = ReadSquare3x3Tensor(thermExp->Get("thermalExpansion")->Get("tensor"),
                                            material,
                                            MECH_THERMAL_EXPANSION_SCALAR,
                                            orthoProps,
                                            MECH_THERMAL_EXPANSION_TENSOR,
                                            Global::REAL);
        material->SetCoefFct(MECH_THERMAL_EXPANSION_TENSOR, te);
      }
    }*/
  }
  
  //**********************************************************************
  //*************  READ MAGNETOSTRICTIVE  ********************************
  //**********************************************************************
  void XMLMaterialHandler::ReadMagStrict(BaseMaterial *material, PtrParamNode pn) {
    //read real magmech coupling tensor
    if (pn->Has("magnetoStrictionTensor_h"))
    {
      PtrCoefFct pctCoef = ReadTensor(pn->Get("magnetoStrictionTensor_h"), Global::COMPLEX);
      material->SetCoefFct( MAGNETOSTRICTION_TENSOR_h, pctCoef);
    }
  }
  
  //**********************************************************************
  //*************  READ ELECTRIC CONDUCTIVITY ******************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadElectricConduction(BaseMaterial *material, PtrParamNode elec)
  {
    // read electric conductivity
    if (elec->Has("electricConductivity")) {
      PtrParamNode cond = elec->Get("electricConductivity");

      if (cond->Has("linear")) {
        MaterialType orthoProps[3] = {
            ELEC_CONDUCTIVITY_1, ELEC_CONDUCTIVITY_2, ELEC_CONDUCTIVITY_3
        };
        ReadSquare3x3Tensor(cond->Get("linear"), material, ELEC_CONDUCTIVITY_SCALAR,
            orthoProps, ELEC_CONDUCTIVITY_TENSOR, Global::COMPLEX);
      }
      
      // we know only nonlinear isotropic material
      if (cond->Has("nonlinear") && cond->Get("nonlinear")->Has("isotropic")) {
        PtrParamNode iso = cond->Get("nonlinear")->Get("isotropic");
        BaseMaterial::MatDescriptorNl info = ReadNonlinDescriptor(iso, material);
        material->SetNonLinMatIso(ELEC_CONDUCTIVITY_SCALAR, info);
      } // nonlinear isotropic material
    }
  }
  
  //**********************************************************************
  //*************  READ HYSTERESIS****************************************
  //**********************************************************************
  void XMLMaterialHandler::ReadHysteresis(BaseMaterial *material, PtrParamNode hystNode, bool isMagnetic){
    PtrParamNode operatorNode = NULL;

    if (hystNode->Has("elecPolarization")){
      operatorNode = hystNode->Get("elecPolarization");
      ReadHystOperator(material, operatorNode, false, isMagnetic);
    } else if(hystNode->Has("magPolarization")){
      operatorNode = hystNode->Get("magPolarization");
      ReadHystOperator(material, operatorNode, false, isMagnetic);
    } else {
      EXCEPTION("Either elecPolarization or magPolarization has to be defined for hysteresis node.")
    }
    
    /*
     * tracing related flags (tracing of hyst operator is done in CoefFunctionHyst.cc)
     * > must be done via mat.xml; as tracing is done for each material separately
     */
    Double trace_JacResolution = 1e-5;
    bool trace_forceCentral = false;
    bool trace_forceRetracing = false;
    
    PtrParamNode tracingNode = NULL;
    if(hystNode->Has("AdaptTracingOfHystOperator")){
      tracingNode = hystNode->Get("AdaptTracingOfHystOperator");
      if(tracingNode->Has("JacResolution")){
        trace_JacResolution = tracingNode->Get("JacResolution")->As<Double>();
      }
      if(tracingNode->Has("forceCentral")){
        trace_forceCentral = tracingNode->Get("forceCentral")->As<bool>();       
      }   
      if(tracingNode->Has("forceRetracing")){
        trace_forceRetracing = tracingNode->Get("forceRetracing")->As<bool>();
      } 
    } 
    
    material->SetScalar(trace_JacResolution, TRACE_JAC_RESOLUTION, Global::REAL );
    if (trace_forceCentral){
      material->SetScalar(1, TRACE_FORCE_CENTRALDIFF);
    } else {
      material->SetScalar(0, TRACE_FORCE_CENTRALDIFF);
    }
    if (trace_forceRetracing){
      material->SetScalar(1, TRACE_FORCE_RETRACING);
    } else {
      material->SetScalar(0, TRACE_FORCE_RETRACING);
    }
    
    PtrParamNode couplingNode = NULL;
    int usePolarization = 0;
    if (hystNode->Has("piezoCouplingAndStrains")){
      couplingNode = hystNode->Get("piezoCouplingAndStrains");
    } else if (hystNode->Has("magstrictCouplingAndStrains")){
      couplingNode = hystNode->Get("magstrictCouplingAndStrains");
    }

    int couplingDefined = 0;
    if (couplingNode != NULL){
      couplingDefined = 1;
      PtrParamNode strainModelingNode = couplingNode->Get("strainModeling");
      int irrStrainImplementation = -1;
      Double strainSat = 0.0;
      Double irrStrains_c1 = 0.0;
      Double irrStrains_c2 = 0.0;
      Double irrStrains_c3 = 0.0;
      Double irrStrains_d0 = 0.0;
      Double irrStrains_d1 = 0.0;
      bool scaleToStrainSat = false;
      bool paramsDefinedForHalfRange = false;
      int coefdim = 1;
      Matrix<Double> ciCoefs = Matrix<Double>(1,1);

      if (strainModelingNode != NULL){
        PtrParamNode innerNode = NULL;
        if (strainModelingNode->Has("muDatWolf")){
          irrStrainImplementation = 0;
          innerNode = strainModelingNode->Get("muDatWolf");
          strainSat = innerNode->Get("strainSat")->As<Double>();
          irrStrains_c1 = innerNode->Get("c1")->As<Double>();
          irrStrains_c2 = innerNode->Get("c2")->As<Double>();
          irrStrains_c3 = innerNode->Get("c3")->As<Double>();
          scaleToStrainSat = innerNode->Get("scaleToStrainSat")->As<bool>();
          paramsDefinedForHalfRange = innerNode->Get("forHalfRange")->As<bool>();
        } else if (strainModelingNode->Has("muDatLoeffler")){
          irrStrainImplementation = 1;
          innerNode = strainModelingNode->Get("muDatLoeffler");
          strainSat = innerNode->Get("strainSat")->As<Double>();
          coefdim = innerNode->Get("dim_ci")->As<Integer>();
          ciCoefs.Resize(1,coefdim);
          ParamTools::AsTensor<double>(innerNode->Get("ci"),1, coefdim, ciCoefs);
          irrStrains_d0 = innerNode->Get("d0")->As<Double>();
          irrStrains_d1 = innerNode->Get("d1")->As<Double>();
          scaleToStrainSat = innerNode->Get("scaleToStrainSat")->As<bool>();
          paramsDefinedForHalfRange = innerNode->Get("forHalfRange")->As<bool>();
        } else if (strainModelingNode->Has("higherOrderPolynomial")){
          irrStrainImplementation = 2;
          innerNode = strainModelingNode->Get("higherOrderPolynomial");
          strainSat = innerNode->Get("strainSat")->As<Double>();
          coefdim = innerNode->Get("dim_ci")->As<Integer>();
          ciCoefs.Resize(1,coefdim);
          ParamTools::AsTensor<double>(innerNode->Get("ci"),1, coefdim, ciCoefs);
          scaleToStrainSat = innerNode->Get("scaleToStrainSat")->As<bool>();
        }
      }
      material->SetScalar(strainSat, S_SATURATION, Global::REAL );
      material->SetScalar(irrStrainImplementation, HYST_IRRSTRAINS);
      material->SetScalar(irrStrains_c1, HYST_IRRSTRAIN_C1, Global::REAL);
      material->SetScalar(irrStrains_c2, HYST_IRRSTRAIN_C2, Global::REAL);
      material->SetScalar(irrStrains_c3, HYST_IRRSTRAIN_C3, Global::REAL);
      material->SetScalar(irrStrains_d0, HYST_IRRSTRAIN_D0, Global::REAL);
      material->SetScalar(irrStrains_d1, HYST_IRRSTRAIN_D1, Global::REAL);
      material->SetScalar(coefdim, HYST_IRRSTRAIN_CI_SIZE);
      material->SetTensor(ciCoefs, HYST_IRRSTRAIN_CI, Global::REAL);
      if (scaleToStrainSat){
        material->SetScalar(1, HYST_IRRSTRAIN_SCALETOSAT);
      } else {
        material->SetScalar(0, HYST_IRRSTRAIN_SCALETOSAT);
      }
      if (paramsDefinedForHalfRange){
        material->SetScalar(1, HYST_IRRSTRAIN_PARAMSFORHALFRANGE);
      } else {
        material->SetScalar(0, HYST_IRRSTRAIN_PARAMSFORHALFRANGE);
      }

      PtrParamNode strainOperatorNode = couplingNode->Get("hystOperatorForStrains");
      if (strainOperatorNode->Has("usePolarization")){
        usePolarization = 1;
      } else if (strainOperatorNode->Has("separateHystOperator")){
        operatorNode = strainOperatorNode->Get("separateHystOperator");
        ReadHystOperator(material, operatorNode, true, isMagnetic);
      }
      material->SetScalar(usePolarization, IRRSTRAIN_REUSE_P);

      PtrParamNode smallSignalNode = couplingNode->Get("smallSignalForm");
      /*
       * -2 : no coupling at all; no coupling node
       * -1 : coupled but without small signal tensors
       *  0 : coupled e-form/h-form (piezo/magstrict)
       *  1 : coupled d-form (piezo)
       *  2 : coupled g-form (magstrict)
       */
      int strainForm = -2;
      if (smallSignalNode != NULL){

        if (smallSignalNode->Has("noSmallSignalCoupling")){
          strainForm = -1;
        }
        if (smallSignalNode->Has("piezo_eform")){
          strainForm = 0;
        }
        if (smallSignalNode->Has("piezo_dform")){
          strainForm = 1;
        }
        if (smallSignalNode->Has("magstrict_hform")){
          strainForm = 0;
        }
        if (smallSignalNode->Has("magstrict_gform")){
          strainForm = 2;
        }
      }
      material->SetScalar(strainForm, HYST_STRAIN_FORM);
    } else {
      // no coupling defined
      material->SetScalar(-2, HYST_STRAIN_FORM);
    }
    material->SetScalar(couplingDefined, HYST_COUPLING_DEFINED );

  }

  void XMLMaterialHandler::ReadHystOperator(BaseMaterial *material, PtrParamNode operatorNode, bool setStrains, bool isMagnetic){
//    std::cout << "ReadHystOperator" << std::endl;
    PtrParamNode model;
    PtrParamNode pWeight = NULL;
    PtrParamNode pAnhyst = NULL;
//    bool anhystTEAM32 = false;
    PtrParamNode pInversion = NULL;
    PtrParamNode pStrainForm = NULL;
    PtrParamNode initialState = NULL;
    PtrParamNode irrStrainNode = NULL;
    bool isVector = true;
    bool isPreisachType = true;
    Double pSat = 1.0;
    Double hSat = 1.0;         

    /*
     * Instead of setting each parameter twice in the form
     *    if (setStrains){
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
    if (setStrains){
      enumOffset = 100;
    }

    if (operatorNode->Has("scalarPreisach")){
      isVector = false;
      model = operatorNode->Get("scalarPreisach");
      if (setStrains){
        material->SetString("scalarPreisach", HYST_MODEL_STRAIN);
        material->SetString("SCALAR", HYSTERESIS_DIM_STRAIN);
      } else {
        material->SetString("scalarPreisach", HYST_MODEL);
        material->SetString("SCALAR", HYSTERESIS_DIM);
      }

      // read input/output saturation of Preisach hysterese model
      if(model->Has("inputSat")){
        hSat = model->Get("inputSat")->As<Double>();
      } else {
        EXCEPTION("Hysteresis model must have parameter inputSat; something went wrong here!");
      }
      
      if(model->Has("outputSat")){
        pSat = model->Get("outputSat")->As<Double>();
      } else {
        EXCEPTION("Hysteresis model must have parameter outputSat; something went wrong here!");
      }
      material->SetScalar(hSat, MaterialType(X_SATURATION+enumOffset), Global::REAL );
      material->SetScalar(pSat, MaterialType(Y_SATURATION+enumOffset), Global::REAL );

      Matrix<Double> directionVector;
      if (model->Has("dirPolarization"))
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

      if (model->Has("weights")){
        pWeight = model->Get("weights");
        // Read in weights separately below as they require the same steps for VectorHysteresis, too
      }

      // new: scalar model might use the vector inversion, too (but not recommended)
      if (model->Has("hystInversion")){
        pInversion = model->Get("hystInversion");
      }

      if (model->Has("initialState")){
        initialState = model->Get("initialState");
      }
    }
    else if (operatorNode->Has("JilesAtherton")){
      isVector = false;

      isPreisachType = false;
      model = operatorNode->Get("JilesAtherton");
      if (setStrains){
        material->SetString("JilesAtherton", HYST_MODEL_STRAIN);
        material->SetString("SCALAR", HYSTERESIS_DIM_STRAIN);
      } else {
        material->SetString("JilesAtherton", HYST_MODEL);
        material->SetString("SCALAR", HYSTERESIS_DIM);
      }

      if (model->Has("jiles_test")){
        material->SetScalar(model->Get("jiles_test")->As<Double>(), MaterialType(JILES_TEST), Global::REAL );
      }
    }
    else if (operatorNode->Has("vectorPreisach_Sutor")){
//      std::cout << "vectorPreisach_Sutor" << std::endl;
      model = operatorNode->Get("vectorPreisach_Sutor");

      material->SetString("vectorPreisach_Sutor", MaterialType(HYST_MODEL+enumOffset));
      material->SetString("VECTOR", MaterialType(HYSTERESIS_DIM+enumOffset));

      // read input/output saturation of Preisach hysterese model      
      if(model->Has("inputSat")){
        hSat = model->Get("inputSat")->As<Double>();
      } else {
        EXCEPTION("Hysteresis model must have parameter inputSat; something went wrong here!");
      }
      
      if(model->Has("outputSat")){
        pSat = model->Get("outputSat")->As<Double>();
      } else {
        EXCEPTION("Hysteresis model must have parameter outputSat; something went wrong here!");
      }
      material->SetScalar(hSat, MaterialType(X_SATURATION+enumOffset), Global::REAL );
      material->SetScalar(pSat, MaterialType(Y_SATURATION+enumOffset), Global::REAL );

      /*
       * new numbering: 1 -> classical vector model (sutor2012)
       *                2 -> revised model (sutor2015)
       *                10 -> classical vector model, matrix based
       *                20 -> revised model, matrix based
       */
      int evalVersion = 2;
      if (model->Has("evalVersion")){
//        std::cout << "evalVersion" << std::endl;
        if (model->Get("evalVersion")->Has("Classical_Version_2012__list_implementation")){
          evalVersion = 1;
        }
        if (model->Get("evalVersion")->Has("Revised_Version_2015__list_implementation")){
          evalVersion = 2;
        }
        if (model->Get("evalVersion")->Has("Classical_Version_2012__matrix_implementation")){
          evalVersion = 10;
        }
        if (model->Get("evalVersion")->Has("Revised_Version_2015__matrix_implementation")){
          evalVersion = 20;
        }
      }

      material->SetScalar(evalVersion, MaterialType(EVAL_VERSION+enumOffset));

      Double rotResistance = 1.0;
      Double angularDistance = 0.0;
      bool enforceSatOutputAtSatInput = true;
      if (model->Has("rotResistance")){
//        std::cout << "rotResistance" << std::endl;
        rotResistance = model->Get("rotResistance")->As<Double>();
      }

      material->SetScalar(rotResistance, MaterialType(ROT_RESISTANCE+enumOffset), Global::REAL);

      if (model->Has("angularDistance")){
//                std::cout << "angularDistance" << std::endl;
        angularDistance = model->Get("angularDistance")->As<Double>();
      }

      material->SetScalar(angularDistance, MaterialType(ANG_DISTANCE+enumOffset), Global::REAL);

      if (model->Has("enforceSatOutputAtSatInput")){
//        std::cout << "enforceSatOutputAtSatInput" << std::endl;
        enforceSatOutputAtSatInput = model->Get("enforceSatOutputAtSatInput")->As<bool>();
      }

      if (enforceSatOutputAtSatInput){
        material->SetScalar(1, MaterialType(SCALETOSAT+enumOffset));
      } else {
        material->SetScalar(0, MaterialType(SCALETOSAT+enumOffset));
      }

      if (model->Has("weights")){
        pWeight = model->Get("weights");
        // Read in weights separately below as they require the same steps for VectorHysteresis, too
      }

      if (model->Has("hystInversion")){
        pInversion = model->Get("hystInversion");
      }

      if (model->Has("initialState")){
        initialState = model->Get("initialState");
      }

      PtrParamNode debugNode = NULL;
      Double angClip = 0.0;
      Double angRes = 1e-16;
      Double ampRes = 1e-16;
      int printOut = 0;
      int bmpRes = 1000;

      if (model->Has("debuggingParameter")){
        debugNode = model->Get("debuggingParameter");
        if (debugNode->Has("angularClipping")){
          angClip = debugNode->Get("angularClipping")->As<Double>();
        }

        if (debugNode->Has("angularResolution")){
          angRes = debugNode->Get("angularResolution")->As<Double>();
        }

        if (debugNode->Has("amplitudeResolution")){
          ampRes = debugNode->Get("amplitudeResolution")->As<Double>();
        }

        /*
         * if printOut > 0, the overlaid rotation and switching state of each printOut timestep will be
         * written to a bmp file of resolution bmpResolution
         */
        if(debugNode->Has("printOut")){
          printOut = debugNode->Get("printOut")->As<Integer>();
        }

        if(debugNode->Has("bmpResolution")){
          bmpRes = debugNode->Get("bmpResolution")->As<Integer>();
        }
      }

      material->SetScalar(angClip, MaterialType(ANG_CLIPPING+enumOffset), Global::REAL);
      material->SetScalar(angRes, MaterialType(ANG_RESOLUTION+enumOffset), Global::REAL);
      material->SetScalar(ampRes, MaterialType(AMP_RESOLUTION+enumOffset), Global::REAL);
      // extra parameter that will not be set for strains
      if (!setStrains){
        material->SetScalar(printOut, PRINT_PREISACH);
        material->SetScalar(bmpRes, PRINT_PREISACH_RESOLUTION);
      }
    }

    else if (operatorNode->Has("vectorPreisach_Mayergoyz")){
      model = operatorNode->Get("vectorPreisach_Mayergoyz");

      material->SetString("vectorPreisach_Mayergoyz", MaterialType(HYST_MODEL+enumOffset));
      material->SetString("VECTOR", MaterialType(HYSTERESIS_DIM+enumOffset));

      PtrParamNode innerModel = NULL;
      PtrParamNode singleModel = NULL;
      if (model->Has("isotropic")){
        // read isotropic Mayergoyz vector model
        innerModel = model->Get("isotropic");

        material->SetScalar(1, MaterialType(PREISACH_MAYERGOYZ_ISOTROPIC+enumOffset) );

        Double lossParam_a = 0;
        Double lossParam_b = 0;
        // for fine-tuning; explanetion can be found in CFS_MAThysteresis.xsd
        bool useAbsoluteValueOfdPhi = false;
        bool normalizeXInExponentOfG = true;
        int restrictionOfPsi = 4;
        Double fixScalingOfXInExponentOfG = 1.0;
        
        if(innerModel->Has("rotLossCorrectionFactors")){
          if(innerModel->Get("rotLossCorrectionFactors")->Has("lossParam_a")){
            lossParam_a = innerModel->Get("rotLossCorrectionFactors")->Get("lossParam_a")->As<Double>();
          }
          if(innerModel->Get("rotLossCorrectionFactors")->Has("lossParam_b")){
            lossParam_b = innerModel->Get("rotLossCorrectionFactors")->Get("lossParam_b")->As<Double>();
          }
          
          if(innerModel->Get("rotLossCorrectionFactors")->Has("useAbsoluteValueOfdPhi")){
            useAbsoluteValueOfdPhi = innerModel->Get("rotLossCorrectionFactors")->Get("useAbsoluteValueOfdPhi")->As<bool>();
          }
          if(innerModel->Get("rotLossCorrectionFactors")->Has("normalizeXInExponentOfG")){
            normalizeXInExponentOfG = innerModel->Get("rotLossCorrectionFactors")->Get("normalizeXInExponentOfG")->As<bool>();
          }
          if(innerModel->Get("rotLossCorrectionFactors")->Has("restrictionOfPsi")){
            restrictionOfPsi = innerModel->Get("rotLossCorrectionFactors")->Get("restrictionOfPsi")->As<Integer>();
          }
          if(innerModel->Get("rotLossCorrectionFactors")->Has("fixScalingOfXInExponentOfG")){
            fixScalingOfXInExponentOfG = innerModel->Get("rotLossCorrectionFactors")->Get("fixScalingOfXInExponentOfG")->As<Double>();
          }
        }
        
        material->SetScalar( lossParam_a, MaterialType(MAYERGOYZ_LOSSPARAM_A+enumOffset), Global::REAL);
        material->SetScalar( lossParam_b, MaterialType(MAYERGOYZ_LOSSPARAM_B+enumOffset), Global::REAL);
        material->SetScalar( fixScalingOfXInExponentOfG, MaterialType(MAYERGOYZ_SCALINGOFXINEXP+enumOffset), Global::REAL);
        
        if(useAbsoluteValueOfdPhi==true){
          material->SetScalar( 1, MaterialType(MAYERGOYZ_USEABSDPHI+enumOffset));
        } else {
          material->SetScalar( 0, MaterialType(MAYERGOYZ_USEABSDPHI+enumOffset));
        }
        if(normalizeXInExponentOfG==true){
          material->SetScalar( 1, MaterialType(MAYERGOYZ_NORMALIZEXINEXP+enumOffset));
        } else {
          material->SetScalar( 0, MaterialType(MAYERGOYZ_NORMALIZEXINEXP+enumOffset));
        }
        material->SetScalar( restrictionOfPsi, MaterialType(MAYERGOYZ_RESTRICTIONOFPSI+enumOffset));
                
        int numDir = 11;
        if (innerModel->Has("numDirections")){
          numDir = innerModel->Get("numDirections")->As<Integer>();
        }
        material->SetScalar(numDir, MaterialType(PREISACH_MAYERGOYZ_NUM_DIR+enumOffset) );

        if (innerModel->Has("ScalarModel")){
          singleModel = innerModel->Get("ScalarModel");
        } else {
          EXCEPTION("Single scalar Preisach model required for isotropic vector model");
        }

        Matrix<Double> startAxis = Matrix<Double>(1,3);
        if (innerModel->Has("startAxis")){
          ParamTools::AsTensor<double>(innerModel->Get("startAxis"),1, 3, startAxis);
          // normalize
          Double startAxisNorm = startAxis.NormL2();
          if (startAxisNorm != 0){
            startAxis[0][0] /= startAxisNorm;
            startAxis[0][1] /= startAxisNorm;
            startAxis[0][2] /= startAxisNorm;
          }
          // random vector will be generated for each element > not done here
        } else {
          startAxis.Init();
          // default case > x axis
          startAxis[0][0] = 1.0;
        }

        material->SetScalar( startAxis[0][0], MaterialType(MAYERGOYZ_STARTAXIS_X+enumOffset), Global::REAL);
        material->SetScalar( startAxis[0][1], MaterialType(MAYERGOYZ_STARTAXIS_Y+enumOffset), Global::REAL);
        material->SetScalar( startAxis[0][2], MaterialType(MAYERGOYZ_STARTAXIS_Z+enumOffset), Global::REAL);

        // read input/output saturation of Preisach hysterese model
        if(singleModel->Has("inputSat")){
          hSat = singleModel->Get("inputSat")->As<Double>();
        } else {
          EXCEPTION("Hysteresis model must have parameter inputSat; something went wrong here!");
        }
        
        if(singleModel->Has("outputSat")){
          pSat = singleModel->Get("outputSat")->As<Double>();
        } else {
          EXCEPTION("Hysteresis model must have parameter outputSat; something went wrong here!");
        }
        
        material->SetScalar(hSat, MaterialType(X_SATURATION+enumOffset), Global::REAL );
        material->SetScalar(pSat, MaterialType(Y_SATURATION+enumOffset), Global::REAL );

        if (singleModel->Has("weights")){
          pWeight = singleModel->Get("weights");
          // Read in weights separately below as they require the same steps for VectorHysteresis, too
        }

        int adaptedToVectorCase = 0;
        if (singleModel->Has("weightsAdaptedToMayergoyzVectorModel")){
          bool adapted = singleModel->Get("weightsAdaptedToMayergoyzVectorModel")->As<bool>();
          if (adapted){
            adaptedToVectorCase = 1;
          }
        }
        material->SetScalar(adaptedToVectorCase, MaterialType(PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR+enumOffset));

      } else if (model->Has("anIsotropic")){
        EXCEPTION("Anisotropic Mayergoyz vector hysteresis model not yet supported");
      }

      int clipOutput = 2;
      if (model->Has("clipOutputToSat")){
        if (model->Get("clipOutputToSat")->Has("noClipping")){
          clipOutput = 0;
        }
        if (model->Get("clipOutputToSat")->Has("clipAmplitude")){
          clipOutput = 1;
        }
        if (model->Get("clipOutputToSat")->Has("clipComponentParallelToInput")){
          // leads to most reasonable results
          clipOutput = 2;
        }
      }
      material->SetScalar(clipOutput, MaterialType(PREISACH_MAYERGOYZ_CLIPOUTPUT+enumOffset));

      if (model->Has("hystInversion")){
        pInversion = model->Get("hystInversion");
      }

      if (model->Has("initialState")){
        initialState = model->Get("initialState");
      }
    }
    else {
      material->SetString("none", MaterialType(HYST_MODEL+enumOffset));

      return;
    }

    UInt isPreisachType_int = 0;
    if (isPreisachType){
      isPreisachType_int = 1;
    }
    // set if selected model is of Preisach-Type or not (e.g., in case of Jiles-Atherton)
    material->SetScalar(isPreisachType_int, MaterialType(HYST_TYPE_IS_PREISACH+enumOffset));

    if (isPreisachType){
      if (pWeight != NULL){
        // Read in weights
        int dim = -1;
        if (pWeight->Has("dim_weights")) dim = pWeight->Get("dim_weights")->As<Integer>();

        material->SetScalar( dim, MaterialType(PREISACH_WEIGHTS_DIM+enumOffset));

        int weightType = 0; // 0 = const, 1 = muDat, 2 = muDatExtended, 3 = givenTensor
        PtrParamNode pWeightInner;
        if (pWeight->Has("weightType")){
          pWeightInner = pWeight->Get("weightType");
        }

        if (pWeightInner->Has("const")){
          weightType = 0;
          Double constValue = pWeightInner->Get("const")->As<Double>();

          material->SetScalar(constValue, MaterialType(PREISACH_WEIGHTS_CONSTVALUE+enumOffset), Global::REAL );
          material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));

        } else if(pWeightInner->Has("muLorentz")){
//          Parameter for Lorentzian weight function as used and defined by
//          Wolf in \n
//          "Generalisiertes Preisach-Modell für die Simulation und Kompensation der Hysterese piezokeramischer Wandler"
//          Mu_Lorentzian(alpha,beta) = \n
//           \t     A/(1 + ((alpha-h)/(h)*sigma)*((alpha-h)/(h)*sigma) ) * \n
//           \t     A/(1 + ((beta+h)/(h)*sigma)*((beta+h)/(h)*sigma) );
//          <xsd:element name="A" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="h1" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="h2" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="sigma1" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="sigma2" type="xsd:double" minOccurs="1"/>
          
          /*
           * Note: to save enums, store A,h1,h2,sigma1 and sigma2 under the same ENUMs as MUDAT-parameter
           * to distinguish betwee MUDAT and Lorentzian weights, we have the weightType
           */
          weightType = 4;
          PtrParamNode lorentzWolf = pWeightInner->Get("muLorentzExtended");
          Double A = lorentzWolf->Get("A")->As<Double>();
          Double h = lorentzWolf->Get("h2")->As<Double>();
          Double sigma = lorentzWolf->Get("sigma2")->As<Double>();
          bool paramsDefinedForHalfRangeBool = lorentzWolf->Get("forHalfRange")->As<bool>();
          int paramsDefinedForHalfRange = 0;
          if(paramsDefinedForHalfRangeBool){
            paramsDefinedForHalfRange = 1;
          }
          material->SetScalar(A, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL );
          material->SetScalar(h, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL );
          material->SetScalar(sigma, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL );
          material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));
          material->SetScalar(paramsDefinedForHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));

        } else if(pWeightInner->Has("muLorentzExtended")){
//          Parameter for Lorentzian weight function as used and defined by
//          Wolf in \n
//          "Generalisiertes Preisach-Modell f��r die Simulation und Kompensation der Hysterese piezokeramischer Wandler"
//          Mu_Lorentzian(alpha,beta) = \n
//           \t     A/(1 + ((alpha-h2)/(h2)*sigma2)*((alpha-h2)/(h2)*sigma2) ) * \n
//           \t     A/(1 + ((beta+h1)/(h1)*sigma1)*((beta+h1)/(h1)*sigma1) );
//          <xsd:element name="A" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="h1" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="h2" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="sigma1" type="xsd:double" minOccurs="1"/>
//          <xsd:element name="sigma2" type="xsd:double" minOccurs="1"/>        
          /*
           * Note: to save enums, store A,h1,h2,sigma1 and sigma2 under the same ENUMs as MUDAT-parameter
           * to distinguish betwee MUDAT and Lorentzian weights, we have the weightType
           */
          weightType = 5;
          PtrParamNode lorentzWolf = pWeightInner->Get("muLorentzExtended");
          Double A = lorentzWolf->Get("A")->As<Double>();
          Double h1 = lorentzWolf->Get("h1")->As<Double>();
          Double h2 = lorentzWolf->Get("h2")->As<Double>();
          Double sigma1 = lorentzWolf->Get("sigma1")->As<Double>();
          Double sigma2 = lorentzWolf->Get("sigma2")->As<Double>();
          bool paramsDefinedForHalfRangeBool = lorentzWolf->Get("forHalfRange")->As<bool>();
          int paramsDefinedForHalfRange = 0;
          if(paramsDefinedForHalfRangeBool){
            paramsDefinedForHalfRange = 1;
          }
          material->SetScalar(A, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL );
          material->SetScalar(h1, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL );
          material->SetScalar(h2, MaterialType(PREISACH_WEIGHTS_MUDAT_H2+enumOffset), Global::REAL );
          material->SetScalar(sigma1, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL );
          material->SetScalar(sigma2, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA2+enumOffset), Global::REAL );
          material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));
          material->SetScalar(paramsDefinedForHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));
          
        } else if(pWeightInner->Has("muDat")){
          weightType = 1;
          PtrParamNode muDat = pWeightInner->Get("muDat");
          Double A = muDat->Get("A")->As<Double>();
          Double h = muDat->Get("h2")->As<Double>();
          Double sigma = muDat->Get("sigma2")->As<Double>();
          Double eta = muDat->Get("eta")->As<Double>();

          int paramsDefinedForHalfRange = 0;
          bool paramsDefinedForHalfRangeBool = muDat->Get("forHalfRange")->As<bool>();
          if (paramsDefinedForHalfRangeBool){
            paramsDefinedForHalfRange = 1;
          }

          material->SetScalar(A, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL );
          material->SetScalar(h, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL );
          material->SetScalar(sigma, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL );
          material->SetScalar(eta, MaterialType(PREISACH_WEIGHTS_MUDAT_ETA+enumOffset), Global::REAL );
          material->SetScalar(paramsDefinedForHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));
          material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));

        } else if (pWeightInner->Has("muDatExtended")){
          weightType = 2;
          PtrParamNode muDatExt = pWeightInner->Get("muDatExtended");
          Double A = muDatExt->Get("A")->As<Double>();
          Double h1 = muDatExt->Get("h1")->As<Double>();
          Double h2 = muDatExt->Get("h2")->As<Double>();
          Double sigma1 = muDatExt->Get("sigma1")->As<Double>();
          Double sigma2 = muDatExt->Get("sigma2")->As<Double>();
          Double eta = muDatExt->Get("eta")->As<Double>();

          int paramsDefinedForHalfRange = 0;
          bool paramsDefinedForHalfRangeBool = muDatExt->Get("forHalfRange")->As<bool>();
          if (paramsDefinedForHalfRangeBool){
            paramsDefinedForHalfRange = 1;
          }

          material->SetScalar(A, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL );
          material->SetScalar(h1, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL );
          material->SetScalar(h2, MaterialType(PREISACH_WEIGHTS_MUDAT_H2+enumOffset), Global::REAL );
          material->SetScalar(sigma1, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL );
          material->SetScalar(sigma2, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA2+enumOffset), Global::REAL );
          material->SetScalar(eta, MaterialType(PREISACH_WEIGHTS_MUDAT_ETA+enumOffset), Global::REAL );
          material->SetScalar(paramsDefinedForHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));
          material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));

        } else if (pWeightInner->Has("weightTensor")){
          weightType = 3;
          Matrix<Double> preisachWeightTensor(dim,dim);
          ParamTools::AsTensor<double>(pWeightInner->Get("weightTensor"), dim, dim, preisachWeightTensor);

          material->SetTensor( preisachWeightTensor, MaterialType(PREISACH_WEIGHTS_TENSOR+enumOffset), Global::REAL);
          material->SetScalar(weightType, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));

        } else {
          EXCEPTION("No valid Preisach weights found");
        }
      } else {
        EXCEPTION("No valid Preisach weights found");
      }

      if (pWeight->Has("anhystereticParameter")){
        pAnhyst = pWeight->Get("anhystereticParameter");
      }
            
    } else {
      std::cout << "Non-Preisach Hysteresis Model found!" << std::endl;
    }

    Double a,b,c,d;
    bool onlyanhyst;
    bool anhystForHalfRange;
//    bool cInAtan; // no longer used; own paramater d instroduced instead
    bool anhystOutputSat;
    if(pAnhyst != NULL){
//        Add anhysteretic part to output of hyst operator to obtain overall POLARIZATION
//            p_total(e) = p(e) + a*atan(b*e + d) + c*e
//        Sutor version (cInAtan=false): p_total(e) = p(e) + a*atan(b*e) + c*e > d = 0
//        Loeffler version (DEFAULT,cInAtan=true): p_total(e) = p(e) + a*atan(b*(e+d_star)) > c = 0; d = b*d_star
//        note1: a,b,c and d are meant for NORMALIZED p and e \n
//        note2: in their works, Sutor and Loeffler model p,e to be normalized to range [-0.5,0.5]
//        in CFS p,e are normalized to [-1,1] such that a,b,c have to be adapted which is done by CFS if
//        flag forHalfRange = true
//        
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
        if(pAnhyst->Has("d")){
          d = pAnhyst->Get("d")->As<Double>();
        } else {
          d = 0.0;
        }
        if(pAnhyst->Has("onlyAnhyst")){
          onlyanhyst = pAnhyst->Get("onlyAnhyst")->As<bool>();
        } else {
          onlyanhyst = false;
        }

        if(pAnhyst->Has("forHalfRange")){
          anhystForHalfRange = pAnhyst->Get("forHalfRange")->As<bool>();
        } else {
          anhystForHalfRange = false;
        }
        if(pAnhyst->Has("anhystPartCountsTowardsOutputSat")){
          anhystOutputSat = pAnhyst->Get("anhystPartCountsTowardsOutputSat")->As<bool>();
        } else {
          anhystOutputSat = true;
        }

    } else {
      a = 0;
      b = 0;
      c = 0;
      d = 0;
      onlyanhyst = false;
      anhystForHalfRange = false;
      anhystOutputSat = true;
    }

    material->SetScalar(a, MaterialType(PREISACH_WEIGHTS_ANHYST_A+enumOffset), Global::REAL);
    material->SetScalar(b, MaterialType(PREISACH_WEIGHTS_ANHYST_B+enumOffset), Global::REAL);
    material->SetScalar(c, MaterialType(PREISACH_WEIGHTS_ANHYST_C+enumOffset), Global::REAL);
    material->SetScalar(d, MaterialType(PREISACH_WEIGHTS_ANHYST_D+enumOffset), Global::REAL);
    if(onlyanhyst){
      material->SetScalar(1, MaterialType(PREISACH_WEIGHTS_ANHYST_ONLY+enumOffset));
    } else {
      material->SetScalar(0, MaterialType(PREISACH_WEIGHTS_ANHYST_ONLY+enumOffset));
    }
    if (anhystForHalfRange){
      material->SetScalar(1, MaterialType(PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE+enumOffset));
    } else {
      material->SetScalar(0, MaterialType(PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE+enumOffset));
    }
    if(anhystOutputSat){
      material->SetScalar(1, MaterialType(PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT+enumOffset));
    } else {
      material->SetScalar(0, MaterialType(PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT+enumOffset));
    }

    Matrix<Double> initialStateTensor = Matrix<Double>(1,3);
    initialStateTensor.Init();
    bool prescribeOutput = false;
    bool scaleBySaturation = false;

    if (initialState != NULL){
      PtrParamNode InOutState = NULL;
      if (initialState->Has("initialOutput")){
        InOutState = initialState->Get("initialOutput");
        prescribeOutput = true;
      } else if (initialState->Has("initialInput")){
        InOutState = initialState->Get("initialInput");
        prescribeOutput = false;
      } else {
        EXCEPTION("Either input or output must be given at this point");
      }
      if (InOutState != NULL){
        if (InOutState->Has("Vector")){
          //std::cout << "InitialState found" << std::endl;
          ParamTools::AsTensor<double>(InOutState->Get("Vector"),1, 3, initialStateTensor);
          //std::cout << "IntialState: " << initialStateTensor.ToString() << std::endl;
        }
        if (InOutState->Has("scaleVectorBySaturation")){
          scaleBySaturation = InOutState->Get("scaleVectorBySaturation")->As<bool>();
        }
      }
    }

    // strain hyst operator will automatically get the same initial input as the polarization
    // if output is prescribed for polarization, the resulting input will be passed to strain operator
    if (!setStrains){
      material->SetScalar( initialStateTensor[0][0], INITIAL_STATE_X, Global::REAL);
      material->SetScalar( initialStateTensor[0][1], INITIAL_STATE_Y, Global::REAL);
      material->SetScalar( initialStateTensor[0][2], INITIAL_STATE_Z, Global::REAL);
      if (scaleBySaturation){
        material->SetScalar(1, PREISACH_SCALEINITIALSTATE);
      } else {
        material->SetScalar(0, PREISACH_SCALEINITIALSTATE);
      }

      if (prescribeOutput){
        material->SetScalar(1, PREISACH_PRESCRIBEOUTPUT);
      } else {
        material->SetScalar(0, PREISACH_PRESCRIBEOUTPUT);
      }
    }

    // common parameter for all three inversion methods
    int maxNumOuterIts = 50;
    double tolH = 1e-12;
    double tolB = 1e-12;
    int tolHrel = 0;
    int tolBrel = 0;
//     typedef enum {LEVENBERGMARQUARDT=0, NEWTON=1, JACOBIFREENEWTON=2, PROJECTEDLM=3, EVERETTBASED=4, FIXPOINT=5, NONE=6} localInversionFlag; 
    // 0 = LM, 1 = Newton, 2 = JacobiFreeNewtonKrylov, 3 = projected LM, 4 = Everett based, 5 = FP

    localInversionFlag inversionMethod;
    if(isVector){
      inversionMethod = LOCAL_NEWTON;
    } else {
      inversionMethod = LOCAL_EVERETTBASED;
    }

    // LM parameter
    int maxNumberRegularizationIterations = 50;
    double alphaRegStart = 0.25;
    double alphaRegMin = 0.001953125;
    double alphaRegMax = 8192.0;
    double trustRegionLow = 0.15;
    double trustRegionMid = 0.35;
    double trustRegionHigh = 0.85;
    double jacRes = 1e-7;
    // -1 = no jacobian > for JacobiFreeNewtonKrylov
    //  0 = forward/backward differences
    //  1 = central differences
    //  2 = forward/backward differences with scaling of diagonal
    int jacImplementation = 2;

    // Newton parameter
    int maxNumberLinesearchIterations = 50;
    double alphaLSMin = 0.001;
    double alphaLSMax = 1.0;
    bool stopLSAtLocalMin = false;
    // jacRes and jacImplementation as for LM

    // JacobiFreeKrylovNewton
    // > no extra parameter compared with Newton

    // parameter for projected Levenberg-Marquardt
    // > "Levenberg-Marquardt methods for constrained nonlinear equations with strong local convergence properties"
    //     -Kanzow,Yamashita,Fukushima
    // > values used by authors in section 4
    Double projLM_rho = 1e-8;
    Double projLM_beta = 0.9;
    Double projLM_sigma = 1e-4;
    Double projLM_gamma = 0.99995;
    Double projLM_p = 2.1;

    // parameters not given clearly in article;
    // mu, tau, c
    // > tau and c are from Armijo type linsearch taken from
    // > https://en.wikipedia.org/wiki/Backtracking_line_search
    // only ranges are known:
    //  mu > 0
    //  tau, c in (0,1)
    //
    Double projLM_mu = 1.0;
    Double projLM_tau = 0.8;
    Double projLM_c = 0.8;

    // FP
    Double convergenceFactor = 1.0;

    // important: for electrostatics, we need no inversion and should not set the
    // parameter above
    bool setInversion = false;
    bool printWarnings = false;

    if (pInversion != NULL){

      setInversion = true;
      if (pInversion->Has("InversionMethod"))
      {
        if(pInversion->Get("InversionMethod")->Has("Fixpoint")){
          inversionMethod = LOCAL_FIXPOINT;
          PtrParamNode invMethod = pInversion->Get("InversionMethod")->Get("Fixpoint");
          if (invMethod->Has("convergenceFactor")){
            convergenceFactor = invMethod->Get("convergenceFactor")->As<double>();
          }
        } else if (pInversion->Get("InversionMethod")->Has("LevenbergMarquardtWithTrustregion")){
//          std::cout << "LevenbergMarquardtWithTrustregion" << std::endl;
          inversionMethod = LOCAL_LEVENBERGMARQUARDT;
          PtrParamNode invMethod = pInversion->Get("InversionMethod")->Get("LevenbergMarquardtWithTrustregion");

          if (invMethod->Has("maxNumberRegularizationIterations")){
            maxNumberRegularizationIterations = invMethod->Get("maxNumberRegularizationIterations")->As<Integer>();
          }

          if (invMethod->Has("alphaRegStart")){
            alphaRegStart = invMethod->Get("alphaRegStart")->As<double>();
          }
          if (invMethod->Has("alphaRegMin")){
            alphaRegMin = invMethod->Get("alphaRegMin")->As<double>();
          }
          if (invMethod->Has("alphaRegMax")){
            alphaRegMax = invMethod->Get("alphaRegMax")->As<double>();
          }

          if (invMethod->Has("trustRegionLow")){
            trustRegionLow = invMethod->Get("trustRegionLow")->As<double>();
          }
          if (invMethod->Has("trustRegionMid")){
            trustRegionMid = invMethod->Get("trustRegionMid")->As<double>();
          }
          if (invMethod->Has("trustRegionHigh")){
            trustRegionHigh = invMethod->Get("trustRegionHigh")->As<double>();
          }
          if (invMethod->Has("jacobiResolution")){
            jacRes = invMethod->Get("jacobiResolution")->As<double>();
          }
          if (invMethod->Has("jacobiImplementation")){
            if (invMethod->Get("jacobiImplementation")->Has("ForwardBackwardDifferences")){
              jacImplementation = 0;
            }
            if (invMethod->Get("jacobiImplementation")->Has("CentralDifferences")){
              jacImplementation = 1;
            }
            if (invMethod->Get("jacobiImplementation")->Has("ForwardBackwardWithScaledDiagonal")){
              jacImplementation = 2;
            }
          }

        } else if (pInversion->Get("InversionMethod")->Has("NewtonWithLinesearch")){
//          std::cout << "NewtonWithLinesearch" << std::endl;
          inversionMethod = LOCAL_NEWTON;
          PtrParamNode invMethod = pInversion->Get("InversionMethod")->Get("NewtonWithLinesearch");

          if (invMethod->Has("numberOfLinesearchSteps")){
            maxNumberLinesearchIterations = invMethod->Get("numberOfLinesearchSteps")->As<Integer>();
          }

          if (invMethod->Has("alphaLSMin")){
            alphaLSMin = invMethod->Get("alphaLSMin")->As<double>();
          }
          if (invMethod->Has("alphaLSMax")){
            alphaLSMax = invMethod->Get("alphaLSMax")->As<double>();
          }

          if (invMethod->Has("stopLinesearchAtLocalMin")){
            stopLSAtLocalMin = invMethod->Get("stopLinesearchAtLocalMin")->As<bool>();
          }

          if (invMethod->Has("jacobiResolution")){
            jacRes = invMethod->Get("jacobiResolution")->As<double>();
          }
          if (invMethod->Has("jacobiImplementation")){
            if (invMethod->Get("jacobiImplementation")->Has("ForwardBackwardDifferences")){
              jacImplementation = 0;
            }
            if (invMethod->Get("jacobiImplementation")->Has("CentralDifferences")){
              jacImplementation = 1;
            }
            if (invMethod->Get("jacobiImplementation")->Has("ForwardBackwardWithScaledDiagonal")){
              jacImplementation = 2;
            }
          }

        } else if (pInversion->Get("InversionMethod")->Has("JacobianFreeNewtonKrylovWithLinesearch")){
//          std::cout << "JacobianFreeNewtonKrylovWithLinesearch" << std::endl;
          inversionMethod = LOCAL_JACOBIFREENEWTON;
          PtrParamNode invMethod = pInversion->Get("InversionMethod")->Get("JacobianFreeNewtonKrylovWithLinesearch");

          if (invMethod->Has("numberOfLinesearchSteps")){
            maxNumberLinesearchIterations = invMethod->Get("numberOfLinesearchSteps")->As<Integer>();
          }

          if (invMethod->Has("alphaLSMin")){
            alphaLSMin = invMethod->Get("alphaLSMin")->As<double>();
          }
          if (invMethod->Has("alphaLSMax")){
            alphaLSMax = invMethod->Get("alphaLSMax")->As<double>();
          }

          if (invMethod->Has("stopLinesearchAtLocalMin")){
            stopLSAtLocalMin = invMethod->Get("stopLinesearchAtLocalMin")->As<bool>();
          }

          jacImplementation = -1;

        } else if (pInversion->Get("InversionMethod")->Has("ProjectedLMWithLinesearch")){
//          std::cout << "ProjectedLMWithLinesearch" << std::endl;
          inversionMethod = LOCAL_PROJECTEDLM;
          PtrParamNode invMethod = pInversion->Get("InversionMethod")->Get("ProjectedLMWithLinesearch");

          if (invMethod->Has("numberOfLinesearchSteps")){
            maxNumberLinesearchIterations = invMethod->Get("numberOfLinesearchSteps")->As<Integer>();
          }

          if (invMethod->Has("alphaLSMin")){
            alphaLSMin = invMethod->Get("alphaLSMin")->As<double>();
          }
          if (invMethod->Has("alphaLSMax")){
            alphaLSMax = invMethod->Get("alphaLSMax")->As<double>();
          }

          if (invMethod->Has("jacobiResolution")){
            jacRes = invMethod->Get("jacobiResolution")->As<double>();
          }
          if (invMethod->Has("jacobiImplementation")){
            if (invMethod->Get("jacobiImplementation")->Has("ForwardBackwardDifferences")){
              jacImplementation = 0;
            }
            if (invMethod->Get("jacobiImplementation")->Has("CentralDifferences")){
              jacImplementation = 1;
            }
            if (invMethod->Get("jacobiImplementation")->Has("ForwardBackwardWithScaledDiagonal")){
              jacImplementation = 2;
            }
          }

          if (invMethod->Has("rho")){
            projLM_rho = invMethod->Get("rho")->As<double>();
          }
          if (invMethod->Has("beta")){
            projLM_beta = invMethod->Get("beta")->As<double>();
          }
          if (invMethod->Has("sigma")){
            projLM_sigma = invMethod->Get("sigma")->As<double>();
          }
          if (invMethod->Has("gamma")){
            projLM_gamma = invMethod->Get("gamma")->As<double>();
          }
          if (invMethod->Has("p")){
            projLM_p = invMethod->Get("p")->As<double>();
          }
          if (invMethod->Has("mu")){
            projLM_mu = invMethod->Get("mu")->As<double>();
          }
          if (invMethod->Has("tau")){
            projLM_tau = invMethod->Get("tau")->As<double>();
          }
          if (invMethod->Has("c")){
            projLM_c = invMethod->Get("c")->As<double>();
          }
        } else if(pInversion->Get("InversionMethod")->Has("EverettBasedInversion")){          
          inversionMethod = LOCAL_EVERETTBASED;
          if(isVector){
            EXCEPTION("EverettBasedInversion only supported for scalar model");
          }
        } else {
          inversionMethod = LOCAL_NOINVERSION;
          EXCEPTION("No valid method selected");
        }

        if (pInversion->Has("maxNumberOuterIterations"))
        {
          maxNumOuterIts = pInversion->Get("maxNumberOuterIterations")->As<Integer>();
        }

        if (pInversion->Has("residualTolH"))
        {
          PtrParamNode tolHNode = pInversion->Get("residualTolH");
          tolH = tolHNode->Get("value")->As<double>();
          bool isRel = tolHNode->Get("isRelative")->As<bool>();
          if(isRel){
            tolHrel = 1;
          } else {
            tolHrel = 0;
          }
        }
        if(pInversion->Has("residualTolB"))
        {
          PtrParamNode tolBNode = pInversion->Get("residualTolB");
          tolB = tolBNode->Get("value")->As<double>();
          bool isRel = tolBNode->Get("isRelative")->As<bool>();
          if(isRel){
            tolBrel = 1;
          } else {
            tolBrel = 0;
          }
        }

        if (pInversion->Has("printWarnings"))
        {
          printWarnings = pInversion->Get("printWarnings")->As<bool>();
        }

      }
    }
    //Hyst operator for strains does not use inversion! Only forward mode used
    if ((setInversion==true) && (setStrains==false)){
      material->SetScalar(maxNumOuterIts, MAX_NUM_IT_HYST_INV);
      material->SetScalar(tolH, RES_TOL_H_HYST_INV, Global::REAL);
      material->SetScalar(tolB, RES_TOL_B_HYST_INV, Global::REAL);
      material->SetScalar(tolHrel, RES_TOL_H_HYST_INV_ISREL);
      material->SetScalar(tolBrel, RES_TOL_B_HYST_INV_ISREL);

      material->SetScalar(inversionMethod, VEC_HYST_INV_METHOD);

      material->SetScalar(maxNumberRegularizationIterations, MAX_NUM_REG_IT_HYST_INV);
      material->SetScalar(alphaRegStart, ALPHA_REG_HYST_INV, Global::REAL);
      material->SetScalar(alphaRegMin, ALPHA_REG_MIN_HYST_INV, Global::REAL);
      material->SetScalar(alphaRegMax, ALPHA_REG_MAX_HYST_INV, Global::REAL);
      material->SetScalar(trustRegionLow, TRUST_LOW_HYST_INV, Global::REAL);
      material->SetScalar(trustRegionMid, TRUST_MID_HYST_INV, Global::REAL);
      material->SetScalar(trustRegionHigh, TRUST_HIGH_HYST_INV, Global::REAL);

      material->SetScalar(jacRes, JAC_RESOLUTION_HYST_INV, Global::REAL);
      material->SetScalar(jacImplementation, JAC_IMPLEMENTATION_HYST_INV);

      material->SetScalar(maxNumberLinesearchIterations, MAX_NUM_LS_IT_HYST_INV);
      material->SetScalar(alphaLSMin, ALPHA_LS_MIN_HYST_INV, Global::REAL);
      material->SetScalar(alphaLSMax, ALPHA_LS_MAX_HYST_INV, Global::REAL);

      if (stopLSAtLocalMin == true){
        material->SetScalar(1, STOP_INV_LS_AT_LOCAL_MIN);
      } else {
        material->SetScalar(0, STOP_INV_LS_AT_LOCAL_MIN);
      }

      if (printWarnings == true){
        material->SetScalar(1, HYST_LOCAL_INVERSION_PRINT_WARNINGS);
      } else {
        material->SetScalar(0, HYST_LOCAL_INVERSION_PRINT_WARNINGS);
      }

      material->SetScalar(projLM_mu, HYST_INV_PROJLM_MU, Global::REAL);
      material->SetScalar(projLM_rho, HYST_INV_PROJLM_RHO, Global::REAL);
      material->SetScalar(projLM_beta, HYST_INV_PROJLM_BETA, Global::REAL);
      material->SetScalar(projLM_sigma, HYST_INV_PROJLM_SIGMA, Global::REAL);
      material->SetScalar(projLM_gamma, HYST_INV_PROJLM_GAMMA, Global::REAL);
      material->SetScalar(projLM_tau, HYST_INV_PROJLM_TAU, Global::REAL);
      material->SetScalar(projLM_c, HYST_INV_PROJLM_C, Global::REAL);
      material->SetScalar(projLM_p, HYST_INV_PROJLM_P, Global::REAL);

      // FP
      material->SetScalar(convergenceFactor, HYST_INV_FP_SAFETYFACTOR, Global::REAL);


    }
//    std::cout << "ReadHystOperator - done" << std::endl;
  }


  //!======================== helper methods ==================================

  bool XMLMaterialHandler::ReadScalar(PtrParamNode ptrNode, std::string& val,
                                      std::string str1, std::string str2) {

    bool isAvailable = false;
    if ( ptrNode->Has(str1) ) {
      PtrParamNode scalar = ptrNode->Get(str1);

      // read the real part
      if ( scalar->Has(str2) ) {
        val = scalar->Get(str2)->As<std::string>();
        isAvailable = true;
      }
    }
    return isAvailable;
  }

  PtrCoefFct XMLMaterialHandler::ReadScalar( PtrParamNode ptrNode,
                                             std::string str,
                                             Global::ComplexPart type )
  {
    std::string realVal = "0.0";
    std::string imagVal = "0.0";

    if ( type == Global::REAL || type == Global::COMPLEX ) {
       ReadScalar( ptrNode, realVal, str, "real" );
    }

    if ( type == Global::IMAG || type == Global::COMPLEX ) {
      ReadScalar( ptrNode, imagVal, str, "imag" );
    }
    return CoefFunction::Generate(mp_, type, realVal, imagVal );
  }

  bool XMLMaterialHandler::ReadScalarLin(PtrParamNode ptrNode, std::string& val,
                                         std::string str1, std::string str2) {
    bool isAvailable = false;
    if ( ptrNode->Has(str1) ) {
      PtrParamNode linear = ptrNode->Get(str1)->Get("linear");
      // read the real part
      if ( linear->Has(str2) ) {
        val = linear->Get(str2)->As<std::string>();
        isAvailable = true;
      }
    }
    return isAvailable;
  }

  PtrCoefFct XMLMaterialHandler::ReadScalarLin( PtrParamNode ptrNode,
                                                std::string str,
                                                Global::ComplexPart type )
  {
    std::string realVal = "0.0";
    std::string imagVal = "0.0";

    if ( type == Global::REAL || type == Global::COMPLEX ) {
      ReadScalarLin( ptrNode, realVal, str, "real" );
    }

    if ( type == Global::IMAG || type == Global::COMPLEX ) {
      ReadScalarLin( ptrNode, imagVal, str, "imag" );
    }
    return CoefFunction::Generate(mp_, type, realVal, imagVal );
  }


  PtrCoefFct XMLMaterialHandler::ReadTensor( PtrParamNode ptrNode,
                                             Global::ComplexPart type )
  {
    // Obtain dimensions
    UInt dim1 = ptrNode->Get("dim1")->As<UInt>();
    UInt dim2;
    if ( !ptrNode->Has("dim2") ) {
      dim2 = dim1;
    } else {
      dim2 = ptrNode->Get("dim2")->As<UInt>();
    }
    UInt numEntries = dim1 * dim2;

    StdVector<std::string> rVals (numEntries), iVals(numEntries);
    rVals.Init("0.0");
    iVals.Init("0.0");
    if (ptrNode->Has("real")) {
      PtrParamNode tensor =  ptrNode->Get("real");
      ParamTools::AsStringTensor( tensor, numEntries, rVals );
    }

    if (ptrNode->Has("imag") ) {
      PtrParamNode tensor =  ptrNode->Get("imag");
      ParamTools::AsStringTensor( tensor, numEntries, iVals );
    }

    return CoefFunction::Generate( mp_, Global::COMPLEX, dim1, dim2, rVals, iVals );
  }


  void XMLMaterialHandler::
  ReadSquare3x3Tensor( PtrParamNode p, BaseMaterial *mat,
                       MaterialType isoPProp,
                       MaterialType* orthoProp,
                       MaterialType tensorProp,
                       Global::ComplexPart part)
  {
    BM::SymmetryType symType = BM::NOSYMMETRY;    
    if ( p->Has("tensor") ) {
      PtrParamNode pTensor = p->Get("tensor");
      PtrCoefFct tensorCoef = ReadTensor( pTensor, part );
      mat->SetCoefFct(tensorProp, tensorCoef);
      symType = BM::GENERAL;
    } //tensor

    else if ( p->Has("isotropic") ) {
      PtrCoefFct isoCoef = ReadScalar( p, "isotropic", part );
      mat->SetCoefFct(isoPProp, isoCoef);
      symType = BM::ISOTROPIC;
    } //isotropic

    else if ( p->Has("orthotropic") ) {
      PtrParamNode pOrtho = p->Get("orthotropic");
      PtrCoefFct orthoCoef;
      orthoCoef = ReadScalar( pOrtho, "value_1", part );
      mat->SetCoefFct(orthoProp[0], orthoCoef);
      orthoCoef = ReadScalar( pOrtho, "value_2", part );
      mat->SetCoefFct(orthoProp[1], orthoCoef);
      orthoCoef = ReadScalar( pOrtho, "value_3", part );
      mat->SetCoefFct(orthoProp[2], orthoCoef);
      symType = BM::ORTHOTROPIC;
    } //orthotropic

    else if ( p->Has("transversalIsotropic") ) {
      PtrParamNode pTransIso = p->Get("transversalIsotropic");
      PtrCoefFct transCoef;
      transCoef = ReadScalar( pTransIso, "value", part );
      mat->SetCoefFct(isoPProp, transCoef);
      transCoef = ReadScalar( pTransIso, "value_3", part );
      mat->SetCoefFct(orthoProp[2], transCoef);
      symType = BM::TRANS_ISOTROPIC;
    } //transversalIsotropic

    // Set symmetry type
    mat->SetSymmetryType( tensorProp, symType );

  }

  // Read a mechanical stiffness tensor
  BM::SymmetryType XMLMaterialHandler::ReadStiffnessTensor(PtrParamNode ptrNode,
                                                           Global::ComplexPart type,
                                                           BM::CoefMap &coefMap)
  {
    PtrCoefFct elastCoef;

    if (ptrNode->HasByVal("tensor", "dim1", "6") || ptrNode->HasByVal("tensor", "dim1", "3")) {
      PtrParamNode elastTensor = ptrNode->Get("tensor");
      elastCoef = ReadTensor( elastTensor, type );
      coefMap[MECH_STIFFNESS_TENSOR] = elastCoef;
      return BM::GENERAL;
    }

    if (ptrNode->Has("isotropic")) {
      PtrParamNode elastIso = ptrNode->Get("isotropic");

      if (elastIso->Has("elasticityModulus") && elastIso->Has("poissonNumber")) {
        elastCoef = ReadScalar(elastIso, "elasticityModulus", type );
        coefMap[MECH_EMODULUS] = elastCoef;
        elastCoef = ReadScalar(elastIso, "poissonNumber", type );
        coefMap[MECH_POISSON] = elastCoef;
      }
      else if (elastIso->Has("shearModulus") && elastIso->Has("compressionModulus")) {
        elastCoef = ReadScalar(elastIso, "shearModulus", type );
        coefMap[MECH_GMODULUS] = elastCoef;
        elastCoef = ReadScalar(elastIso, "compressionModulus", type );
        coefMap[MECH_KMODULUS] = elastCoef;
      }
      if (elastIso->Has("lameParameterMu") && elastIso->Has("lameParameterLambda")) {
        elastCoef = ReadScalar(elastIso, "lameParameterMu", type );
        coefMap[MECH_LAME_MU] = elastCoef;
        elastCoef = ReadScalar(elastIso, "lameParameterLambda", type );
        coefMap[MECH_LAME_LAMBDA] = elastCoef;
      }

      return BM::ISOTROPIC;
    }

    if (ptrNode->Has("transversalIsotropic")) {
      PtrParamNode transIso = ptrNode->Get("transversalIsotropic");

      elastCoef = ReadScalar(transIso, "elasticityModulus", type);
      coefMap[MECH_EMODULUS] = elastCoef;
      elastCoef = ReadScalar(transIso, "elasticityModulus_3", type);
      coefMap[MECH_EMODULUS_3] = elastCoef;

      elastCoef = ReadScalar(transIso, "poissonNumber", type);
      coefMap[MECH_POISSON] = elastCoef;
      elastCoef = ReadScalar(transIso, "poissonNumber_3", type);
      coefMap[MECH_POISSON_3] = elastCoef;

      elastCoef = ReadScalar(transIso, "shearModulus", type);
      coefMap[MECH_GMODULUS] = elastCoef;
      elastCoef = ReadScalar(transIso, "shearModulus_3", type);
      coefMap[MECH_GMODULUS_3] = elastCoef;

      return BM::TRANS_ISOTROPIC;
    }

    if (ptrNode->Has("orthotropic")) {
      PtrParamNode elastOrth = ptrNode->Get("orthotropic");

      elastCoef = ReadScalar(elastOrth, "elasticityModulus_1", type );
      coefMap[MECH_EMODULUS_1] = elastCoef;
      elastCoef = ReadScalar(elastOrth, "elasticityModulus_2", type );
      coefMap[MECH_EMODULUS_2] = elastCoef;
      elastCoef = ReadScalar(elastOrth, "elasticityModulus_3", type );
      coefMap[MECH_EMODULUS_3] = elastCoef;

      elastCoef = ReadScalar(elastOrth, "poissonNumber_12", type );
      coefMap[MECH_POISSON_12] = elastCoef;
      elastCoef = ReadScalar(elastOrth, "poissonNumber_13", type );
      coefMap[MECH_POISSON_13] = elastCoef;
      elastCoef = ReadScalar(elastOrth, "poissonNumber_23", type );
      coefMap[MECH_POISSON_23] = elastCoef;

      elastCoef = ReadScalar(elastOrth, "shearModulus_12", type );
      coefMap[MECH_GMODULUS_12] = elastCoef;
      elastCoef = ReadScalar(elastOrth, "shearModulus_13", type );
      coefMap[MECH_GMODULUS_13] = elastCoef;
      elastCoef = ReadScalar(elastOrth, "shearModulus_23", type );
      coefMap[MECH_GMODULUS_23] = elastCoef;

      return BM::ORTHOTROPIC;
    }

    return BM::NOSYMMETRY;
  }

  // Read the Kelvin-Voigt viscous tensor
  BM::SymmetryType XMLMaterialHandler::ReadViscousTensorKV(PtrParamNode ptrNode, Global::ComplexPart type, BM::CoefMap &coefMap)
  {
    PtrCoefFct coef;
    if (ptrNode->HasByVal("tensor", "dim1", "6") || ptrNode->HasByVal("tensor", "dim1", "3")) {
      PtrParamNode tensor = ptrNode->Get("tensor");
      coef = ReadTensor(tensor, type);
      coefMap[MECH_KV_VISCOUS_TENSOR] = coef;
      return BM::GENERAL;
    }

    if (ptrNode->Has("isotropic")) {
      EXCEPTION("'Isotropic' declaration of viscous tensor is currently not supported!");
    }

    if (ptrNode->Has("transversalIsotropic")) {
      EXCEPTION("'transversalIsotropic' declaration of viscous tensor is currently not supported!");
    }

    if (ptrNode->Has("orthotropic")) {
      EXCEPTION("'orthotropic' declaration of viscous tensor is currently not supported!");
    }

    return BM::NOSYMMETRY;
  }

  // Read a smooth stiffness tensor
  BM::SymmetryType XMLMaterialHandler::ReadStiffnessTensorSmooth(PtrParamNode ptrNode,
                                                                  Global::ComplexPart type,
                                                                  BM::CoefMap &coefMap)
  {
    PtrCoefFct elastCoef;

    if (ptrNode->HasByVal("tensor", "dim1", "6") || ptrNode->HasByVal("tensor", "dim1", "3")) {
      PtrParamNode elastTensor = ptrNode->Get("tensor");
      elastCoef = ReadTensor( elastTensor, type );
      coefMap[SMOOTH_STIFFNESS_TENSOR] = elastCoef;
      return BM::GENERAL;
    }

    if (ptrNode->Has("isotropic")) {
      PtrParamNode elastIso = ptrNode->Get("isotropic");

      if (elastIso->Has("elasticityModulus") && elastIso->Has("poissonNumber")) {
        elastCoef = ReadScalar(elastIso, "elasticityModulus", type );
        coefMap[SMOOTH_EMODULUS] = elastCoef;
        elastCoef = ReadScalar(elastIso, "poissonNumber", type );
        coefMap[SMOOTH_POISSON] = elastCoef;
      }

      return BM::ISOTROPIC;
    }

    return BM::NOSYMMETRY;
  }


  void XMLMaterialHandler::ReadXYValues( PtrParamNode paramNode,
                                         const std::string &xName,
                                         const std::string &yName,
                                         Vector<Double>& xValues,
                                         Vector<Double>& yValues ) {
    // H values
    UInt dimH, dimB;
    PtrParamNode hVal = paramNode->Get(xName);
    dimH = hVal->Get("dim1")->As<UInt>();

    Matrix<Double> hMatrix;
    ParamTools::AsTensor( hVal->Get("real"), dimH, 1, hMatrix );

    //B values
    PtrParamNode bVal = paramNode->Get(yName);
    dimB = bVal->Get("dim1")->As<UInt>();

    if ( dimH != dimB ) {
      EXCEPTION("Dimensions of " << xName << " and " << yName << " do not match");
    }

    Matrix<Double> bMatrix;
    ParamTools::AsTensor( bVal->Get("real"), dimB, 1, bMatrix );

    xValues.Resize(dimH);
    yValues.Resize(dimB);
    for (UInt i=0; i<dimH; i++) {
      xValues[i] = hMatrix[0][i];
      yValues[i] = bMatrix[0][i];
    }
  }

  void XMLMaterialHandler::ReadRayleighDamping(PtrParamNode paramNode, BaseMaterial *material)
  {
    // get the type of Rayleigh damping
    PtrParamNode rayleighTypeNode = paramNode->GetChild();
    std::string rayleighTypeName = rayleighTypeNode->GetName();
    // call the corresponding function to compute the Rayleigh coefficients
    if (rayleighTypeName == "coefficients") {
      // If both alpha and beta are set, we can use them directly
      material->SetCoefFct(RAYLEIGH_ALPHA, CoefFunction::Generate(mp_, Global::REAL, rayleighTypeNode->Get("alpha")->As<std::string>()));
      material->SetCoefFct(RAYLEIGH_BETA, CoefFunction::Generate(mp_, Global::REAL, rayleighTypeNode->Get("beta")->As<std::string>()));
    }
    else if (rayleighTypeName == "computeFromTwoPoints") {
      // compute alpha and beta from two points with specified loss tangens and frequency, where
      // 1/w * alpha + w * beta = tanDelta
      StdVector<double> frequencies(2);
      StdVector<double> lossTanDeltas(2);
      ParamNodeList twoPointNodes = rayleighTypeNode->GetChildren();
      double tanDelta1 = twoPointNodes[0]->Get("lossTangensDelta")->As<double>();
      double tanDelta2 = twoPointNodes[1]->Get("lossTangensDelta")->As<double>();
      double w1 = twoPointNodes[0]->Get("frequency")->As<double>();
      double w2 = twoPointNodes[1]->Get("frequency")->As<double>();
      w1 *= 2.0 * M_PI;
      w2 *= 2.0 * M_PI;
      double jDet = (w2 / w1 - w1 / w2);
      double a = (w2 * tanDelta1 - w1 * tanDelta2) / jDet;
      double b = (-tanDelta1 / w2 + tanDelta2 / w1) / jDet;
      std::string alpha = lexical_cast<std::string>(a);
      std::string beta = lexical_cast<std::string>(b);
      // set the coefFunctions for later use
      material->SetCoefFct(RAYLEIGH_ALPHA, CoefFunction::Generate(mp_, Global::REAL, alpha));
      material->SetCoefFct(RAYLEIGH_BETA, CoefFunction::Generate(mp_, Global::REAL, beta));
    }
    else if (rayleighTypeName == "computeFromTanDeltaAtFrequency") {
      // compute alpha and beta from a single point with a center frequency and a delta frequency
      // and a loss tangens delta that is specified on [f_c*(1+df), f_c*(1-df)]
      // 1/w * alpha + w * beta = tanDelta
      double tanDelta = rayleighTypeNode->Get("lossTangensDelta")->As<double>();
      double ratioDeltaF = rayleighTypeNode->Get("ratioDeltaF")->As<double>();
      double wc = rayleighTypeNode->Get("centerFreq")->As<double>();
      wc *= 2.0 * M_PI;
      double w1 = wc * (1 - ratioDeltaF);
      double w2 = wc * (1 + ratioDeltaF);
      double jDet = (w2 / w1 - w1 / w2);
      double a = tanDelta * (w2 - w1) / jDet;
      double b = tanDelta * (-1 / w2 + 1 / w1) / jDet;
      std::string alpha = lexical_cast<std::string>(a);
      std::string beta = lexical_cast<std::string>(b);
      // set the coefFunctions for later use
      material->SetCoefFct(RAYLEIGH_ALPHA, CoefFunction::Generate(mp_, Global::REAL, alpha));
      material->SetCoefFct(RAYLEIGH_BETA, CoefFunction::Generate(mp_, Global::REAL, beta));
    }
    else {
      EXCEPTION("Unknown Rayleigh damping type '" << rayleighTypeName << "'. Possible types are: 'coefficients', 'computeFromTwoPoints', 'computeFromTanDeltaAtFrequency'.");
    }
    // set the type of Rayleigh damping
    material->SetRayleighType(BaseMaterial::ALPHA_BETA);
  }

  void XMLMaterialHandler::ReadLossTanDeltaDamping(PtrParamNode paramNode, BaseMaterial *material)
  {
    material->SetString(paramNode->Get("value")->As<std::string>(), LOSS_TANGENS_DELTA);
    material->SetRayleighType(BaseMaterial::ADAPTED_LOSS_TANGENS);
  }

  void XMLMaterialHandler::ReadKelvinVoigtDamping(PtrParamNode paramNode, BaseMaterial *material)
  {
    PtrParamNode linearNode = paramNode->Get("linear");
    if (linearNode) {
      PtrParamNode tensorNode = linearNode->Get("tensor");
      if (tensorNode) {
        BM::SymmetryType symType = BM::NOSYMMETRY;
        BM::CoefMap coefMap;
        // read the linear viscous tensor
        symType = ReadViscousTensorKV(linearNode, Global::COMPLEX, coefMap);
        for (const auto &coefIt : coefMap) {
          material->SetCoefFct(coefIt.first, coefIt.second);
        }
        material->SetSymmetryType(MECH_KV_VISCOUS_TENSOR, symType);
      }
      else {
        EXCEPTION("Kelvin-Voigt viscous tensor must be given as 'tensor'");
      }
    }
    else {
      EXCEPTION("Kelvin-Voigt viscous tensor must be 'linear'");
    }
  }

  BaseMaterial::MatDescriptorNl
  XMLMaterialHandler::ReadNonlinDescriptor(PtrParamNode paramNode,
                                           BaseMaterial *material)
  {
    BaseMaterial::MatDescriptorNl info;
    info.approxType = NO_APPROX_TYPE;
    info.measAccuracy = 0.01;
    info.maxVal = 1000;
    info.fileName = "";
    info.factor = 1.0;

    // read dependency
    if (paramNode->Has("dependency")) {
      std::string dep = paramNode->Get("dependency")->As<std::string>();
      material->SetString(dep, NONLIN_DEPENDENCY);
    }

    // read approximation type
    if (paramNode->Has("approxType")) {
      std::string type = paramNode->Get("approxType")->As<std::string>();
      info.approxType = ApproxCurveTypeEnum.Parse(type);
    }

    // read measurement accuracy
    if (paramNode->Has("measAccuracy")) {
      info.measAccuracy = paramNode->Get("measAccuracy")->As<Double>();
    }

    // read maximum value for approximation
    if (paramNode->Has("maxApproxVal")) {
      info.maxVal = paramNode->Get("maxApproxVal")->As<Double>();
    }

    // read name of function file
    if (paramNode->Has("dataName")) {
      info.fileName = paramNode->Get("dataName")->As<std::string>().c_str();
    }

    // read factor
    if (paramNode->Has("factor")) {
      info.factor = paramNode->Get("factor")->As<Double>();
    }

    return info;
  }

  } // end of namespace
