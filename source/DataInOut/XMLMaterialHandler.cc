// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_xerces.hh>
#include "XMLMaterialHandler.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/programOptions.hh"

// header for materials
#include "Materials/electroMagneticMaterial.hh"
#include "Materials/electrostaticMaterial.hh"
#include "Materials/heatMaterial.hh"
#include "Materials/acousticMaterial.hh"
#include "Materials/mechanicMaterial.hh"
#include "Materials/piezoMaterial.hh"
#include "Materials/flowMaterial.hh"
#include "Materials/thermoelasticMaterial.hh"
#include "Materials/pyroelectricMaterial.hh"

// Note, that the methods ComputeIso/OrthoMechStiffnesTensor were commented out
// in revision 7562 and are not in the code -> check the repository!

namespace CoupledField {

  XMLMaterialHandler::XMLMaterialHandler( const std::string & fileName )
    : MaterialHandler( fileName) {

    parser_ = NULL;

    // Create a ParamNode and parse the material file
    std::string schema = progOpts->GetSchemaPathStr();
    schema += "/CFS-Material/CFS_Material.xsd";
  
    // Initialize our xerces dom parser to handle the  xml file
    Xerces* xerces = new Xerces(fileName, schema);

    parser_ = xerces->CreateParamNodeInstance();

    // release the xerces ressources, the parser_ is not affected
    delete xerces;
  }
  
  XMLMaterialHandler::~XMLMaterialHandler() {
      
    delete parser_;
  }
  
  BaseMaterial * XMLMaterialHandler::
  LoadMaterial( const std::string matName, 
               const MaterialClass matClass ) {

    BaseMaterial * material = NULL;
    
    std::string strMatClass;

    Enum2String(matClass,strMatClass);
    
    if(!parser_->Has("material", "name", matName))
      EXCEPTION("Cannot find material '" << matName << "'");
    
    // first get the material element:  <material name="iron">
    ParamNode* pn = NULL;
    pn = parser_->Get("material", "name", matName);
    
    if( !pn ) {
      EXCEPTION( "Material with name '" << matName 
                 << "' could not be found in material file!" );
    }
    // the the requested material class: <mechanical>
    pn = pn->Get(strMatClass);  
   
    if ( matClass == PIEZO ) {
      material = new PiezoMaterial();
      ReadPiezo( material, pn);
    }
    else if ( matClass == MECHANIC ) {
      material = new MechanicMaterial();
      ReadMechanic( material, pn );
    }    
    else if ( matClass == FLUID ) {
      material = new AcousticMaterial();
      ReadAcoustic( material, pn );
    }
    else if ( matClass == ELECTROMAGNETIC ) {
      material = new ElectroMagneticMaterial();
      ReadMagnetic( material, pn );
    }
    else if ( matClass == ELECTROSTATIC ) {
      material = new ElectroStaticMaterial();
      ReadElectrostatic( material, pn );
    }
    else if ( matClass == THERMIC ) {
      material = new HeatMaterial();
      ReadThermic( material, pn );
    }
    else if ( matClass == FLOW ) {
      material = new FlowMaterial();
      ReadFlow( material, pn );
    }
    else if ( matClass == PYROELECTRIC ) {
      material = new PyroelectricMaterial();
      ReadPyroelectric( material,pn );
    }
    else if ( matClass == THERMOELASTIC ) {
      material = new ThermoelasticMaterial();
      ReadThermoelastic( material, pn );
    }
    else {
      EXCEPTION( "material type:" << matClass << " not defined" );
    }
    // Finalize setup of material
    material->Finalize();

    return material;
  }

  //**********************************************************************
    //*************  READ PIEZO ********************************************
      //**********************************************************************

  void XMLMaterialHandler::ReadPiezo(BaseMaterial *material, ParamNode* pn) 
  {

    //read real piezo coupling tensor
    if(pn->Has("piezoCouplingTensor"))
    {
      Matrix<Double> couplingTensor(3,6);

      ParamNode* pct = pn->Get("piezoCouplingTensor");
      if(pct->Has("real"))
      {
        ParamTools::AsTensor<double>(pct->Get("real"), 3, 6, couplingTensor);
        material->SetTensor( couplingTensor, PIEZO_TENSOR, REAL );
      }
      if(pct->Has("imag"))
      {
        ParamTools::AsTensor<double>(pct->Get("imag"), 3, 6, couplingTensor);
        material->SetTensor( couplingTensor, PIEZO_TENSOR, IMAG );
      }
    } 

    //read nonlinearity of a coupling coefficient
    if(pn->Has("piezoCouplingCoefficient", "nonlinear", "function"))
    {
      ParamNode* pcc = pn->Get("piezoCouplingCoefficient", 
                               std::string("nonlinear"), 
                               "function");
      if(pcc->Has("entry"))
        material->SetScalar(pcc->Get("entry")->AsInt(), NONLIN_COEFFICIENT);
      
      if(pcc->Has("dependency"))
        material->SetScalar(pcc->Get("dependency")->AsString(), NONLIN_DEPENDENCY);
      
      if(pcc->Has("approxType"))
        material->SetScalar(pcc->Get("approxType")->AsString(), NONLIN_APPROXIMATION_TYPE);

      if(pcc->Has("dataName"))
        material->SetScalar(pcc->Get("dataName")->AsString(), NONLIN_DATA_NAME);
    }

    // Print material information to .info-file
    Info->PrintMaterial(material );
  }

//**********************************************************************
//*************  READ MECHANICS ****************************************
//**********************************************************************
  void XMLMaterialHandler::ReadMechanic(BaseMaterial *material, ParamNode* mech) 
  {
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

    

    //read material density
    if(mech->Has("density"))
      material->SetScalar(mech->Get("density")->AsDouble(), DENSITY, REAL);

    // quite a lot is elasitcity
    if(mech->Has("elasticity"))
    {
      ParamNode* elast = mech->Get("elasticity");

      if(elast->Has("tensor", "dim1", "6"))
      {
        ParamNode* tens = elast->Get("tensor", std::string("dim1"), "6");
        Matrix<Double> elasticityTensor(6,6);

        //read real elasticity tensor   
        if(tens->Has("real"))
        {
          ParamTools::AsTensor<double>(tens->Get("real"),6,6,elasticityTensor); 
          material->SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, REAL); 
          flagElastTensorReal = true;
        }
        if(tens->Has("imag"))
        {
          ParamTools::AsTensor<double>(tens->Get("imag"),6,6,elasticityTensor); 
          material->SetTensor( elasticityTensor, MECH_STIFFNESS_TENSOR, IMAG ); 
          flagElastTensorImag = true;
        }
      } // end tensor  
 
      // check values for isotropic      
      if(elast->Has("isotropic"))
      {
        // read the real part
        if(elast->Get("isotropic")->Has("real"))
        {
          ParamNode* real = elast->Get("isotropic")->Get("real");

          // read real elasticity modulus
          if(real->Has("elasticityModulus"))
          {
            material->SetScalar(real->Get("elasticityModulus")->AsDouble(), MECH_EMODULUS, REAL ); 
            flagEModulReal = true;
          }
          
          // read real Poisson number
          if(real->Has("poissonNumber"))
          {
            material->SetScalar(real->Get("poissonNumber")->AsDouble(), MECH_POISSON, REAL ); 
            flagPoissonReal = true;
          }
        }
        // read the imaginary part
        if(elast->Get("isotropic")->Has("imag"))
        {
          ParamNode* imag = elast->Get("isotropic")->Get("imag");
          
          //read imaginary elasticity modulus
          if(imag->Has("elasticityModulus"))
          {
            material->SetScalar(imag->Get("elasticityModulus")->AsDouble(), MECH_EMODULUS, IMAG ); 
            flagEModulImag = true;
          }

          // read imaginary Poisson number
          if(imag->Has("poissonNumber"))
          {
            material->SetScalar(imag->Get("poissonNumber")->AsDouble(), MECH_POISSON, IMAG ); 
            flagPoissonImag = true;
          }
        }
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
          
        ParamNode* real = elast->Get("orthotropic")->Get("real");
        
        //read orthotropic elasticity modulus
        if(real->Has("elasticityModulus_1"))
        {
          material->SetScalar(real->Get("elasticityModulus_1")->AsDouble(), MECH_EMODULUS_X, REAL ); 
          flagEModulXReal = true;
        }

        if(real->Has("elasticityModulus_2"))
        {
          material->SetScalar(real->Get("elasticityModulus_2")->AsDouble(), MECH_EMODULUS_Y, REAL ); 
          flagEModulYReal = true;
        }

        if(real->Has("elasticityModulus_3"))
        {
          material->SetScalar(real->Get("elasticityModulus_3")->AsDouble(), MECH_EMODULUS_Z, REAL ); 
          flagEModulZReal = true;
        }

        // read orthotropic Poisson numbers
        if(real->Has("poissonNumber_12"))
        {
          material->SetScalar(real->Get("poissonNumber_12")->AsDouble(), MECH_POISSON_XY, REAL ); 
          flagPoissonXYReal = true;
        }

        if(real->Has("poissonNumber_23"))
        {
          material->SetScalar(real->Get("poissonNumber_23")->AsDouble(), MECH_POISSON_YZ, REAL ); 
          flagPoissonYZReal = true;
        }

        if(real->Has("poissonNumber_13"))
        {
          material->SetScalar(real->Get("poissonNumber_13")->AsDouble(), MECH_POISSON_XZ, REAL ); 
          flagPoissonXZReal = true;
        }
    
        // read orthotropic shear modulus
        if(real->Has("shearModulus_23"))
        {
          material->SetScalar(real->Get("shearModulus_23")->AsDouble(), MECH_GMODULUS_YZ, REAL ); 
          flagShearModulYZReal = true;
        }

        if(real->Has("shearModulus_31"))
        {
          material->SetScalar(real->Get("shearModulus_31")->AsDouble(), MECH_GMODULUS_ZX, REAL ); 
          flagShearModulZXReal = true;
        }

        if(real->Has("shearModulus_12"))
        {
          material->SetScalar(real->Get("shearModulus_12")->AsDouble(), MECH_GMODULUS_XY, REAL ); 
          flagShearModulXYReal = true;
        }
      }  // orthotropic      
    } // end of elasticity


    if (flagEModulReal==true && 
        flagPoissonReal==true && 
        flagElastTensorReal==false) {
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
      EXCEPTION( "mechanical stiffness tensor is over determined.\n"
                 << " You specified the tensor as well as E-Modul "
                 << "and Poisson number" );
    }
    else if (flagEModulReal==false && 
             flagEModulXReal==false && 
             flagElastTensorReal==false) {
      EXCEPTION ("mechanical stiffness must be specified somehow." );
    }
    else {
      EXCEPTION( "mechanical stiffness tensor can not be computed." );
    }

    // elasticityCoefficient of type <elasticityCoefficient nonlinear="function">
    if(mech->Has("elasticityCoefficient", "nonlinear", "function"))
    {
      ParamNode* ec = mech->Get("elasticityCoefficient", 
                                std::string("nonlinear"),
                                "function");
      if(ec->Has("entry")) 
        material->SetScalar(ec->Get("entry")->AsInt(), NONLIN_COEFFICIENT);

     if(ec->Has("dependency"))
       material->SetScalar(ec->Get("dependency")->AsString(), NONLIN_DEPENDENCY );

     if(ec->Has("approxType"))
       material->SetScalar(ec->Get("approxType")->AsString(), NONLIN_APPROXIMATION_TYPE );

     if(ec->Has("dataName"))
       material->SetScalar(ec->Get("dataName")->AsString(), NONLIN_DATA_NAME );
    }; // end of <elasticityCoefficient nonlinear="function">


    //read coefficients for irreversible mechanical strain
    if(mech->Has("irreversibleStrainCoefficient"))
    {
      ParamNode* isc = mech->Get("irreversibleStrainCoefficient");
      // the dimension is only printed in the old param handler version 7562
      //if(isc->Has("dim")) std::cout << "dim=" << isc->Get("dim")->AsInt() << std::endl;
      
      if(isc->Has("coeffs"))
      {
        // read matrix
        Matrix<Double> matrixCoeffs(5,1);
        ParamTools::AsTensor<double>(isc->Get("coeffs"),1,5,matrixCoeffs); 

        // transform to vector
        Vector<Double> coeffs;
        coeffs.Resize( matrixCoeffs.GetSizeCol());
        for( UInt i=0; i<matrixCoeffs.GetSizeCol(); i++)
          coeffs[i] = matrixCoeffs[0][i];
          
        material->SetVector( coeffs, COEFF_STRAIN_IRREVERSIBLE, REAL ); 
      }
    } // end of irreversibleStrainCoefficient

    // read mechanical damping
    if(mech->Has("mechanicalDamping"))
    {
      // first rayleigh damping
      if(mech->Get("mechanicalDamping")->Has("rayleigh"))
      {
        ParamNode* r = mech->Get("mechanicalDamping")->Get("rayleigh");

        if(r->Has("alpha"))
         material->SetScalar(r->Get("alpha")->AsDouble(), RAYLEIGH_ALPHA, REAL);
         
        if(r->Has("beta"))
         material->SetScalar(r->Get("beta")->AsDouble(), RAYLEIGH_BETA, REAL);

        if(r->Has("lossTangensDelta"))
         material->SetScalar(r->Get("lossTangensDelta")->AsDouble(), LOSS_TANGENS_DELTA, REAL);

        if(r->Has("measuredFreq"))
         material->SetScalar(r->Get("measuredFreq")->AsDouble(), RAYLEIGH_FREQUENCY, REAL);
      }
      if(mech->Get("mechanicalDamping")->Has("fractional"))
      {
        ParamNode* f = mech->Get("mechanicalDamping")->Get("fractional");
        
        if(f->Has("alg"))        
          material->SetScalar(f->Get("alg")->AsString(), FRACTIONAL_ALG );
          
        if(f->Has("memory"))        
          material->SetScalar(f->Get("memory")->AsInt(), FRACTIONAL_MEMORY );
          
        if(f->Has("interpolation"))        
          material->SetScalar(f->Get("interpolation")->AsString(), FRACTIONAL_INTERPOL );
      }
    }

    // Print information to info file
    Info->PrintMaterial( material);
  }


//**********************************************************************
//*************  READ ACOUSTICS ****************************************
//**********************************************************************
  void XMLMaterialHandler::ReadAcoustic(BaseMaterial *material, ParamNode* acou)
  {
    //read density
    if(acou->Has("density"))
      material->SetScalar(acou->Get("density")->AsDouble(), DENSITY, REAL ); 
      
    //read compression modulus
    if(acou->Has("compressionModulus"))
      material->SetScalar(acou->Get("compressionModulus")->AsDouble(), ACOU_BULK_MODULUS, REAL );

    // check for acousticDamping
    if(acou->Has("acousticDamping"))
    {
      ParamNode* ad = acou->Get("acousticDamping");
      
      // check rayleigh
      if(ad->Has("rayleigh"))
      {
        ParamNode* r = ad->Get("rayleigh");
        
        if(r->Has("alpha"))
          material->SetScalar(r->Get("alpha")->AsDouble(), RAYLEIGH_ALPHA, REAL );
          
        if(r->Has("beta"))
          material->SetScalar(r->Get("beta")->AsDouble(), RAYLEIGH_BETA, REAL );

        if(r->Has("lossTangensDelta"))
          material->SetScalar(r->Get("lossTangensDelta")->AsDouble(), LOSS_TANGENS_DELTA, REAL );
       
        if(r->Has("measuredFreq"))
          material->SetScalar(r->Get("measuredFreq")->AsDouble(), RAYLEIGH_FREQUENCY, REAL );
      } // end of acousticDamping:rayleigh
      
      // read alpha0 of thermo viscous damping
      if(ad->Has("thermoViscous"))
      {
        if(ad->Get("thermoViscous")->Has("alpha0"))
          material->SetScalar(ad->Get("thermoViscous")->Get("alpha0")->AsDouble(), ACOU_ALPHA, REAL );
      }

      // read fractional damping
      if(ad->Has("fractional"))
      {
        ParamNode* f = ad->Get("fractional");
        
        if(f->Has("alpha0")) 
          material->SetScalar(f->Get("alpha0")->AsDouble(), ACOU_ALPHA, REAL );

        // read exponent of fractional damping      
        if(f->Has("y")) 
          material->SetScalar(f->Get("y")->AsDouble(), FRACTIONAL_EXPONENT, REAL );
      }
    } // end of acousticDamping

    //read acoustic non linearity
    if(acou->Has("acousticNonlinear"))
    {
      if(acou->Get("acousticNonlinear")->Has("bOverA"))
        material->SetScalar(acou->Get("acousticNonlinear")->Get("bOverA")->AsDouble(), BOVERA, REAL );
    }  

    // Print material information to info-file
    Info->PrintMaterial( material );
  }

//**********************************************************************
//*************  READ ELECTROSTATICS ************************************
//**********************************************************************
  void XMLMaterialHandler::ReadElectrostatic(BaseMaterial *material, ParamNode* elec)
  {
    // check for permittivity
    if(elec->Has("permittivity"))
    {
      ParamNode* p = elec->Get("permittivity");
      
      // check for tensor with dim1 = 3 <tensor dim1="3">
      if(p->Has("tensor", "dim1", "3"))
      {
        Matrix<Double> permittivityTensor(3,3);

        // read real permittivity tensor 
        if(p->Get("tensor", std::string("dim1"), "3")->Has("real"))
        {
          ParamNode* tensor =  p->Get("tensor",
                                      std::string("dim1"),
                                      "3")->Get("real");        
          ParamTools::AsTensor<double>(tensor, 3, 3, permittivityTensor);
          material->SetTensor(permittivityTensor, ELEC_PERMITTIVITY, REAL);
        }

        // read imaginary permittivity tensor
        if(p->Get("tensor", std::string("dim1"), "3")->Has("imag"))
        {
          ParamNode* tensor =  p->Get("tensor",
                                      std::string("dim1"),
                                      "3")->Get("imag");
          ParamTools::AsTensor<double>(tensor, 3, 3, permittivityTensor);
          material->SetTensor(permittivityTensor, ELEC_PERMITTIVITY, IMAG);
        }
      } // end of <tensor dim1="3">
      
      // check for isotropic permittivity
      // KILLME isotropic permittivity is NOT set in r7562!!
    } // end of permittivity
    
    // check for <permittivityCoefficient nonlinear="function">
    if(elec->Has("permittivityCoefficient", "nonlinear", "function"))
    {
      ParamNode* pc = elec->Get("permittivityCoefficient", 
                                std::string("nonlinear"),
                                "function");
      
      // read nonlinearity of a permittivity coefficient
      if(pc->Has("entry"))
        material->SetScalar(pc->Get("entry")->AsInt(), NONLIN_COEFFICIENT);
        
      // read non linear dependency of a permittivity coefficient
      if(pc->Has("dependency"))
        material->SetScalar(pc->Get("dependency")->AsString(), NONLIN_DEPENDENCY);

      // read non linear approxType of a permittivity coefficient
      if(pc->Has("approxType"))
        material->SetScalar(pc->Get("approxType")->AsString(), NONLIN_APPROXIMATION_TYPE);

      // read non linear data name of a permittivity coefficient        
      if(pc->Has("dataName"))
        material->SetScalar(pc->Get("dataName")->AsString(), NONLIN_DATA_NAME);
    } // end of permittivityCoefficient

    //read Preisach hysterese model
    if(elec->Has("hystModel"))
    {
      if(elec->Get("hystModel")->Has("preisach"))
      {
        ParamNode* p = elec->Get("hystModel")->Get("preisach");
        
        // force name
        material->SetScalar("preisach", HYST_MODEL);

        // read E saturation of Preisach hysterese model
        if(p->Has("eSat"))
          material->SetScalar(p->Get("eSat")->AsDouble(), X_SATURATION, REAL ); 
 
        // read P saturation of Preisach hysterese model
        if(p->Has("pSat"))
          material->SetScalar(p->Get("pSat")->AsDouble(), Y_SATURATION, REAL ); 

        // read P saturation of Preisach hysterese model
        if(p->Has("Pr"))
          material->SetScalar(p->Get("Pr")->AsDouble(), Y_REMANENCE, REAL ); 

        // read direction of polarization
        if(p->Has("dirP"))
        {
          int dir = p->Get("dirP")->AsInt();
          
          if(dir == 1) material->SetScalar("X", P_DIRECTION );
          if(dir == 2) material->SetScalar("Y", P_DIRECTION );
          if(dir == 3) material->SetScalar("Z", P_DIRECTION );
          
          if(dir != 1 && dir != 2 && dir != 3)
            EXCEPTION(dir << " is valid coordinate direction for electric preisach "
                      << " hysteresis model polarization");
        }
        
        // read weight dimension of Preisach hysterese model for weights
        int dim = -1;
        if(p->Has("dim")) dim = p->Get("dim")->AsInt();
    
        // read real permittivity tensor    
        if(p->Has("weights"))
        {
          Matrix<Double> preisachWeightTensor(dim,dim);
          ParamTools::AsTensor<double>(p->Get("weights"), dim, dim, preisachWeightTensor);
          material->SetTensor( preisachWeightTensor, PREISACH_WEIGHTS, REAL);
        }
      }
    }

    // Print information to info file
    Info->PrintMaterial( material );

  }

//**********************************************************************
//*************  READ MAGNETIC *****************************************
//**********************************************************************
  void XMLMaterialHandler::ReadMagnetic(BaseMaterial *material, ParamNode* mag)
  {
    // read electric conductivity
    if(mag->Has("electricConductivity"))
      material->SetScalar(mag->Get("electricConductivity")->AsDouble(), MAG_CONDUCTIVITY, REAL);
    
    // read magnetic permeability
    if(mag->Has("magneticPermeability"))
    {
      if(mag->Get("magneticPermeability")->Has("linear"))
      {
        ParamNode* lin = mag->Get("magneticPermeability")->Get("linear");
        double eps = 1e-10;
        
        if(lin->Has("isotropic"))
        {
          if(lin->Get("isotropic")->AsDouble() < eps)
            EXCEPTION("Magnetic permeability is near zero. Check material database");
          material->SetScalar(lin->Get("isotropic")->AsDouble(), MAG_PERMEABILITY, REAL );
        }
  
        if(lin->Has("orthotropic"))
        {
          ParamNode* ortho = lin->Get("orthotropic");
          bool permOrtho_1=false, permOrtho_2=false, permOrtho_3=false;
          
          if(ortho->Has("permeability_1"))
          {
            if(ortho->Get("permeability_1")->AsDouble() < eps)
              EXCEPTION("Magnetic permeability is near zero; Check material database");
            material->SetScalar(ortho->Get("permeability_1")->AsDouble(), MAG_PERMEABILITY_1, REAL); 
            permOrtho_1 = true;  
          }
          
          if(ortho->Has("permeability_2"))
          {
            if(ortho->Get("permeability_2")->AsDouble() < eps)
              EXCEPTION("Magnetic permeability is near zero; Check material database");
            material->SetScalar(ortho->Get("permeability_2")->AsDouble(), MAG_PERMEABILITY_2, REAL); 
            permOrtho_2 = true;  
          }
  
          if(ortho->Has("permeability_3"))
          {
            if(ortho->Get("permeability_3")->AsDouble() < eps)
              EXCEPTION("Magnetic permeability is near zero; Check material database");
            material->SetScalar(ortho->Get("permeability_3")->AsDouble(), MAG_PERMEABILITY_3, REAL); 
            permOrtho_3 = true;  
          }
        
          // check, if there is an orthotropic permeability!!
          if (permOrtho_1 == true && permOrtho_2 == true && permOrtho_3 == true)
            material->SetSymmetryType(BaseMaterial::ORTHOTROPIC);
        } // end of linear orthotropic
      } // end of linear

      // we know only nonlinear isotropic material
      if(mag->Get("magneticPermeability")->Has("nonlinear") && 
         mag->Get("magneticPermeability")->Get("nonlinear")->Has("isotropic"))
      {
        ParamNode* iso = mag->Get("magneticPermeability")->Get("nonlinear")->Get("isotropic");
        // In r7562  dependency and  approxType are not set in Material
        
        // read nonlinear approxType of magnetic permeability
        if(iso->Has("measAccuracy"))
          material->SetScalar(iso->Get("measAccuracy")->AsDouble(), DATA_ACCURACY, REAL );
                  
        // read nonlinear approxType of magnetic permeability
        if(iso->Has("maxApproxVal"))
          material->SetScalar(iso->Get("maxApproxVal")->AsDouble(), MAX_APPROX_VAL, REAL );

        // read nonlinear dataName of magnetic permeability
        if(iso->Has("dataName"))
          material->SetNonlinFileName(iso->Get("dataName")->AsString().c_str());
      } // nonlinear isotropic material   
    } // end of magneticPermeability  


    //read Preisach hysterese model
    if(mag->Has("hystModel"))
    {
      if(mag->Get("hystModel")->Has("preisach"))
      {
        ParamNode* p = mag->Get("hystModel")->Get("preisach");
        
        // force name
        material->SetScalar("preisach", HYST_MODEL);

        // read E saturation of Preisach hysterese model
        if(p->Has("eSat"))
          material->SetScalar(p->Get("eSat")->AsDouble(), X_SATURATION, REAL ); 
 
        // read P saturation of Preisach hysterese model
        if(p->Has("pSat"))
          material->SetScalar(p->Get("pSat")->AsDouble(), Y_SATURATION, REAL ); 

        // read direction of polarization
        if(p->Has("dirP"))
        {
          int dir = p->Get("dirP")->AsInt();
          
          if(dir == 1) material->SetScalar("X", P_DIRECTION );
          if(dir == 2) material->SetScalar("Y", P_DIRECTION );
          if(dir == 3) material->SetScalar("Z", P_DIRECTION );
          
          if(dir != 1 && dir != 2 && dir != 3)
            EXCEPTION(dir << " is valid coordinate direction for electric preisach "
                      << " hysteresis model polarization");
        }
        
        // read weight dimension of Preisach hysterese model for weights
        int dim = -1;
        if(p->Has("dim")) dim = p->Get("dim")->AsInt();
    
        // read real permittivity tensor    
        if(p->Has("weights"))
        {
          Matrix<Double> preisachWeightTensor(dim,dim);
          ParamTools::AsTensor<double>(p->Get("weights"), dim, dim, preisachWeightTensor);
          material->SetTensor( preisachWeightTensor, PREISACH_WEIGHTS, REAL);
        }
      }
    }


    // Print information to info file
    Info->PrintMaterial( material ); 
  }

//**********************************************************************
//*************  READ THERMIC ******************************************
//**********************************************************************
  void XMLMaterialHandler::ReadThermic(BaseMaterial *material, ParamNode* therm)
  {
    // read density
    if(therm->Has("density"))
      material->SetScalar(therm->Get("density")->AsDouble(), DENSITY, REAL);

    // read heat capacity
    if(therm->Has("heatCapacity"))
      material->SetScalar(therm->Get("heatCapacity")->AsDouble(), HEAT_CAPACITY, REAL);

    // read thermal conductivity
    if(therm->Has("thermalConductivity"))
      {
        ParamNode* thc = therm->Get("thermalConductivity");
        if(thc->Has("isotropic"))

          {
            material->SetScalar(thc->Get("isotropic")->AsDouble(), HEAT_CONDUCTIVITY, REAL);

          }
        else if(thc->Has("tensor"))
          {
            // can only be a real 3x3 tensor
            Matrix<double> tensor(3,3);
            ParamNode* tens_pn = thc->Get("tensor",
                                          std::string("dim1"),
                                          "3")->Get("real");

            ParamTools::AsTensor<double>(tens_pn, 3, 3, tensor);
            material->SetTensor(tensor, HEAT_CONDUCTIVITY_TENSOR, REAL);

          }
      }
    
    // Print information to info file
    Info->PrintMaterial( material );
  }
  
  //**********************************************************************
    //*************  READ FLOW *********************************************
      //**********************************************************************
    void XMLMaterialHandler::ReadFlow(BaseMaterial *material, ParamNode* flow)
    {    
      // read density
      if(flow->Has("density"))
        material->SetScalar(flow->Get("density")->AsDouble(), DENSITY, REAL);
      
      // dynamicViscosity is NOT set in r7562 
      
      // Print information to info file
      Info->PrintMaterial( material );
    }

  //**********************************************************************
    //*************  READ PYROELECTRIC *************************************
      //**********************************************************************
    void XMLMaterialHandler::ReadPyroelectric(BaseMaterial *material, 
                                              ParamNode* pyro){
      
      if (pyro->Has("pyrocoefficient")){
        ParamNode* py = pyro->Get("pyrocoefficient");
        if(py->Has("tensor"))
          {
            // can only be a real 3x3 tensor
            Matrix<double> tensor(3,3);
            ParamNode* tens_pn = py->Get("tensor",
                                         std::string("dim1"),
                                         "3")->Get("real");
            ParamTools::AsTensor<double>(tens_pn, 3, 3, tensor);
            material->SetTensor(tensor,PYROCOEFFICIENT_TENSOR,REAL);
          }
      }
      Info->PrintMaterial( material );
    }

  //**********************************************************************
    //*************  READ THERMOELASTIC ************************************
      //**********************************************************************
    void XMLMaterialHandler::ReadThermoelastic(BaseMaterial *material,
                                               ParamNode* thermExp) {

      if(thermExp->Has("thermalExpansion")){
        ParamNode* te = thermExp->Get("thermalExpansion");
        if(te->Has("tensor"))
          {
            // can only be a real 3x3 tensor
            Matrix<double> tensor(3,3);
            ParamNode* tens_pn = te->Get("tensor",
                                         std::string("dim1"),
                                         "3")->Get("real");
            ParamTools::AsTensor<double>(tens_pn, 3, 3, tensor);
            material->SetTensor(tensor,THERMAL_EXPANSION_TENSOR,REAL);
          }
      }
      Info->PrintMaterial( material );
    }
}

  //   //Double      doubValue;
//     //Integer     inteValue, dim;
//     std::string striValue;
//     Matrix<Double> Pyrocoefficient_Tensor(3,3);

//     // Construct vectors for restricted search parameter
//     StdVector<std::string> keyVec;
//     StdVector<std::string> attrVec;
//     StdVector<std::string> valVec;

//     //read real permittivity tensor
//     const unsigned int dim1=3, dim2=3;
//     keyVec = "material","pyroelectric","pyrocoefficient","tensor","real";
//     attrVec= "name"    ,""        ,""            ,"dim1";
//     valVec =  matName  ,""        ,""            ,"3";
//     if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
//       parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
//                                    dim1, dim2, Pyrocoefficient_Tensor );
//       material->SetTensor( Pyrocoefficient_Tensor, PYROCOEFFICIENT_TENSOR, REAL ); 
//            // std::cerr << "real Pyrocoefficient_Tensor=" << std::endl << Pyrocoefficient_Tensor << std::endl;
//     }

//     //read imaginary permittivity tensor
//     keyVec = "material","pyroelectric","pyrocoefficient","tensor","imag";
//     attrVec= "name"    ,""        ,""            ,"dim1";
//     valVec =  matName  ,""        ,""            ,"3";
//     if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
//       parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
//                                    dim1, dim2, Pyrocoefficient_Tensor );
//       material->SetTensor( Pyrocoefficient_Tensor, PYROCOEFFICIENT_TENSOR, IMAG ); 
//       // std::cerr << "imaginary permittivityTensor=" << std::endl << permittivityTensor << std::endl;
//     }
 
//     // Print information to info file
//     Info->PrintMaterial( material );
//  }

//**********************************************************************
//*************  READ THERMOELASTIC ************************************
//**********************************************************************
//   void XMLMaterialHandler::ReadThermoelastic(BaseMaterial *material,
//                                     const std::string matName) {
//     Double      doubValue;
// 
//     // Construct vectors for restricted search parameter
//     StdVector<std::string> keyVec;
//     StdVector<std::string> attrVec;
//     StdVector<std::string> valVec;
// 
//     //read thermal expansion
//     keyVec = "material","thermoelastic","thermalExpansion";
//     attrVec= "name"    ,"";
//     valVec =  matName  ,"";
//     if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
//       parser_->Get( keyVec, attrVec, valVec, doubValue );
//       material->SetScalar( doubValue, THERMAL_EXPANSION, REAL );
//        //std::cerr << "thermalExpansion=" << doubValue << std::endl;
//     }



//     std::string striValue;
//     Matrix<Double> thermalExpansion_Tensor(3,3);

//     // Construct vectors for restricted search parameter
//     StdVector<std::string> keyVec;
//     StdVector<std::string> attrVec;
//     StdVector<std::string> valVec;

//     //read real permittivity tensor
//     const unsigned int dim1=3, dim2=3;
//     keyVec = "material","thermoelastic","thermalExpansion","tensor","real";
//     attrVec= "name"    ,""        ,""            ,"dim1";
//     valVec =  matName  ,""        ,""            ,"3";
//     if (parser_->ContainElem( keyVec, attrVec, valVec ) ) {
//       parser_->GetDim1xDim2Tensor( keyVec, attrVec, valVec, 
//                                    dim1, dim2, thermalExpansion_Tensor );
//       material->SetTensor( thermalExpansion_Tensor, THERMAL_EXPANSION_TENSOR, REAL ); 
//             //std::cerr << "real thermalExpansion_Tensor=" << std::endl << thermalExpansion_Tensor << std::endl;
//     }


//     // Print information to info file
//     Info->PrintMaterial( material );
//   }

  //}
