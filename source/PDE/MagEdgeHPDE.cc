#include <fstream>

#include "MagEdgeHPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionHarmBalance.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"



// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/GradientOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/MultiHarmonicDriver.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(magEdgeHPde, "magEdgeHPde")


  // **************
  //  Constructor
// **************
   MagEdgeHPDE::MagEdgeHPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticHEdge";
    pdematerialclass_ = ELECTROMAGNETIC;
    formulation_ = MagBasePDE::EDGE;
    nonLin_ = false;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    // default is false
    useGradFields_ = paramNode->Get("useGradientFields")->As<bool>();


    if ( !is3d_ )
      EXCEPTION("MagEdgeHPDE is just implemented for 3D setups!");
  }


  // *************
  //  Destructor
  // *************
  MagEdgeHPDE::~MagEdgeHPDE() {
  }

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagEdgeHPDE::InitNonLin() {
    SinglePDE::InitNonLin();
  }

  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgeHPDE::DefineIntegrators() {
    this->DefineStandardIntegrators();
    // in MagBasePDE
    DefineCoilIntegrators();
  }


  LinearForm* MagEdgeHPDE::GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId)
  {
    LinearForm * ret = NULL;
    ret = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(), factor, coef, updatedGeo_);
    return ret;
  }


  void MagEdgeHPDE::DefineCoilIntegrators()
  {

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    std::map<Coil::IdType, shared_ptr<Coil>>::iterator coilIt;
    coilIt = coils_.begin();
    for (; coilIt != coils_.end(); coilIt++)
    {
      Coil &actCoil = *(coilIt->second);
      // run over all parts
      std::map<RegionIdType, shared_ptr<Coil::Part>>::iterator partIt;
      partIt = actCoil.parts_.begin();

      if (actCoil.sourceType_ != Coil::CURRENT)
      {
        EXCEPTION("Only current excitation is implemented in H-formulation!!");
      }
      else
      {

        for (partIt = actCoil.parts_.begin();
             partIt != actCoil.parts_.end();
             partIt++)
        {
          Coil::Part &actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
          actSDList->SetRegion(actRegion);
          LinearForm* curInt = NULL;


          std::map<UInt, PtrCoefFct> jFct;
          CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_, CoefXpr::OP_MULT);
          PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

          // ----------- IMPORTANT (start) -------------------------------------
          // Here the value for the current density vector j is determined.
          // The linearform of the formulation is as follows
          //         (rho_art j, curl(h')) ,
          // where rho_art is an artificial resistivity that occurs due to the
          // penalty method to end up in this formulation.
          // A common way to calculate this coefficient is
          //         rho_art = mu*10^(n/2) ,
          // where mu is the permeability in this region and n is parameter
          // that should be chosen rather small. However, since coils are usually 
          // in regions where mu ~ mu0 this factor is harcoded and only n can be 
          // altered.  Anyway, this should work also for the nonlinear hysteretic case,
          // since the coil region is usually linear material.
          // (Honestly, I did not know how to implement it differently)
          // get the penaltyParameter from the .xml file and set the default value
          Double penaltyParameter = 5;
          myParam_->GetValue("penaltyFunctionParameter", penaltyParameter, ParamNode::PASS);
          CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect*(1/(1.256637061e-06*std::pow(10,penaltyParameter/2)))), CoefXpr::OP_DIV);
          jFct[0] = CoefFunction::Generate(mp_, part, jVec);
          // ----------- IMPORTANT (end) -------------------------------------


          coilCurrentDens_[actRegion] = jFct[0];
          curInt = GetCurrentDensityInt( 1.0, jFct[0] );
          curInt->SetName("CoilIntegrator");
          LinearFormContext * coilContext = new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );

        } // loop: parts
      }
    } // loop: coils
  }


  void MagEdgeHPDE::DefineStandardIntegrators(){

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // get the penaltyParameter from the .xml file and set the default value
    Double penaltyParameter = 5;
    myParam_->GetValue("penaltyFunctionParameter", penaltyParameter, ParamNode::PASS);

    std::cout << penaltyParameter << std::endl;

    // get FEFunction and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY]; // MAG_POTENTIAL
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
        // set current region and materials
        actRegion = regions_[iRegion];
        actMat    = materials_[actRegion];
        // get current region name
        std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

        // create new entity list
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( actRegion );

        // Set the FE-Ansatz for the current region
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
        std::string polyId = curRegNode->Get("polyId")->As<std::string>();
        std::string integId = curRegNode->Get("integId")->As<std::string>();
        feSpace->SetRegionApproximation(actRegion,polyId,integId);

        // pass entitylist ot fespace / fefunction
        feFunc->AddEntityList( actSDList );

        // ====================================================================
        // STIFFNESS INTEGRATOR [START] (curlN, curlN)
        // ====================================================================
        // get material coefficient (in this case: artificial conductivity that regularizes the PDE and is defined as 1/sigma_artificial = mu*10^(n/2)
        // where n is an integer number)
        PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "100");
        //mu_regularize = actMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
        Double mu_regularize;
        materials_[actRegion]->GetScalar( mu_regularize, MAG_PERMEABILITY_SCALAR, Global::REAL );
        PtrCoefFct rho;
        rho = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*std::pow(10,penaltyParameter/2)));

        BaseBDBInt* curlcurl = NULL;
        curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), rho,1.0, updatedGeo_);
        curlcurl->SetName("CurlCurlIntegrator");

        BiLinFormContext * stiffContext = new BiLinFormContext(curlcurl, STIFFNESS );
        stiffContext->SetEntities( actSDList, actSDList );
        stiffContext->SetFeFunctions( feFunc, feFunc );
        assemble_->AddBiLinearForm( stiffContext );

        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,curlcurl) );
        // ====================================================================
        // STIFFNESS INTEGRATOR [END]
        // ====================================================================
        // ====================================================================
        // MASS INTEGRATOR [START] (N, N)
        // ====================================================================
        PtrCoefFct mu = NULL;
        mu = actMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);

        BaseBDBInt *massInt = NULL;
        //massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),mu,1.0);
        massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,1.0,updatedGeo_);
        massInt->SetName("MassIntegrator");

        BiLinFormContext * massContext = new BiLinFormContext(massInt, STIFFNESS );
        massContext->SetEntities( actSDList, actSDList );
        massContext->SetFeFunctions( feFunc, feFunc );
        assemble_->AddBiLinearForm( massContext );
        // insert mass integrator to list of defined mass integrators
        massInts_[actRegion] = massInt;
        // ====================================================================
        // MASS INTEGRATOR [END]
        // ====================================================================
        }// end for regions
    }


    void MagEdgeHPDE::DefinePrimaryResults() {

      StdVector<std::string> vecComponents;
      vecComponents = "x", "y", "z";

      // === MAGNETIC RESULT FIELD INTENSITY ===
      shared_ptr<ResultInfo> potInfo(new ResultInfo);
      potInfo->resultType = MAG_FIELD_INTENSITY;
      potInfo->dofNames = vecComponents;
      potInfo->unit = "A/m";
      potInfo->definedOn = ResultInfo::ELEMENT;
      potInfo->entryType = ResultInfo::VECTOR;

      feFunctions_[MAG_FIELD_INTENSITY]->SetResultInfo(potInfo);
      DefineFieldResult( feFunctions_[MAG_FIELD_INTENSITY], potInfo );

    }

    void MagEdgeHPDE::DefinePostProcResults() {

      StdVector<std::string> vecComponents;
      vecComponents = "x", "y", "z";

      // === MAGNETIC FLUX DENSITY ===
      shared_ptr<ResultInfo> magFlux(new ResultInfo);
      magFlux->resultType = MAG_FLUX_DENSITY;
      magFlux->SetVectorDOFs(dim_, isaxi_);
      magFlux->dofNames = vecComponents;
      magFlux->unit = "Vs/m^2";
      magFlux->definedOn = ResultInfo::ELEMENT;
      magFlux->entryType = ResultInfo::VECTOR;

    }

    void MagEdgeHPDE::FinalizePostProcResults() {
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();
    }


    std::map<SolutionType, shared_ptr<FeSpace> >
    MagEdgeHPDE::CreateFeSpaces(const std::string& formulation,
                                PtrParamNode infoNode ) {
        //ok default case so we create grid based approximation H1 elements
        //and standard Gauss integration
        std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
        if(formulation == "default" || formulation == "H_CURL"){
        PtrParamNode potSpaceNode = infoNode->Get("magFieldIntensity");
        crSpaces[MAG_FIELD_INTENSITY] =
            FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_ );
        crSpaces[MAG_FIELD_INTENSITY]->Init(solStrat_);
        }else{
        EXCEPTION("The formulation " << formulation
                    << "of magnetic edge PDE is not known!");
        }
        return crSpaces;
    }


} // end of namespace
