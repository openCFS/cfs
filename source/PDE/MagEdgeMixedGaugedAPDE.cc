#include <fstream>

#include "MagEdgeMixedGaugedAPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "FeBasis/H1/FeSpaceH1Hi.hh"
#include "FeBasis/H1/H1Elems.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"


// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(MagEdgeMixedGaugedAPDE, "MagEdgeMixedGaugedAPDE")


  // **************
  //  Constructor
  // **************
  MagEdgeMixedGaugedAPDE::MagEdgeMixedGaugedAPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdgeMixedSFG";
    pdematerialclass_ = ELECTROMAGNETIC;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgeMixedGaugedAPDE is just implemented for 3D setups!");
  }


  // *************
  //  Destructor
  // *************
  MagEdgeMixedGaugedAPDE::~MagEdgeMixedGaugedAPDE() {
  }

  // ********************************************************************************
  //  METHOD THAT CALLS ALL METHODS FOR THE DEFINITION OF INTEGRATORS
  // ********************************************************************************
  void MagEdgeMixedGaugedAPDE::DefineIntegrators() {
    this->DefineStandardIntegrators();
    DefineCoilIntegrators();
  }




  // ********************************************************************************
  //  DEFINITION OF STIFFNESS INTEGRATORS
  // ********************************************************************************
  void MagEdgeMixedGaugedAPDE::DefineStandardIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // get FEFunction and space
    shared_ptr<BaseFeFunction> VecFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> ScaFeFunc = feFunctions_[LAGRANGE_MULT];
    shared_ptr<FeSpace> VecFeSpace = VecFeFunc->GetFeSpace();
    shared_ptr<FeSpace> ScaFeSpace = ScaFeFunc->GetFeSpace();
    PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      // set current region and materials
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];

      // get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // Create new entity list, based on the region
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // Set the FE-Ansatz for the current region for the vector quantitiy
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string VecPolyId = curRegNode->Get("VecPolyId")->As<std::string>();
      std::string VecIntegId = curRegNode->Get("VecIntegId")->As<std::string>();
      VecFeSpace->SetRegionApproximation(actRegion, VecPolyId, VecIntegId);

      // Set the FE-Ansatz for the current region for the scalar quantitiy
      std::string ScaPolyId = curRegNode->Get("ScaPolyId")->As<std::string>();
      std::string ScaIntegId = curRegNode->Get("ScaIntegId")->As<std::string>();
      ScaFeSpace->SetRegionApproximation(actRegion, ScaPolyId, ScaIntegId);
      
      // Pass entitylist to fespace / fefunction for magnetic vector and electric scalar potential
      VecFeFunc->AddEntityList( actSDList );
      ScaFeFunc->AddEntityList( actSDList );


      // ========================================================================
      // STIFFNESS INTEGRATOR (start)
      // (I) : (nu curla, curla') + (gradu, a') = (j,a')    for all a' in H(curl)
      // (II): (a, gradu')                      = 0         for all u' in H1
      // ========================================================================

      // upper left matrix: (nu curla, curla')
      // ==============================
      PtrCoefFct nu = NULL;
      nu = actMat->GetScalCoefFnc( MAG_RELUCTIVITY_SCALAR, Global::REAL); // get magnetic reluctivity
      BaseBDBInt* stiffUpperLeft = NULL;
      stiffUpperLeft = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), nu, 1.0, updatedGeo_) ;
      stiffUpperLeft->SetName("CurlNCurlNIntegrator");
      BiLinFormContext * stiffUpperLeftContext = new BiLinFormContext(stiffUpperLeft, STIFFNESS );
      stiffUpperLeftContext->SetEntities( actSDList, actSDList );
      stiffUpperLeftContext->SetFeFunctions( VecFeFunc, VecFeFunc );
      assemble_->AddBiLinearForm( stiffUpperLeftContext );

      // upper right matrix: (gradu, a')
      // ==============================
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      BiLinearForm* stiffUpperRight = NULL;
      stiffUpperRight = new ABInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                                    new GradientOperator<FeH1,3,1,Double>(),
                                    constOne, 1.0, updatedGeo_);
      stiffUpperRight->SetName("GradVNIntegrator");
      BiLinFormContext * stiffUpperRightContext = new BiLinFormContext(stiffUpperRight, STIFFNESS );
      stiffUpperRightContext->SetEntities( actSDList, actSDList );
      stiffUpperRightContext->SetFeFunctions( VecFeFunc, ScaFeFunc );
      assemble_->AddBiLinearForm( stiffUpperRightContext );

      // lower left matrix: (a, gradu')
      // ==============================
      BiLinearForm* stiffLowerLeft = NULL;
      stiffLowerLeft = new ABInt<>(new GradientOperator<FeH1,3,1,Double>(),
                                    new IdentityOperator<FeHCurl,3,1,Double>(),
                                    constOne, 1.0, updatedGeo_);
      stiffLowerLeft->SetName("GradVIdentityAIntegratorLowerLeft");
      BiLinFormContext * stiffLowerLeftContext = new BiLinFormContext(stiffLowerLeft, STIFFNESS );
      stiffLowerLeftContext->SetEntities( actSDList, actSDList );
      stiffLowerLeftContext->SetFeFunctions(ScaFeFunc, VecFeFunc  );
      assemble_->AddBiLinearForm( stiffLowerLeftContext );
      // ====================================================================
      // STIFFNESS INTEGRATOR (end)
      // ====================================================================
    } // end for regions
  } // end DefineIntegrators


  LinearForm* MagEdgeMixedGaugedAPDE::GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId)
  {
    LinearForm * ret = NULL;
    ret = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3,1, Double>(), factor, coef, updatedGeo_);
    return ret;
  }

  // ********************************************************************************
  //  DEFINITION OF RHS INTEGRATORS (that are coming solely from coils in this PDE)
  // ********************************************************************************
  void MagEdgeMixedGaugedAPDE::DefineCoilIntegrators()
  {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
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
        EXCEPTION("Only current excitation is implemented in mixedSFG formulation");
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

          CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
          jFct[0] = CoefFunction::Generate(mp_, part, jVec);

          // ====================================================================
          // RHS INTEGRATOR (start)
          // (I) : (j, a')
          // (II): 0
          // ====================================================================
          coilCurrentDens_[actRegion] = jFct[0];
          curInt = GetCurrentDensityInt( 1.0, jFct[0] );
          curInt->SetName("CoilIntegrator");
          LinearFormContext * coilContext = new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );
          // ====================================================================
          // RHS INTEGRATOR (end)
          // ====================================================================
        } // loop: parts
      }
    } // loop: coils
  }





  // ********************************************************************************
  // TIME-STEPPING SECTION
  // ********************************************************************************
  void MagEdgeMixedGaugedAPDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    // Important: Create a new time scheme just for the elec potential unknowns, as otherwise the
    // size of the vectors does not match!
    GLMScheme * scheme2 = new Trapezoidal(gamma);
    shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(scheme2, 0, nlType) );
    feFunctions_[LAGRANGE_MULT]->SetTimeScheme(myScheme2);
  }

  // ********************************************************************************
  // PRIMARY RESULTS:
  // 1) h on the edges of the mesh
  // 1) u lagrange multiplier on the nodes of the mesh
  // ********************************************************************************
  void MagEdgeMixedGaugedAPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;

    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL], potInfo );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";


    // === LAGRANGE MULTIPLIER ===
    shared_ptr<BaseFeFunction> scalFct = feFunctions_[LAGRANGE_MULT];
    shared_ptr<ResultInfo> res2(new ResultInfo);
    res2->resultType = LAGRANGE_MULT;
    res2->dofNames = "";
    res2->unit = "V";
    res2->definedOn = ResultInfo::NODE;
    res2->entryType = ResultInfo::SCALAR;
    res2->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
    results_.Push_back( res2 );
    availResults_.insert( res2 );
    scalFct->SetResultInfo(res2);
    DefineFieldResult( scalFct, res2 );
  }

  // ********************************************************************************
  // POSTPROCESSING RESULTS:
  // - MAG_FLUX_DENSTIY
  // - MAG_FIELD_INTENSITY
  // ********************************************************************************
  void MagEdgeMixedGaugedAPDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    if(!fluxDensityDefined_){
      DefineMagFluxDensity();
    }
    PtrCoefFct bFunc = this->GetCoefFct(MAG_FLUX_DENSITY);

/*     // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> magIntens(new ResultInfo);
    magIntens->resultType = MAG_FIELD_INTENSITY;
    magIntens->SetVectorDOFs(dim_, isaxi_);
    magIntens->dofNames = vecComponents;
    magIntens->unit = "A/m";
    magIntens->definedOn = ResultInfo::ELEMENT;
    magIntens->entryType = ResultInfo::VECTOR;
    DefineFieldResult( fieldIntensity_, magIntens );
    availResults_.insert( magIntens ); */
  }

  // ********************************************************************************
  // FINALIZE POSTPROCESSING RESULTS
  // ********************************************************************************
  void MagEdgeMixedGaugedAPDE::FinalizePostProcResults() {
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();
    }

  // ********************************************************************************
  // CREATE FE-SPACES
  // ********************************************************************************
  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeMixedGaugedAPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    PtrParamNode VecSpaceNode = infoNode->Get("magFieldIntensity");
    crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, VecSpaceNode, FeSpace::HCURL, ptGrid_ );
    crSpaces[MAG_POTENTIAL]->Init(solStrat_);

    PtrParamNode ScaSpaceNode = infoNode->Get("lagrangeMultiplier");
    crSpaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, ScaSpaceNode, FeSpace::H1, ptGrid_);
    crSpaces[LAGRANGE_MULT]->Init(solStrat_);
    return crSpaces;
  }
} // end of namespace
