#include <fstream>

#include "MagEdgeMixedAHPDE.hh"

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
DEFINE_LOG(magEdgeMixedAHPde, "magEdgeMixedAHPde")


  // **************
  //  Constructor
  // **************
  MagEdgeMixedAHPDE::MagEdgeMixedAHPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdgeMixedAH";
    pdematerialclass_ = ELECTROMAGNETIC;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgeMixedAHPDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));

  }


  // *************
  //  Destructor
  // *************
  MagEdgeMixedAHPDE::~MagEdgeMixedAHPDE() {
  }

  // ********************************************************************************
  //  METHOD THAT CALLS ALL METHODS FOR THE DEFINITION OF INTEGRATORS
  // ********************************************************************************
  void MagEdgeMixedAHPDE::DefineIntegrators() {
    this->DefineStandardIntegrators();
    DefineCoilIntegrators();
  }

  // ********************************************************************************
  //  DEFINITION OF STIFFNESS INTEGRATORS
  // ********************************************************************************
  void MagEdgeMixedAHPDE::DefineStandardIntegrators() {


    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // get FEFunction and space
    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> magIntensFeFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> magVecPotFeSpace = magVecPotFeFunc->GetFeSpace();
    shared_ptr<FeSpace> magIntensFeSpace = magIntensFeFunc->GetFeSpace();
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

      // Set the FE-Ansatz for the current region for the vector quantities
      // (magnetic vector potential and magnetic field intensity)
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string VecPolyId = curRegNode->Get("polyId")->As<std::string>();
      std::string VecIntegId = curRegNode->Get("integId")->As<std::string>();
      magVecPotFeSpace->SetRegionApproximation(actRegion, VecPolyId, VecIntegId);
      magIntensFeSpace->SetRegionApproximation(actRegion, VecPolyId, VecIntegId);

      // Pass entitylist to fespace / fefunction for magnetic vector potential and magnetic field intensity
      magVecPotFeFunc->AddEntityList( actSDList );
      magIntensFeFunc->AddEntityList( actSDList );


      // ====================================================================
      // STIFFNESS INTEGRATOR (start)
      // (I) : (mu h, h') - (a, curlh') = 0          for all h' in H(curl)
      // (II): (curlh, a')              = (j,a')     for all a' in H(curl)
      // ====================================================================

      // upper left matrix: (mu N, N')
      // ==============================
      PtrCoefFct mu = NULL;
      mu = actMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);

      BaseBDBInt *stiffUpperLeft = NULL;
      stiffUpperLeft = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,1.0,updatedGeo_);
      stiffUpperLeft->SetName("NNIntegrator");
      BiLinFormContext * stiffUpperLeftContext = new BiLinFormContext(stiffUpperLeft, STIFFNESS );
      stiffUpperLeftContext->SetEntities( actSDList, actSDList );
      stiffUpperLeftContext->SetFeFunctions( magIntensFeFunc, magIntensFeFunc );
      assemble_->AddBiLinearForm( stiffUpperLeftContext );

      // upper right matrix: (N, curlN')
      // ==============================
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      PtrCoefFct minus_constOne = CoefFunction::Generate( mp_, Global::REAL, "-1.0");
      BiLinearForm* stiffUpperRight = NULL;
      stiffUpperRight = new ABInt<>(new CurlOperator<FeHCurl,3,Double>(),
                                    new IdentityOperator<FeHCurl,3,1,Double>(),
                                    minus_constOne, 1.0, updatedGeo_);
      stiffUpperRight->SetName("NCurlNIntegrator");
      BiLinFormContext * stiffUpperRightContext = new BiLinFormContext(stiffUpperRight, STIFFNESS );
      stiffUpperRightContext->SetEntities( actSDList, actSDList );
      stiffUpperRightContext->SetFeFunctions( magVecPotFeFunc, magIntensFeFunc );
      assemble_->AddBiLinearForm( stiffUpperRightContext );

      // lower left matrix: (curlN, N')
      // ==============================
      BiLinearForm* stiffLowerLeft = NULL;
      stiffLowerLeft = new ABInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                                    new CurlOperator<FeHCurl,3,Double>(),
                                    constOne, 1.0, updatedGeo_);
      stiffLowerLeft->SetName("CurlNNIntegrator");
      BiLinFormContext * stiffLowerLeftContext = new BiLinFormContext(stiffLowerLeft, STIFFNESS );
      stiffLowerLeftContext->SetEntities( actSDList, actSDList );
      stiffLowerLeftContext->SetFeFunctions( magIntensFeFunc, magVecPotFeFunc );
      assemble_->AddBiLinearForm( stiffLowerLeftContext );
      // ====================================================================
      // STIFFNESS INTEGRATOR (end)
      // ====================================================================
    } // end for regions
  } // end DefineIntegrators


  LinearForm* MagEdgeMixedAHPDE::GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId)
  {
    LinearForm * ret = NULL;
    ret = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3,1, Double>(), factor, coef, updatedGeo_);
    return ret;
  }


  // ********************************************************************************
  //  DEFINITION OF RHS INTEGRATORS (that are coming solely from coils in this PDE)
  // ********************************************************************************
  void MagEdgeMixedAHPDE::DefineCoilIntegrators()
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
          // (I) : 0
          // (II): (j, a')
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
          // (I) : 0
          // (II): (j, a')
          // ====================================================================
        } // loop: parts
      }
    } // loop: coils
  }



  // ********************************************************************************
  // TIME-STEPPING SECTION
  // ********************************************************************************
  void MagEdgeMixedAHPDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    GLMScheme * scheme2 = new Trapezoidal(gamma);
    shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(scheme2, 0, nlType) );
    feFunctions_[MAG_FIELD_INTENSITY]->SetTimeScheme(myScheme2);
  }

  // ********************************************************************************
  // PRIMARY RESULTS:
  // 1) h on the edges of the mesh
  // 1) a (mag. vector poterntial) on the edges of the mesh
  // ********************************************************************************

  void MagEdgeMixedAHPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;
    potInfo->SetFeFunction(feFunctions_[MAG_POTENTIAL]);

    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL], potInfo );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";

    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<BaseFeFunction> scalFct = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<ResultInfo> res2(new ResultInfo);
    res2->resultType = MAG_FIELD_INTENSITY;
    res2->dofNames = vecComponents;
    res2->unit = "A/m";
    res2->definedOn = ResultInfo::ELEMENT;
    res2->entryType = ResultInfo::VECTOR;
    res2->SetFeFunction(feFunctions_[MAG_FIELD_INTENSITY]);

    feFunctions_[MAG_FIELD_INTENSITY]->SetResultInfo(res2);
    DefineFieldResult( feFunctions_[MAG_FIELD_INTENSITY], res2 );

  }

  // ********************************************************************************
  // POSTPROCESSING RESULTS: not needed, since h is already the source field hs
  //                         and curla is not necessary yet.
  // ********************************************************************************
  void MagEdgeMixedAHPDE::DefinePostProcResults() {
  }

  // ********************************************************************************
  // FINALIZE POSTPROCESSING RESULTS
  // ********************************************************************************
  void MagEdgeMixedAHPDE::FinalizePostProcResults() {
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();
    }


  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeMixedAHPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    PtrParamNode magVecPpotSpaceNode = infoNode->Get("magPotential");
    crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, magVecPpotSpaceNode, FeSpace::HCURL, ptGrid_ );
    crSpaces[MAG_POTENTIAL]->Init(solStrat_);

    PtrParamNode VecSpaceNode = infoNode->Get("magFieldIntensity");
    crSpaces[MAG_FIELD_INTENSITY] = FeSpace::CreateInstance(myParam_, VecSpaceNode, FeSpace::HCURL, ptGrid_ );
    crSpaces[MAG_FIELD_INTENSITY]->Init(solStrat_);

    return crSpaces;
  }

} // end of namespace

