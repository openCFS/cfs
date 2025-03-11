#include <fstream>

#include "MagEdgeMixedSFGPDE.hh"

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
DEFINE_LOG(MagEdgeMixedSFGPDE, "MagEdgeMixedSFGPDE")


  // **************
  //  Constructor
  // **************
  MagEdgeMixedSFGPDE::MagEdgeMixedSFGPDE( Grid * aptgrid, PtrParamNode paramNode,
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
      EXCEPTION("MagEdgeMixedSFGPDE is just implemented for 3D setups!");
  }


  // *************
  //  Destructor
  // *************
  MagEdgeMixedSFGPDE::~MagEdgeMixedSFGPDE() {
  }

  // ********************************************************************************
  //  METHOD THAT CALLS ALL METHODS FOR THE DEFINITION OF INTEGRATORS
  // ********************************************************************************
  void MagEdgeMixedSFGPDE::DefineIntegrators() {
    this->DefineStandardIntegrators();
    DefineCoilIntegrators();
  }




  // ********************************************************************************
  //  DEFINITION OF STIFFNESS INTEGRATORS
  // ********************************************************************************
  void MagEdgeMixedSFGPDE::DefineStandardIntegrators() {

    RegionIdType actRegion;

    // get FEFunction and space
    shared_ptr<BaseFeFunction> VecFeFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<BaseFeFunction> ScaFeFunc = feFunctions_[LAGRANGE_MULT];
    shared_ptr<FeSpace> VecFeSpace = VecFeFunc->GetFeSpace();
    shared_ptr<FeSpace> ScaFeSpace = ScaFeFunc->GetFeSpace();

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      // set current region and materials
      actRegion = regions_[iRegion];

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


      // ====================================================================
      // STIFFNESS INTEGRATOR (start)
      // (I) : curlh curlh' + gradu h'
      // (II): h gradu'
      // ====================================================================

      // upper left matrix: curlN curlN
      // ==============================
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      BaseBDBInt* stiffUpperLeft = NULL;
      stiffUpperLeft = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), constOne, 1.0, updatedGeo_) ;
      stiffUpperLeft->SetName("CurlNCurlNIntegrator");
      BiLinFormContext * stiffUpperLeftContext = new BiLinFormContext(stiffUpperLeft, STIFFNESS );
      stiffUpperLeftContext->SetEntities( actSDList, actSDList );
      stiffUpperLeftContext->SetFeFunctions( VecFeFunc, VecFeFunc );
      assemble_->AddBiLinearForm( stiffUpperLeftContext );

      // upper right matrix: gradV N
      // ==============================
      BiLinearForm* stiffUpperRight = NULL;
      stiffUpperRight = new ABInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                                    new GradientOperator<FeH1,3,1,Double>(),
                                    constOne, 1.0, updatedGeo_);
      stiffUpperRight->SetName("GradVNIntegrator");
      BiLinFormContext * stiffUpperRightContext = new BiLinFormContext(stiffUpperRight, STIFFNESS );
      stiffUpperRightContext->SetEntities( actSDList, actSDList );
      stiffUpperRightContext->SetFeFunctions( VecFeFunc, ScaFeFunc );
      assemble_->AddBiLinearForm( stiffUpperRightContext );

      // lower left matrix: N gradV
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


  LinearForm* MagEdgeMixedSFGPDE::GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId)
  {
    LinearForm * ret = NULL;
    ret = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(), factor, coef, updatedGeo_);
    return ret;
  }

  // ********************************************************************************
  //  DEFINITION OF RHS INTEGRATORS (that are coming solely from coils in this PDE)
  // ********************************************************************************
  void MagEdgeMixedSFGPDE::DefineCoilIntegrators()
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
          // (I) : j curlh'
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
          // (I) : j curlh'
          // (II): 0
          // ====================================================================
        } // loop: parts
      }
    } // loop: coils
  }





  // ********************************************************************************
  // TIME-STEPPING SECTION
  // ********************************************************************************
  void MagEdgeMixedSFGPDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_FIELD_INTENSITY]->SetTimeScheme(myScheme);

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
  void MagEdgeMixedSFGPDE::DefinePrimaryResults() {

      StdVector<std::string> vecComponents;
      vecComponents = "x", "y", "z";

      // === MAGNETIC RESULT FIELD INTENSITY ===
      shared_ptr<ResultInfo> potInfo(new ResultInfo);
      potInfo->resultType = MAG_FIELD_INTENSITY;
      potInfo->SetFeFunction(feFunctions_[MAG_FIELD_INTENSITY]);
      potInfo->dofNames = vecComponents;
      potInfo->unit = "A/m";
      potInfo->definedOn = ResultInfo::ELEMENT;
      potInfo->entryType = ResultInfo::VECTOR;
      feFunctions_[MAG_FIELD_INTENSITY]->SetResultInfo(potInfo);
      DefineFieldResult( feFunctions_[MAG_FIELD_INTENSITY], potInfo );
      hdbcSolNameMap_[MAG_FIELD_INTENSITY] = "fluxParallel";


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
      // -----------------------------------
      //  Define xml-names of Dirichlet BCs
      // -----------------------------------
      hdbcSolNameMap_[LAGRANGE_MULT] = "fluxParallel";

  }

  // ********************************************************************************
  // POSTPROCESSING RESULTS: not needed, since h is already the source field hs
  // ********************************************************************************
  void MagEdgeMixedSFGPDE::DefinePostProcResults() {
  }

  // ********************************************************************************
  // FINALIZE POSTPROCESSING RESULTS
  // ********************************************************************************
  void MagEdgeMixedSFGPDE::FinalizePostProcResults() {
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();
    }

  // ********************************************************************************
  // CREATE FE-SPACES
  // ********************************************************************************
  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeMixedSFGPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    PtrParamNode VecSpaceNode = infoNode->Get("magFieldIntensity");
    crSpaces[MAG_FIELD_INTENSITY] = FeSpace::CreateInstance(myParam_, VecSpaceNode, FeSpace::HCURL, ptGrid_ );
    crSpaces[MAG_FIELD_INTENSITY]->Init(solStrat_);

    PtrParamNode ScaSpaceNode = infoNode->Get("lagrangeMultiplier");
    crSpaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, ScaSpaceNode, FeSpace::H1, ptGrid_);
    crSpaces[LAGRANGE_MULT]->Init(solStrat_);
    return crSpaces;
  }
} // end of namespace

