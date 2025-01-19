#include <fstream>

#include "MagEdgeRegulSFGPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"

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
DEFINE_LOG(MagEdgeRegulSFGPDE, "MagEdgeRegulSFGPDE")


  // **************
  //  Constructor
  // **************
  MagEdgeRegulSFGPDE::MagEdgeRegulSFGPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdgeRegulSFG";
    pdematerialclass_ = ELECTROMAGNETIC;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    // check if we have a 2d setup
    // bool is2d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "plane";
    // if ( !is2d )
    //   EXCEPTION("MagEdgeRegulSFGPDE is just implemented for 2D setups!");

    mu0_ = 4.0 * M_PI * 1e-7;
    rho0_ = mu0_ * std::pow(10.0, (5.0/2.0));
  }


  // *************
  //  Destructor
  // *************
  MagEdgeRegulSFGPDE::~MagEdgeRegulSFGPDE() {
  }

  // ********************************************************************************
  //  METHOD THAT CALLS ALL METHODS FOR THE DEFINITION OF INTEGRATORS
  // ********************************************************************************
  void MagEdgeRegulSFGPDE::DefineIntegrators() {
    this->DefineStandardIntegrators();
    DefineCoilIntegrators();
  }




  // ********************************************************************************
  //  DEFINITION OF STIFFNESS INTEGRATORS
  // ********************************************************************************
  void MagEdgeRegulSFGPDE::DefineStandardIntegrators() {

    RegionIdType actRegion;

    // get FEFunction and space
    shared_ptr<BaseFeFunction> VecFeFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> VecFeSpace = VecFeFunc->GetFeSpace();

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

      // Pass entitylist to fespace / fefunction for magnetic vector and electric scalar potential
      VecFeFunc->AddEntityList( actSDList );


      // left matrix part: mu0 N N
      // ==============================
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      BaseBDBInt* stiffUpperLeft = NULL;
      if(is3d_){
        stiffUpperLeft = new BBInt<>(new  IdentityOperator<FeHCurl,3, 1, Double>(), constOne, mu0_, updatedGeo_) ;
      }else{
        stiffUpperLeft = new BBInt<>(new  IdentityOperator<FeHCurl,2, 1, Double>(), constOne, mu0_, updatedGeo_) ;
      }  
      stiffUpperLeft->SetName("NNIntegrator");
      BiLinFormContext * stiffUpperLeftContext = new BiLinFormContext(stiffUpperLeft, STIFFNESS );
      stiffUpperLeftContext->SetEntities( actSDList, actSDList );
      stiffUpperLeftContext->SetFeFunctions( VecFeFunc, VecFeFunc );
      assemble_->AddBiLinearForm( stiffUpperLeftContext );

      // right matrix part: rho curlN curlN
      // ==============================
      BiLinearForm* stiffUpperRight = NULL;
      if(is3d_){
        stiffUpperRight = new BBInt<>(new CurlOperator<FeHCurl,3,Double>(),
                                    constOne, rho0_, updatedGeo_);
      }else{
        stiffUpperRight = new BBInt<>(new CurlOperator<FeHCurl,2,Double>(),
                                    constOne, rho0_, updatedGeo_);
      } 
      
      stiffUpperRight->SetName("CurlNCurlNIntegrator");
      BiLinFormContext * stiffUpperRightContext = new BiLinFormContext(stiffUpperRight, STIFFNESS );
      stiffUpperRightContext->SetEntities( actSDList, actSDList );
      stiffUpperRightContext->SetFeFunctions( VecFeFunc, VecFeFunc );
      assemble_->AddBiLinearForm( stiffUpperRightContext );

      // ====================================================================
      // STIFFNESS INTEGRATOR (end)
      // ====================================================================
    } // end for regions
  } // end DefineIntegrators


  LinearForm* MagEdgeRegulSFGPDE::GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId)
  {
    LinearForm * ret = NULL;
    if(is3d_){
      ret = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(), factor, coef, updatedGeo_);
    }else{
      ret = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 2, Double>(), factor, coef, updatedGeo_);
    }
    
    return ret;
  }

  // ********************************************************************************
  //  DEFINITION OF RHS INTEGRATORS (that are coming solely from coils in this PDE)
  // ********************************************************************************
  void MagEdgeRegulSFGPDE::DefineCoilIntegrators()
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
          // ====================================================================
          coilCurrentDens_[actRegion] = jFct[0];
          curInt = GetCurrentDensityInt( rho0_, jFct[0] );
          curInt->SetName("CoilIntegrator");
          LinearFormContext * coilContext = new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );
          // ====================================================================
          // RHS INTEGRATOR (end)
          // (I) : j curlh'
          // ====================================================================
        } // loop: parts
      }
    } // loop: coils
  }





  // ********************************************************************************
  // TIME-STEPPING SECTION
  // ********************************************************************************
  void MagEdgeRegulSFGPDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_FIELD_INTENSITY]->SetTimeScheme(myScheme);
  }

  // ********************************************************************************
  // PRIMARY RESULTS:
  // 1) h on the edges of the mesh
  // 1) u lagrange multiplier on the nodes of the mesh
  // ********************************************************************************
  void MagEdgeRegulSFGPDE::DefinePrimaryResults() {

      StdVector<std::string> vecComponents;
      if(is3d_){
        vecComponents = "x", "y", "z";
      }else{
        vecComponents = "x", "y";
      }
      

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
      hdbcSolNameMap_[MAG_FIELD_INTENSITY] = "fluxParallel"; // actually flux-normal!!!!
  }

  // ********************************************************************************
  // POSTPROCESSING RESULTS: not needed, since h is already the source field hs
  // ********************************************************************************
  void MagEdgeRegulSFGPDE::DefinePostProcResults() {
  }

  // ********************************************************************************
  // FINALIZE POSTPROCESSING RESULTS
  // ********************************************************************************
  void MagEdgeRegulSFGPDE::FinalizePostProcResults() {
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();
    }

  // ********************************************************************************
  // CREATE FE-SPACES
  // ********************************************************************************
  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeRegulSFGPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    PtrParamNode VecSpaceNode = infoNode->Get("magFieldIntensity");
    crSpaces[MAG_FIELD_INTENSITY] = FeSpace::CreateInstance(myParam_, VecSpaceNode, FeSpace::HCURL, ptGrid_ );
    crSpaces[MAG_FIELD_INTENSITY]->Init(solStrat_);
    return crSpaces;
  }
} // end of namespace

