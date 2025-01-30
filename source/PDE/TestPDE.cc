
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "TestPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Utils/StdVector.hh"

#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
//#include "Materials/Models/VectorPreisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"


#include "Driver/Assemble.hh"

//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/LaplOp.hh"
#include "Forms/Operators/IdentityOperator.hh"

//new postprocessing concept
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {

  DEFINE_LOG(testpde, "pde.test")

  // ======================================================
  // SET SOLUTION INFORMATION
  // ======================================================
  TestPDE::TestPDE(Grid * aptgrid, PtrParamNode paramNode,
                   PtrParamNode infoNode,
                   shared_ptr<SimState> simState, Domain* domain)
  :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    pdename_           = "testPDE";
    pdematerialclass_  = TESTMAT;
    nonLin_            = false;
  
    //! Always use updated Lagrangian formulation 
    updatedGeo_        = true;
  }



  std::map<SolutionType, shared_ptr<FeSpace> >
  TestPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("testDof");
      crSpaces[TEST_DOF] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[TEST_DOF]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of TEST PDE is not known!");
    }
    return crSpaces;
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void TestPDE::InitNonLin() {

    SinglePDE::InitNonLin();

  }


  void TestPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;  

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;

    //get FEFunction and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[TEST_DOF];
    shared_ptr<FeSpace> mySpace = feFunc->GetFeSpace();
  
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];
    
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
    
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
    
      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
    
      // pass entitylist of fespace / fefunction
      feFunc->AddEntityList( actSDList );

      // ====================================================================
      // stiffness integrator
      // ====================================================================
      PtrCoefFct beta = actSDMat->GetScalCoefFnc( TEST_BETA, Global::REAL );

      BaseBDBInt* stiffInt = NULL;
      if( dim_ == 2 ) {
        stiffInt = new BBInt<>(new GradientOperator<FeH1,2>(), beta,1.0, updatedGeo_ );
      } else {
        stiffInt = new BBInt<>(new GradientOperator<FeH1,3>(), beta,1.0, updatedGeo_ );
      }
      stiffInt->SetName("StiffnessIntegrator");

      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunc,feFunc);
      stiffInt->SetFeSpace( mySpace );

      assemble_->AddBiLinearForm( stiffIntDescr );
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

      // ====================================================================
      // mass integrator
      // ====================================================================
      PtrCoefFct alpha = actSDMat->GetScalCoefFnc( TEST_ALPHA, Global::REAL );
      BiLinearForm *massInt = NULL;
      if(dim_==2)
        massInt = new BBInt<>(new IdentityOperator<FeH1,2,1,Double>(), alpha,1.0, updatedGeo_ );
      else
        massInt = new BBInt<>(new IdentityOperator<FeH1,3,1,Double>(), alpha,1.0, updatedGeo_ );

      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace( mySpace );

      BiLinFormContext *massContext =  new BiLinFormContext(massInt, DAMPING );
    
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunc,feFunc);
      assemble_->AddBiLinearForm( massContext );
    }

  }

  void TestPDE::DefineRhsLoadIntegrators() {
  
    LOG_TRACE(testpde) << "Defining rhs load integrators for test PDE";
  
    // Get FESpace and FeFunction of electric potential
    shared_ptr<BaseFeFunction> myFct = feFunctions_[TEST_DOF];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
  
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dofNames;
  
  
    bool coefUpdateGeo = true;
    // =====================
    //  TEST SOURCE DENSITY
    // =====================
    LOG_DBG(testpde) << "Reading source density of test PDE";
  
    ReadRhsExcitation( "testSourceDensity", dofNames,
                       ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Test source density must be defined on elements")
          }
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();
    
      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                         Complex(1.0), coef[i], coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                                        1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("TestSourceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for
  }


  void TestPDE::DefineSolveStep() {

    solveStep_ = new StdSolveStep(*this);

  }

  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void TestPDE::InitTimeStepping()
  {
    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    GLMScheme* scheme = nullptr;
    if (integrationScheme)
    {
      scheme = GetXmlDefinedScheme(integrationScheme);
    }
    else
    {
      // Until now no effective mass formulation in the trapezoidal
      // integration scheme is implemented!
      scheme = new Trapezoidal(1.0);
    }

    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0));
    feFunctions_[TEST_DOF]->SetTimeScheme(myScheme);
  }


  void TestPDE::DefinePrimaryResults() {

    // === TEMPERATURE ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = TEST_DOF;
    
    res1->dofNames = "";
    res1->unit = "?";
    res1->definedOn = ResultInfo::MapSolTypeToDefinedOn(TEST_DOF);
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[TEST_DOF]->SetResultInfo(res1);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    res1->SetFeFunction(feFunctions_[TEST_DOF]);
    DefineFieldResult( feFunctions_[TEST_DOF], res1 );



    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[TEST_DOF] = "testGround";
    idbcSolNameMap_[TEST_DOF] = "testPotential";

    // === TEST RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = TEST_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    availResults_.insert( rhs );
    rhsFeFunctions_[TEST_DOF]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[TEST_DOF], rhs );
  }

  void TestPDE::DefinePostProcResults() {

    // === TEST FIELD ===
    StdVector<std::string> vecComponents;
    if( dim_ == 3 ) {
      vecComponents = "x", "y", "z";
    }
    else if( isaxi_ ) {
      vecComponents = "r", "z";
    }
    else {
      vecComponents = "x", "y";
    }
    shared_ptr<BaseFeFunction> feFct = feFunctions_[TEST_DOF];
    shared_ptr<ResultInfo> field(new ResultInfo);
    field->resultType = TEST_FIELD;
    field->dofNames = vecComponents;
    field->unit = "??";
    field->definedOn = ResultInfo::MapSolTypeToDefinedOn(TEST_FIELD);
    field->entryType = ResultInfo::VECTOR;
    availResults_.insert( field );
    shared_ptr<CoefFunctionFormBased> fieldFunc;
    if( isComplex_ ) {
      fieldFunc.reset(new CoefFunctionBOp<Complex>(feFct, field));
    } else {
      fieldFunc.reset(new CoefFunctionBOp<Double>(feFct, field));
    }
    DefineFieldResult( fieldFunc, field );
    stiffFormCoefs_.insert(fieldFunc);
  }
  
} // end of namespace CoupledField
