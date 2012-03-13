
#include "AcousticPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/ConvectivePierceOperator.hh"

#include "FeBasis/FeFunctions.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

#include "Domain/Results/ResultFunctor.hh"
#include "Domain/Results/ExternalFieldFunctors.hh"

#include "Driver/Assemble.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

  DECLARE_LOG(acousticpde)
   DEFINE_LOG(acousticpde, "pde.acoustic")


  AcousticPDE::AcousticPDE( Grid* aGrid, PtrParamNode paramNode )
              : SinglePDE( aGrid, paramNode ){

    pdename_           = "acoustic";
    pdematerialclass_  = FLUID;
    maxTimeDerivOrder_ = 2;
    nonLin_            = false;
    InitialCondition_  = 0.0;

    //check for pressure or potential formulation
    std::string pdeFormulation = myParam_->Get("formulation")->As<std::string>();
    if(pdeFormulation == "default"){
      formulation_ = ACOU_PRESSURE;
    }else{
      formulation_ = SolutionTypeEnum.Parse(pdeFormulation);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticPDE::CreateFeSpaces( const std::string&  formulation,
                  PtrParamNode infoNode ){

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      std::string form = SolutionTypeEnum.ToString(formulation_);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[formulation_] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptgrid_);
      crSpaces[formulation_]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of acoustic PDE is not known!");
    }
    return crSpaces;
  }

  void AcousticPDE::DefineIntegrators(){

    RegionIdType actRegion;
    // BaseMaterial * actSDMat = NULL;

    //type of geometry
    std::string geometryType;
    param->Get("domain")->GetValue("geometryType", geometryType );

    // convert to tensor type
    // SubTensorType tensorType = FULL;
    if (geometryType == "plane") {
      // tensorType = PLANE_STRAIN;
    } else if (geometryType == "axi") {
      // tensorType = AXI;
      isaxi_ = true;
    }

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      // actSDMat = it->second;

      // Get current region name
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);

      Double density=0.0,compressibility=0.0,c0=0.0;

      materials_[actRegion]->GetScalar(density,DENSITY,Global::REAL);
      materials_[actRegion]->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);
      c0 = std::sqrt(compressibility/density);

      // ====================================================================
      // standard stiffness integrator
      // ====================================================================
      shared_ptr<CoefFunction> coeffK
                = CoefFunction::Generate(Global::REAL, "1.0");

      BiLinearForm * stiffInt = NULL;
      if( dim_ == 2 ) {
        shared_ptr<CoefFunction> coeffK
          = CoefFunction::Generate(Global::REAL, "1.0");
        stiffInt = new BBInt<GradientOperator<FeH1,2>, Double >(coeffK,1.0 );
      }
      else{
        shared_ptr<CoefFunction> coeffK
          = CoefFunction::Generate(Global::REAL, "1.0");
        stiffInt = new BBInt<GradientOperator<FeH1,3>, Double >(coeffK,1.0 );
      }
      stiffInt->SetName("LaplaceIntegrator");

      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[formulation_]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
      stiffInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );

      // ====================================================================
      // standard mass integrator
      // ====================================================================

      BiLinearForm *massInt = NULL;

      shared_ptr<CoefFunction> coeffM
        = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0/(c0*c0)));
      if(dim_==2)
        massInt = new BBInt<IdentityOperator<FeH1,2,1,Double>, Double  >(coeffM,1.0 );
      else
        massInt = new BBInt<IdentityOperator<FeH1,3,1,Double>, Double  >(coeffM,1.0 );

      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace() );

      BiLinFormContext *massContext =  new BiLinFormContext(massInt, MASS );

      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
      assemble_->AddBiLinearForm( massContext );

      // ====================================================================
      // check for flow (Pierce equation)
      // ====================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if(flowId != "") {
        if ( formulation_ != ACOU_POTENTIAL )
          EXCEPTION("Pierce-Equation just possible in velocity potential formulation" );

        //Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow","name",flowId.c_str());
        if(isComplex_){
          shared_ptr<FieldFunctor<Complex> > fct = dynamic_pointer_cast<FieldFunctor<Complex> >(meanFlowFunctor_);
          fct->AddRegion(actRegion,flowNode);
        }else{
          shared_ptr<FieldFunctor<Double> > fct = dynamic_pointer_cast<FieldFunctor<Double> >(meanFlowFunctor_);
          fct->AddRegion(actRegion,flowNode);
        }
        std::cout << "Do Flow" << std::endl;

        //now create the integrators
        BiLinearForm *convectiveStiff = NULL;
        BiLinearForm *convectiveDamp = NULL;
        if( dim_ == 2 ) {
          convectiveDamp  = new ABInt<IdentityOperator<FeH1,2,1>,ConvectiveOperator<FeH1,2,1> >(coeffM, 2.0);
          convectiveStiff = new ABInt<ConvectivePierceOperator<FeH1,2,1>,ConvectiveOperator<FeH1,2,1> >(coeffM, -1.0);
        } else {
          convectiveDamp  = new ABInt<IdentityOperator<FeH1,3,1>,ConvectiveOperator<FeH1,3,1> >(coeffM, 2.0);
          convectiveStiff = new ABInt<ConvectivePierceOperator<FeH1,3,1>,ConvectiveOperator<FeH1,3,1> >(coeffM, -1.0);
        }
        convectiveDamp->SetBCoefFunctionOpB(meanFlowFunctor_);
        convectiveDamp->SetName("convectiveDampPierce");
        convectiveStiff->SetBCoefFunctionOpB(meanFlowFunctor_);
        convectiveStiff->SetName("convectiveStiffPierce");

        BiLinFormContext *convectiveContextDamp  =  new BiLinFormContext(convectiveDamp, DAMPING );
        BiLinFormContext *convectiveContextStiff =  new BiLinFormContext(convectiveDamp, STIFFNESS );

        convectiveContextDamp->SetEntities( actSDList, actSDList );
        convectiveContextDamp->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
        convectiveContextStiff->SetEntities( actSDList, actSDList );
        convectiveContextStiff->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);

        assemble_->AddBiLinearForm( convectiveContextDamp );
        assemble_->AddBiLinearForm( convectiveContextStiff );
      }
    }
  }

  void AcousticPDE::DefineSurfaceIntegrators( ){
    //========================================================================================
    // ABC boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );

      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType sRegId = ptgrid_->GetRegion().Parse(regionName);

        // --- Set the FE ansatz for the current region ---
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());
        std::string polyId = curRegNode->Get("polyId")->As<std::string>();
        std::string integId = curRegNode->Get("integId")->As<std::string>();
        feFunctions_[formulation_]->GetFeSpace()->SetRegionApproximation(sRegId, polyId,integId);

        RegionIdType aRegion = ptgrid_->GetRegion().Parse(volRegName);
        //
        Double density=0.0,compressibility=0.0,c0;

        //
        materials_[aRegion]->GetScalar(density,DENSITY,Global::REAL);
        materials_[aRegion]->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);
        c0 = std::sqrt(compressibility/density);

        shared_ptr<CoefFunction> coeffDamp
                  = CoefFunction::Generate(Global::REAL, "1.0");

        BiLinearForm * abcInt = NULL;

        if( dim_ == 2 ) {
          abcInt = new BBInt<IdentityOperator<FeH1,2,1> >(coeffDamp,1.0/c0 );
        } else {
          abcInt = new BBInt<IdentityOperator<FeH1,3,1> >(coeffDamp,1.0/c0 );
        }

        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, DAMPING );

        abcContext->SetEntities( actSDList, actSDList );
        abcContext->SetFeFunctions( feFunctions_[formulation_] , feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }
    }


  }


  void AcousticPDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticPDE::SetInitialCondition(){

  }

  void AcousticPDE::CalcResults( shared_ptr<BaseResult> result ){
    switch (result->GetResultInfo()->resultType ) {
    case ACOU_PRESSURE:
    case ACOU_POTENTIAL:
      feFunctions_[formulation_]->ExtractResult( result );
      break;

    default:
      WARN( "Resulttype not computable by acoustic PDE" );
      break;
    }
  }

  //!  Define available postprocessing results
  void AcousticPDE::DefinePrimaryResults(){

    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    if ( formulation_ ==  ACOU_PRESSURE) {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "";
      res1->unit = "Pa";
    } else {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "";
      res1->unit = "m^2/s";
    }
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[formulation_]->SetResultInfo(res1);
    results_.Push_back( res1 );
    availResults_.insert( res1 );

    res1->SetFeFunction(feFunctions_[formulation_]);


    // === ACOUSTIC RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ACOU_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    availResults_.insert( rhs );
    postProcResults_[ACOU_RHS_LOAD] = ACOU_RHS_LOAD;

    //creates the mean flow
    StdVector<std::string> velDofNames;
    std::string geometryType;
    param->Get("domain")->GetValue("geometryType", geometryType );

    if( geometryType == "3d" ) {
      velDofNames = "x", "y", "z";
    } else if( geometryType == "plane" ) {
      velDofNames = "x", "y";
    } else if( geometryType == "axi" ) {
      velDofNames = "r", "z";
    }

    CreateMeanFlowFunction(velDofNames);
  }

   void AcousticPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
     //// === MEAN FLUIDMECH VELOCITY ===
     shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
     flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
     flowvelocity->dofNames = dofNames;
     flowvelocity->unit = "m/s";

     flowvelocity->definedOn = ResultInfo::NODE;
     flowvelocity->entryType = ResultInfo::VECTOR;



     shared_ptr<BaseFeFunction> meanFunction;
     std::string form = SolutionTypeEnum.ToString(MEAN_FLUIDMECH_VELOCITY);
     PtrParamNode feSpaceNode = infoNode_->Get("feSpaces");
     PtrParamNode potSpaceNode = feSpaceNode->Get(form);
     shared_ptr<FeSpace> tmpSpace = FeSpace::CreateInstance(myParam_,potSpaceNode,
                                                            FeSpace::H1, ptgrid_);

     if(isComplex_){
       meanFunction.reset(new FeFunction<Complex>());
       meanFunction->SetFeSpace(tmpSpace);
       meanFunction->SetResultInfo(flowvelocity);
       meanFunction->SetGrid(ptgrid_);
       meanFunction->SetPDE(this);
       flowvelocity->SetFeFunction(meanFunction);
       if(dim_==2)
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,2,2,Complex>,Complex >(meanFunction,flowvelocity));
       else
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,3,3,Complex>,Complex >(meanFunction,flowvelocity));
     }else{
       meanFunction.reset(new FeFunction<Double>());
       meanFunction->SetFeSpace(tmpSpace);
       meanFunction->SetResultInfo(flowvelocity);
       meanFunction->SetGrid(ptgrid_);
       meanFunction->SetPDE(this);
       flowvelocity->SetFeFunction(meanFunction);
       if(dim_==2)
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,2,2,Double>,Double >(meanFunction,flowvelocity));
       else
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,3,3,Double>,Double >(meanFunction,flowvelocity));
     }

     results_.Push_back( flowvelocity );
     availResults_.insert( flowvelocity );

   }

  //! Init the time stepping
  void AcousticPDE::InitTimeStepping(){
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::NEWMARK, 0) );

    feFunctions_[formulation_]->SetTimeScheme(myScheme);

  }
}
