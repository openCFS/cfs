
#include "AcousticPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"


#include "FeBasis/FeFunctions.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

#include "Domain/Results/ResultFunctor.hh"

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
    std::string pdeFormulation;
    paramNode->GetValue("pdeFormulation",pdeFormulation,ParamNode::EX);
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
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[formulation_]->Init();
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
      // stiffness integrator
      // ====================================================================

      shared_ptr<CoefFunction> coeffK
                = CoefFunction::Generate(Global::REAL, "1.0");

//      BBInt< GradientOperator,FeH1,Double > * stiffInt;
      BiLinearForm * stiffInt = NULL;
      if( dim_ == 2 ) {
         // if( isComplex_ ){
         //   shared_ptr<CoefFunction> coeffK
         //             = CoefFunction::Generate(Global::COMPLEX, "1.0","0.0");
         //   stiffInt = new BBInt<GradientOperator<FeH1,2,Complex>, Complex >(coeffK,1.0 );
         // }else{
            shared_ptr<CoefFunction> coeffK
                      = CoefFunction::Generate(Global::REAL, "1.0");
            stiffInt = new BBInt<GradientOperator<FeH1,2>, Double >(coeffK,1.0 );
          //}
      }else{
        //if( isComplex_ ){
        //  shared_ptr<CoefFunction> coeffK
        //            = CoefFunction::Generate(Global::COMPLEX, "1.0","0.0");
        //  stiffInt = new BBInt<GradientOperator<FeH1,3,Complex>, Complex >(coeffK,1.0 );
        //}else{
          shared_ptr<CoefFunction> coeffK
                    = CoefFunction::Generate(Global::REAL, "1.0");
          stiffInt = new BBInt<GradientOperator<FeH1,3>, Double >(coeffK,1.0 );
        //}
      }
      stiffInt->SetName("LaplaceIntegrator");
      //stiffInt = new BBInt<GradientOperator,FeH1,Double >(coeffK,1.0 );

      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[formulation_]->AddEntityList( actSDList );

      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
      stiffInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );

      // ====================================================================
      // mass integrator
      // ====================================================================

      //shared_ptr<CoefFunction> coeffM
      //    = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0/(c0*c0)));

      BiLinearForm *massInt = NULL;

      //if( isComplex_ ){
      //  shared_ptr<CoefFunction> coeffM
      //      = CoefFunction::Generate(Global::COMPLEX, lexical_cast<std::string>(1.0/(c0*c0)),"0.0");
      //  if(dim_==2)
      //    massInt = new BBInt<IdentityOperator<FeH1,2,1,Complex>, Complex >(coeffK,1.0 );
      //  else
      //    massInt = new BBInt<IdentityOperator<FeH1,3,1,Complex>, Complex >(coeffK,1.0 );
      //}else{
        shared_ptr<CoefFunction> coeffM
                  = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0/(c0*c0)));
        if(dim_==2)
          massInt = new BBInt<IdentityOperator<FeH1,2,1,Double>, Double  >(coeffK,1.0 );
        else
          massInt = new BBInt<IdentityOperator<FeH1,3,1,Double>, Double  >(coeffK,1.0 );
      //}


      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace() );

      BiLinFormContext *massContext =  new BiLinFormContext(massInt, MASS );

      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
      assemble_->AddBiLinearForm( massContext );
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

 }

  //! Init the time stepping
  void AcousticPDE::InitTimeStepping(){
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::NEWMARK, 0) );

    feFunctions_[formulation_]->SetTimeScheme(myScheme);

  }
}
