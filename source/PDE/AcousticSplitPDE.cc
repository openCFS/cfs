#include "AcousticSplitPDE.hh" /// ADAPT

#include "General/defs.hh"	/// ADAPT ?? der ganzen header files

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


//new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/SurfaceOperators.hh"

#include "FeBasis/FeFunctions.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <def_expl_templ_inst.hh>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"	/// ADAPT pseudo?

namespace CoupledField{

  DECLARE_LOG(acousticsplitpde)	/// ADAPT
   DEFINE_LOG(acousticsplitpde, "pde.acousticsplit") /// ADAPT


  AcousticSplitPDE::AcousticSplitPDE( Grid* aGrid, PtrParamNode paramNode, /// ADAPT
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "acousticSplit"; /// ADAPT
    pdematerialclass_  = FLUID; /// ADAPT
    nonLin_            = false; ///
    
    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;


    //check for pressure or potential formulation
    std::string pdeFormulation = myParam_->Get("formulation")->As<std::string>();
    if(pdeFormulation == "default"){
      formulation_ = SPLIT_SCALAR;
    }else{
      formulation_ = SolutionTypeEnum.Parse(pdeFormulation);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticSplitPDE::CreateFeSpaces( const std::string&  formulation, PtrParamNode infoNode ){

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      std::string form = SolutionTypeEnum.ToString(formulation_);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[formulation_] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[formulation_]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of split PDE is not known!");
    }

    return crSpaces;
  }
  
  

  void AcousticSplitPDE::DefineIntegrators(){

    RegionIdType actRegion;
    // BaseMaterial * actSDMat = NULL;

    //type of geometry
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    //flag indicating frequency PML formulation
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      // actSDMat = it->second;

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
      
      //=======================================================================
      // Generate coefficient functions
      //=======================================================================
      PtrCoefFct val1 = CoefFunction::Generate( mp_, Global::REAL, "1.0");;

      // ====================================================================
      // Take account for pml (frequency domain only)
      // ====================================================================
      shared_ptr<CoefFunction> coeffPMLScal, coeffPMLVec;
      bool isMapping = false;
      if( dampingList_[actRegion] == MAPPING ) {
        std::string dampId;
        curRegNode->GetValue("dampingId",dampId);
        if(analysistype_ == HARMONIC){
          EXCEPTION("Harmonic analysis not allowed!");
        }else{
          PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("mapping","id",dampId.c_str());
          coeffPMLVec.reset(new CoefFunctionMapping<Double>(pmlNode,val1,actSDList,regions_,true));
          coeffPMLScal.reset(new CoefFunctionMapping<Double>(pmlNode,val1,actSDList,regions_,false));
          // store pml factor
          //matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLVec);
          isMapping = true;
        }
      }


      // ====================================================================
      // standard stiffness integrator
      // ====================================================================
      // factor for damping matrix: factor / c0
      BaseBDBInt * stiffInt = NULL;
      if( dim_ == 2 ) {
        if(isMapping){
          stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,2,Double>(),
                                        coeffPMLScal, 1.0, updatedGeo_ );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), val1 , 1.0, updatedGeo_ );
        }
      }
      else{
        if(isMapping){
          stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,3,Double>(),
              coeffPMLScal, 1.0, updatedGeo_ );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), val1, 1.0, updatedGeo_ );
        }
      }
      
      stiffInt->SetName("LaplaceIntegrator");

      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[formulation_]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
      stiffInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_[actRegion] = stiffInt;


      }
    }

  void AcousticSplitPDE::DefineNcIntegrators() { //as it is
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        if (dim_ == 2)
          DefineMortarCoupling<2,1>(formulation_, *ncIt);
        else
          DefineMortarCoupling<3,1>(formulation_, *ncIt);
        break;
      case NC_NITSCHE:
        if (dim_ == 2)
          DefineNitscheCoupling<2,1>(formulation_, *ncIt );
        else
          DefineNitscheCoupling<3,1>(formulation_, *ncIt );
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }
  
  void AcousticSplitPDE::DefineSurfaceIntegrators( ){

  }

  void AcousticSplitPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(acousticsplitpde) << "Defining rhs load integrators for splitting PDE";
    
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> empty;

    // =====================================
    //  rhsValues for e.g. for splitting
    // =====================================
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, updatedGeo_ );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      coef[i]->SetConservative(true);
      this->rhsFeFunctions_[formulation_]->AddLoadCoefFunction(coef[i], ent[i]);
    }

  }

  void AcousticSplitPDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticSplitPDE::DefinePrimaryResults(){	/// ADAPT ?? Formulation
    StdVector<std::string> vecComponents;
    if( dim_ == 3 ) {
      vecComponents = "x", "y", "z";
    }
    else if( isaxi_ ) {
      vecComponents = "phi";
    }
    else {
      vecComponents = "z";
    }
    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    if ( formulation_ ==  SPLIT_SCALAR) {
      res1->resultType = SPLIT_SCALAR;
      res1->dofNames = "";
      res1->unit = "m^2/s";
      res1->definedOn = ResultInfo::NODE;
      res1->entryType = ResultInfo::SCALAR;
    } else {
      res1->resultType = SPLIT_VECTOR;
      res1->dofNames = vecComponents;
      res1->unit = "m^2/s";
      res1->definedOn = ResultInfo::NODE;
      res1->entryType = ResultInfo::VECTOR;
    }

    feFunctions_[formulation_]->SetResultInfo(res1);
    results_.Push_back( res1 );
    res1->SetFeFunction(feFunctions_[formulation_]);
    DefineFieldResult( feFunctions_[formulation_], res1 );
    

    // === SPLIT RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = SPLIT_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "m^3/s";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    this->rhsFeFunctions_[formulation_]->SetResultInfo(rhs);
    DefineFieldResult( this->rhsFeFunctions_[formulation_], rhs );
    results_.Push_back( rhs );
    availResults_.insert( rhs );

  }
  
  void AcousticSplitPDE::FinalizePostProcResults(){
    //first call base class method
    SinglePDE::FinalizePostProcResults();

  }

  void AcousticSplitPDE::DefinePostProcResults(){

  }

  //! Init the time stepping
  void AcousticSplitPDE::InitTimeStepping(){
    // Use complete implicit scheme
    Double gamma = 1;
    GLMScheme * scheme = new Trapezoidal(gamma);
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0) );
    if(formulation_==SPLIT_SCALAR)
      feFunctions_[SPLIT_SCALAR]->SetTimeScheme(myScheme);
    else
      feFunctions_[SPLIT_VECTOR]->SetTimeScheme(myScheme);

  }

}
