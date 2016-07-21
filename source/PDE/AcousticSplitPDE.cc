#include "AcousticSplitPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


//new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
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
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <def_expl_templ_inst.hh>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

namespace CoupledField{

  DECLARE_LOG(acousticsplitpde)
   DEFINE_LOG(acousticsplitpde, "pde.acousticsplit")


  AcousticSplitPDE::AcousticSplitPDE( Grid* aGrid, PtrParamNode paramNode,
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "acousticSplit";
    pdematerialclass_  = FLUID;
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
  AcousticSplitPDE::CreateFeSpaces( const std::string&  formulation, PtrParamNode infoNode ) {
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
  

  void AcousticSplitPDE::ReadDampingInformation() {
    std::map<std::string, DampingType> idDampType;

    // try to get dampingList
    PtrParamNode dampListNode = myParam_->Get( "dampingList", ParamNode::PASS );
    if( dampListNode ) {

      // get specific damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();

      for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->As<std::string>();

        // determine type of damping
        DampingType actType;
        String2Enum( dampString, actType );

        // store damping type string
        idDampType[actId] = actType;

      }
    }

    // Run over all region
    ParamNodeList regionNodes =
        myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptGrid_->GetRegion().Parse( actRegionName );

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }
      dampingList_[actRegionId] = idDampType[actDampingId];
    }
  }

  void AcousticSplitPDE::DefineIntegrators(){
    RegionIdType actRegion;

    //type of geometry
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;

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
      PtrCoefFct val1 = CoefFunction::Generate( mp_, Global::REAL, "1.0");

      // ====================================================================
      // check for vorticity field, which is just piped through the simulation
      // ====================================================================
      std::string vorticityId = curRegNode->Get("vorticityId")->As<std::string>();
      PtrCoefFct regionVortic;

      if ( vorticityId != "" ) {
        std::set<UInt> definedDofs;
        //just real valued vorticity allowed
        shared_ptr<ResultInfo> vorticInfo = GetResultInfo(FLUIDMECH_VORTICITY);
        //Add the region information
        PtrParamNode vorticNode = myParam_->Get("flowList")->GetByVal("vorticity","name",vorticityId.c_str());

        ReadUserFieldValues( actSDList, vorticNode, vorticInfo->dofNames, vorticInfo->entryType,
            isComplex_, regionVortic, definedDofs, updatedGeo_ );
        vorticityCoef_->AddRegion( actRegion, regionVortic );
      }

      // ====================================================================
      // check for density field, which is just piped through the simulation
      // ====================================================================
      std::string densityId = curRegNode->Get("densityId")->As<std::string>();
      PtrCoefFct regionDensity;

      if ( densityId != "" ) {
        std::set<UInt> definedDofs;
        //just real valued density allowed
        shared_ptr<ResultInfo> densityInfo = GetResultInfo(FLUIDMECH_DENSITY);
        //Add the region information
        PtrParamNode densityNode = myParam_->Get("flowList")->GetByVal("density","name",densityId.c_str());

        ReadUserFieldValues( actSDList, densityNode, densityInfo->dofNames, densityInfo->entryType,
            isComplex_, regionDensity, definedDofs, updatedGeo_ );
        densityCoef_->AddRegion( actRegion, regionDensity );
      }

      // ====================================================================
      // check for velocity field, which is just piped through the simulation
      // ====================================================================
      std::string velocityId = curRegNode->Get("flowId")->As<std::string>();
      PtrCoefFct regionVelocity; // ??

      if ( velocityId != "" ) {
        std::set<UInt> definedDofs;
        //just real valued velocity allowed
        shared_ptr<ResultInfo> velocityInfo = GetResultInfo(FLUIDMECH_VELOCITY);
        //Add the region information
        PtrParamNode velocityNode = myParam_->Get("flowList")->GetByVal("flow","name",velocityId.c_str());

        ReadUserFieldValues( actSDList, velocityNode, velocityInfo->dofNames, velocityInfo->entryType,
            isComplex_, regionVelocity, definedDofs, updatedGeo_ );
        velocityCoef_->AddRegion( actRegion, regionVelocity );
      }

      BaseBDBInt * stiffInt = NULL;
      if(formulation_ == SPLIT_SCALAR){
        // ====================================================================
        // Take account for mapping
        // ====================================================================
        shared_ptr<CoefFunction> coeffMAPScal, coeffMAPVec;
        bool isMapping = false;
        if( dampingList_[actRegion] == MAPPING ) {
          std::string dampId;
          curRegNode->GetValue("dampingId",dampId);
          if(analysistype_ == HARMONIC){
            EXCEPTION("Harmonic analysis not allowed!");
          }else{
            PtrParamNode mapNode = myParam_->Get("dampingList")->GetByVal("mapping","id",dampId.c_str());
            coeffMAPVec.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,true));
            coeffMAPScal.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,false));
            isMapping = true;
          }
        }


      // ====================================================================
      // standard stiffness integrator Scalar
      // ====================================================================

        if( dim_ == 2 ) {
          if(isMapping){
            stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,2,Double>(),
                                          coeffMAPScal, 1.0, updatedGeo_ );
            stiffInt->SetBCoefFunctionOpB(coeffMAPVec);
          }else{
            stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), val1 , 1.0, updatedGeo_ );
          }
        }
        else{
          if(isMapping){
            stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,3,Double>(),
                coeffMAPScal, 1.0, updatedGeo_ );
            stiffInt->SetBCoefFunctionOpB(coeffMAPVec);
          }else{
            stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), val1, 1.0, updatedGeo_ );
          }
        }
      }else{
        // ====================================================================
        // standard stiffness integrator Vector
        // ====================================================================

          if( dim_ == 2 ) {
              stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), val1 , 1.0, updatedGeo_ );
          }
          else{
            // TODO edge elements
            EXCEPTION("Not implemented yet!");
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

      if(formulation_ == SPLIT_SCALAR){
        // initially, check for regularization factor
        Double regularizationFactor = 1e-6;
        // ============================
        // Standard Mass Regularization Matrix
        // ============================
        bool regulize = true;
        if ( regulize ) {
          // do not use gradients for non-conductive regions (for regularization
          // only the lowest order mass term is used)
          //useGrad = false;

          Matrix<Double> reluc;
          // add region to set of "regularized" regions
          regularizedRegions_.insert(actRegion);
        }

        BaseBDBInt *massInt;

        BiLinFormContext * massContext;
        if ( regulize) {
          // we have to guarantee, that we add some mass to curl-curl integrator.
          // Additionally, the integrator gets scaled by the edge size for a uniform
          // conditioning
          massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeH1,3,Double>(),
              val1,regularizationFactor);
          massInt->SetName("MassIntegrator");
          massContext =  new BiLinFormContext(massInt, STIFFNESS );

          massContext->SetEntities( actSDList, actSDList );
          massContext->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
          massInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

          assemble_->AddBiLinearForm( massContext );
        }
      }

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

    // TODO formulation_ onlz if this ...

  }

  void AcousticSplitPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(acousticsplitpde) << "Defining rhs load integrators for splitting PDE";
    
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> empty;

    // Get FESpace and FeFunction of formulation
    shared_ptr<BaseFeFunction> myFct = feFunctions_[formulation_];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    LinearForm * lin = NULL;
    Double factor = 1.0;

    // =====================================
    //  rhsValues for e.g. for splitting
    // =====================================
    if ( formulation_ ==  SPLIT_SCALAR) {
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, updatedGeo_ );
    }else{
      ReadRhsExcitation( "rhsValues", empty, ResultInfo::VECTOR, isComplex_,
                                ent, coef, updatedGeo_ );
    }
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      coef[i]->SetConservative(true);
      this->rhsFeFunctions_[formulation_]->AddLoadCoefFunction(coef[i], ent[i]);
    } //for

    // ================
    //  RHS DENSITY
    // ================
    if ( formulation_ ==  SPLIT_SCALAR) {
      LOG_DBG(acousticsplitpde) << "Reading rhs densities";
      ReadRhsExcitation( "rhsDensity", empty,
                         ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("Rhs density must be defined on elements")
        }

        if(isComplex_) {
          lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                           Complex(factor), coef[i], updatedGeo_);
        } else  {
          lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
              factor, coef[i], updatedGeo_);
        }

        lin->SetName("RhsDensityInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
      } // for
    }else{
      LOG_DBG(acousticsplitpde) << "Reading rhs densities";
      ReadRhsExcitation( "rhsDensity", empty,
                         ResultInfo::VECTOR, isComplex_, ent, coef, updatedGeo_ );
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("Rhs density must be defined on elements")
        }
        if( dim_ == 2 ) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                             Complex(factor), coef[i], updatedGeo_);
          } else  {
            lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                factor, coef[i], updatedGeo_);
          }
        }else{
          // TODO edge elements
        }
        lin->SetName("RhsDensityInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
      } // for
    }





  }

  void AcousticSplitPDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticSplitPDE::DefinePrimaryResults(){
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

    StdVector<std::string> vecDofNames;
    if( ptGrid_->GetDim() == 3 ) {
      vecDofNames = "x", "y", "z";
    } else {
      if( ptGrid_->IsAxi() ) {
        vecDofNames = "r", "z";
      } else {
        vecDofNames = "x", "y";
      }
    }

    // === Primary result according to PDE definition ===
    shared_ptr<ResultInfo> potent( new ResultInfo);
    if ( formulation_ ==  SPLIT_SCALAR) {
      potent->resultType = SPLIT_SCALAR;
      potent->dofNames = "";
      potent->unit = "m^2/s";
      potent->definedOn = ResultInfo::NODE;
      potent->entryType = ResultInfo::SCALAR;
    } else {
      // TODO if 2D
      potent->resultType = SPLIT_VECTOR;
      potent->dofNames = "";
      potent->unit = "m^2/s";
      potent->definedOn = ResultInfo::NODE;
      potent->entryType = ResultInfo::SCALAR;
    }

    feFunctions_[formulation_]->SetResultInfo(potent);
    results_.Push_back( potent );
    potent->SetFeFunction(feFunctions_[formulation_]);
    DefineFieldResult( feFunctions_[formulation_], potent );
    
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    if ( formulation_ ==  SPLIT_SCALAR) {
      hdbcSolNameMap_[SPLIT_SCALAR] = "homDir";
      idbcSolNameMap_[SPLIT_SCALAR] = "inhomDir";
    } else {
      hdbcSolNameMap_[SPLIT_VECTOR] = "homDir";
      idbcSolNameMap_[SPLIT_VECTOR] = "inhomDir";

    }

    // === SPLIT RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    if ( formulation_ ==  SPLIT_SCALAR) {
    rhs->resultType = SPLIT_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "m^3/s";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    } else {
      rhs->resultType = SPLIT_RHS_LOAD;
      rhs->dofNames = vecComponents;
      rhs->unit = "m^3/s";
      rhs->definedOn = results_[0]->definedOn;
      rhs->entryType = ResultInfo::VECTOR;
    }
    this->rhsFeFunctions_[formulation_]->SetResultInfo(rhs);
    DefineFieldResult( this->rhsFeFunctions_[formulation_], rhs );
    results_.Push_back( rhs );
    availResults_.insert( rhs );



  }
  
//  void AcousticSplitPDE::FinalizePostProcResults(){
//    //first call base class method
//    SinglePDE::FinalizePostProcResults();
//
//  }

  void AcousticSplitPDE::DefinePostProcResults(){
    shared_ptr<BaseFeFunction> feFct = feFunctions_[formulation_];

    StdVector<std::string> vecDofNames;
    if( ptGrid_->GetDim() == 3 ) {
      vecDofNames = "x", "y", "z";
    } else {
      if( ptGrid_->IsAxi() ) {
        vecDofNames = "r", "z";
      } else {
        vecDofNames = "x", "y";
      }
    }

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

    shared_ptr<ResultInfo> vel;
    PtrCoefFct velFct;
    shared_ptr<CoefFunctionFormBased>  velFctPot;

    if ( formulation_ ==  SPLIT_SCALAR) {
      // === Acoustic_VELOCITY pertubation===
      vel.reset(new ResultInfo);
      vel->resultType = SPLIT_SCALAR_VELOCITY;
      vel->dofNames = vecDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::NODE;
      // Velocity v = grad Potential
      if( isComplex_ ) {
        velFctPot.reset(new CoefFunctionBOp<Complex>(feFct, vel, 1.0));
      } else {
        velFctPot.reset(new CoefFunctionBOp<Double>(feFct, vel, 1.0));
      }
      velFct = velFctPot;
      stiffFormCoefs_.insert(velFctPot);
      DefineFieldResult( velFct, vel );

      // === POTENTIAL ENERGY ===
      shared_ptr<ResultFunctor> keFuncPot;
      shared_ptr<ResultInfo> potEnergy(new ResultInfo);
      potEnergy->resultType = SPLIT_POT_ENERGY;
      potEnergy->dofNames = "";
      potEnergy->unit = "Ws";
      potEnergy->entryType = ResultInfo::SCALAR;
      potEnergy->definedOn = ResultInfo::REGION;
      availResults_.insert ( potEnergy );
      if( isComplex_ ) {
        keFuncPot.reset(new EnergyResultFunctor<Complex>(feFct, potEnergy, 0.5));
      } else {
        keFuncPot.reset(new EnergyResultFunctor<Double>(feFct, potEnergy, 0.5));
      }
      resultFunctors_[SPLIT_POT_ENERGY] = keFuncPot;
      stiffFormFunctors_.insert(keFuncPot);


      // TODO adapt for 3D
      // === FLOW FIELD ===
      // vorticity
      shared_ptr<ResultInfo> vorticity( new ResultInfo);
      vorticity->resultType = FLUIDMECH_VORTICITY;
      vorticity->dofNames = vecComponents;
      vorticity->unit = "1/s";

      vorticity->definedOn = ResultInfo::NODE;
      vorticity->entryType = ResultInfo::VECTOR;

      vorticityCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, 1,1,isComplex_));
      DefineFieldResult( vorticityCoef_, vorticity );

      //density
      shared_ptr<ResultInfo> density( new ResultInfo);
      density->resultType = FLUIDMECH_DENSITY;
      density->dofNames = "";
      density->unit = "kg/m^3";

      density->definedOn = ResultInfo::NODE;
      density->entryType = ResultInfo::SCALAR;

      densityCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,isComplex_));
      DefineFieldResult( densityCoef_, density );

      //velocity
      shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
      flowvelocity->resultType = FLUIDMECH_VELOCITY;
      flowvelocity->dofNames = vecDofNames;
      flowvelocity->unit = "m/s";

      flowvelocity->definedOn = ResultInfo::NODE;
      flowvelocity->entryType = ResultInfo::VECTOR;

      velocityCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      DefineFieldResult( velocityCoef_, flowvelocity );
  //    results_.Push_back( flowvelocity );
  //    availResults_.insert( flowvelocity );

      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      // === ACOUSTIC Source FIELD LAMB===
      // LAMB
      shared_ptr<ResultInfo> lamb( new ResultInfo);
      lamb->resultType = SPLIT_LAMB;
      lamb->dofNames = vecDofNames;
      lamb->unit = "kg/(ms)^2";

      lamb->definedOn = ResultInfo::NODE;
      lamb->entryType = ResultInfo::VECTOR;

      CoefXprBinOp minusFct = CoefXprBinOp(mp_, velocityCoef_, velFct,
                                             CoefXpr::OP_SUB);
      PtrCoefFct uMinusGradScalar = CoefFunction::Generate(mp_, part, minusFct);

      CoefXprBinOp crossFct = CoefXprBinOp(mp_, vorticityCoef_, uMinusGradScalar,
                                                 CoefXpr::OP_CROSS);
      PtrCoefFct vorticityCrossU = CoefFunction::Generate(mp_, part, crossFct);

      CoefXprVecScalOp multiFct = CoefXprVecScalOp(mp_, vorticityCrossU, densityCoef_,
                                                 CoefXpr::OP_MULT);
      PtrCoefFct lambVec = CoefFunction::Generate(mp_, part, crossFct);
  //    DefineFieldResult( lambVec, lamb );

  // COPY constructor problem
      // #1 möchte gerne Lamb vector als ergebnis
      // #2 muss ich nach set derivative noch was amchen?
      //divLamb
      shared_ptr<ResultInfo> divLamb( new ResultInfo);
      divLamb->resultType = SPLIT_DIVLAMB;
      divLamb->dofNames = "";
      divLamb->unit = "kg/(m^3s^2)";

      divLamb->definedOn = ResultInfo::NODE;
      divLamb->entryType = ResultInfo::SCALAR;

  //    PtrCoefFct divLambVec;
  //    divLambVec = lambVec;  // TODO ... copy constructor coeffunction
      UInt gDim = this->domain_->GetGrid()->GetDim();
      UInt dDim = lamb->dofNames.GetSize();
      lambVec->SetDerivativeOperation(CoefFunction::VECTOR_DIVERGENCE,gDim,dDim);
      DefineFieldResult( lambVec, divLamb );

    }else{
      // === incompressible flow_VELOCITY pertubation===
      vel.reset(new ResultInfo);
      vel->resultType = SPLIT_VECTOR_VELOCITY;
      vel->dofNames = vecDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::NODE;
      // Velocity v = curl Potential
      // TODO formulation_ onlz if this ... (FE function space)?
      if( isComplex_ ) {
        velFctPot.reset(new CoefFunctionBOp<Complex>(feFct, vel, 1.0));
      } else {
        velFctPot.reset(new CoefFunctionBOp<Double>(feFct, vel, 1.0));
      }
      velFct = velFctPot;
      stiffFormCoefs_.insert(velFctPot);
      DefineFieldResult( velFct, vel );


      // TODO formulation_ onlz if this ... (Results and post processing)

    }

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


//  void AcousticSplitPDE::FinilizeBeforTimeStep() {
//    RegionIdType actRegion;
//
//    // Define integrators for "standard" materials
//    std::map<RegionIdType, BaseMaterial*>::iterator it;
//    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();
//
//    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
//      // Set current region and material
//      actRegion = it->first;
//
//      // Get current region name
//      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
//
//      // create new entity list
//      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
//      actSDList->SetRegion( actRegion );
//
//      // --- Set the FE ansatz for the current region ---
//      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
//
//
//
//    }
//
// }
}
