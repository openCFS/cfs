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
#include "Domain/CoefFunction/CoefFunctionMulti.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <cmath>


#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

namespace CoupledField{

   DEFINE_LOG(acousticsplitpde, "pde.acousticsplit")


  AcousticSplitPDE::AcousticSplitPDE( Grid* aGrid, PtrParamNode paramNode,
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "acousticSplit";
    pdematerialclass_  = ACOUSTIC;
    nonLin_            = false;
    updatedGeo_        = false;
    formulation_ = SPLIT_SCALAR;

    //check for split scalar or vector formulation
    std::string formulation = myParam_->Get("formulation")->As<std::string>();
    if(formulation != "default"){
      formulation_ = SolutionTypeEnum.Parse(formulation);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticSplitPDE::CreateFeSpaces( const std::string&  formulation, PtrParamNode infoNode ) {
    if(!(formulation == "default" || formulation == "H1" || formulation == "HCurl")){
      EXCEPTION("The formulation " << formulation << " of split PDE is not known!");
    }
    if(analysistype_ == HARMONIC){
      EXCEPTION("Harmonic analysis not allowed!");
    }

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    std::string form = SolutionTypeEnum.ToString(formulation_);
    PtrParamNode potSpaceNode = infoNode->Get(form);
    crSpaces[formulation_] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,
            ((formulation == "HCurl" && dim_== 3) ? FeSpace::HCURL : FeSpace::H1),
            ptGrid_);

    crSpaces[formulation_]->Init(solStrat_);
    return crSpaces;
  }
  

  void AcousticSplitPDE::ReadDampingInformation() {
    std::map<std::string, DampingType> idDampType;

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // In this case the "damping" layer is an infinite mapping of the solution,
    // since it is a coordinate transformation it can be set up using the PML-like
    // framework.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++) {
      actRegion = regions_[iRegion];

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
      // Generate coefficient functions ("MATERIAL")
      //=======================================================================
      PtrCoefFct val1 = CoefFunction::Generate( mp_, Global::REAL, "1.0");

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

      // ====================================================================
      // Take account for mapping of an infinite domain
      // ====================================================================
      shared_ptr<CoefFunction> coeffMAPScal, coeffMAPVec;
      bool isMapping = false;
      if( dampingList_[actRegion] == MAPPING ) {
        std::string dampId;
        curRegNode->GetValue("dampingId",dampId);
        PtrParamNode mapNode = myParam_->Get("dampingList")->GetByVal("mapping","id",dampId.c_str());
        coeffMAPVec.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,true));
        coeffMAPScal.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,false));
        isMapping = true;
      }

      // ============================
      // Standard Mass Regularization Matrix K \approx K + regFac*M
      // ============================
      // initially, check for regularization factor
      Double regularizationFactor = 1e-6;
      regularizedRegions_.insert(actRegion);

      BaseBDBInt * stiffInt = NULL;
      BaseBDBInt *massInt = NULL;
      switch(formulation_) {
         case SPLIT_SCALAR  :
           // ====================================================================
           // standard stiffness integrator Scalar (includes infinite mapping)
           // ====================================================================
           if( dim_ == 2 ) {
             if(isMapping){
               stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,2,Double>(),
                                             coeffMAPScal, 1.0, updatedGeo_ );
               stiffInt->SetBCoefFunctionOpB(coeffMAPVec);
               massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeH1,3,Double>(),
                   val1,regularizationFactor);
             }else{
               stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), val1 , 1.0, updatedGeo_ );
               massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeH1,3,Double>(),
                   val1,regularizationFactor);
             }
           }
           else{
             if(isMapping){
               stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,3,Double>(),
                   coeffMAPScal, 1.0, updatedGeo_ );
               stiffInt->SetBCoefFunctionOpB(coeffMAPVec);
               massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeH1,3,Double>(),
                   val1,regularizationFactor);
             }else{
               stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), val1, 1.0, updatedGeo_ );
               massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeH1,3,Double>(),
                   val1,regularizationFactor);
             }
           }
           stiffInt->SetName("LaplaceIntegrator");
           break; //optional

         case SPLIT_VECTOR  :
           // ====================================================================
           // standard stiffness integrator Vector TODO include mapping for the problem
           // ====================================================================
           if( dim_ == 2 ) {
             if(isMapping){

               stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,2,Double>(),
                                             coeffMAPScal, 1.0, updatedGeo_ );
               stiffInt->SetBCoefFunctionOpB(coeffMAPVec);
               massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeH1,3,Double>(),
                   val1,regularizationFactor);
             }else{
               stiffInt = new BBInt<Double>(new CurlOperator<FeH1,2,Double>(), val1 , 1.0, updatedGeo_ );
               massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeH1,3,Double>(),
                   val1,regularizationFactor);
             }
           }
           else{
             if(isMapping){
               //TODO implement mapping Vector Potential e.g. implementation vector potential
               stiffInt = new BBInt<Double>(new ScaledCurlOperator<3,Double>(),
                   coeffMAPScal, 1.0, updatedGeo_ );
               stiffInt->SetBCoefFunctionOpB(coeffMAPVec);
             }
             else
             {
               stiffInt = new BBInt<Double>(new CurlOperator<FeHCurl,3,Double>(), val1 , 1.0, updatedGeo_ );
               massInt = new BBInt<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),
                 val1,regularizationFactor);
             }
           }
           stiffInt->SetName("CurlCurlIntegrator");
            break; //optional

         default : //Optional
            EXCEPTION("Unknown PDE fromulation.")
      }
      
      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[formulation_]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
      stiffInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

      massInt->SetName("MassIntegrator");
      BiLinFormContext *massContext =  new BiLinFormContext(massInt, STIFFNESS );

      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
      massInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm( massContext );

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

    // Get FESpace and FeFunction of formulation
    shared_ptr<BaseFeFunction> myFct = feFunctions_[formulation_];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    StdVector<std::string> vecDofNames = myFct->GetResultInfo()->dofNames;

    LinearForm * lin = NULL;
    Double factor = 1.0;

    // ===============================================
    //  rhsValues int(test div(u)) or int(test rot(u))
    // ===============================================
    if ( formulation_ ==  SPLIT_SCALAR || dim_==2) {
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, updatedGeo_ );
    }

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      coef[i]->SetConservative(true);
      this->rhsFeFunctions_[formulation_]->AddLoadCoefFunction(coef[i], ent[i]);
    } //for

    // ==============================
    //  RHS DENSITY div(u) or rot(u)
    // ==============================
    //  RHS DENSITY rot(u)
    LOG_DBG(acousticsplitpde) << "Reading RHS densities" << ((formulation_== SPLIT_VECTOR) ?  " CurlCurl": " Poisson") ;
    ReadRhsExcitation( (dim_==2) ? "rhsDensity" : "rhsDensity3d",
         (dim_==2) ? empty : vecDofNames,
             (formulation_== SPLIT_VECTOR) ? ResultInfo::VECTOR : ResultInfo::SCALAR,
                 isComplex_, ent, coef, updatedGeo_ );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if (ent[i]->GetType() == EntityList::NODE_LIST)
        EXCEPTION("Rhs density must be defined on elements")

      switch(formulation_) {
         case SPLIT_VECTOR  :
            if( dim_ == 2 ) {
              if(isComplex_) {
                lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                                 Complex(factor), coef[i], updatedGeo_);
              } else  {
                lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                    factor, coef[i], updatedGeo_);
              }
            }else{
                lin = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                    factor, coef[i], updatedGeo_);
            }
            break; //optional

         default : //Optional
           if(isComplex_) {
             lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                              Complex(factor), coef[i], updatedGeo_);
           } else  {
             lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                 factor, coef[i], updatedGeo_);
           }
      } // end switch case

      lin->SetName("RhsDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // end for

    // ===========================
    //  NORMAL or TANGENTIAL VELOCITY (surface)
    // We read the total velocity vector (in the case of tangential velocity)
    // and its components are projected during the integration
    // by tangential coordinates spanned by the vectors of the jacobian of the shape functions
    // In this sense we have a gradient basis on the surface element.
    // For the normal velocity component we only use the normal component as input.
    // ===========================
    LOG_DBG(acousticsplitpde) << "Reading neumann bc";
    std::set<RegionIdType> regions (regions_.Begin(), regions_.End() );

    switch(formulation_) {
       case SPLIT_VECTOR  :
           WARN("Boundary values must be defined on nodes!")
           if(myParam_->Get("bcsAndLoads")->Has("normalVelocity"))
             EXCEPTION("Normal Velocity is not used in the scalar potential formulation! - uncomment tag in XML input ");

           if ( dim_ == 2){
             ReadRhsExcitation( "tangentialVelocity", vecDofNames, ResultInfo::SCALAR, isComplex_,
                                   ent, coef, updatedGeo_ );
             if(myParam_->Get("bcsAndLoads")->Has("tangentialVelocity3D"))
               EXCEPTION("Tangential Velocity 3D is not used in the scalar potential formulation! - uncomment tag in XML input ");
           }else{
             ReadRhsExcitation( "tangentialVelocity3D", vecDofNames, ResultInfo::VECTOR, isComplex_,
                                   ent, coef, updatedGeo_ );
             if(myParam_->Get("bcsAndLoads")->Has("tangentialVelocity"))
               EXCEPTION("Tangential Velocity 2D is not used in the scalar potential formulation! - uncomment tag in XML input ");
           }
          break; //optional

       default : //Optional
         ReadRhsExcitation( "normalVelocity", empty, ResultInfo::SCALAR, isComplex_,
                               ent, coef, updatedGeo_ );
         if(myParam_->Get("bcsAndLoads")->Has("tangentialVelocity") || myParam_->Get("bcsAndLoads")->Has("tangentialVelocity3D"))
           EXCEPTION("Tangential Velocity is not used in the scalar potential formulation! - uncomment tag in XML input ");
    } // end switch case

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if(ent[i]->GetType() != EntityList::SURF_ELEM_LIST)
        EXCEPTION("Neumann boundary can only be defined on surface elements");

      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
              Complex(factor), coef[i], regions, updatedGeo_);
        } else {
          lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,2,1>(),
              factor, coef[i], regions, updatedGeo_);
        }
      } else  {
        switch(formulation_) {
           case SPLIT_VECTOR  :
             if(isComplex_) {
               lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                   Complex(factor), coef[i],regions, updatedGeo_);
             }else{
               lin = new BUIntegrator<Double,true>( new IdentityOperator<FeHCurl,3,1,Double>(),
                   factor, coef[i],regions, updatedGeo_);
             }
             break; //optional

           default : //Optional
             if(isComplex_) {
               lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
                   Complex(factor), coef[i], regions, updatedGeo_);
             } else {
               lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,3,1>(),
                   factor, coef[i] , regions, updatedGeo_);
             }
        } // end switch case
      }

      lin->SetName("inhomNeumannBC");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
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

    // === Primary result according to PDE definition ===
    shared_ptr<ResultInfo> potent( new ResultInfo);
    potent->resultType = formulation_;
    potent->unit = "m^2/s";
    if ( formulation_ ==  SPLIT_SCALAR) {
      potent->dofNames = "";
      potent->definedOn = ResultInfo::NODE;
      potent->entryType = ResultInfo::SCALAR;
    }
    else
    {
      potent->dofNames = vecComponents;
      potent->definedOn = ( ptGrid_->GetDim() == 3 ) ? ResultInfo::ELEMENT : ResultInfo::NODE;
      potent->entryType = ResultInfo::VECTOR;
    }

    feFunctions_[formulation_]->SetResultInfo(potent);
    results_.Push_back( potent );
    potent->SetFeFunction(feFunctions_[formulation_]);
    DefineFieldResult( feFunctions_[formulation_], potent );
    
    // === SPLIT RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = SPLIT_RHS_LOAD;
    rhs->unit = "m^3/s";
    rhs->definedOn = results_[0]->definedOn;
    if ( formulation_ ==  SPLIT_SCALAR) {
      rhs->dofNames = "";
      rhs->entryType = ResultInfo::SCALAR;
    }
    else
    {
      rhs->dofNames = vecComponents;
      rhs->entryType = ResultInfo::VECTOR;
    }

    this->rhsFeFunctions_[formulation_]->SetResultInfo(rhs);
    DefineFieldResult( this->rhsFeFunctions_[formulation_], rhs );
    results_.Push_back( rhs );
    availResults_.insert( rhs );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[formulation_] = "homDir";
    idbcSolNameMap_[formulation_] = "inhomDir";

  }
  
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

    // 0. //=== FLOW FIELD ===
    //velocity
    shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
    flowvelocity->resultType = FLUIDMECH_VELOCITY;
    flowvelocity->dofNames = vecDofNames;
    flowvelocity->unit = "m/s";

    flowvelocity->definedOn = ResultInfo::NODE;
    flowvelocity->entryType = ResultInfo::VECTOR;

    velocityCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    DefineFieldResult( velocityCoef_, flowvelocity );

    // I. //Definitions for velocity results
    shared_ptr<ResultInfo> velocityResult;
    PtrCoefFct velocityResultFct;
    shared_ptr<CoefFunctionFormBased>  velocityResultFctPot;
    velocityResult.reset(new ResultInfo);
    velocityResult->dofNames = vecDofNames;
    velocityResult->unit = "m/s";
    velocityResult->entryType = ResultInfo::VECTOR;
    velocityResult->definedOn = ResultInfo::ELEMENT;
    if( isComplex_ ) {
      velocityResultFctPot.reset(new CoefFunctionBOp<Complex>(feFct, velocityResult, 1.0));
    } else {
      velocityResultFctPot.reset(new CoefFunctionBOp<Double>(feFct, velocityResult, 1.0));
    }
    velocityResultFct = velocityResultFctPot;
    stiffFormCoefs_.insert(velocityResultFctPot);

    // II. //Definitions for the complimentary velocity results
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    // === Corrected Velocity scalar split equation===
    shared_ptr<ResultInfo> velCorr( new ResultInfo);
    velCorr->dofNames = vecDofNames;
    velCorr->unit = "m/s";
    velCorr->definedOn = ResultInfo::NODE;
    velCorr->entryType = ResultInfo::VECTOR;

    // === complementary flow_VELOCITY perturbation===
    velocityResult->resultType = ( formulation_ ==  SPLIT_SCALAR) ? SPLIT_SCALAR_VELOCITY : SPLIT_VECTOR_VELOCITY;
    DefineFieldResult( velocityResultFct, velocityResult );
    // this computes the complement of the Helmholtz decomposition
    velCorr->resultType = ( formulation_ ==  SPLIT_SCALAR) ? SPLIT_VECTOR_VELOCITY : SPLIT_SCALAR_VELOCITY;
    CoefXprBinOp minusFct = CoefXprBinOp(mp_, velocityCoef_, velocityResultFct,
                                           CoefXpr::OP_SUB);
    PtrCoefFct velCorrVec = CoefFunction::Generate(mp_, part, minusFct);

    DefineFieldResult( velCorrVec, velCorr );

    //  III. //=== POTENTIAL ENERGY ===
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
  }

  //! Init the time stepping
  void AcousticSplitPDE::InitTimeStepping()
  {
    if (GetDomain()->GetAdaptiveData())
        EXCEPTION("Adaptive timestepping is not supported for AcousticSplitPDE: variable-step BDF2 is incompatible with acoustic wave simulation. Use fixed deltaT.");

    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    // Create scheme from XML or use default Trapezoidal
    GLMScheme* scheme = nullptr;
    if (integrationScheme)
    {
      scheme = GetXmlDefinedScheme(integrationScheme);
    }
    else
    {
      scheme = new Trapezoidal(1.0);
    }

    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0));
    feFunctions_[formulation_]->SetTimeScheme(myScheme);
  }
}
