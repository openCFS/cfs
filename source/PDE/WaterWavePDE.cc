#include "WaterWavePDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/DivOperator.hh"

#include "FeBasis/FeFunctions.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionDiagTensorFromScalar.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <cmath>


#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

  DEFINE_LOG(waterWavepde, "pde.waterWave")


  WaterWavePDE::WaterWavePDE( Grid* aGrid, PtrParamNode paramNode,
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "waterWave";
    pdematerialclass_  = ACOUSTIC;
    nonLin_            = false;
    isMechCoupled_     = false;
    
    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;

    isTimeDomPML_      = false;

    isAPML_ = false;

    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);

    //type of geometry
    isaxi_ = ptGrid_->IsAxi();

    //check, if subtype is surface gravity waves
    isSurfaceGravityWave_ = false;
    std::string subType = myParam_->Get("subType")->As<std::string>();
    if( subType == "surfaceGravityWave")
      isSurfaceGravityWave_ = true;

    // compute surface wave velocity
    g_ = CoefFunction::Generate( mp_, Global::REAL, "9.81");
    PtrCoefFct omega = CoefFunction::Generate( mp_, Global::REAL, "2*pi*f");
    //      PtrCoefFct gravityC0 =
    Cdeep_ = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp( mp_, g_, omega, CoefXpr::OP_DIV));
    PtrCoefFct omega2 = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_,omega,omega,CoefXpr::OP_MULT));
    kdeep_ = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp( mp_, omega2 , g_, CoefXpr::OP_DIV)); // omega^2/g
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  WaterWavePDE::CreateFeSpaces( const std::string&  formulation,
                  PtrParamNode infoNode ){

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      std::string form = SolutionTypeEnum.ToString(WATER_PRESSURE);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[WATER_PRESSURE] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[WATER_PRESSURE]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of water wave PDE is not known!");
    }

    // ===================================
    // Check for transient PML
    // ===================================
    if(this->analysistype_ == TRANSIENT && isTimeDomPML_){
      //now define the additional uknowns
      EXCEPTION("Trasient PML not implemented yet!");
    }
    return crSpaces;
  }

  void WaterWavePDE::ReadDampingInformation()
  {
    // get regions of current PDE
    ParamNodeList regionParamNodes = myParam_->Get("regionList")->GetChildren();
    // corresponding region id
    RegionIdType actRegionId;
    // corresponding region name and damping id
    std::string actRegionName, actDampingId;
    // try to get the dampingList and return if it is not specified
    PtrParamNode dampListNode = myParam_->Get("dampingList", ParamNode::PASS);
    if (dampListNode) {
      // get the individual damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();
      // map of damping ids from the xml and corresponding damping types
      std::map<std::string, DampingType> idDampType;

      // Run over all region param nodes and assign the required damping information
      for (UInt iRegion = 0; iRegion < regionParamNodes.GetSize(); ++iRegion) {
        regionParamNodes[iRegion]->GetValue("name", actRegionName);
        regionParamNodes[iRegion]->GetValue("dampingId", actDampingId);

        // pass if no damping is defined
        if (actDampingId == "")
          continue;

        // parse region id from region name
        actRegionId = ptGrid_->GetRegion().Parse(actRegionName);

        // now, read the damping information from the dampingList
        for (UInt iChild = 0; iChild < dampNodes.GetSize(); ++iChild) {
          std::string dampString = dampNodes[iChild]->GetName();
          std::string actId = dampNodes[iChild]->Get("id")->As<std::string>();
          // only consider damping information for the current damping id
          if (actId != actDampingId)
            continue;

          // determine type of damping
          DampingType actType;
          String2Enum(dampString, actType);

          // store damping type string
          idDampType[actId] = actType;
          // break after the information is set as only one damping ID per region is possible
          break;
        }

        // check actDampingId was indeed registerd above
        if (idDampType.count(actDampingId) == 0)
          EXCEPTION("Damping with id '" << actDampingId << "' of region '" << actRegionName << "' was not found. Is it defined in the 'dampingList'?");

        // assign damping type to the region
        dampingList_[actRegionId] = idDampType[actDampingId];

        // if Rayleigh damping is specified, parse and store the additional damping information
        if (dampingList_[actRegionId] == RAYLEIGH) {
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }
        else if (dampingList_[actRegionId] == ADAPTED_LOSS_TANGENS_DELTA) {
          if (!(analysistype_ == BasePDE::HARMONIC))
            EXCEPTION("Adapted loss tangent delta damping is only allowed for harmonic analysis.");
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetFreqAdaptedRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }
        else if (dampingList_[actRegionId] == GLOBAL_RAYLEIGH) {
          EXCEPTION("Global Rayleigh damping is not yet implemented.");
          if (dampNodes.GetSize() != 1)
            EXCEPTION("Global Rayleigh damping does not allow for additional damping nodes defined.");
          RaylDampingData actRaylCoeffs;
          actRaylCoeffs.alpha = dampNodes[0]->Get("alpha")->As<std::string>();
          actRaylCoeffs.beta = dampNodes[0]->Get("beta")->As<std::string>();
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
        }

        // set flag to compute extra integrators for transient PML
        if (dampingList_[actRegionId] == PML && analysistype_ == BasePDE::TRANSIENT) {
          isTimeDomPML_ = true;
        }
      }
    }

    // read the transform list and store the transform types
    std::map<std::string, DampingType> transformType;
    if (myParam_->Has("transformList")) {
        ParamNodeList transformNodes = myParam_->Get("transformList")->GetChildren();
        std::string strType, strId;
        for (UInt k = 0; k < transformNodes.GetSize(); k++) {
            strType = transformNodes[k]->GetName();
            strId = transformNodes[k]->Get("id")->As<std::string>();
            // determine type of damping
            DampingType actType;
            String2Enum( strType, actType );
            if( transformType.count( strId ) > 0 ) {
                EXCEPTION( "Transform id '" << strId << "' not unique");
            }
            else {
                transformType[strId] = actType;
            }
        }
    }

    // loop over all regions and determine transform
    RegionIdType actRegion;
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
        // Set current region and material
        actRegion = it->first;
        // actSDMat = it->second;

        // Get current region name
        std::string regionName = ptGrid_->GetRegion().ToString(actRegion);      // functions for coordinate transformations (PML or infinite mapping)

        // create new entity list
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( actRegion );

        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());

        shared_ptr<CoefFunction> coeffTransScal, coeffTransVec;

        // transform list
        bool isTransform = false;
        std::string actDampingId, actTransformId;
        curRegNode->GetValue( "dampingId", actDampingId ); // dampingId defaults to empty if not set
        if ( !(actDampingId == "") && curRegNode->Has("transforms")) { // dampingId and transform list
            EXCEPTION("dampingId and transformList cannot be used at the same time - check region '" <<regionName<<"'" );
        }
        else if ( !(actDampingId == "") && !(curRegNode->Has("transforms"))
                  && !(dampingList_[actRegion] == RAYLEIGH)) {// only (non-rayleigh) damping id set
            LOG_TRACE(waterWavepde) << regionName<< ": id =" <<actDampingId<<"\n";
            //EXCEPTION("please use the new <transfromList>");
            std::cout << "please use the new <transfromList> to define a PML or mapping layer\n";
            PtrParamNode transNode = myParam_->Get("dampingList")->GetByVal("pml","id",actDampingId);
            coeffTransVec.reset(new CoefFunctionPML<Complex>(transNode,Cdeep_,actSDList,regions_,true));
            coeffTransScal.reset(new CoefFunctionPML<Complex>(transNode,Cdeep_,actSDList,regions_,false));
            isTransform = true;
        }
        else if (curRegNode->Has("transforms") && (actDampingId == "") ) { // only transform list
            LOG_TRACE(waterWavepde) << "list\n";
            // read damping ids and multiply transform
            // check for PML and analysis type
            ParamNodeList transformNodes = curRegNode->Get("transforms")->GetChildren();
            for (UInt k = 0; k < transformNodes.GetSize(); k++) {
                transformNodes[k]->GetValue( "name", actTransformId );
                // get PML/mapping definition for Id
                PtrParamNode transNode;
                PtrCoefFct vec,scal;
                if ( transformType.count(actTransformId) > 0 ) {
                    if (transformType.at(actTransformId)==PML) {
                        transNode = myParam_->Get("transformList")->GetByVal("pml","id",actTransformId);
                        vec.reset(new CoefFunctionPML<Complex>(transNode,Cdeep_,actSDList,regions_,true));
                        scal.reset(new CoefFunctionPML<Complex>(transNode,Cdeep_,actSDList,regions_,false));
                    }
                    else if (transformType.at(actTransformId)==MAPPING) {
                        transNode = myParam_->Get("transformList")->GetByVal("mapping","id",actTransformId);
                        PtrCoefFct sos = CoefFunction::Generate( mp_, Global::REAL, "1.0");
                        if(analysistype_ == HARMONIC){
                            PtrCoefFct one = CoefFunction::Generate( mp_, Global::REAL, "1.0");
                            sos = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp( mp_, one, kdeep_, CoefXpr::OP_DIV));
                            vec.reset(new CoefFunctionMapping<Complex>(transNode,sos,actSDList,regions_,true));
                            scal.reset(new CoefFunctionMapping<Complex>(transNode,sos,actSDList,regions_,false));
                        }
                        else {
                            vec.reset(new CoefFunctionMapping<Double>(transNode,sos,actSDList,regions_,true));
                            scal.reset(new CoefFunctionMapping<Double>(transNode,sos,actSDList,regions_,false));
                        }
                    }
                    else {
                        EXCEPTION("this should not happen")
                    }
                }
                else {
                    EXCEPTION("transform id '" <<actTransformId<<"' not found in transformList");
                }
                if (k==0) { // initialize
                    coeffTransVec = vec;//.reset(new CoefFunctionPML<Complex>(transNode,Cdeep_,actSDList,regions_,true));
                    coeffTransScal = scal; //.reset(new CoefFunctionPML<Complex>(transNode,Cdeep_,actSDList,regions_,false));
                }
                else { // multiply
                    shared_ptr<CoefFunction> prodVec,prodScal;
                    prodVec = CoefFunction::Generate( mp_, Global::COMPLEX,
                            CoefXprBinOp(mp_, coeffTransVec, vec, CoefXpr::OP_MULT_COMP) );
                    prodScal = CoefFunction::Generate( mp_, Global::COMPLEX,
                            CoefXprBinOp(mp_, coeffTransScal, scal, CoefXpr::OP_MULT) );
                    coeffTransVec = prodVec;
                    coeffTransScal = prodScal;
                }
                isTransform = true;
            }
        }
        else { // no mapping or damping
            LOG_TRACE(waterWavepde) << "normal\n";
        }
        // save to transformList
        shared_ptr< std::pair<PtrCoefFct,PtrCoefFct> > Fcts ; // Initialize as NULL pointer
        if (isTransform) { // set if we have transform
            // first scalar, second vector
            Fcts.reset( new std::pair<PtrCoefFct,PtrCoefFct>(coeffTransScal,coeffTransVec ) );
        }
        transformFctList_[actRegion] = Fcts;
    }
  }

  void WaterWavePDE::DefineIntegrators(){

    RegionIdType actRegion;
    // BaseMaterial * actSDMat = NULL;

    //type of geometry
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[WATER_PRESSURE]->GetFeSpace();
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      // Set current region and material
      actRegion = regions_[iRegion];
      // actSDMat = it->second;
      // save the density
      PtrCoefFct dens = materials_[actRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
      matCoefs_[ELEM_DENSITY]->AddRegion(actRegion, materials_[actRegion]->GetScalCoefFnc( DENSITY, Global::REAL ) );

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
      PtrCoefFct factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");

      // ====================================================================
      // standard stiffness integrator
      // ====================================================================
      BaseBDBInt * stiffInt = NULL;
      // if we have a transform defined on the region we need to use scaled operators
      if(transformFctList_[actRegion]){
          LOG_TRACE(waterWavepde) << "transform on region '"<<regionName<<"'\n";
          PtrCoefFct coeffTransScal = transformFctList_[actRegion]->first;
          PtrCoefFct coeffTransVec = transformFctList_[actRegion]->second;
          if(coeffTransVec->IsComplex()) {
              LOG_TRACE(waterWavepde) << "  -> should be PML\n";
              if ( dim_ == 2) {
                  stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,2,Complex>(), coeffTransScal, 1.0, updatedGeo_ );
                  stiffInt->SetBCoefFunctionOpB(coeffTransVec);
              } else {
                  stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,3,Complex>(), coeffTransScal, 1.0, updatedGeo_ );
                  stiffInt->SetBCoefFunctionOpB(coeffTransVec);
              }
          }
          else {
              LOG_TRACE(waterWavepde) << "  -> should be only mapping\n";
              if (dim_ == 2) {
                  stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,2,Double>(),coeffTransScal, 1.0, updatedGeo_ );
                  stiffInt->SetBCoefFunctionOpB(coeffTransVec);
              } else {
                  stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,3,Double>(),coeffTransScal, 1.0, updatedGeo_ );
                  stiffInt->SetBCoefFunctionOpB(coeffTransVec);
              }
          }
      }
      // use the standard operators
      else {
          if (dim_==2) {
              stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), factor, 1.0, updatedGeo_ );
          } else {
              stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), factor, 1.0, updatedGeo_ );
          }
      }

      stiffInt->SetName("LaplaceIntegrator");

      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );

      // check for damping
      if (dampingList_[actRegion] == RAYLEIGH || dampingList_[actRegion] == ADAPTED_LOSS_TANGENS_DELTA || dampingList_[actRegion] == GLOBAL_RAYLEIGH) {
        RaylDampingData &actDamp = (regionRaylDamping_[actRegion]);
        stiffIntDescr->SetSecDestMat(DAMPING, actDamp.beta);
      }

      feFunctions_[WATER_PRESSURE]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[WATER_PRESSURE],feFunctions_[WATER_PRESSURE]);
      stiffInt->SetFeSpace( feFunctions_[WATER_PRESSURE]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

    }
  }

  template<UInt DIM>
  void WaterWavePDE::DefineTransientPMLInts(shared_ptr<ElemList> eList, std::string id){

    EXCEPTION("REFACTOR");
  }

  
  void WaterWavePDE::DefineSurfaceIntegrators( ){
    //========================================================================================
    // boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
        //free surface condition for gravity waves
        ParamNodeList freeSurfaceNodes = bcNode->GetList( "freeSurfaceCondition" );
        for( UInt i = 0; i < freeSurfaceNodes.GetSize(); i++ ) {
            std::string regionName = freeSurfaceNodes[i]->Get("name")->As<std::string>();
            shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
            std::string volRegName = freeSurfaceNodes[i]->Get("volumeRegion")->As<std::string>();
            LOG_TRACE(waterWavepde) << "free surf of "<< volRegName << "\n";

            // define necessary factors
            PtrCoefFct factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");
            // factor for mass matrix: 1 / gravity
            PtrCoefFct coeffMass = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, factor, g_, CoefXpr::OP_DIV ) );

            // setup integrator
            BiLinearForm * gravityInt = NULL;
            // check if volume region has mapping or PML
            RegionIdType volRegion = ptGrid_->GetRegion().Parse(volRegName);
            if (transformFctList_[volRegion]) { // we have mapping / PML in the region
                // read mapping list from volume
                PtrCoefFct coeffTransScal = transformFctList_[volRegion]->first;
                PtrCoefFct prod;
                if (coeffTransScal->IsComplex()) { // PML
                    LOG_TRACE(waterWavepde) << "  -> has PML\n";
                    prod = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, coeffMass, coeffTransScal, CoefXpr::OP_MULT ) );
                } else { // MAPPING
                    LOG_TRACE(waterWavepde) << "  -> has only mapping\n";
                    prod = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, coeffMass, coeffTransScal, CoefXpr::OP_MULT ) );
                }
                coeffMass = prod;
            }
            if (coeffMass->IsComplex()) {
                if( dim_ == 2 ) {
                    gravityInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1>(), coeffMass, 1.0, updatedGeo_ );
                } else {
                    gravityInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1>(), coeffMass, 1.0, updatedGeo_ );
                }
            } else {
                if( dim_ == 2 ) {
                    gravityInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), coeffMass, 1.0, updatedGeo_ );
                } else {
                    gravityInt = new BBInt<>(new IdentityOperator<FeH1,3,1>(), coeffMass, 1.0, updatedGeo_ );
                }
            }
            gravityInt->SetName("gravityWaveIntegrator");
            BiLinFormContext *gravityContext = new BiLinFormContext(gravityInt, MASS);

            // check for damping
            if (dampingList_[volRegion] == RAYLEIGH || dampingList_[volRegion] == ADAPTED_LOSS_TANGENS_DELTA || dampingList_[volRegion] == GLOBAL_RAYLEIGH) {
              RaylDampingData &actDamp = (regionRaylDamping_[volRegion]);
              gravityContext->SetSecDestMat(DAMPING, actDamp.alpha);
            }

            gravityContext->SetEntities( actSDList, actSDList );
            gravityContext->SetFeFunctions( feFunctions_[WATER_PRESSURE] , feFunctions_[WATER_PRESSURE]);
            feFunctions_[WATER_PRESSURE]->AddEntityList( actSDList );
            assemble_->AddBiLinearForm( gravityContext );
        } //free surface condition for gravity waves
      }
    }
  }


  void WaterWavePDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(waterWavepde) << "Defining rhs load integrators for WaterWave PDE";
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[WATER_PRESSURE];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> empty;
    LinearForm * lin = NULL;

    // obtain density
    shared_ptr<CoefFunctionMulti> densFct = matCoefs_[ELEM_DENSITY];
    shared_ptr<CoefFunctionSurf> surfDens(new CoefFunctionSurf(false));
    surfDens->SetVolumeCoefs( densFct->GetRegionCoefs() );

    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

    bool coefUpdateGeo;

    // ===========================
    //  general surface load
    // ===========================
    StdVector<std::string> volumeRegions;
    ReadRhsExcitation( "surfaceLoad", empty, ResultInfo::SCALAR, isComplex_, ent, coef,coefUpdateGeo,volumeRegions);
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      LOG_TRACE(waterWavepde) << "  surface region number "<< i+1;
      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++)  {
          // check dimension
          UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
          if( elemDim != (dim_-1) ) {
              EXCEPTION("surfaceLoad can only be defined on surface elements");
          }
          //TODO: one should get the bounding volume element for each surface element and get the transform from there!
      }
      PtrCoefFct exValue;
      // check for volume region, if defined get transorm
      std::string volRegName = volumeRegions[i];
      PtrCoefFct mapFact;
      if (!(volRegName=="")) {
          LOG_TRACE(waterWavepde) << "  -> volume region: '"<< volRegName <<"'";
          RegionIdType actRegion = ptGrid_->GetRegionId(volRegName);
          // if we have a transform defined on the region we need to use scaled operators
          if(transformFctList_[actRegion]){
              PtrCoefFct mapFact = transformFctList_[actRegion]->first;
              LOG_TRACE(waterWavepde) << "  -> volumeRegion has transform";
              if (mapFact->IsComplex() && coef[i]->IsComplex()) { // PML
                  LOG_TRACE(waterWavepde) << "  -> complex";
                  exValue = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, coef[i], mapFact, CoefXpr::OP_MULT ) );
              } else { // MAPPING
                  LOG_TRACE(waterWavepde) << "  -> real";
                  exValue = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, coef[i], mapFact, CoefXpr::OP_MULT ) );
              }
          }
          else {
              EXCEPTION("volumeRegion has no transform defined");
          }
      }
      else {
          exValue = coef[i]; //mapFact = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      }
      // define integrators
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
                                                1.0, exValue, volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,2,1>(),
                                               1.0, exValue, volRegions, coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
                                                1.0, exValue, volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,3,1>(),
                                               1.0, exValue , volRegions, coefUpdateGeo);
        }
      }
      lin->SetName("SurfaceLoadIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

    // ===========================
    //  DISPLACEMENT (surface)
    // ===========================
    LOG_DBG(waterWavepde) << "Reading total velocity";
    StdVector<std::string> vecDofNames;
    if(dim_ == 3)
      vecDofNames = "x", "y", "z";
    if(dim_ == 2 && !isaxi_)
      vecDofNames = "x", "y";
    if(dim_ == 2 && isaxi_)
      EXCEPTION("Axi is not implemented!");
      // vecDofNames = "r", "z";
    ReadRhsExcitation( "displacement", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef,coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("displacement can only be defined on surface elements");
      }
      // in this case the pressure can be related to the normal velocity as
      // p_n = - Omega^2*v_n*rho
      PtrCoefFct omega2Rho = CoefFunction::Generate( mp_, Global::COMPLEX,
          CoefXprBinOp(mp_,CoefFunction::Generate( mp_, Global::COMPLEX,"4*pi*pi*f*f","0.0"),surfDens, CoefXpr::OP_MULT)) ;
      PtrCoefFct exValue = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_,omega2Rho, coef[i], CoefXpr::OP_MULT) );

      if( dim_ == 2) {
        lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(), 1.0, exValue, volRegions, coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(), 1.0, exValue, volRegions, coefUpdateGeo);
      }

      lin->SetName("VelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }


    // =====================================
    //  rhsValues
    // =====================================
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      coef[i]->SetConservative(true);
      this->rhsFeFunctions_[WATER_PRESSURE]->AddLoadCoefFunction(coef[i], ent[i]);
    }

    // ===========================
    //  TOTAL ACCELERATION (surface)
    // ===========================
    LOG_DBG(waterWavepde) << "Reading acceleration";
    ReadRhsExcitation( "acceleration", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef,coefUpdateGeo );
    // complex part for expressions
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("accelertation can only be defined on surface elements");
      }
      PtrCoefFct rhoAcc = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, surfDens, coef[i], CoefXpr::OP_MULT) );
      if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(), -1.0, rhoAcc, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,2>(), -1.0, rhoAcc, volRegions, coefUpdateGeo);
          }
      } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(), -1.0, rhoAcc, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,3>(), -1.0, rhoAcc, volRegions, coefUpdateGeo);
          }
      }
      lin->SetName("AccelerationIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

  }

  void WaterWavePDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void WaterWavePDE::DefinePrimaryResults(){

    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1( new ResultInfo);

    res1->resultType = WATER_PRESSURE;
    res1->dofNames = "";
    res1->unit = "Pa";

    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[WATER_PRESSURE]->SetResultInfo(res1);
    results_.Push_back( res1 );
    res1->SetFeFunction(feFunctions_[WATER_PRESSURE]);
    DefineFieldResult( feFunctions_[WATER_PRESSURE], res1 );
    
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[WATER_PRESSURE] = "soundSoft";
    idbcSolNameMap_[WATER_PRESSURE] = "pressure";
    
    // === ACOUSTIC RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = WATER_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = ResultInfo::NODE;
    rhs->entryType = ResultInfo::SCALAR;
    this->rhsFeFunctions_[WATER_PRESSURE]->SetResultInfo(rhs);
    DefineFieldResult( this->rhsFeFunctions_[WATER_PRESSURE], rhs );
    results_.Push_back( rhs );
    availResults_.insert( rhs );

    //creates vector dofs
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

    // === PML DAMPING FACTORS ===
    //if( matCoefs_.find(PML_DAMP_FACTOR) != matCoefs_.end() ) {
    shared_ptr<ResultInfo> pml ( new ResultInfo );
    pml->resultType = PML_DAMP_FACTOR;
    pml->dofNames = vecDofNames;
    //pml->dofNames = "";
    pml->unit = "";
    pml->definedOn = ResultInfo::ELEMENT;
    pml->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> pmlFct(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, 
                                                               isComplex_));
    //matCoefs_[PML_DAMP_FACTOR] = pmlFct;
    DefineFieldResult(pmlFct, pml);
    //}

    // === PML AUX Variables ===
    if(this->isTimeDomPML_){
      if(!this->isAPML_ && dim_ == 3){
        shared_ptr<ResultInfo> pmlScal ( new ResultInfo );
        pmlScal->resultType = WATER_PMLAUXSCALAR;
        pmlScal->dofNames = "";
        pmlScal->unit = "-";
        pmlScal->definedOn = ResultInfo::NODE;
        pmlScal->entryType = ResultInfo::SCALAR;
        feFunctions_[WATER_PMLAUXSCALAR]->SetResultInfo(pmlScal);
        results_.Push_back( pmlScal );
        pmlScal->SetFeFunction(feFunctions_[WATER_PMLAUXSCALAR]);
        DefineFieldResult( feFunctions_[WATER_PMLAUXSCALAR], pmlScal );
      }

      shared_ptr<ResultInfo> pmlVec ( new ResultInfo );
      pmlVec->resultType = WATER_PMLAUXVEC;
      pmlVec->dofNames = vecDofNames;
      pmlVec->unit = "-";
      pmlVec->definedOn = ResultInfo::NODE;
      pmlVec->entryType = ResultInfo::VECTOR;
      feFunctions_[WATER_PMLAUXVEC]->SetResultInfo(pmlVec);
      results_.Push_back( pmlVec );
      pmlVec->SetFeFunction(feFunctions_[WATER_PMLAUXVEC]);
      DefineFieldResult( feFunctions_[WATER_PMLAUXVEC], pmlVec );
    }

  }
  
  void WaterWavePDE::FinalizePostProcResults(){
    //first call base class method
    SinglePDE::FinalizePostProcResults();

    // complex part for expressions
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // TorqueDensityTensor
    // compute r x sigma (cross product of location vector and stress tensor)
    // The cross product can be written as Matrix R times vector (inner product), thus generalizes to R in sigma
    // since sigma is diagonal with p (pressure) as entries we compute R p (Matrix R times scalar p)
    PtrCoefFct presCoef = GetCoefFct(WATER_PRESSURE);
    assert(presCoef);
    shared_ptr<CoefFunctionMulti> tdCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[WATER_TDT]);
    // setup the matrix
    StdVector<std::string> crossProductMatrixElements;
    // set the nonzero entries
    if ( ptGrid_->GetDim() == 3 ) { // for general 3D tensor (unsymmatric)
      crossProductMatrixElements = "0","-z","y" , "z","0","-x" , "-y","x","0";
    } else { // in 2d we have only a result in z-direction (Mz = - P*nx*ry + P*ny*rx )
      crossProductMatrixElements = "-y","x";
    }
    // generate the matrix R and multiply by p
    PtrCoefFct crossProductMat;
    if(!isComplex_){
      crossProductMat = CoefFunction::Generate(mp_,Global::REAL,crossProductMatrixElements);
    } else {
      StdVector<std::string> zeroElements; zeroElements.Resize(crossProductMatrixElements.GetSize(),"0");
      crossProductMat = CoefFunction::Generate(mp_,Global::COMPLEX,crossProductMatrixElements,zeroElements);
    }
    PtrCoefFct tdtXpr = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, presCoef, crossProductMat, CoefXpr::OP_MULT));
    // now set the expression for all regions
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ) {
      RegionIdType actRegion = *regIt;
      tdCoef->AddRegion( actRegion, tdtXpr );
    }

  }

  void WaterWavePDE::DefinePostProcResults(){

    shared_ptr<BaseFeFunction> feFct = feFunctions_[WATER_PRESSURE];

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

    // === DENSITY ===
    shared_ptr<ResultInfo> density ( new ResultInfo );
    density->resultType = ELEM_DENSITY;
    density->dofNames = "";
    density->unit = "kg/m^3";
    density->definedOn = ResultInfo::ELEMENT;
    density->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> densFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false )); // we do not have complex density
    matCoefs_[ELEM_DENSITY] = densFct;
    DefineFieldResult(densFct, density);

    // === PARTICLE_POSITION ===
    shared_ptr<ResultInfo> pos(new ResultInfo);
    pos->resultType = WATER_POSITION;
    pos->dofNames = vecDofNames;
    pos->unit = "m";
    pos->entryType = ResultInfo::VECTOR;
    pos->definedOn = ResultInfo::ELEMENT;
    // pressure gradient
    shared_ptr<CoefFunctionFormBased> presGradFct;
    shared_ptr<BaseFeFunction> presFct = feFunctions_[WATER_PRESSURE];
    if( isComplex_ ) {
      presGradFct.reset(new CoefFunctionBOp<Complex>(presFct, pos, 1.0));
    } else {
      presGradFct.reset(new CoefFunctionBOp<Double>(presFct, pos, 1.0));
    }
    stiffFormCoefs_.insert(presGradFct);
    // u = 1/(rho*omega^2) * grad(p)
    PtrCoefFct oneOverOmega2rho = CoefFunction::Generate( mp_, Global::REAL,
              CoefXprBinOp( mp_, CoefFunction::Generate( mp_, Global::REAL,"1.0"),
                CoefXprBinOp(mp_,CoefFunction::Generate( mp_, Global::REAL, "4*pi*pi*f*f"), densFct, CoefXpr::OP_MULT ),
              CoefXpr::OP_DIV ));
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    PtrCoefFct posFct = CoefFunction::Generate( mp_, part, CoefXprBinOp( mp_, oneOverOmega2rho, presGradFct, CoefXpr::OP_MULT ) );
    DefineFieldResult( posFct, pos );

    // === WATER PRESSURE TENSOR ===
    // p I
    shared_ptr<ResultInfo> presTens(new ResultInfo);
    presTens->resultType = WATER_PRES_TENS;
    StdVector<std::string> tensorComponentNames;
    if(dim_==3){
      tensorComponentNames = "xx", "yy", "zz", "yz", "xz", "xy";
    } else { // 2d
      tensorComponentNames = "xx", "yy", "xy";
    }
    presTens->dofNames = tensorComponentNames;
    presTens->unit = MapSolTypeToUnit(WATER_PRES_TENS);
    presTens->entryType = ResultInfo::TENSOR;
    presTens->definedOn = ResultInfo::ELEMENT;
    presTens->SetFeFunction(feFunctions_[WATER_PRESSURE]);
    availResults_.insert( presTens );
    StdVector<PtrCoefFct> presTensDiagValues;
    presTensDiagValues = StdVector<PtrCoefFct>(dim_);
    for(UInt i = 0; i < dim_; i++){
      presTensDiagValues[i] = this->GetCoefFct( WATER_PRESSURE );
    }
    std::string presTensSubType;
    if (dim_==2) {
      presTensSubType = "plane";
    } else { // should be 3D
      presTensSubType = "3d";
    }
    shared_ptr<CoefFunction> presTensCoef (new CoefFunctionDiagTensorFromScalar(presTensDiagValues,presTensSubType));
    DefineFieldResult( presTensCoef, presTens );

    // === WATER_SURFACE_TRACTION ===
    shared_ptr<ResultInfo> surfaceTractionInfo;
    shared_ptr<CoefFunctionSurf> surfaceTractionFct;
    surfaceTractionInfo.reset(new ResultInfo);
    surfaceTractionInfo->resultType = WATER_SURFACE_TRACTION;
    surfaceTractionInfo->dofNames = vecDofNames;
    surfaceTractionInfo->unit = "Pa";
    surfaceTractionInfo->entryType = ResultInfo::VECTOR;
    surfaceTractionInfo->definedOn = ResultInfo::SURF_ELEM;
    surfaceTractionFct.reset(new CoefFunctionSurf(true, 1.0, surfaceTractionInfo));
    DefineFieldResult(surfaceTractionFct, surfaceTractionInfo);
    surfCoefFcts_[surfaceTractionFct] = feFunctions_[WATER_PRESSURE];

    // === FLUID-MECHANIC REACTION FORCE (= integral of surface traction over the surface region ) ===
    shared_ptr<ResultInfo> reactionForceInfo;
    reactionForceInfo.reset(new ResultInfo);
    reactionForceInfo->resultType = WATER_SURFACE_FORCE;
    reactionForceInfo->dofNames = vecDofNames;
    reactionForceInfo->unit = MapSolTypeToUnit(WATER_SURFACE_FORCE);
    reactionForceInfo->entryType = ResultInfo::VECTOR;
    reactionForceInfo->definedOn = ResultInfo::SURF_REGION;
    // Integrate surface traction
    shared_ptr<ResultFunctor> reactionForceFct;
    if (isComplex_)
      reactionForceFct.reset(new ResultFunctorIntegrate<Complex>(surfaceTractionFct, feFct, reactionForceInfo));
    else
      reactionForceFct.reset(new ResultFunctorIntegrate<Double>(surfaceTractionFct, feFct, reactionForceInfo));
    resultFunctors_[WATER_SURFACE_FORCE] = reactionForceFct;
    availResults_.insert(reactionForceInfo);
    
    // === TorqueDensityTensor TDT (= cross product of location vector with stress tensor ) ===
    // implemented as a vector, such that the normal mapping with CoefFunctionSurf works
    shared_ptr<ResultInfo> torqueDensityTensorInfo;
    torqueDensityTensorInfo.reset(new ResultInfo);
    torqueDensityTensorInfo->resultType = WATER_TDT;
    StdVector<std::string> tensorFullComponentNames;
    if(dim_==3){
      tensorFullComponentNames = "xx","xy","yz","yx","yy","yz","zx","zy","zz";
    } else { // 2d
      tensorFullComponentNames = "x","y";
    }
    torqueDensityTensorInfo->dofNames = tensorFullComponentNames;
    torqueDensityTensorInfo->unit = MapSolTypeToUnit(WATER_TDT);
    torqueDensityTensorInfo->entryType = ResultInfo::TENSOR;
    torqueDensityTensorInfo->definedOn = ResultInfo::ELEMENT;
    // actual computation (expression) is defined in finalize, 
    // Tensor valued results are (can be?) evaluated from vector-type coefFunctions
    shared_ptr<CoefFunctionMulti> torqueDensityTensorFct(new CoefFunctionMulti(CoefFunction::VECTOR, tensorFullComponentNames.GetSize(), 1, isComplex_));
    DefineFieldResult(torqueDensityTensorFct, torqueDensityTensorInfo);

    // === FLUID-MECHANIC REACTION TORQUE DENSITY (= cross product of location vector with surface traction ) ===
    // since the surface traction cannot be used in CoefExpressions we map a "torque density tensor" in normal direction
    shared_ptr<ResultInfo> surfaceTorqueDensityInfo;
    shared_ptr<CoefFunctionSurf> surfaceTorqueDensityFct;
    surfaceTorqueDensityInfo.reset(new ResultInfo);
    surfaceTorqueDensityInfo->resultType = WATER_SURFACE_TORQUE_DENSITY;
    surfaceTorqueDensityInfo->dofNames = vecDofNames;
    surfaceTorqueDensityInfo->unit = MapSolTypeToUnit(WATER_SURFACE_TORQUE_DENSITY);
    surfaceTorqueDensityInfo->entryType = ResultInfo::VECTOR;
    surfaceTorqueDensityInfo->definedOn = ResultInfo::SURF_ELEM;
    surfaceTorqueDensityFct.reset(new CoefFunctionSurf(true, 1.0, surfaceTorqueDensityInfo));
    DefineFieldResult(surfaceTorqueDensityFct, surfaceTorqueDensityInfo);
    surfCoefFcts_[surfaceTorqueDensityFct] = torqueDensityTensorFct;

    // === FLUID-MECHANIC REACTION TORQUE (= integral of surface torque density over the surface region ) ===
    shared_ptr<ResultInfo> reactionTorqueInfo;
    reactionTorqueInfo.reset(new ResultInfo);
    reactionTorqueInfo->resultType = WATER_SURFACE_TORQUE;
    StdVector<std::string> torqueDofNames;
    if(dim_==3)
      torqueDofNames = vecDofNames;
    else // 2D, we only have a z-direction, but since Vector-results are expected to have size=dim_ we add a dummy component
      torqueDofNames = "z","dummy";
    reactionTorqueInfo->dofNames = torqueDofNames;
    reactionTorqueInfo->unit = MapSolTypeToUnit(WATER_SURFACE_TORQUE);
    reactionTorqueInfo->entryType = ResultInfo::VECTOR;
    reactionTorqueInfo->definedOn = ResultInfo::SURF_REGION;
    // Integrate surface torque
    shared_ptr<ResultFunctor> reactionTorqueFct;
    if (isComplex_)
        reactionTorqueFct.reset(new ResultFunctorIntegrate<Complex>(surfaceTorqueDensityFct, feFct, reactionTorqueInfo));
    else
        reactionTorqueFct.reset(new ResultFunctorIntegrate<Double>(surfaceTorqueDensityFct, feFct, reactionTorqueInfo));
    resultFunctors_[WATER_SURFACE_TORQUE] = reactionTorqueFct;
    availResults_.insert(reactionTorqueInfo);
  }
  
  //! Init the time stepping
  void WaterWavePDE::InitTimeStepping()
  {
    if (GetDomain()->GetAdaptiveData())
        EXCEPTION("Adaptive timestepping is not supported for WaterWavePDE: variable-step BDF2 is incompatible with wave simulation. Use fixed deltaT.");

    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    PtrParamNode timeStepAlphaNode = this->myParam_->Get("timeStepAlpha", ParamNode::PASS);
    if (integrationScheme && timeStepAlphaNode)
      throw Exception("Both 'integrationScheme' and 'timeStepAlpha' are specified for the water wave PDE. "
                      "Please use 'integrationScheme' only, as it provides more flexibility and "
                      "supersedes the legacy 'timeStepAlpha' parameter.");

    auto makeScheme = [&]() -> GLMScheme* {
      if (integrationScheme)
        return GetXmlDefinedScheme(integrationScheme);
      else {
        Double alpha = this->myParam_->Get("timeStepAlpha")->As<Double>();
        return new Newmark(0.5, 0.25, alpha);
      }
    };

    if (this->isTimeDomPML_)
    {
      // Basically the choice for alpha scheme needs to be done everytime we have
      // a damping matrix, not just for PML
      shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(makeScheme(), 0));
      shared_ptr<BaseTimeScheme> vecScheme(new TimeSchemeGLM(makeScheme(), 0));
      feFunctions_[WATER_PRESSURE]->SetTimeScheme(acouScheme);
      feFunctions_[WATER_PMLAUXVEC]->SetTimeScheme(vecScheme);

      if (!this->isAPML_ && dim_ == 3)
      {
        shared_ptr<BaseTimeScheme> scalScheme(new TimeSchemeGLM(makeScheme(), 0));
        feFunctions_[WATER_PMLAUXSCALAR]->SetTimeScheme(scalScheme);
      }
    }
    else
    {
      shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(makeScheme(), 0));
      feFunctions_[WATER_PRESSURE]->SetTimeScheme(acouScheme);
    }
  }


  template void WaterWavePDE::DefineTransientPMLInts<2>(shared_ptr<ElemList>, std::string);
  template void WaterWavePDE::DefineTransientPMLInts<3>(shared_ptr<ElemList>, std::string);
