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
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <def_expl_templ_inst.hh>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

  DECLARE_LOG(waterWavepde)
  DEFINE_LOG(waterWavepde, "pde.waterWave")


  WaterWavePDE::WaterWavePDE( Grid* aGrid, PtrParamNode paramNode,
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "waterWave";
    pdematerialclass_  = FLUID;
    nonLin_            = false;
    isMechCoupled_     = false;
    
    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;

    isTimeDomPML_      = false;

    isAPML_ = false;

    //check, if subtype is surface gravity waves
    isSurfaceGravityWave_ = false;
    std::string subType = myParam_->Get("subType")->As<std::string>();
    if( subType == "surfaceGravityWave")
      isSurfaceGravityWave_ = true;

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

//      if(dim_==3 && !this->isAPML_){
//        PtrParamNode scalarpml = infoNode->Get("TransientPMLScalarAuxVar");
//        crSpaces[WATER_PMLAUXSCALAR] =
//            FeSpace::CreateInstance(myParam_,scalarpml,FeSpace::H1, ptGrid_);
//        crSpaces[WATER_PMLAUXSCALAR]->Init(solStrat_);
//      }
//
//      PtrParamNode vectorPML = infoNode->Get("TransientPMLVectorAuxVar");
//      crSpaces[WATER_PMLAUXVEC] =
//          FeSpace::CreateInstance(myParam_,vectorPML,FeSpace::H1, ptGrid_);
//      crSpaces[WATER_PMLAUXVEC]->Init(solStrat_);

    }
    return crSpaces;
  }
  
  
  void WaterWavePDE::ReadDampingInformation() {
    std::map<std::string, DampingType> idDampType;
    std::map<std::string, shared_ptr<RaylDampingData> > idRaylData;

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

    // Run over all region and set entry in "regionNonLinId"
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
      if(dampingList_[actRegionId] == PML &&
          analysistype_ == BasePDE::TRANSIENT ) {
        isTimeDomPML_ = true;
      }
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

    //flag indicating frequency PML formulation
    bool harmonicPML = false;
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

      // compute surface wave velocity
      PtrCoefFct gravity = CoefFunction::Generate( mp_, Global::REAL, "9.81");
      PtrCoefFct omega = CoefFunction::Generate( mp_, Global::REAL, "2*pi*f");
      PtrCoefFct gravityC0 =
          CoefFunction::Generate( mp_,  Global::REAL,
                    CoefXprBinOp( mp_, gravity, omega, CoefXpr::OP_DIV));

      PtrCoefFct factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");

      // ====================================================================
      // Take account for pml (frequency domain only)
      // ====================================================================
      shared_ptr<CoefFunction> coeffPMLScal, coeffPMLVec;
      //just coeffunctions for the transformation of jacobians
      shared_ptr<CoefFunction> coeffPMLfactor;

      if( dampingList_[actRegion] == PML ) {
        std::string dampId;
        curRegNode->GetValue("dampingId",dampId);
        if(analysistype_ == HARMONIC){
          PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",dampId.c_str());
          coeffPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode,gravityC0,actSDList,regions_,true));
          coeffPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode,gravityC0,actSDList,regions_,false));

          // store pml factor
          //matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLVec);
          coeffPMLfactor = CoefFunction::Generate( mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, factor,coeffPMLScal,CoefXpr::OP_MULT));
          harmonicPML = true;
        }else{
          if(dim_==2)
            DefineTransientPMLInts<2>(actSDList,dampId);
          else
            DefineTransientPMLInts<3>(actSDList,dampId);
        }
      }else{
        harmonicPML = false;
      } // damping == PML

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
          coeffMAPVec.reset(new CoefFunctionMapping<Double>(mapNode,factor,actSDList,regions_,true));
          coeffMAPScal.reset(new CoefFunctionMapping<Double>(mapNode,factor,actSDList,regions_,false));
          isMapping = true;
        }
      }

      // ====================================================================
      // standard stiffness integrator
      // ====================================================================
      BaseBDBInt * stiffInt = NULL;
      if( dim_ == 2 ) {
        if(isMapping){
          stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,2,Double>(),
              coeffPMLScal, 1.0, updatedGeo_ );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }
        else if(harmonicPML){
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,2,Complex>(),
                                        coeffPMLfactor, 1.0, updatedGeo_ );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), factor, 
                                       1.0, updatedGeo_ );
        }
      }
      else{
        if(isMapping){
                    stiffInt = new BBInt<Double>(new ScaledGradientOperator<FeH1,3,Double>(),
                        coeffMAPScal, 1.0, updatedGeo_ );
                    stiffInt->SetBCoefFunctionOpB(coeffMAPVec);
                  }
        else if(harmonicPML){
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,3,Complex>(),
                                        coeffPMLfactor, 1.0, updatedGeo_ );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), factor,
                                       1.0, updatedGeo_ );
        }
      }
      
      stiffInt->SetName("LaplaceIntegrator");

      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[WATER_PRESSURE]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[WATER_PRESSURE],feFunctions_[WATER_PRESSURE]);
      stiffInt->SetFeSpace( feFunctions_[WATER_PRESSURE]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_[actRegion] = stiffInt;

    }
  }

  template<UInt DIM>
  void WaterWavePDE::DefineTransientPMLInts(shared_ptr<ElemList> eList, std::string id){

    EXCEPTION("REFACTOR");
  }

  
  void WaterWavePDE::DefineSurfaceIntegrators( ){
    //========================================================================================
    // ABC boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      //free surface condition for gravity waves
      ParamNodeList freeSurfaceNodes = bcNode->GetList( "freeSurfaceCondition" );
      for( UInt i = 0; i < freeSurfaceNodes.GetSize(); i++ ) {
        std::string regionName = freeSurfaceNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = freeSurfaceNodes[i]->Get("volumeRegion")->As<std::string>();

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        // gravity
        PtrCoefFct gravity = CoefFunction::Generate( mp_, Global::REAL, "9.81");
        PtrCoefFct dens = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );

        PtrCoefFct factor;
        factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");

        // factor for damping matrix: factor / gravity
        PtrCoefFct coeffMass =
        CoefFunction::Generate( mp_, Global::REAL,
                               CoefXprBinOp(mp_, factor, gravity, CoefXpr::OP_DIV ) );
        BiLinearForm * gravityInt = NULL;
        if( dim_ == 2 ) {
          gravityInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), coeffMass, 1.0, updatedGeo_ );
        } else {
          gravityInt = new BBInt<>(new IdentityOperator<FeH1,3,1>(), coeffMass, 1.0, updatedGeo_ );
        }

        gravityInt->SetName("gravityWaveIntegrator");
        BiLinFormContext *gravityContext = new BiLinFormContext(gravityInt, MASS);

        gravityContext->SetEntities( actSDList, actSDList );
        gravityContext->SetFeFunctions( feFunctions_[WATER_PRESSURE] , feFunctions_[WATER_PRESSURE]);
        feFunctions_[WATER_PRESSURE]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( gravityContext );
      } //free surface condition for gravity waves

      //PML for free surface condition for gravity waves
      ParamNodeList pmlfreeSurfaceNodes = bcNode->GetList( "pml4FreeSurfaceCondition" );
      for( UInt i = 0; i < pmlfreeSurfaceNodes.GetSize(); i++ ) {
        if(analysistype_ != HARMONIC)
          EXCEPTION("pml4FreeSurfaceCondition just working for harmonic anlaysis");

        std::string regionName = pmlfreeSurfaceNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = pmlfreeSurfaceNodes[i]->Get("volumeRegion")->As<std::string>();

        RegionIdType volRegion = ptGrid_->GetRegion().Parse(volRegName);

        // create new entity list
        shared_ptr<ElemList> volSDList( new ElemList(ptGrid_ ) );
        volSDList->SetRegion( volRegion );

        PtrCoefFct factor;
        factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");

        // compute surface wave velocity
        PtrCoefFct gravity = CoefFunction::Generate( mp_, Global::REAL, "9.81");
        PtrCoefFct omega = CoefFunction::Generate( mp_, Global::REAL, "2*pi*f");
        PtrCoefFct gravityC0 = CoefFunction::Generate( mp_,  Global::REAL,
                                   CoefXprBinOp( mp_, gravity, omega, CoefXpr::OP_DIV));

        //just coef-functions for the transformation of jacobians
        shared_ptr<CoefFunction> coeffPMLScal;
        shared_ptr<CoefFunction> coeffPMLfactor;
        shared_ptr<CoefFunction> coeffPMLMass;
        if( dampingList_[volRegion] == PML ) {
          std::string dampId = pmlfreeSurfaceNodes[i]->Get("dampingId")->As<std::string>();
          if(analysistype_ == HARMONIC){
            PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",dampId.c_str());
            coeffPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode,gravityC0,volSDList,regions_,false));
            coeffPMLfactor = CoefFunction::Generate( mp_, Global::COMPLEX,
                                              CoefXprBinOp(mp_, factor, coeffPMLScal,CoefXpr::OP_MULT));
          }else
            EXCEPTION("PML for gravity waves just for harmonic case");
        }
        else
          EXCEPTION("pml4FreeSurfaceCondition needs definition of dampingId");

        // factor for damping matrix: factor / gravity
        PtrCoefFct coeffMassPML =
            CoefFunction::Generate( mp_, Global::COMPLEX,
                               CoefXprBinOp(mp_, coeffPMLfactor, gravity, CoefXpr::OP_DIV ) );
        BiLinearForm * gravityIntPML = NULL;
        if( dim_ == 2 ) {
          gravityIntPML = new BBInt<Complex>(new IdentityOperator<FeH1,2,1>(), coeffMassPML, 1.0, updatedGeo_ );
        } else {
          gravityIntPML = new BBInt<Complex>(new IdentityOperator<FeH1,3,1>(), coeffMassPML, 1.0, updatedGeo_ );
        }

        gravityIntPML->SetName("gravityWaveIntegratorPML");
        BiLinFormContext *gravityContextPML = new BiLinFormContext(gravityIntPML, MASS);

        gravityContextPML->SetEntities( actSDList, actSDList );
        gravityContextPML->SetFeFunctions( feFunctions_[WATER_PRESSURE] , feFunctions_[WATER_PRESSURE]);
        feFunctions_[WATER_PRESSURE]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( gravityContextPML );
        //std::cout << "Have added gravityWaveIntegratorPML\n" << std::endl;
      }
    }
  }


  void WaterWavePDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(waterWavepde) << "Defining rhs load integrators for acoustic PDE";
    
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[WATER_PRESSURE];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> empty;
    LinearForm * lin = NULL;

    // obtain density
//    shared_ptr<CoefFunctionMulti> densFct  = matCoefs_[ELEM_DENSITY];
//    shared_ptr<CoefFunctionSurf> surfDens(new CoefFunctionSurf(false));
//    surfDens->SetVolumeCoefs( densFct->GetRegionCoefs() );
    PtrCoefFct surfDens = CoefFunction::Generate( mp_, Global::REAL, "1.0");
           
    // In the case of acou-mech coupling we have to multiply the
    // integrators by -densiy
    Double scalFactor = 1.0;
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

    bool coefUpdateGeo;
    

    // ===========================
    //  general surface load
    // ===========================
    ReadRhsExcitation( "surfaceLoad", empty, ResultInfo::SCALAR, isComplex_,
                       ent, coef,coefUpdateGeo );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("surfaceLoad can only be defined on surface elements");
      }

      PtrCoefFct exValue;
      if ( isMechCoupled_ == true ) {
        scalFactor = -1.0;
        exValue =
            CoefFunction::Generate( mp_, part,
                                   CoefXprBinOp(mp_, coef[i],surfDens, CoefXpr::OP_MULT) );
      } else {
        exValue = coef[i];
      }
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
                                                scalFactor, exValue, volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,2,1>(),
                                               scalFactor, exValue, volRegions, coefUpdateGeo);
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
    } // general surface load

    // =====================================
    //  rhsValues for e.g. for aeroacoustics
    // =====================================
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      coef[i]->SetConservative(true);
      this->rhsFeFunctions_[WATER_PRESSURE]->AddLoadCoefFunction(coef[i], ent[i]);
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
    rhs->definedOn = results_[0]->definedOn;
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

  }

  void WaterWavePDE::DefinePostProcResults(){
  }

  //! Init the time stepping
  void WaterWavePDE::InitTimeStepping(){

    Double alpha = this->myParam_->Get("timeStepAlpha")->As<Double>();

    if(this->isTimeDomPML_){
      //basically the choice for alpha scheme needs to be done everytime we have
      //a damping matrix not just for PML

      //scheme for main unknown
      GLMScheme * scheme1 = new Newmark(0.5,0.25,alpha);
      GLMScheme * scheme2 = new Newmark(0.5,0.25,alpha);
      shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(scheme1,0));
      shared_ptr<BaseTimeScheme> vecScheme(new TimeSchemeGLM(scheme2,0));

      feFunctions_[WATER_PMLAUXVEC]->SetTimeScheme(vecScheme);
      feFunctions_[WATER_PRESSURE]->SetTimeScheme(acouScheme);

      if(!this->isAPML_ && dim_ == 3){
        GLMScheme * scheme3 = new Newmark(0.5,0.25,alpha);
        shared_ptr<BaseTimeScheme> scalScheme(new TimeSchemeGLM(scheme3,0));
        feFunctions_[WATER_PMLAUXSCALAR]->SetTimeScheme(scalScheme);
      }
    }else{
      //GLMScheme * scheme1 = new Newmark(0.8,0.4225,-0.3);
      //GLMScheme * scheme1 = new Newmark(0.6,0.3025,alpha);
      GLMScheme * scheme1 = new Newmark(0.5,0.25,alpha);
      shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(scheme1,0));
      feFunctions_[WATER_PRESSURE]->SetTimeScheme(acouScheme);
    }
  }

}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template void WaterWavePDE::DefineTransientPMLInts<2>(shared_ptr<ElemList>, std::string);
  template void WaterWavePDE::DefineTransientPMLInts<3>(shared_ptr<ElemList>, std::string);
#endif
