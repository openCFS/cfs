// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include <boost/serialization_hdf5/hdf5_oarchive.hpp>
#include <boost/serialization_hdf5/hdf5_iarchive.hpp>

#include "PDE/SinglePDE.hh"

// for coordinate handling
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"

#include "OLAS/algsys/basesystem.hh"
#include "OLAS/algsys/standardsys.hh"

// header for scripting
#include <def_use_scripting.hh>
#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

// header for logging
#include "DataInOut/Logging/cfslog.hh"

// header for Materialhandling
#include "DataInOut/MaterialHandler.hh"

// header for Solvestep and assemble
#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"
#include "Driver/singleDriver.hh"
#include "Driver/transientdriver.hh"
#include "Driver/harmonicDriver.hh"

// header for iterative coupling
#include "CoupledPDE/pdecoupling.hh"

// header for memento/restart handling
#include "MatVec/vectorSerialization.hh"
#include "pdememento.hh"

// header for resultHandling
#include "DataInOut/resultHandler.hh"
#include "DataInOut/ResultCache.hh"

#include "DataInOut/postProc.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "Forms/linearForm.hh"

using std::string;

namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(pde)
  DEFINE_LOG(pde, "pde")

  SinglePDE::SinglePDE( Grid *aptgrid, ParamNode* paramNode ) :
    StdPDE( aptgrid, paramNode ),
    ptError_(NULL),
    tolSpaceErr_(0.0),
    pdeId_(NO_PDE_ID),
    isDirectCoupled_(false),
    isInitialized_(false),
    usePenalty_(true),
    maxTimeDerivOrder_(0),
    mHandle_(domain->GetMathParser()->GetNewHandle()) // Obtain mathParser handler
  {
    // Register functions for scripting
    RegisterFunctions();
  }


  // **************
  //   Destructor
  // **************
  SinglePDE::~SinglePDE() {


    // Delete algebraic system only if
    // PDE is not direct coupled
    if ( isDirectCoupled_ == false ) {
      if ( needsAlgsys_ && (algsys_ != NULL))
        delete algsys_;
      if (solveStep_) delete solveStep_;
      if (TS_alg_) delete TS_alg_;
      if (assemble_) delete assemble_;
    }


    if (sol_) delete sol_;
    if( !isDirectCoupled_ ) {
      if (solVec_) delete solVec_;
      if (rhsVec_) delete rhsVec_;
    }

    if ( needSolPrev_ && (solPrev_ != NULL))
      delete solPrev_;


    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      if (it->second != NULL) delete it->second;
    }
    materials_.clear();

  }


  // ********
  //   Init
  // ********
  void SinglePDE::Init( UInt sequenceStep, InfoNode* base ) {

    sequenceStep_ = sequenceStep;

    infoNode_ = base == NULL ? info->Get("PDE")->Get(pdename_) : base->Get(pdename_);
    infoNode_->Get(InfoNode::HEADER)->Get("sequenceStep")->SetValue(sequenceStep);

    StdVector<RegionIdType> allIDs;
    LOG_TRACE(pde) << pdename_ << ": Starting Initialization";


    // =====================================================================
    // Get type of analysis
    // =====================================================================
    LOG_TRACE(pde) << pdename_ << ": Obtaining analysistye";
    analysistype_ = domain->GetSingleDriver()->GetAnalysisType();

    // NOTE: The concept of isAlwaysStatic bites with Direct Coupling
    //       and must be re-designed
    if ( isAlwaysStatic_ == true &&
         analysistype_ == TRANSIENT ) {
      analysistype_ = STATIC;
    }

    isComplex_ = IsComplex(analysistype_);


    // =====================================================================
    // get regions/subdomains for PDE
    // =====================================================================

    LOG_TRACE(pde) << pdename_ << ": Obtaining regions";


    // Obtain regions the pde is defined on
    StdVector<ParamNode*> regionNodes =
      myParam_->Get("regionList")->GetList("region");

    // output to info-file
    InfoNode* list = infoNode_->Get(InfoNode::HEADER);

    // output and set subdoms_
    for( UInt i = 0; i < regionNodes.GetSize(); i++ )
    {
      InfoNode* in_ = list->Get("region");
      in_->Get("name")->SetValue(regionNodes[i]->Get("name")->AsString());
      bool complexMat = false;
      regionNodes[i]->Get("complexMaterial",  complexMat, false );
      
      RegionIdType actRegionId = ptgrid_->GetRegion().Parse(regionNodes[i]->Get("name")->AsString());
      subdoms_.Push_back( actRegionId );
      complexMatData_[actRegionId] = complexMat;
    }


    // Generate a fitting algebraic system only if PDE is NOT
    // direct coupled
    if( needsAlgsys_ == true ) {
      if ( isDirectCoupled_ == false) {
        algsys_ = new StandardSystem(FindLinearSystem(pdename_));
      }

      // Get parameter and report object of OLAS
      olasParams_ = algsys_->GetOLASParams();
      olasReport_ = algsys_->GetOLASReport();

      // Obtain unique pde identifier
      pdeId_ = algsys_->ObtainPDEId( pdename_ );


      // Determine, if this is a parallel run
      // and pass this information to OLAS
      bool parallel = false;

      olasParams_->SetValue( "Parallel", parallel );
    }

    // =====================================================================
    // create assemble class
    // =====================================================================

    // Create new assemble class with according analysistype
    if( isDirectCoupled_ == false && needsAlgsys_ == true) {
      assemble_ = new Assemble( algsys_, analysistype_, maxTimeDerivOrder_ );
    }

    // Determine if solution is of complex type or not
    if ( analysistype_ == HARMONIC ) {
      sol_ = new NodeStoreSol<Complex>;
      //if(!isDirectCoupled_ ) {
        solVec_ = new Vector<Complex>;
        rhsVec_ = new Vector<Complex>;
      //}
    }
    else {
      sol_ = new NodeStoreSol<Double>;

     // if(!isDirectCoupled_ ) {
        solVec_ = new Vector<Double>;
        rhsVec_ = new Vector<Double>;

      //}

      if ( needSolPrev_ ) {
        solPrev_ = new NodeStoreSol<Double>;
      }
    }

    // =====================================================================
    // read in material data
    // =====================================================================
    LOG_TRACE(pde) << pdename_ << ": Reading material information";
    ReadMaterialData();

    // =====================================================================
    // read in damping information
    // =====================================================================
    LOG_TRACE(pde) << pdename_ << ": Reading damping information";
    ReadDampingInformation( );

    // =====================================================================
    // read in NonLinearities
    // =====================================================================
    LOG_TRACE(pde) << pdename_ << ": Initializing non-linearities";
    InitNonLin();

    // =====================================================================
    // trigger definition of available results
    // =====================================================================
    DefineAvailResults();

    // =====================================================================
    // initialize EQN-object and Storeresults class
    // =====================================================================

    // Name of linear system depends of coupling type
    std::string systemName = pdename_;
    if ( isDirectCoupled_ )
      systemName  = "direct";

    // How do we want to treat inhomogeneous Dirichlet boundary conditions?
    {
      std::string aux = "penalty";
      ParamNode * systemsNode =
        param->Get("sequenceStep", std::string("index"), sequenceStep_)
        ->Get("linearSystems",false);
      if( systemsNode ) {
        ParamNode * mySystemNode = systemsNode->
          Get("system", "name", systemName, false);
        if( mySystemNode) {
          ParamNode * setupNode = mySystemNode->Get("setup", false);
          if( setupNode ) {
            setupNode->Get("idbcHandling", aux, false );
          }
        }
      }
      usePenalty_ = aux == "penalty" ? true : false;
      infoNode_->Get(InfoNode::HEADER)->Get("idbc", "handling of IDBCs")->SetValue(aux);
    }

    //determine which kind of equationmap will be needed
    if(results_.GetSize() == 1){
      if(results_[0]->fctType->IsDiscontinuous()){
        eqnMap_ = shared_ptr<DiscontinuousEqnMap>(new DiscontinuousEqnMap( ptgrid_, pdeId_, usePenalty_ ));
      }else{
        eqnMap_ = shared_ptr<EqnMap>(new EqnMap( ptgrid_, pdeId_, usePenalty_ ));
      }
    }else if(results_.GetSize() == 2){
      if(!results_[0]->fctType->IsDiscontinuous() && results_[1]->fctType->IsDiscontinuous()){
        eqnMap_ = shared_ptr<MixedEqnMap>(new MixedEqnMap( ptgrid_, pdeId_, usePenalty_ ));
      }else if(results_[0]->fctType->IsDiscontinuous() && results_[1]->fctType->IsDiscontinuous()){
        eqnMap_ = shared_ptr<DiscontinuousEqnMap>(new DiscontinuousEqnMap( ptgrid_, pdeId_, usePenalty_ ));
      }else{
        eqnMap_ = shared_ptr<EqnMap>(new EqnMap( ptgrid_, pdeId_, usePenalty_ ));
      }
    }else{
      //this is the defulat case
      //more cases can be implemented as needed
      eqnMap_ = shared_ptr<EqnMap>(new EqnMap( ptgrid_, pdeId_, usePenalty_ ));
    }




    // Create a new equation map


    // =====================================================================
    // read in boundary conditions
    // =====================================================================

    // Incorporate values of memento here, if the values are used
    // as "dirichlet"-values from a previous simulation run
    // (e.g. obtained from a previous static run)
    if( memento_ != NULL
        && mementoAsDirichlet_ == true ) {
      IncorporateMemento();
    }

    LOG_TRACE(pde) << pdename_ << ": Reading boundary conditions";
    ReadBCs();
    ReadSpecialBCs();
    //
#ifdef USE_SCRIPTING
    // Trigger event for scripting
    StdVector<std::string> args;
    args.Push_back( pdename_ );
    messenger->TriggerEvent( CFSMessenger::CFS_ReadBCs, args );
    args.Clear();
#endif

    // =====================================================================
    // Create time stepping algorithm
    // =====================================================================
    if ( analysistype_ == TRANSIENT &&
         isDirectCoupled_ == false) {
      InitTimeStepping();
      if ( TS_alg_ != NULL ) {
        Double dt;
        dt = dynamic_cast<TransientDriver*>(domain->GetSingleDriver())
          ->GetDeltaT();
      }
    }



    // =====================================================================
    // define the integrators for PDE and initialize eqn object
    // =====================================================================

    // Call initialization of (bi)linear integrators
    LOG_TRACE(pde) << pdename_ << ": Defining integrators";
    DefineIntegrators();

    // Print information about defined integrators
    if( !isDirectCoupled_ && needsAlgsys_ == true)
      assemble_->ToInfo(infoNode_->Get(InfoNode::HEADER)->Get("integrators"));

    // now we know about nonlinearities and we can trigger the
    // material objects to perform the approximations of the nonlinear
    // sampled data
    std::map<RegionIdType, BaseMaterial*>::iterator itMat;
    for ( itMat = materials_.begin(); itMat != materials_.end(); itMat++ ) {
      itMat->second->InitApproxCurves();
    }


    // Finish equation mapping
    LOG_TRACE(pde) << pdename_ << ": Mapping Equations";
    eqnMap_->SetHomoDirichletBCs ( hdBcs_ );
    eqnMap_->SetInhomDirichletBCs( idBcs_ );
    eqnMap_->SetConstraints( constraints_ );
    eqnMap_->Finalize();

    infoNode_->Get(InfoNode::HEADER)->Get("equations", "number of equations in linear system")->SetValue(eqnMap_->GetNumEqns());

    // writes numerous lists about mapping -> might be huge!
    if(progOpts->DoListMapping())
      eqnMap_->ToInfo(infoNode_->Get(InfoNode::HEADER)->Get("mapping"));

    numPDENodes_ = eqnMap_->GetNumLocalNodes();
    numElems_ = eqnMap_->GetNumLocalElems();

    // Initialize Storesolution class
    sol_->SetNumSolutions(results_.GetSize());
    sol_->SetNumNodes(numPDENodes_);
    for (UInt iSol=0; iSol<results_.GetSize(); iSol++) {
      sol_->SetSolutionType(results_[iSol]->resultType,iSol);
      sol_->SetNumDofs( results_[iSol]->dofNames.GetSize(),
                        results_[iSol]->resultType );
    }

    // Note: this is only a temporary solution
    sol_->SetPtrEQNData(eqnMap_.get(), ptgrid_);
    sol_->SetRegions( subdoms_ );
    sol_->Init();

    if ( solPrev_ != NULL ) {
      solPrev_->SetNumSolutions(results_.GetSize());
      solPrev_->SetNumNodes(numPDENodes_);
      for (UInt iSol=0; iSol<results_.GetSize(); iSol++) {
        solPrev_->SetSolutionType(results_[iSol]->resultType,iSol);
        solPrev_->SetNumDofs( results_[iSol]->dofNames.GetSize(),
                              results_[iSol]->resultType );
      }

      // Note: this is only a temporary solution
      solPrev_->SetPtrEQNData(eqnMap_.get(), ptgrid_);
      solPrev_->SetRegions( subdoms_ );
      solPrev_->Init();
    }

    if(!isDirectCoupled_ ) {
      solVec_->Resize(eqnMap_->GetNumEqns());
      solVec_->Init();
      rhsVec_->Resize(eqnMap_->GetNumEqns());
      rhsVec_->Init();

      if ( analysistype_ == HARMONIC ) {
        sol_->SetAlgSysDataPointer(solVec_->GetSize(),
                                   dynamic_cast<Vector<Complex>&>(*solVec_).GetPointer() );
      } else {
        sol_->SetAlgSysDataPointer(solVec_->GetSize(),
                                   dynamic_cast<Vector<Double>&>(*solVec_).GetPointer() );

      }
    }


    if ( analysistype_ == TRANSIENT &&
         isDirectCoupled_ == false) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain->GetSingleDriver())
        ->GetDeltaT();
      TS_alg_->Init( dt, eqnMap_->GetNumEqns() );
    }

//    // =====================================================================
//    // Set the initial conditions
//    // =====================================================================
//    if ( analysistype_ == TRANSIENT){
//    	SetInitialCondition();
//    }

    // =====================================================================
    // define which solution types have to be saved
    // =====================================================================
    LOG_TRACE(pde) << pdename_ << ": Reading store results";

    ReadStoreResults();
    ReadSpecialResults();

    // Set information at sol_ object
    sol_->SetResults( results_ );


    if ( solPrev_ != NULL )
      solPrev_->SetResults( results_ );

    PreparePDE4Computation();

    //! Define step solution driver
    if ( isDirectCoupled_ == false ) {
      LOG_TRACE(pde) << pdename_ << ": Defining solveStep class";
      DefineSolveStep();
    }

    // Call event procedure for scripting
#ifdef USE_SCRIPTING
        // Trigger event for scripting
        args.Clear();
        args.Push_back( pdename_ );
        messenger->TriggerEvent( CFSMessenger::CFS_PdeInit, args );
#endif

    // Finally set the initialization flag to true
    isInitialized_ = true;
    LOG_TRACE(pde) << pdename_ << ": Finished initializaton";
  }


  void SinglePDE::SaveSolution( const Double * ptSol, UInt size ) {

    Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVec_);
    solHelp.Resize(size);
    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }

    sol_->SetAlgSysDataPointer( size, solHelp.GetPointer() );

  }

  void SinglePDE::SaveSolution( const Complex * ptSol, UInt size ) {

    Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*solVec_);
    solHelp.Resize(size);
    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }

    sol_->SetAlgSysDataPointer( size, solHelp.GetPointer() );

  }

  void SinglePDE::SavePrevSolution( const Double * ptSolPrev, UInt size ) {

    Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVecPrev_);

    solHelp.Resize(size);

    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSolPrev[i];
    }

    solPrev_->SetAlgSysDataPointer( size, solHelp.GetPointer() );

  }

  void SinglePDE::SaveRHS( const Double * ptSol, UInt size ) {

     Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*rhsVec_);
     solHelp.Resize(size);
     for ( UInt i = 0; i < size; i++ ) {
       solHelp[i] = ptSol[i];
     }
   }

   void SinglePDE::SaveRHS( const Complex * ptSol, UInt size ) {

     Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*rhsVec_);
     solHelp.Resize(size);
     for ( UInt i = 0; i < size; i++ ) {
       solHelp[i] = ptSol[i];
     }
   }


   /** can generally be called multiple times. We overwrite old values! Brute force but keeps data size */
   void SinglePDE::WriteGeneralPDEdefines()
   {
    // loads
    InfoNode* base = infoNode_->Get(InfoNode::HEADER)->Get("loads");

    for(unsigned int i = 0, in = loads_.GetSize(); i < in; i++)
    {
      LoadBc const & actBc = *loads_[i];
      EntityList const & actList = *actBc.entities;

      InfoNode* in = base->Get(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
      in->Get("value")->SetValue(actBc.value);
      in->Get("phase")->SetValue(actBc.phase);
    }

    // Homogeneous Dirichlet BC
    base = infoNode_->Get(InfoNode::HEADER)->Get("homDirichletBC");

    for(unsigned int i = 0, in = hdBcs_.GetSize(); i < in; i++)
    {
      HomDirichletBc const & actBc = *hdBcs_[i];
      EntityList const & actList = *actBc.entities;

      InfoNode* in = base->Get(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
    }

    // Inhomogeneous Dirichlet BC
    base = infoNode_->Get(InfoNode::HEADER)->Get("inhomDirichletBC");

    for(unsigned int i = 0, in = idBcs_.GetSize(); i < in; i++)
    {
      InhomDirichletBc const & actBc = *idBcs_[i];
      EntityList const & actList = *actBc.entities;

      InfoNode* in = base->Get(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
      in->Get("value")->SetValue(actBc.value);
      in->Get("phase")->SetValue(actBc.phase);
    }

    // Inhomogeneous Neumann BC
    base = infoNode_->Get(InfoNode::HEADER)->Get("inhomNeumannBC");

    for(unsigned int i = 0, in = inBcs_.GetSize(); i < in; i++)
    {
      InhomNeumannBc const & actBc = *inBcs_[i];
      EntityList const & actList = *actBc.entities;

      InfoNode* in = base->Get(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
      in->Get("value")->SetValue(actBc.value);
      in->Get("phase")->SetValue(actBc.phase);
    }

    // constraints
    base = infoNode_->Get(InfoNode::HEADER)->Get("constraints");
    // periodic boundary conditions blow this up.
    if(constraints_.GetSize() <= 5 || progOpts->DoListMapping())
    {
      for(unsigned int i = 0, in = constraints_.GetSize(); i < in; i++)
      {
        Constraint const & actBc = *constraints_[i];
        EntityList const & masterList = *actBc.masterEntities;
        EntityList const & slaveList = *actBc.slaveEntities;

        InfoNode* in = base->Get("pair", "master", masterList.GetName());
        in->Get("slave")->SetValue(slaveList.GetName());
        // the names are repeated for the different dofs
        std::string dof = actBc.result->GetDofName(actBc.masterDof);
        if(!in->Has("dof", dof))
          in->Get("dof", InfoNode::APPEND)->SetValue(dof);

        in->Get("periodic")->SetValue(actBc.periodic);
      }
    }
    else
    {
      if(constraints_.GetSize() > 5)
      {
        base->Get("number")->SetValue(constraints_.GetSize());
        base->Get(InfoNode::COMMENT)->SetValue("run cfs with -l for list");
      }
    }
  }

  void SinglePDE::ReadStoreResults() {

    ResultSet::iterator it;

    // fetch result node and leave, if none is present
    ParamNode * resultNode = myParam_->Get("storeResults", false);
    if( !resultNode )
      return;

    // Iterate over all availabe results
    for (it = availResults_.begin(); it != availResults_.end(); it++ ) {
       CheckStoreResult(*it);
    }
    Info->PrintF( pdename_, "\n");

  }


  bool SinglePDE::CheckStoreResult(shared_ptr<ResultInfo> candidate)
  {

    StdVector<std::string> regionNames, nodeNames, writeResults, actOutDest;
    StdVector<std::string> postProcNames, outDestNames, neighborRegions;
    UInt saveBegin, saveEnd, saveInc;
    std::string quantity, complexFormatString, listElemName, entityName;
    ComplexFormat complexFormat;
    shared_ptr<EntityList> actList;

    EntityList::ListType entityType;
    EntityList::DefineType defineType;
    ResultHandler * resHandler = domain->GetResultHandler();

    // initialize map for relating EntityUnknownType and name of xml-element
    using std::make_pair;
    std::map<ResultInfo::EntityUnknownType, std::string> elemNames;
    std::map<ResultInfo::EntityUnknownType, bool> isHistory;
    elemNames.insert(make_pair(ResultInfo::NODE, "nodeResult"));
    elemNames.insert(make_pair(ResultInfo::PFEM, "nodeResult"));
    elemNames.insert(make_pair(ResultInfo::ELEMENT, "elemResult"));
    elemNames.insert(make_pair(ResultInfo::SURF_ELEM, "surfElemResult"));
    elemNames.insert(make_pair(ResultInfo::REGION, "regionResult"));
    elemNames.insert(make_pair(ResultInfo::SURF_REGION, "surfRegionResult"));

    isHistory.insert(make_pair(ResultInfo::NODE, false));
    isHistory.insert(make_pair(ResultInfo::PFEM, false));
    isHistory.insert(make_pair(ResultInfo::ELEMENT, false));
    isHistory.insert(make_pair(ResultInfo::SURF_ELEM, false));
    isHistory.insert(make_pair(ResultInfo::REGION, true));
    isHistory.insert(make_pair(ResultInfo::SURF_REGION, true));
    

    // fetch result node and leave, if none is present
    ParamNode * resultNode = myParam_->Get("storeResults", false);
    if( !resultNode )
      return false;

      // Convert enum
      quantity = SolutionTypeEnum.ToString(candidate->resultType);
      LOG_DBG(pde) << "Searching for storeResults of quantity '"
                   << quantity << "'";

      // try to catch possible errors
      try {

        // Get type of result
        std::string xmlElemName = elemNames[candidate->definedOn];
        if( xmlElemName == "" ){
          return false;
        }

        // Remember current result node
        ParamNode* actResultNode =
          resultNode->Get(xmlElemName, "type", quantity, false );

        // Check on which entity type the result is defined on
        switch(candidate->definedOn)
        {
        case ResultInfo::NODE:
        case ResultInfo::PFEM:
          entityType = EntityList::NODE_LIST;
          break;
        case ResultInfo::REGION:
        case ResultInfo::SURF_REGION:
          entityType = EntityList::REGION_LIST;
          break;
        case ResultInfo::SURF_ELEM:
          entityType = EntityList::SURF_ELEM_LIST;
          break;
        case ResultInfo::ELEMENT:
          entityType = EntityList::ELEM_LIST;
          break;
        default:
          EXCEPTION("Type of 'definedOn' was not found");
        }

        // intialize variables
        neighborRegions.Clear();
        regionNames.Clear();
        //outDestNames.Clear();

        // ========== Look for defineType 'REGION' ==========
        // 1a) Look if result is defined on 'allRegions'
        defineType = EntityList::REGION;



        // if no node was found, continue with next result
        if( !actResultNode) {
          return false;
        }

        // determine complexFormat
        complexFormatString = actResultNode->Get("complexFormat")->AsString();
        String2Enum( complexFormatString, complexFormat );

        // otherwise check, if result is to be saved on "allRegions"
        if( actResultNode->Has("allRegions" ) ) {
          ptgrid_->GetRegion().ToString(subdoms_,regionNames);

          ParamNode * allRegionsNode = actResultNode->Get("allRegions");

          std::string allPostProcName, allOutDestName;

          // fetch postProcNames
          allRegionsNode->Get("postProcId", allPostProcName );
          postProcNames.Resize( regionNames.GetSize() );
          postProcNames.Init( allPostProcName );

          //fetch outDestName
          allRegionsNode->Get("outputIds", allOutDestName );
          outDestNames.Resize( regionNames.GetSize() );
          outDestNames.Init( allOutDestName );

          // fetch saveBegin, saveEnd and saveInc
          allRegionsNode->Get("saveBegin", saveBegin );
          allRegionsNode->Get("saveEnd", saveEnd );
          allRegionsNode->Get("saveInc", saveInc );

          // fetch writeResult flag
          std::string writeResult;
          allRegionsNode->Get("writeResult", writeResult );
          writeResults.Resize( regionNames.GetSize() );
          writeResults.Init( writeResult );

        } else {

          StdVector<ParamNode*> regionNodes;
          ParamNode * listNode = NULL;
          // 1b) Look for regions the result is defined on
          if(candidate->definedOn == ResultInfo::NODE ||
             candidate->definedOn == ResultInfo::PFEM ||
             candidate->definedOn == ResultInfo::ELEMENT ||
             candidate->definedOn == ResultInfo::REGION ) {
            listNode = actResultNode->Get("regionList", false);
            if( listNode )
              regionNodes = listNode->GetList("region");
          } else if(candidate->definedOn == ResultInfo::SURF_ELEM ||
                    candidate->definedOn == ResultInfo::SURF_REGION ) {
            listNode = actResultNode->Get("surfRegionList", false);
            if( listNode )
              regionNodes = listNode->GetList("surfRegion");

            // fetch entry with neighboring regions
            for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
              neighborRegions.Push_back( regionNodes[i]->
                                         Get("neighborRegion")->AsString() );
            }
          }

          // only enter, at least one region is present
          if( listNode ) {
            // fetch saveBegin, saveEnd and saveInc
            listNode->Get( "saveBegin", saveBegin );
            listNode->Get( "saveEnd", saveEnd );
            listNode->Get( "saveInc", saveInc );

            // iterate over all regions
            regionNames.Clear();
            postProcNames.Clear();
            outDestNames.Clear();
            writeResults.Clear();
            for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
              regionNames.Push_back( regionNodes[i]->Get("name")->AsString() );
              postProcNames.Push_back( regionNodes[i]->Get("postProcId")->AsString() );
              outDestNames.Push_back( regionNodes[i]->Get("outputIds")->AsString() );
              writeResults.Push_back( regionNodes[i]->Get("writeResult")->AsString() );
            }
          }
        }

        // Check, if any region was found for this result type
        if( regionNames.GetSize() != 0 ) {
          candidate->complexFormat = complexFormat;

          // iterate over all regions
          for( UInt iRegion = 0; iRegion < regionNames.GetSize(); iRegion++ )
          {
            actList = ptgrid_->GetEntityList( entityType, regionNames[iRegion],
                                              defineType );
            shared_ptr<BaseResult> actSol;
            if( isComplex_ ) {
              actSol = shared_ptr<BaseResult>(new Result<Complex>());
            } else {
              actSol = shared_ptr<BaseResult>(new Result<Double>());
            }

            // intialize result object
            actSol->SetResultInfo(candidate);
            actSol->SetEntityList( actList );
            resultLists_[candidate].Push_back( actSol );

            // extract all output destinations and determine bool flag for writeResult
            SplitStringList( outDestNames[iRegion], actOutDest, ',' );
            bool writeResult = writeResults[iRegion] == "yes"  ? true : false ;

            // pass result to resulthandler
            resHandler->RegisterResult( actSol, saveBegin, saveInc, saveEnd,
                                        actOutDest,
                                        postProcNames[iRegion], writeResult,
                                        isHistory[candidate->definedOn] );

            // if neighboring region is present, store related volume
            // neighbor region
            if( neighborRegions.GetSize() > 0 ) {
              if( !neighborRegions[iRegion].empty() ) {
                surfNeighborRegions_[actSol] =
                  ptgrid_->GetRegion().Parse( neighborRegions[iRegion] );
              }
            }

          }
        }


        // ========== Look for defineType  node/elemList (history)  ==========

        std::string entityTypeName;
        StdVector<std::string> histNames;
        neighborRegions.Clear();

        ParamNode * histNode = NULL;
        StdVector<ParamNode*> histEntities;

        if(candidate->definedOn == ResultInfo::NODE
            || candidate->definedOn == ResultInfo::PFEM ) {
          defineType = EntityList::NAMED_NODES;
          histNode = actResultNode->Get("nodeList",false);
          if( histNode )
            histEntities = histNode->GetList("nodes");
          entityTypeName = "nodes";

        } else if(candidate->definedOn == ResultInfo::ELEMENT ) {
          defineType = EntityList::NAMED_ELEMS;
          histNode = actResultNode->Get("elemList",false);
          if( histNode )
            histEntities = histNode->GetList("elems");
          entityTypeName = "elements";

        } else if(candidate->definedOn == ResultInfo::SURF_ELEM ) {
          defineType = EntityList::NAMED_ELEMS;
          histNode = actResultNode->Get("surfElemList",false);
          if( histNode)
            histEntities = histNode->GetList("surfElems");
          entityTypeName = "surfElems";

          // fetch entry with neighboring regions
          // fetch entry with neighboring regions
          for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
            neighborRegions.Push_back( histEntities[i]->
                                       Get("neighborRegion")->AsString() );
          }
        }

        // only proceed, if any history result is defined+
        if( histNode ) {

          // fetch saveBegin, saveEnd and saveInc
          histNode->Get("saveBegin", saveBegin );
          histNode->Get("saveEnd", saveEnd );
          histNode->Get("saveInc", saveInc );

          // iterate over all regions
          histNames.Clear();
          postProcNames.Clear();
          outDestNames.Clear();
          writeResults.Clear();
          for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
            histNames.Push_back( histEntities[i]->Get("name")->AsString() );
            postProcNames.Push_back( histEntities[i]->Get("postProcId")->AsString() );
            outDestNames.Push_back( histEntities[i]->Get("outputIds")->AsString() );
            writeResults.Push_back( histEntities[i]->Get("writeResult")->AsString() );
          }
        }

        if( histNames.GetSize() > 0 )
        {
          candidate->complexFormat = complexFormat;

          // iterate over all entityNames
          for( UInt i = 0; i < histNames.GetSize(); i++ )
          {
            actList = ptgrid_->GetEntityList( entityType,
                                              histNames[i], defineType );
            shared_ptr<BaseResult> actSol;
            if( isComplex_ ) {
              actSol = shared_ptr<BaseResult>(new Result<Complex>());
            } else {
              actSol = shared_ptr<BaseResult>(new Result<Double>());
            }

            // Set result info and entitylist at the result object
            actSol->SetResultInfo(candidate);
            actSol->SetEntityList( actList );
            resultLists_[candidate].Push_back( actSol );

            // extract all output destinations and determine bool flag for writeResult
            SplitStringList( outDestNames[i], actOutDest, ',' );
            bool writeResult = (writeResults[i] == "yes"  ? true : false );

            resHandler->RegisterResult( actSol, saveBegin, saveInc, saveEnd,
                                        actOutDest,
                                        postProcNames[i], writeResult, true );

            // if neighboring region is present, store related volume
            // neighbor region
            if( neighborRegions.GetSize() > 0 ) {
              if( !neighborRegions[i].empty() ) {
                surfNeighborRegions_[actSol] =
                  ptgrid_->GetRegion().Parse( neighborRegions[i] );
              }
            }

          }
        }
      } catch( Exception &ex ) {
        RETHROW_EXCEPTION(ex, "Could not determine storeResults for quantity '"
                          << quantity << "' within pde '" << pdename_ << "'" );
      }
    return true;
  }

  void SinglePDE::WriteResultsInFile( const UInt kstep,
                                      const Double actTimeFreq ) {
    LOG_DBG(pde) << "WriteResultsInFile() kstep: " <<  kstep
                 << " actTimeFreq: " << actTimeFreq;
    ResultMap::iterator it = resultLists_.begin();
    ResultHandler * resHandler = domain->GetResultHandler();

    // iterate over all results
    for( ; it != resultLists_.end(); it++ ) {
      ResultList & actList = it->second;

      // iterate over all solutions for each result type
      for( UInt i = 0; i < actList.GetSize(); i++ ) {

        // get string representation of quantity and entity list
        std::string listName;
        std::string quantity = SolutionTypeEnum.ToString(actList[i]->GetResultInfo()->resultType);
        listName = actList[i]->GetEntityList()->GetName();

        // Only calculate result, if needed
        if( resHandler->IsResultNeeded( actList[i] ) ) 
        {
          try 
          {
            CalcResults( actList[i] );
          } 
          catch (Exception &ex ) 
          {
            RETHROW_EXCEPTION( ex, "Could not calculate result '" << quantity
                               << "' on '" << listName << "' in pde '" << pdename_ << "'");
          }
          try 
          {
            resHandler->UpdateResult( actList[i] );
          } 
          catch (Exception &ex ) 
          {
            RETHROW_EXCEPTION( ex, "Could not write result '" << quantity
                               << "' on '" << listName
                               << "' to output file(s) in pde '" << pdename_ << "'");
          }
        }
      }
    }

#ifdef USE_SCRIPTING
    StdVector<std::string> context;
    context.Push_back( pdename_ );
    context.Push_back( lexical_cast<std::string>(solveStep_->GetActStep() ) );

    if ( analysistype_ == TRANSIENT ||
         analysistype_ == STATIC ) {
      context.Push_back( lexical_cast<std::string>(solveStep_->GetActTime() ) );
    } else {
      context.Push_back( lexical_cast<std::string>(solveStep_->GetActFreq() ) );
    }
    messenger->TriggerEvent( CFSMessenger::CFS_CalcResults,
                             context );
#endif
  }


  // **********
  //   SetBCs
  // **********
  void SinglePDE::SetBCs() {

     // Trigger setting of BC from script file
 #ifdef USE_SCRIPTING
     StdVector<std::string> context;
     context.Push_back( pdename_ );
     context.Push_back( lexical_cast<std::string>(solveStep_->GetActStep() ) );

     if ( analysistype_ == TRANSIENT ||
          analysistype_ == STATIC ) {
       context.Push_back( lexical_cast<std::string>(solveStep_->GetActTime() ) );
     } else {
       context.Push_back( lexical_cast<std::string>(solveStep_->GetActFreq() ) );
     }
     messenger->TriggerEvent( CFSMessenger::CFS_SetBCs, context );
 #endif

     UInt dof;
     Double val;
     StdVector<UInt> nodes;
     Integer eqnNr;
     Vector<Double> globCoord;

     // get global coordinate system and math parser
     CoordSystem * coosy = domain->GetCoordSystem();
     MathParser * parser = domain->GetMathParser();


     // ---------------------------
     // INHOMOGENEOUS DIRICHLET BC
     // ---------------------------

     Double stepVal, phase = 0.0;
     StdVector<Double>  val_tfunc_vec, dirVal_vec;

     // get step value according to type of analysis
     if (analysistype_ == HARMONIC)
       stepVal = solveStep_->GetActFreq();
     else
       stepVal = solveStep_->GetActTime();

     for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {

       // Get grip of actual idBC
       InhomDirichletBc const & actBc = *(idBcs_[i]);

       dof = actBc.dof;
       // std::cout << "dof=" << dof << std::endl;

       // Get EntityIterator
       EntityIterator it = actBc.entities->GetIterator();

       // set info for input function
       ResultCache::SetInfo(ResultCache::OUT_REAL,
                            dof,
                            actBc.entities->GetName(),
                            actBc.result->resultType,
                            stepVal);

       for ( it.Begin(); !it.IsEnd(); it++ ) {

         try {
           StdVector<Integer> eqns;
           eqnMap_->GetEqns( eqns, *actBc.result, it, dof  );

           // loop over all equations for the dirichlet conditions
           for( UInt iEqn = 0; iEqn < eqns.GetSize(); iEqn++){
             eqnNr = eqns[iEqn];

             // omit all equations, which are homogeneous dirichlet
             // boundary condition (0) or constraint slave eqn (<0)
             if( eqnNr <= 0 ) continue;

             // If iterator points to a node, pass the current coordinate
             // to the parser
             if( it.GetType() == EntityList::NODE_LIST ) {
               // Get node coordinate
               ptgrid_->GetNodeCoordinate( globCoord, it.GetNode() );
               parser->SetCoordinates( mHandle_, *coosy, globCoord );
               ResultCache::SetIndex(it.GetPos());
             }
             else {
               // this case needs to be implemented ...
             }

             // Now evaluate value of IDBC
             ResultCache::SetOutputType(ResultCache::OUT_AMPL);
             parser->SetExpr( mHandle_, actBc.value );
             val = parser->Eval( mHandle_ );

             // Sanity check. This should not happen, but might appear
             // in the case that the same node/dof belongs to a region
             // with hom. and a region with inhom. Dirichlet BCs. This
             // problem was already encountered!
             if (eqnNr == 0) {

               EXCEPTION( "Got eqn number 0 for inhom Dirichlet BC! "
                          << "Probably you have a node/dof that belongs to both "
                          << "a region with hom. and one with inhom. Dirichlet BCs."
                          << " Go check your .mesh file!" );
             }

             // Transform Dirichlet boundary conditions for effmass-formulation
             if (effectiveMass_) {
               val = TS_alg_->DirichletBC4EffMassMatrix(val,eqnNr);
             }

             // Case of complex-valued entries
             if (analysistype_ == HARMONIC ) {

               ResultCache::SetOutputType(ResultCache::OUT_PHASE);
               parser->SetExpr( mHandle_, actBc.phase );
               phase = parser->Eval( mHandle_ );
               Complex complexValue( val * cos( phase / 180 * PI ),
                                     val * sin( phase / 180 * PI ) );
               algsys_->SetDirichlet(  pdeId_, eqnNr, complexValue);
             }
             else {
               //   std::cout << "IHDBC val=" << val << std::endl;
               algsys_->SetDirichlet( pdeId_, eqnNr, val );
             }
           }
         } catch (Exception& ex ) {
           RETHROW_EXCEPTION(ex, pdename_ << ": Could not apply Inhom. Dirichlet boundary condition "
                             << " for nodes '" <<  actBc.entities->GetName()
                             << "' with value '" << actBc.value << "'" );
         }
       }
     }
  }



  void SinglePDE::ReadBCs() {


    // fetch "bcsAndLoads" parameter node, if present.
    // otherwise leave
    ParamNode * bcsNode = myParam_->Get("bcsAndLoads", false );
    if( !bcsNode )
      return;

    std::string name, resultName, dof, entType, value, phase;
    EntityList::DefineType defineType;
    shared_ptr<ResultInfo> actResultInfo;


    // =====================================================================
    // homogeneous Dirichlet BC
    // =====================================================================

    // fetch paramnodes for hdbc
    StdVector<ParamNode*> hdbcNodes = bcsNode->GetList("dirichletHom");

    // iterate over all parameter nodes
    for( UInt i = 0; i < hdbcNodes.GetSize(); i++ ) {

      try {
        // read parameters
        dof.clear();
        hdbcNodes[i]->Get( "name", name );
        hdbcNodes[i]->Get( "quantity", resultName );
        hdbcNodes[i]->Get( "dof", dof, false );
        hdbcNodes[i]->Get( "entityType", entType );

        // fetch related resultInfo object
        actResultInfo = GetResultInfo( SolutionTypeEnum.Parse(resultName) );

        // Create homogeneous boundary condition
        if( entType == "nodeList" ) {
          defineType = EntityList::NAMED_NODES;
        } else {
          defineType = EntityList::REGION;
        }

        shared_ptr<HomDirichletBc> actBc ( new HomDirichletBc );
        EntityList::ListType listType = EntityList::listType.Parse(entType);
        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( listType,
                                  name, defineType );
        actBc->entities = actList;
        actBc->result = actResultInfo;
        actBc->eqnMap = eqnMap_;
        if( dof.empty() ) {
          actBc->dof = 1;
        } else {
          actBc->dof = actResultInfo->GetDofIndex( dof );
        }

        // add definition
        hdBcs_.Push_back( actBc );
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create homogeneous boundary conditions on '"
                           << name << "'" );
      }
    }

    //=====================================================================
    // inhomogeneous Dirichlet BC
    // =====================================================================

    // fetch paramnodes for idbc
    StdVector<ParamNode*> idbcNodes = bcsNode->GetList("dirichletInhom");

    // iterate over all parameter nodes
    for( UInt i = 0; i < idbcNodes.GetSize(); i++ ) {

      try {
        // read parameters
        dof.clear();
        idbcNodes[i]->Get( "name", name );
        idbcNodes[i]->Get( "quantity", resultName );
        idbcNodes[i]->Get( "dof", dof, false );
        idbcNodes[i]->Get( "value", value );
        idbcNodes[i]->Get( "phase", phase );
        idbcNodes[i]->Get( "entityType", entType );

        // fetch related resultInfo object
        actResultInfo = GetResultInfo( SolutionTypeEnum.Parse(resultName) );

        // Create inhomogeneous boundary condition
        shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
        if( entType == "nodeList" ) {
          defineType = EntityList::NAMED_NODES;
        } else {
          defineType = EntityList::REGION;
        }

        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
                                  name, defineType );
        actBc->entities = actList;
        actBc->result = actResultInfo;
        actBc->eqnMap = eqnMap_;
        if( dof.empty() ) {
          actBc->dof = 1;
        } else {
          actBc->dof = actResultInfo->GetDofIndex( dof );
        }
        actBc->value = value;
        actBc->phase = phase;

        // add definition
        idBcs_.Push_back( actBc );
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create inhomogeneous boundary conditions on '"
                           << name << "'" );
      }
     }


    // =====================================================================
    // inhomogeneous Neumann BC
    // =====================================================================

    // fetch paramnodes for inbc
    StdVector<ParamNode*> inbcNodes = bcsNode->GetList("neumannInhom");

    // iterate over all parameter nodes
    for( UInt i = 0; i < inbcNodes.GetSize(); i++ ) {
      try {
        dof.clear();
        inbcNodes[i]->Get( "name", name );
        inbcNodes[i]->Get( "quantity", resultName );
        inbcNodes[i]->Get( "dof", dof, false );
        inbcNodes[i]->Get( "value", value );
        inbcNodes[i]->Get( "phase", phase );
        inbcNodes[i]->Get( "entityType", entType );

        // fetch related resultInfo object
        actResultInfo = GetResultInfo( SolutionTypeEnum.Parse(resultName) );

        // Create inhomogeneous Neumann boundary condition
        shared_ptr<InhomNeumannBc> actBc ( new InhomNeumannBc );
        if( entType == "nodeList" ) {
          defineType = EntityList::NAMED_NODES;
        } else {
          defineType = EntityList::REGION;
        }

        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
                                  name, defineType );

        actBc->entities = actList;
        actBc->result = actResultInfo;
        actBc->eqnMap = eqnMap_;
        if( dof.empty() ) {
          actBc->dof = 1;
        } else {
          actBc->dof = actResultInfo->GetDofIndex( dof );
        }
        actBc->value = value;
        actBc->phase = phase;

        // add definition
        inBcs_.Push_back( actBc );
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create inhomogeneous Neumann conditions on '"
                           << name << "'" );
      }
    }

    // =====================================================================
    // Constraint Conditions
    // =====================================================================

    // fetch paramnodes for constraint
    StdVector<ParamNode*> csNodes = bcsNode->GetList("constraint");
    std::string masterDof, slaveDof;

    // iterate over all parameter nodes
    for( UInt i = 0; i < csNodes.GetSize(); i++ ) {
      try {
        csNodes[i]->Get( "name", name );
        csNodes[i]->Get( "quantity", resultName );
        csNodes[i]->Get( "entityType", entType );
        csNodes[i]->Get( "masterDof", masterDof );
        csNodes[i]->Get( "slaveDof", slaveDof );

        // fetch related resultInfo object
        actResultInfo = GetResultInfo( SolutionTypeEnum.Parse(resultName) );

        // Create constraint condition
        shared_ptr<Constraint> actBc ( new Constraint );
        if( entType == "nodeList" ) {
          defineType = EntityList::NAMED_NODES;
        } else {
          defineType = EntityList::REGION;
        }

        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
                                  name, defineType );

        actBc->masterEntities = actList;
        actBc->slaveEntities = actList;
        if( masterDof.empty() ) {
          actBc->masterDof = 1;
        } else {
          actBc->masterDof = actResultInfo->GetDofIndex( masterDof );
        }
        if( slaveDof.empty() ) {
          actBc->slaveDof = 1;
        } else {
          actBc->slaveDof = actResultInfo->GetDofIndex( masterDof );
        }

        actBc->result = actResultInfo;
        actBc->eqnMap = eqnMap_;

        // add definition
        constraints_.Push_back( actBc );
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create constraints on '"
                           << name << "'" );
      }
    }

    // =====================================================================
    // Periodic boundary conditions
    // =====================================================================

    // fetch paramnodes for constraint
    StdVector<ParamNode*> prNodes = bcsNode->GetList("periodic");
    std::string masterName, slaveName;

    // iterate over all parameter nodes
    for( UInt i = 0; i < prNodes.GetSize(); i++ ) {
      try {
        prNodes[i]->Get( "master", masterName );
        prNodes[i]->Get( "slave", slaveName );
        prNodes[i]->Get( "dof", dof );
        prNodes[i]->Get( "quantity", resultName );

        // fetch related resultInfo object
        actResultInfo = GetResultInfo( SolutionTypeEnum.Parse(resultName) );

        // get entitylists
        NodeList masterList( ptgrid_ ), slaveList( ptgrid_ );
        masterList.SetNamedNodes( masterName );
        slaveList.SetNamedNodes( slaveName );

        // ensure, that both lists have the same length
        if( masterList.GetSize() != slaveList.GetSize() ) {
          EXCEPTION( "Node lists '" << masterName << "' and '"
                     << slaveName << "' have different size" );
        }

        // iterate over all master nodes and try to find "nearest"
        // node in slave list
        Vector<Double> mLoc, sLoc, diff;
        Double minDist, dist;
        StdVector<UInt> nodes(2);
        EntityIterator masterIt = masterList.GetIterator();
        for( masterIt.Begin(); !masterIt.IsEnd(); masterIt++ ) {

          minDist = 1e42;

          // obtain nodal coordinate
          ptgrid_->GetNodeCoordinate( mLoc, masterIt.GetNode() );
          nodes.Init();
          nodes[0] = masterIt.GetNode();

          // iterate over all slave nodes and find the one with minimum
          // distance
          EntityIterator slaveIt = slaveList.GetIterator();
          for( slaveIt.Begin(); !slaveIt.IsEnd(); slaveIt++ ) {
            ptgrid_->GetNodeCoordinate( sLoc, slaveIt.GetNode() );
            diff = mLoc - sLoc;
            dist = diff.NormL2();
            if( dist < minDist) {
              minDist = dist;
              nodes[1] = slaveIt.GetNode();
            }
          }
          shared_ptr<NodeList> nodePair(new NodeList( ptgrid_ ) );
          nodePair->SetNodes( nodes );

          // create new constraint condition
          shared_ptr<Constraint> actBc ( new Constraint );
          actBc->masterEntities = nodePair;
          actBc->slaveEntities = nodePair;
          if( dof.empty() ) {
            actBc->masterDof = 1;
          } else {
            actBc->masterDof = actResultInfo->GetDofIndex( dof );
          }
          actBc->slaveDof = actBc->masterDof;
          actBc->result = actResultInfo;
          actBc->eqnMap = eqnMap_;
          actBc->periodic = true;

          // add definition
          constraints_.Push_back( actBc );
        }
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create periodic boundary on '"
                           << name << "'" );
      }
    }

    // =====================================================================
    // Load definitions
    // =====================================================================

    // fetch paramnodes for loads
    ReadLoads(bcsNode->GetList("load"), loads_);
    assemble_->AddLoads( loads_ );
  }


  void SinglePDE::ReadLoads(StdVector<ParamNode*> loadNodes, LoadList& out_list)
  {

    std::string name, resultName, dof, entType, value, phase, weight;
    EntityList::DefineType defineType;
    shared_ptr<ResultInfo> actResultInfo;


    // iterate over all parameter nodes
    for( UInt i = 0; i < loadNodes.GetSize(); i++ ) {
      try {
        dof.clear();
        loadNodes[i]->Get( "name", name );
        loadNodes[i]->Get( "quantity", resultName );
        loadNodes[i]->Get( "dof", dof, false );
        loadNodes[i]->Get( "value", value );
        loadNodes[i]->Get( "phase", phase );
        loadNodes[i]->Get( "weight", weight, false ); // only mechanical for optimization
        loadNodes[i]->Get( "entityType", entType );


        // fetch related resultInfo object
        actResultInfo = GetResultInfo( SolutionTypeEnum.Parse(resultName) );

        // Create load condition
        shared_ptr<LoadBc> actLoad( new LoadBc );
        if( entType == "nodeList" ) {
          defineType = EntityList::NAMED_NODES;
        } else {
          defineType = EntityList::REGION;
        }

        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
                                  name, defineType );

        actLoad->entities = actList;
        actLoad->result = actResultInfo;
        actLoad->eqnMap = eqnMap_;
        if ( dof.empty() ) {
          actLoad->dof = 1;
        } else {
          actLoad->dof = actResultInfo->GetDofIndex(dof);
        }
        actLoad->value = value;
        actLoad->phase = phase;
        actLoad->weight = weight;
        out_list.Push_back( actLoad);
      }
      catch (Exception & ex )
        RETHROW_EXCEPTION( ex, "Can not create load condition on '"
                           << name << "'" );
    }
  }

  void SinglePDE::ReadRegionLoads(){
    ReadRegionLoadsFromXML(myParam_->Get("bcsAndLoads", false), regionLoads_);
  }
  
  void SinglePDE::ReadRegionLoadsFromXML(ParamNode* bcNode, std::map<RegionIdType, RegionLoad>& regloads) {

    StdVector<std::string> names, dofs, refCoord, type, phase;
    StdVector<std::string> tempNames, tempDofs,  tempPhase;
    StdVector<std::string>  tempRefCoord, tempType;
    StdVector<RegionIdType> regionIds;
    StdVector<UInt> vecComp;
    StdVector<std::string> loadVec, tempLoadVec, tempLoad(dim_);
    UInt locDof = 0;
    Integer index = -1;

    // Check, if function was called by a scripting command
#ifdef USE_SCRIPTING
    if ( messenger->IsEvaluating() == true ) {

      // obtain parameters from messenger object
      // Note: If scripting is used, only one region load
      // can be specified per call
      SCRIPT_GET( std::string,name);
      SCRIPT_GET( std::string, value );
      SCRIPT_GET( std::string, dof );
      SCRIPT_GET( std::string, refCoordSys );
      SCRIPT_GET( std::string, type );

      // Copy single entries into vectors
      tempNames.Push_back( name );
      tempLoadVec.Push_back( value );
      tempDofs.Push_back( dof );
      tempRefCoord.Push_back( refCoordSys );
      tempType.Push_back( type );

    } else {
#endif
      // obtain parameters from ParamHandler
      // Note: Here all region loads are read (in contrast
      // when called by an external script)

      // try to get bcsAndLoads node
      if( !bcNode )
        return;
      StdVector<ParamNode*> loadNodes = bcNode->GetList("regionLoad");


      for( UInt i = 0; i < loadNodes.GetSize(); i++ ) {

        ParamNode * actNode = loadNodes[i];

        tempNames.Push_back( actNode->Get("name")->AsString() );
        tempLoadVec.Push_back( actNode->Get("value")->AsString() );
        tempPhase.Push_back( actNode->Get("phase")->AsString() );
        tempDofs.Push_back( actNode->Get("dof")->AsString() );
        tempRefCoord.Push_back( actNode->Get("coordSysId")->AsString() );
        tempType.Push_back( actNode->Get("type")->AsString() );
      }

#ifdef USE_SCRIPTING
    }
#endif

    // --- Common part for scripting and parameter file ---


    // Now sort the names and remove double entries
    for (UInt i = 0; i < tempNames.GetSize(); i++) {
      index = names.Find(tempNames[i]);
      if ( index == -1) {
        names.Push_back(tempNames[i]);
      }
    }

    // Convert region names to ID - vector
    ptgrid_->GetRegion().Parse(names, regionIds );


    // loop over all regions
    for (UInt i = 0; i < names.GetSize(); i++) {

      // get for each name all related entries for value,
      // dof, refCoordSys and type
      loadVec.Clear();
      phase.Clear();
      dofs.Clear();
      refCoord.Clear();
      type.Clear();
      for (UInt iEntry = 0; iEntry < tempNames.GetSize(); iEntry++ ) {
        if ( names[i] == tempNames[iEntry] ) {
          loadVec.Push_back(tempLoadVec[iEntry]);
          phase.Push_back(tempPhase[iEntry]);
          dofs.Push_back(tempDofs[iEntry]);
          refCoord.Push_back(tempRefCoord[iEntry]);
          type.Push_back(tempType[iEntry]);
        }
      }

      // check if all entries for  refCoord and type
      // are the same
      for (UInt k=0; k<refCoord.GetSize(); k++) {
        if ( refCoord[k] != refCoord[0] ||
             type[k] != type[0] ) {
          EXCEPTION( "MechPDE::DefineRegionLoads: The region load on region "
                     << names[i] << " has not for all dofs the same entry for "
                     << "refCoord or type (total/unit)!" );
        }
      }

      // Check if an entry already exists for this region
      RegionLoad * curLoad;

      std::map<RegionIdType, RegionLoad>::iterator it;
      it = regloads.find( regionIds[i] );

      if ( it == regloads.end() ) {
        regloads.insert( std::map<RegionIdType, RegionLoad>::value_type( regionIds[i],
                                                                             RegionLoad( dim_, isaxi_ ) ) );
      }
      it = regloads.find( regionIds[i] );
      curLoad = & (*it).second;

      // -- Fill in the data we have so far --
      curLoad->name = ptgrid_->GetRegion().ToString( regionIds[i] );
      curLoad->phase = phase[0];

      if ( curLoad->refCoord != std::string()
           && curLoad->refCoord != refCoord[0] ) {
        EXCEPTION( "Inconsistent definition of time data for regionLoads" );
      } else {
        curLoad->refCoord = refCoord[0];
      }

      if ( curLoad->volume < EPS ) {
        curLoad->volume = ptgrid_->CalcVolumeOfRegion(regionIds[i], isaxi_);
      }

      // now create local load vector
      for (UInt iDim=0; iDim < loadVec.GetSize(); iDim++) {
        locDof = domain->GetCoordSystem(refCoord[iDim])->
          GetVecComponent(dofs[iDim]);
        curLoad->value[locDof-1] = loadVec[iDim];
        curLoad->type = type[iDim];
      }
    }
  }


  void SinglePDE::ReadMaterialData() {
    UInt i, numRegions;
    
    // get list of parameter nodes for region definitions
    StdVector<ParamNode*> regionNodes;

    ParamNode * regionListNode = param->
      Get("domain")->Get("regionList", false );
    if( regionListNode)
      regionNodes = regionListNode->GetList("region");

    numRegions = regionNodes.GetSize();

    // obtain pointer to materialHandler
    MaterialHandler * matLoader = NULL;
    matLoader = domain->GetMaterialHandler();


    // -------------------
    // NORMAL MATERIALS
    // -------------------
    std::string region, material, composite, refCoordSys;

    // iterate over all regions
    for( i = 0; i < numRegions; ++i ) {

      try{
        // get data from node
        regionNodes[i]->Get( "name", region );
        regionNodes[i]->Get( "material", material );
        regionNodes[i]->Get( "coordSysId", refCoordSys );

        // get regionId
        RegionIdType actRegionId = ptgrid_->GetRegion().Parse( region );

        // if no material is set, continue with next loop run
        if( material.empty() )
          continue;

        // if region is not contained for current pde, simply continue
        // with next loop
        if( subdoms_.Find( actRegionId) < 0 )
          continue;

        // print logging information
        Info->PrintF( pdename_, "Material '%s' for region '%s' (ID = %d) "
                      "follows\n", material.c_str(), region.c_str(),
                      actRegionId );
        // Read data
        materials_[actRegionId] = matLoader->
          LoadMaterial( material, pdematerialclass_ );

        // Check for local coordinate system
        if( !refCoordSys.empty() ) {
          CoordSystem * actCoosy =
            domain->GetCoordSystem( refCoordSys);
          materials_[actRegionId]->SetCoordSys( actCoosy );
        }

        // Check for material rotation parameters
        ParamNode * rotNode = regionNodes[i]->Get("matRotation", false);

        Vector<Double> rotVec (3);
        rotVec.Init();

        // NOTE: If no rotation is specified and the dimension is
        // 2D, -> material is rotated by
        // alpha = -90 and gamma = -90 degree,
        // so that we pick by default the yz-plane
        if( !rotNode ) {
          if( dim_ == 2) {
            rotVec[0] = -90.0;
            rotVec[2] = -90.0;
            materials_[actRegionId]->
              RotateAllTensorsByRotationAngles( rotVec, true );
          }
          continue;
        } else {
          rotNode->Get( "alpha", rotVec[0] );
          rotNode->Get( "beta", rotVec[1] );
          rotNode->Get( "gamma", rotVec[2] );

          materials_[actRegionId]->
            RotateAllTensorsByRotationAngles( rotVec, true );
        }

      } catch (Exception& ex ) {
        std::string matClassString;
        Enum2String( pdematerialclass_, matClassString );
        RETHROW_EXCEPTION(ex, "Could not assign material '"
                          << material << "' of materialClass '"
                          << matClassString << "' to region '"
                          << region << "' within pde '"
                          << pdename_ << "'");
      }
    }



    // -------------------
    // COMPOSITE MATERIALS
    // -------------------

    // iterate over all regions
    for( i = 0; i < numRegions; ++i ) {

      try{
        // get data from node
        regionNodes[i]->Get( "name", region );
        regionNodes[i]->Get( "composite", composite );

        // get regionId
        RegionIdType actRegionId = ptgrid_->GetRegion().Parse( region );

        // if no composite is set, continue with next loop run
        if( composite == "" )
          continue;

        // print logging information
        std::ostringstream out;
        out << "Composite material '" << composite << "' for "
            << "region '" << region << "' (ID = " << actRegionId
            << ") follows:\n";
        Info->PrintF( pdename_, out.str().c_str());
        out.str("");

        // get composite node
        ParamNode * compNode = NULL;
        try {
          compNode = param->Get("domain")
            ->Get("composite", "name", composite );
        } catch( Exception& ex ) {
          RETHROW_EXCEPTION( ex, "No composite material defined with name '"
                             << composite << "'" );
        }

        // get laminaNodes
        StdVector<ParamNode*> laminaNodes = compNode->GetList("lamina");

        // Create new lamina and fill ine materials and thicknesses
        Composite & myMat = compositeMaterials_[actRegionId];
        myMat.name = composite;

        // iterate over all single laminas
        for( UInt j = 0; j < laminaNodes.GetSize(); j++ ) {

          // fetch data for lamina
          std::string lamMaterial;
          Double lamThickness, lamOrientation;

          laminaNodes[j]->Get( "material", lamMaterial);
          laminaNodes[j]->Get( "thickness", lamThickness);
          laminaNodes[j]->Get( "orientation", lamOrientation);

          // Print information
          out << " --- Lamina " << j+1 << ": thickness = "
              << lamThickness
              << " m, orientation = " << lamOrientation << "---\n";
          Info->PrintF( pdename_, out.str().c_str());
          out.str("");

          myMat.thickness.Push_back( lamThickness );
          myMat.orientation.Push_back( lamOrientation );
          myMat.materials.Push_back( matLoader->
                                     LoadMaterial( lamMaterial,
                                                   pdematerialclass_ ) );

        } // over single laminae
      } catch (Exception& ex ) {
        RETHROW_EXCEPTION( ex, "Could not create composite material '"
                           << composite << "'");
      }
    } // over composite

    // once again: loop over all regions and make sure that there is either
    // a normal material or a composite material
    for( i = 0; i < numRegions; ++i ) {
      regionNodes[i]->Get( "name", region );
      RegionIdType actRegionId = ptgrid_->GetRegion().Parse( region );
      if (subdoms_.Find(actRegionId) < 0)
        continue;
      if ((materials_.find(actRegionId) == materials_.end())
          && (compositeMaterials_.find(actRegionId) 
              == compositeMaterials_.end())) {
        EXCEPTION("Region '" << region << "' has no material assigned.");
      }
    }
  }



  // ======================================================
  // GET /SET  METHODS
  // ======================================================

  //! Activate the direct coupling
  void SinglePDE::SetDirectCoupling () {


    isDirectCoupled_ = true;
  }



  // ======================================================
  // ALGSYS SECTION (SOLVER, ...)
  // ======================================================
  void SinglePDE::DefineAlgSys() 
  {
    // First check if the PDE needs an algebraic system at all
    if( needsAlgsys_ == false ) {
      return;
    }

    // If PDE is not direct coupled then the PDE has to register
    // at the algebraic system and obtain an Id.
    // Afterwards the matrix-graph has to be set up
    if ( isDirectCoupled_ == false ) {

      // forward the complete xml description of the linear system
      // might be NULL!
      // algsys_->SetXML(FindLinearSystem(pdename_));

      // Set linear system parameters for OLAS
      ReadOlasParams( pdename_ );
      olasInfo_ = info->Get("OLAS")->Get(pdename_);
      
      
      // Initialize the matrix graph object
      algsys_->GraphSetupInit(1);

      algsys_->RegisterPDE( pdeId_, eqnMap_->GetNumEqns(),
                            eqnMap_->GetNumLastFreeDof() );

      //assemble_->SetPDEId( pdeId_ );
      solveStep_->SetPDEId( pdeId_ );


      // trigger the creation and assembly of the matrix graph
      algsys_->AssembleInit( pdeId_, pdeId_, false );
      assemble_->SetupMatrixGraph(pdeId_, pdeId_);
      algsys_->AssembleDone( pdeId_, pdeId_, false );

      // finish the assembly of the matrix graph
      algsys_->GraphSetupDone();

      // obtain reordering of the matrix graph and pass it to the EQN-object.
      StdVector<UInt> newOrder;
      algsys_->GetReordering( pdeId_, newOrder );
      eqnMap_->ReorderMapping( newOrder );
    }

    // pass information about dofs, number of dirichlet equations
    // and constraints to the algebraic system
    //algsys_->SetBlockSize( pdeId_, 1 );

    UInt numDir = eqnMap_->GetNumInHomDirichletEqns() +
      numCouplingBcs_;
    algsys_->SetNumDirichletBCs(pdeId_, numDir );

    // create matrices and solver object, if PDE is not direct coupled
    if ( isDirectCoupled_ == false ) {
      CreateMatrices_Solver();

      // =====================================================================
      // Set the initial conditions
      // =====================================================================
      if ( analysistype_ == TRANSIENT){
      	SetInitialCondition();
      }

      // Incorporate values of memento here, if the values are used
      // as "start"-values from a previous simulation run
      // (e.g. obtained from a previous static run)
      if( memento_ != NULL
          && mementoAsDirichlet_ ==  false) {
        IncorporateMemento();
      }
    }

  }



  // ======================================================
  // Adaptivity SECTION
  // ======================================================

#ifdef ADAPTGRID
  void SinglePDE::RefineMesh( const Integer level) {


    Integer          numChilds;
    Integer          numElems;
    StdVector<Elem*> elemssd;
    Double           theta_s;
    Double           coarseFactor;
    Double           tol4Elm;
    Integer          numRefLoops=0;
    Integer          numRefinements;
    Integer          iem;

    // get max num elements for the domain,on which we have the equation
    numElems=ptgrid_->GetMaxnumElem(level,subdoms_);

    // get pointer to array with elements
    ptgrid_->GetElemSD(elemssd,subdoms_[0],level);

    // if element is marked, then value of the array's element is equal 1,
    // else 0.
    markingElems_.Resize(numElems);
    marlingElems_.Init();

    if (!conf->ifget("safety_factor_for_space_adaptivity",theta_s)) {
      theta_s=1.0;
    }

    coarseFactor = 0.7;
    tol4Elm = theta_s*tolSpaceErr_;

    std::cout << " tolerance space error: " << tolSpaceErr_ << std::endl;

    // loop over elements
    for (iem=0; iem<numElems; iem++) {
      elemssd[iem]->refinementFlag = 0;
      if (errorMap_[iem]>tol4Elm) {
        elemssd[iem]->refinementFlag = 1;

        numChilds=elemssd[iem]->ptElem->GetNumChilds();

        numRefinements = defineRefinements(errorMap_[iem],tol4Elm,numChilds);

        if (numRefinements>numRefLoops)
          numRefLoops = numRefinements;

        markingElems_[iem]=numRefinements;
        elemssd[iem]->refinementNumber = numRefinements;
      }
      else {
        if (errorMap_[iem] < coarseFactor*tol4Elm) {
          elemssd[iem]->refinementFlag = 0;
        }
        else {
          elemssd[iem]->refinementFlag = 0;
        }
      }
    }

    // Fuehren die Verfeinerung durch
    std::cout << " number of refinement loops \n" << numRefLoops << std::endl;

    ptgrid_->Refine(numRefLoops);

  }


  bool SinglePDE::TestError(const Integer level) {


    if (!ptError_) ConstructorError();

    // Berechnung der Fehlerkarte
    Double            totalErr;
    Vector<Double>    solVec;

    solVec.Resize(sol_.size());
    int i;
    for (i=0; i<sol_.size(); i++)
      solVec[i]=sol_[0][i];
    //      sol_.toVector(solVec,1);

    ptError_->CalcErrorMap(solVec,subdoms_,ptgrid_,errorMap_,totalErr,level);

    std::cout << " space error: " << totalErr <<
      " tolerance: " << tolSpaceErr_ << std::endl;

    if (totalErr > tolSpaceErr_) return true;
    else return false;
  }


  //In this fnc we delete old pointer to Error-object & create new
  void SinglePDE::ConstructorError() {


    if (ptError_) delete ptError_;

    ptError_=new SpaceErrorEstimator();
    ptError_->Init(this);
  }


  void SinglePDE::WriteErrorInfo(WriteResults * ptmeshes) {
    ptmeshes->WriteElemSolution(errorMap_,0,0,"ERR-errorMap");
    ptmeshes->WriteElemSolution(markingElems_,0,0,"ERR-markedElems");
  }

#endif // AdaptGrid

  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void SinglePDE::CalcInputCoupling() {


    std::string errMsg;
    StdVector<UInt> * nodes;
    SingleVector * val;
    Integer eqnNr;
    UInt couplingDof;

    // Determine maximal allowed equation number for algebraic system
    Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // Reset counter for boundary conditions
    couplingBCsCounter_ = 0;


    // Outer loop over all INPUT coupling terms
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {

      //    ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);
      couplingDof = ptCoupling_->GetInputDof(i);

      // Up to now, Coupling is only possible with
      // Real valued solutions
      Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);


      switch(ptCoupling_->GetInputType(i)) {

        // -------------------
        // COORDINATE COUPLING
        // -------------------
      case COORD:

        ptCoupling_->GetInputNodes(i, nodes);
        ptgrid_->SetNodeOffset( *nodes , help );

        break;

        // -------------------
        // RHS COUPLING
        // -------------------
      case RHS:
        ptCoupling_->GetInputNodes(i, nodes);

        for ( UInt dof = 0; dof < couplingDof; dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
            eqnNr = eqnMap_->GetNodeEqn( (*nodes)[j], dof+1 );
            if ( eqnNr != 0 && eqnNr <= maxAllowedEqn ) {
              algsys_->SetNodeRHS( help[ dof + couplingDof * j ], pdeId_,
                                   eqnNr );
            }

#ifndef NDEBUG
            /*
            else if ( eqnNr > maxAllowedEqn ) {
              (*debug) << "SinglePDE::CalcInputCoupling: "
                       << "(" << pdename_ << ") "
                       << "Refused to pass "
                       << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
                       << "it execeeds numLastFreeDof = " << maxAllowedEqn
                       << std::endl;
            }
            else if ( eqnNr == 0 ) {
              (*debug) << "SinglePDE::CalcInputCoupling: "
                       << "(" << pdename_ << ") "
                       << "Refused to pass node " << (*nodes)[j] << "to SetNodeRHS(), since "
                       << "it is fixed by hom. Dirichlet BC" << std::endl;
            }
            */
#endif

          }
        }
        break;

        // -----------------------
        // InhomDirichlet COUPLING
        // -----------------------
      case ID_BC:

        // How do we want to treat inhomogeneous Dirichlet boundary
        // conditions?
        {
          if ( usePenalty_ == false ) {
            EXCEPTION( "Cannot use inhom. Dirichlet coupling together with "
                       << "IDBC elimination, since the equation numbering "
                       << "object does currently not have the information "
                       << "required to put those values at the end of the "
                       << "equation number interval! Please set idbcHandling="
                       << '"' << "penalty" << '"' << " in your xml-file" );
          }
        }

        // Set flag that the boundary conditions have to be incorporated
        updateCouplingBCs_ = true;

        ptCoupling_->GetInputNodes(i, nodes);

        for ( UInt dof = 0; dof < ptCoupling_->GetInputDof(i); dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize();
                j++, couplingBCsCounter_++) {

            eqnNr = eqnMap_->GetNodeEqn( (*nodes)[j], dof+1 );

            if (eqnNr==0) {
              EXCEPTION( "The specified coupling node has no equation number" );
        //        std::cerr << "node=" << *nodes)[j]
            }

            algsys_->SetDirichlet( pdeId_, eqnNr,
                                   help[dof+j*couplingDof] );
          }
        }
        break;

      case MAT:
        break;

      default:
        EXCEPTION( "Not implemented yet" );
        break;
      }
    }
  }

  void SinglePDE::ResetCoupling() {

    iterCoupledCounter_ = 0;


  }

  void SinglePDE::SetIDBC( const std::string &name,
                           const std::string &dofString,
                           const std::string &value,
                           const std::string &phase ) {

    // Try ro find existing entry in IDBC vector
    Integer index = -1;

    for( UInt i = 0; i < idBcs_.GetSize(); i++ ) {
      if( idBcs_[i]->entities->GetName() == name ) {
        index = i;
        break;
      }
    }

    if( index == -1 ) {
      EXCEPTION( "Entities '" << name << "' not found for "
                 << "Dirichlet values!" );
    } else {
      idBcs_[index]->value = value;
      idBcs_[index]->phase = phase;
    }
  }


 void SinglePDE::WriteRestart( ) {

   if (pdename_=="mechanic" || pdename_=="acoustic"
	   || pdename_=="fluidMech" || pdename_=="smooth") {
     std::string simName = progOpts->GetSimName();

     // create name of output file
     std::string restartFileName = simName+"_"+pdename_+".restart";

     std::ofstream writeTo(restartFileName.c_str(), std::ios::binary);
     
     if( !writeTo.good() ) {
       EXCEPTION( "Could not write to restart file '" << restartFileName
                  << "'!" );
     }

#if 0
     // START OF HDF5 SERIALIZATION DEMO!
     H5::H5File* h5file = new H5::H5File( std::string("SinglePDE.h5").c_str(), H5F_ACC_TRUNC );

     boost::archive::hdf5_oarchive ho(h5file);
     // test compression
     ho.set_compression(6);
     //     ho << BOOST_SERIALIZATION_NVP(temp);
     ho << BOOST_SERIALIZATION_NVP(simName);
     std::vector<UInt> test(3);
     test[0] = 5;
     test[1] = 17;
     test[2] = 894;
     // ho << BOOST_SERIALIZATION_NVP(test);
     delete h5file;
     // END OF HDF5 SERIALIZATION DEMO!
#endif

     boost::archive::binary_oarchive outArchive(writeTo);
     GetMemento(memento_);
     PDEMemento const & temp = *memento_;
     outArchive << temp;
     writeTo.close();
   } else {
     EXCEPTION( "Restart functionality only for "
                << "mechanic and acoustic tested" );
   }


 }

  void SinglePDE::ReadRestart(UInt &startStep ) {

    if (pdename_=="mechanic" || pdename_=="acoustic"
 	   || pdename_=="fluidMech" || pdename_=="smooth") {
      std::string simName = progOpts->GetSimName();
      std::string restartFileName = simName+"_"+pdename_+".restart";

      std::ifstream readRestart(restartFileName.c_str(), std::ios_base::in );
      if( !readRestart.good() ) {
        EXCEPTION( "Could not open restart file '" << restartFileName << "'!" );
      }
      boost::archive::binary_iarchive inArchive( readRestart );

      shared_ptr<PDEMemento> myMemento (new PDEMemento );
      inArchive >> *myMemento;
      readRestart.close();
      SetMemento( myMemento, false );

      startStep = myMemento->GetRestartStep();
    } else  {
      EXCEPTION( "Restart functionality only for "
                 << "mechanic and acoustic tested" );
    }

  }


  void SinglePDE::GetMemento( shared_ptr<PDEMemento>& memento) {

    // create new memento
    shared_ptr<PDEMemento> myMemento (new PDEMemento() );

    // first get memento of coupling object
    if (isIterCoupled_) {
      ptCoupling_->GetMemento(myMemento->couplingMemento_);
      myMemento->isIterCoupled_ = true;
    }

    // then write own data to PDEMemento
    myMemento->analysisType_ = analysistype_;
    myMemento->gridFileName_ = progOpts->GetMeshFileStr();
    myMemento->stepNum_ = domain->GetSingleDriver()->GetActStep(pdename_);

    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {

      // --- Real values --
      Vector<Double> & solReal =
        dynamic_cast<NodeStoreSol<Double>&>(*(sol_)).GetAlgSysVector();
      UInt numDofs = results_[0]->dofNames.GetSize();
      StdVector<Integer> eqns;
      for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ){

        // create new entity list
        shared_ptr<NodeList> actSDList( new NodeList(ptgrid_ ) );
        actSDList->SetNodesOfRegion( subdoms_[iRegion] );

        // create new vector
        Vector<Double> * values = new Vector<Double>;
        Vector<Double> deriv1, deriv2, tn_1;

        values->Resize( actSDList->GetSize() * numDofs );
        values->Init();
        if (analysistype_ == TRANSIENT ) {
          deriv1.Resize( actSDList->GetSize() * numDofs );
          deriv2.Resize( actSDList->GetSize() * numDofs );
            tn_1.Resize( actSDList->GetSize() * numDofs );
        }

        EntityIterator it = actSDList->GetIterator();
        UInt pos = 0;
        for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
          eqnMap_->GetEqns( eqns, *results_[0], it );
          for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
            if ( eqns[iDof] != 0 ) {
              (*values)[numDofs*pos+iDof] = solReal[abs(eqns[iDof])-1];
              if( analysistype_ == TRANSIENT ) {
                deriv1[numDofs*pos+iDof] =
                  TS_alg_->GetDeriv1()[(abs(eqns[iDof])-1)];
                deriv2[numDofs*pos+iDof] =
                    TS_alg_->GetDeriv2()[(abs(eqns[iDof]-1))];
                  tn_1[numDofs*pos+iDof] =
                    TS_alg_->GetOld1()[(abs(eqns[iDof]-1))];
              }
            } else {
              (*values)[numDofs*pos+iDof] = 0.0;
              if ( analysistype_ == TRANSIENT ) {
                deriv1[numDofs*pos+iDof] = 0.0;
                deriv2[numDofs*pos+iDof] = 0.0;
                  tn_1[numDofs*pos+iDof] = 0.0;
              }
            }
          }
        }

        // pass vector to memento object
        std::string regionName = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
        myMemento->solution_[regionName] = values;
        if ( analysistype_ == TRANSIENT ) {
          myMemento->solDeriv1_[regionName] = deriv1;
          myMemento->solDeriv2_[regionName] = deriv2;
          myMemento->sol_tn_1_[regionName] = tn_1;
        }
      }
    } else {

      // --- Complex values --

      // store current frequency
      myMemento->freq_ = dynamic_cast<HarmonicDriver*>(domain->GetSingleDriver())
        ->GetActFreq();

      Vector<Complex> & solComp =
        dynamic_cast<NodeStoreSol<Complex>&>(*(sol_)).GetAlgSysVector();
      UInt numDofs = results_[0]->dofNames.GetSize();
      StdVector<Integer> eqns;
      for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ){

        // create new entity list
        shared_ptr<NodeList> actSDList( new NodeList(ptgrid_ ) );
        actSDList->SetNodesOfRegion( subdoms_[iRegion] );

        // create new vector
        Vector<Complex> * values = new Vector<Complex>;
        values->Resize( actSDList->GetSize() * numDofs );
        values->Init();
        EntityIterator it = actSDList->GetIterator();
        UInt pos = 0;
        for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
          eqnMap_->GetEqns( eqns, *results_[0], it );
          for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
            if ( eqns[iDof] != 0 ) {
              (*values)[numDofs*pos+iDof] = solComp[abs(eqns[iDof])-1];
            } else {
              (*values)[numDofs*pos+iDof] = 0.0;
            }
          }
        }

        // pass vector to memento object
        std::string regionName = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
        myMemento->solution_[regionName] = values;
      }
    }

    // now memento is initialized
    myMemento->isSet_ = true;

    // pass back memento
    memento = myMemento;
  }


  void SinglePDE::SetMemento( shared_ptr<PDEMemento>& memento,
                              bool mementoAsDirichlet ) {

    if( isInitialized_ == true ) {
      EXCEPTION( "SetMemento may only be called, if the method "
                 << "SinglePDE::Init() was not called yet!" );
    }

    memento_ = memento;
    mementoAsDirichlet_ = mementoAsDirichlet;
  }

  void SinglePDE::IncorporateMemento( ) {

    // if there is no memento present -> leave
    if( !memento_ ) {
      return;
    }

   // if there is no information in the memento just leave
    if ( memento_->isSet_ == false ) {
      return;
    }

    // check that memento is bases on the same grid file
    if( memento_->gridFileName_ != progOpts->GetMeshFileStr() ) {
      EXCEPTION( "Error in reading in memento: Memento is based on grid '"
                 << memento_->gridFileName_
                 << ", whereas the current simulation is based on grid '"
                 << progOpts->GetMeshFileStr() );
    }

    UInt numDofs = results_[0]->dofNames.GetSize();
    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {

      // convert solution to transient StoreSolution type
      Vector<Double> & solVec =
        (dynamic_cast<NodeStoreSol<Double> &>(*sol_)).GetAlgSysVector();


      Vector<Double> solDeriv1, solDeriv2, solTn_1;
      if( TS_alg_->GetDeriv1().GetSize() != 0 ) {
        solDeriv1 = TS_alg_->GetDeriv1();
      } else {
        solDeriv1.Resize( solVec.GetSize() );
        solDeriv1.Init();
      }

      if( TS_alg_->GetDeriv2().GetSize() != 0 ) {
        solDeriv2 = TS_alg_->GetDeriv2();
      } else {
        solDeriv2.Resize( solVec.GetSize() );
        solDeriv2.Init();
      }

      if( TS_alg_->GetOld1().GetSize() != 0 ) {
        solTn_1 = TS_alg_->GetOld1();
      } else {
        solTn_1.Resize( solVec.GetSize() );
        solTn_1.Init();
      }

      
      if ( memento_->analysisType_ == HARMONIC ) {

        // -> perform complex-to-real adjustment

        // frequency of memento
        Double freq = memento_->freq_;

        // Iterate over all regions
        for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ) {

          // Check for related region in memento object
          std::string name = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
          if( memento_->solution_.find( name) != memento_->solution_.end() ) {

            // get grip of vector and derivatives of memento
            Vector<Complex> const & sol =
              dynamic_cast<const Vector<Complex>& >(*(memento_->solution_[name]) );

            // create entitylist
            shared_ptr<NodeList> nodes (new NodeList(ptgrid_));
            nodes->SetNodesOfRegion( subdoms_[iRegion] );

            // iterate over all entries
            EntityIterator it = nodes->GetIterator();
            StdVector<Integer> eqns;
            UInt pos = 0;
            for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
              eqnMap_->GetEqns( eqns, *results_[0], it );
              for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
                UInt actPos = numDofs*pos+iDof;
                if ( eqns[iDof] != 0
                     && (unsigned int) (std::abs(eqns[iDof])-1) < solVec.GetSize() ) {
                  solVec[eqns[iDof]-1] = sol[actPos].real();
                  if ( analysistype_ == TRANSIENT ) {
                    solDeriv1[eqns[iDof]-1] = -2*PI*freq* sol[actPos].imag();
                    solDeriv2[eqns[iDof]-1] =
                      -4*PI*PI*freq*freq*sol[actPos].real();
                  }
                }
              }
            }
          }
        }

        // pass derivatives to timestepping algorithm
        if( analysistype_ == TRANSIENT ) {
          TS_alg_->SetDeriv1( solDeriv1 );
          TS_alg_->SetDeriv2( solDeriv2 );
        }

      } else {

        // check if derivatives are also needed
        bool needDeriv = ( analysistype_ == TRANSIENT ) &&
          (memento_->analysisType_ == TRANSIENT);

        // Iterate over all regions
        for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ) {

          // Check for related region in memento object
          std::string name = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
          if( memento_->solution_.find( name) != memento_->solution_.end() ) {

            // get grip of vector and derivatives of memento
            Vector<Double> const & sol =
              dynamic_cast<const Vector<Double>& >(*(memento_->solution_[name]) );
            Vector<Double> const & deriv1 = memento_->solDeriv1_[name];
            Vector<Double> const & deriv2 = memento_->solDeriv2_[name];
            Vector<Double> const & tn_1 = memento_->sol_tn_1_[name];

            // create entitylist
            shared_ptr<NodeList> nodes (new NodeList(ptgrid_));
            nodes->SetNodesOfRegion( subdoms_[iRegion] );

            // iterate over all entries
            EntityIterator it = nodes->GetIterator();
            StdVector<Integer> eqns;
            UInt pos = 0;
            for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
              eqnMap_->GetEqns( eqns, *results_[0], it );
              for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
                if ( eqns[iDof] != 0
                     && (unsigned int) (std::abs(eqns[iDof])-1) < solVec.GetSize() ) {
                  solVec[eqns[iDof]-1] = sol[numDofs*pos+iDof];
                  if ( needDeriv ) {
                    solDeriv1[eqns[iDof]-1] = deriv1[numDofs*pos+iDof];
                    solDeriv2[eqns[iDof]-1] = deriv2[numDofs*pos+iDof];
                      solTn_1[eqns[iDof]-1] = tn_1[numDofs*pos+iDof];
                  }
                }
              }
            }
          }
        }

        // pass derivatives to timestepping algorithm
        if( needDeriv ) {
          TS_alg_->SetDeriv1( solDeriv1 );
          TS_alg_->SetDeriv2( solDeriv2 );
          TS_alg_->SetOld1( solTn_1 );
        }
      }
    } else if ( analysistype_ == HARMONIC ) {

      // check value-usage type
      if( mementoAsDirichlet_ != true) {
        EXCEPTION( "For an harmonic simulation only the usage "
                   << "of a memento as Drichlet values makes sense!" );
      }

          // iterate over all regions of pde
      for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {

        std::string regionName = ptgrid_->GetRegion().ToString( subdoms_[i] );

        // try to find related region in memento object
        if( memento_->solution_.find( regionName ) !=
            memento_->solution_.end() ) {

          Vector<Complex> const & regionSol =
            dynamic_cast<Vector<Complex>&>(*memento_->solution_[regionName]);

          // create entitylist
          shared_ptr<NodeList> nodes (new NodeList(ptgrid_));
          nodes->SetNodesOfRegion( subdoms_[i] );

          // iterate over all entries
          EntityIterator it = nodes->GetIterator();
          UInt pos = 0;
          for( it.Begin(); !it.IsEnd(); it++, pos++ ) {

            for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
              Complex val = regionSol[pos * numDofs + iDof];

              // create idbc-condition and append to class container
              shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
              shared_ptr<NodeList> actList (new NodeList(ptgrid_) );
              StdVector<UInt> nodeList(1);
              nodeList[0] = it.GetNode();
              actList->SetNodes(nodeList);

              actBc->entities = actList;
              actBc->result = results_[0];
              actBc->eqnMap = eqnMap_;
              actBc->dof = iDof+1;
              actBc->value = lexical_cast<std::string>(std::abs( val ) );
              actBc->phase = lexical_cast<std::string>(std::atan2( val.imag(), val.real())
                                    *180/PI );

              // append idbc at end of list
              idBcs_.Push_back( actBc );
            }
          }
        }
      }
    }


  }

  template<class TYPE>
  void SinglePDE::ExtractResult( shared_ptr<BaseResult> res,
                                 BaseNodeStoreSol * ptStoreSol ) {

    StdVector<Integer> eqnNums;
    ResultInfo & actDof = *(res->GetResultInfo() );
    UInt numDofs = actDof.dofNames.GetSize();

    Vector<TYPE> & solHelp =
      dynamic_cast<NodeStoreSol<TYPE>&> ( *ptStoreSol ).GetAlgSysVector();

    EntityIterator it = res->GetEntityList()->GetIterator();
    Vector<TYPE> & actSol = dynamic_cast<Result<TYPE>&>
      (*(res)).GetVector();
    actSol.Resize( res->GetEntityList()->GetSize() *
                   actDof.dofNames.GetSize() );

    for( it.Begin(); !it.IsEnd(); it++ ) {

      // get equation numbes
      eqnMap_->GetEqns( eqnNums, actDof, it );
      for( UInt iDof = 0; iDof < eqnNums.GetSize(); iDof++ ) {
        if( eqnNums[iDof] != 0 ) {
          actSol[it.GetPos()*numDofs+iDof] = solHelp[abs(eqnNums[iDof])-1];
        } else {
          actSol[it.GetPos()*numDofs+iDof] = 0.0;
        }
      }
    }

  }

  void SinglePDE::ExtractDerivResult( shared_ptr<BaseResult> res, UInt deriv ) {

    StdVector<Integer> eqnNums;


    ResultInfo & actDof = *(res->GetResultInfo() );
    UInt numDofs = actDof.dofNames.GetSize();
    EntityIterator it = res->GetEntityList()->GetIterator();

    if( res->GetEntryType() == BaseMatrix::DOUBLE ) {

      // === TRANSIENT CASE ===
      const Vector<Double>& (TimeStepping::*fct) () const;
      switch ( deriv ) {
      case 1:
        fct = &TimeStepping::GetDeriv1;
        break;
      case 2:
        fct = &TimeStepping::GetDeriv2;
        break;
      default :
        EXCEPTION( "Only derivatives up to order 2 possible" );
      }

      const Vector<Double> & solHelp = (TS_alg_->*fct)();
      Vector<Double> & actSol = (dynamic_cast<Result<Double>&>
                                 (*(res))).GetVector();
      actSol.Resize( res->GetEntityList()->GetSize() *
                     actDof.dofNames.GetSize() );

      // iterate over all elements
      for( it.Begin(); !it.IsEnd(); it++ ) {
        eqnMap_->GetEqns( eqnNums, *results_[0], it );

        // iterate over all dofs
        for( UInt iDof = 0; iDof < eqnNums.GetSize(); iDof++ ) {
          if( eqnNums[iDof] != 0 ) {
            actSol[it.GetPos()*numDofs+iDof] = solHelp[abs(eqnNums[iDof])-1];
          } else {
            actSol[it.GetPos()*numDofs+iDof] = 0.0;
          }
        }
      }
    } else {
      // === HARMONIC CASE ===
      Double omega = solveStep_->GetActFreq() * 2 * PI;

      // determine correct factor
      Complex factor = Complex(0.0, 0.0);
      switch( deriv ) {
      case 1:
        factor = Complex( 0.0, omega );
        break;
      case 2:
        factor = Complex( -omega*omega, 0.0 );
        break;
      default :
        EXCEPTION( "Only derivatives up to order 2 possible" );
      }

      Vector<Complex> & solHelp =
        dynamic_cast<Vector<Complex>& > (*solVec_);
      Vector<Complex> & actSol = dynamic_cast<Result<Complex>&>
         (*(res)).GetVector();
      actSol.Resize( res->GetEntityList()->GetSize() *
                     actDof.dofNames.GetSize() );

      // iterate over all elements
      for( it.Begin(); !it.IsEnd(); it++ ) {
        eqnMap_->GetEqns( eqnNums, *results_[0], it );

        // iterate over all dofs
        for( UInt iDof = 0; iDof < eqnNums.GetSize(); iDof++ ) {
          if( eqnNums[iDof] != 0 ) {
            actSol[it.GetPos()*numDofs+iDof] = factor * solHelp[abs(eqnNums[iDof])-1];
          } else {
            actSol[it.GetPos()*numDofs+iDof] = 0.0;
          }
        }
      }

    }
  }

  // Instantiate functions
  template void SinglePDE::ExtractResult<Double>( shared_ptr<BaseResult> res,
                                                  BaseNodeStoreSol * ptStoreSol );
  template void SinglePDE::ExtractResult<Complex>( shared_ptr<BaseResult> res,
                                                   BaseNodeStoreSol * ptStoreSol );
  template<class TYPE>
  void SinglePDE::ExtractRhsResult( shared_ptr<BaseResult> res,
                                    shared_ptr<ResultInfo> eqnResultInfo ) {

    StdVector<Integer> eqnNums;
    ResultInfo & actDof = *(res->GetResultInfo() );
    UInt numDofs = actDof.dofNames.GetSize();

    Vector<TYPE> & solHelp =
      dynamic_cast<Vector<TYPE>&> (*rhsVec_ );

    EntityIterator it = res->GetEntityList()->GetIterator();
    Vector<TYPE> & actSol = dynamic_cast<Result<TYPE>&>
      (*(res)).GetVector();
    actSol.Resize( res->GetEntityList()->GetSize() *
                   actDof.dofNames.GetSize() );

    for( it.Begin(); !it.IsEnd(); it++ ) {

      // get equation numbes
      eqnMap_->GetEqns( eqnNums, *eqnResultInfo, it );
      for( UInt iDof = 0; iDof < eqnNums.GetSize(); iDof++ ) {
        if( eqnNums[iDof] != 0 ) {
          actSol[it.GetPos()*numDofs+iDof] = solHelp[abs(eqnNums[iDof])-1];
        } else {
          actSol[it.GetPos()*numDofs+iDof] = 0.0;
        }
      }
    }

  }
  // Instantiate functions
  template void
  SinglePDE::ExtractRhsResult<Double>( shared_ptr<BaseResult> res,
                                       shared_ptr<ResultInfo> eqnResultInfo );
  template void
  SinglePDE::ExtractRhsResult<Complex>( shared_ptr<BaseResult> res,
                                        shared_ptr<ResultInfo> eqnResultInfo );

  // ======================================================
  // SCRIPTING SECTION
  // ======================================================

  void SinglePDE::RegisterFunctions() {
    typedef FctPointer<SinglePDE> FCPT;
    StdVector<ArgList> a;
    StdVector<FCPT*> pt;
    StdVector<std::string> name;

    // --- IDBC ---
    a.Push_back();
    a.Last().RegisterParam("name", ArgList::STRING );
    a.Last().RegisterParam("dof", ArgList::STRING );
    a.Last().RegisterParam("value", ArgList::STRING );
    a.Last().RegisterParam("phase", ArgList::STRING );
    pt.Push_back( new FCPT( this, &SinglePDE::Wrap_IDBC) );
    name.Push_back( "setIdbc" );


    // --- GetValue ---
    a.Push_back();
    a.Last().RegisterParam( "name", ArgList::STRING );
    pt.Push_back( new FCPT( this, &SinglePDE::Wrap_GetValue) );
    name.Push_back( "getValue" );

    // Now register all functions with scripting
    for (UInt i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }
  }

  void SinglePDE::Wrap_IDBC() {
    SCRIPT_GET( std::string, name );
    SCRIPT_GET( std::string, dof );
    SCRIPT_GET( std::string, value );
    SCRIPT_GET( std::string, phase );
    SetIDBC(name, dof, value, phase );
  }

  void SinglePDE::Wrap_GetValue() {
    SCRIPT_GET( std::string, name );
    StdVector<UInt> nodeNrs;
    ptgrid_->GetNodesByName( nodeNrs, name );
    for (UInt i=0; i<nodeNrs.GetSize(); i++) {
      Double val;
      sol_->Get(nodeNrs[i]-1,0,val);
      SCRIPT_RETVAL.Push_back(lexical_cast<std::string>(val));
    }
  }


  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void SinglePDE::ReadDataPML(std::string& dampingTypePML,
                              Matrix<Double>& inner,
                              Double& dampPML,
                              ParamNode * actNode ) {


    // Check, if pml node has a child "propRegion"
    ParamNode * propRegionNode = actNode->Get( "propRegion", false );

    // If no propagation region is defined explicitly, we
    // let the method GetPMLLayerData() extract the geometric information
    // for the propagation region
    if( propRegionNode ) {

      //resize data for ptopagation region
      inner.Resize(2,dim_);
      inner.Init();

      //xMin
      propRegionNode->Get( "xMin", inner[0][0] );

      //xMax
      propRegionNode->Get( "xMax", inner[1][0] );

      //yMin
      propRegionNode->Get( "yMin", inner[0][1] );

      //yMax
      propRegionNode->Get( "yMax", inner[1][1] );

      if ( dim_ == 3 ) {
        //zMin
        propRegionNode->Get( "zMin", inner[0][2] );

        //zMax
        propRegionNode->Get( "zMax", inner[1][2] );
      }
    }

    //get type of damping function
    actNode ->Get( "type", dampingTypePML );

    //get factor for damping function
    actNode->Get( "dampFactor", dampPML );

  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void SinglePDE::GetPMLLayerData(Matrix<Double>& inner,
                                  Matrix<Double>& outer,
                                  RegionIdType actRegion )  {


    // outstream for info-File
    std::ostringstream out;
    out.clear();
    out << "PML for region '" << ptgrid_->GetRegion().ToString(actRegion) << "':" << std::endl;

    if ( inner.GetNumCols() != dim_ ) {

      //we have to compute it, since the user has not specified it
      inner.Resize(2,dim_);
      inner.Init();

      for (UInt isd = 0; isd < subdoms_.GetSize(); isd++) {
        if ( subdoms_[isd] != actRegion ) {
          StdVector<Elem*> elemssd;
          ptgrid_->GetElems(elemssd, subdoms_[isd] );

          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            StdVector<UInt> & connecth = elemssd[actEl]->connect;

            Matrix<Double> ptCoord;
            ptgrid_->GetElemNodesCoord(ptCoord, connecth,  false );
            for (UInt i=0; i< ptCoord.GetNumCols(); i++) {
              //minInnerX
              if ( ptCoord[0][i] < inner[0][0] )
                inner[0][0] = ptCoord[0][i];

              //minInnerY
              if ( ptCoord[1][i] < inner[0][1] )
                inner[0][1] = ptCoord[1][i];

              if ( dim_ > 2 ) {
                //minInnerZ
                if ( ptCoord[2][i] < inner[0][2] )
                  inner[0][2] = ptCoord[2][i];
              }

              //maxInnerX
              if ( ptCoord[0][i] > inner[1][0] )
                inner[1][0] = ptCoord[0][i];

              //maxInnerY
              if ( ptCoord[1][i] > inner[1][1] )
                inner[1][1] = ptCoord[1][i];

              if ( dim_ > 2 ) {
                //maxInnerZ
                if ( ptCoord[2][i] > inner[1][2] )
                  inner[1][2] = ptCoord[2][i];
              }
            }
          }
        }
      }
    }

    out << "Acoustic propagation coordinates:\n"
        << "   xmin = " << inner[0][0] << std::endl
        << "   xmax = " << inner[1][0] << std::endl
        << "   ymin = " << inner[0][1] << std::endl
        << "   ymax = " << inner[1][1] << std::endl;
    if ( dim_ == 3) {
      out << "   zmin = " << inner[0][2] << std::endl
          << "   zmax = " << inner[1][2] << std::endl;
    }


    outer.Resize(inner.GetNumRows(),inner.GetNumCols());
    // set outer boundary values to max-value of acoustic propagation region
    outer[0][0] = outer[1][0] = inner[1][0];
    outer[0][1] = outer[1][1] = inner[1][1];
    if (inner.GetNumCols() > 2 ) {
      outer[0][2] = outer[1][2] = inner[1][2];
    }

    StdVector<Elem*> elemssd;
    ptgrid_->GetElems(elemssd, actRegion );

    for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
      StdVector<UInt> & connecth = elemssd[actEl]->connect;

      Matrix<Double> ptCoord;
      ptgrid_->GetElemNodesCoord(ptCoord, connecth,  false );
      for (UInt i=0; i< ptCoord.GetNumCols(); i++) {
        //minXPML
        if ( ptCoord[0][i] < outer[0][0] )
          outer[0][0] = ptCoord[0][i];

        //minYPML
        if ( ptCoord[1][i] < outer[0][1] )
          outer[0][1] = ptCoord[1][i];

        if (inner.GetNumCols() > 2 ) {
          //minZPML
          if ( ptCoord[2][i] < outer[0][2] )
            outer[0][2] = ptCoord[2][i];
        }

        //maxXPML
        if ( ptCoord[0][i] > outer[1][0] )
          outer[1][0] = ptCoord[0][i];

        //maxYPML
        if ( ptCoord[1][i] > outer[1][1] )
          outer[1][1] = ptCoord[1][i];

        if (inner.GetNumCols() > 2 ) {
          //maxZPML
          if ( ptCoord[2][i] > outer[1][2] )
            outer[1][2] = ptCoord[2][i];
        }
      }
    }
//     outer[0][0] = 0;
//     outer[1][0] = 1.5;
//     outer[0][1] = 0;
//     outer[1][1] = 2;


    out << "PML layer coordinates:\n"
        << "   xmin = " << outer[0][0] << std::endl
        << "   xmax = " << outer[1][0] << std::endl
        << "   ymin = " << outer[0][1] << std::endl
        << "   ymax = " << outer[1][1] << std::endl;
    if ( dim_ == 3) {
      out << "   zmin = " << outer[0][2] << std::endl
          << "   zmax = " << outer[1][2] << std::endl;
    }

    out << std::endl;
    Info->PrintF( pdename_, out.str().c_str() );

  }



  bool SinglePDE::IsRegionPiezoHyst( std::string regionName ) {

    bool isPiezo = false;
    bool isHyst = false;

    // get direct coupled PDE
    DirectCoupledPDE* ptCoupledPDE = domain->GetDirectCoupledPDE();

    //! Get couplings object
    StdVector<BasePairCoupling*>* couplings = ptCoupledPDE->GetCouplingsObject();

    UInt idx = 0;
    for ( UInt i=0; i<couplings->GetSize(); i++ ) {
      if ( (*couplings)[i]->GetName() == "piezoDirect" ) {
        isPiezo = true;
        idx = i;
      }
    }

    if ( isPiezo ) {

      //      found = true;
      ParamNode * nonLinNode = NULL;
      nonLinNode = (*couplings)[idx]->GetParamNode()->Get("nonLinList", false);
      if ( !nonLinNode )
        return false;

      ParamNode* regionNode = NULL;
      if (  (*couplings)[idx]->GetParamNode()->
            Get("regionList")->Has("region", "name", regionName) ) {
        regionNode = (*couplings)[idx]->GetParamNode()->
          Get("regionList")->Get("region", "name", regionName);
      }
      if( !regionNode)
        return false;

      // check for nonLin Id
      std::string nonLinId;
      regionNode->Get("nonLinId", nonLinId, false);

      if (nonLinId == "" )
        return false;

      // check for hystersis nonlinearity
      isHyst = nonLinNode->Has("hysteresis", "id", nonLinId);
      return isHyst;
    }

    return false;

  }

  bool SinglePDE::BelongsPDE2PiezoHyst( ) {

  	bool isHyst = false;
  	RegionIdType actRegion;

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Set current region and material
      actRegion = it->first;

      // Get current region node
      std::string regionName = ptgrid_->GetRegion().ToString( actRegion );

      if ( IsRegionPiezoHyst( regionName ) )
      	isHyst = true;
    }

    return isHyst;
  }
  

  bool SinglePDE::IsRegionMicroPiezo( std::string regionName ) {

    bool isPiezo = false;
    bool isMicroPiezo = false;

    // get direct coupled PDE
    DirectCoupledPDE* ptCoupledPDE = NULL;
    
    ptCoupledPDE = domain->GetDirectCoupledPDE();

    if( !ptCoupledPDE )
      return false;

    //! Get couplings object
    StdVector<BasePairCoupling*>* couplings = ptCoupledPDE->GetCouplingsObject();

    UInt idx = 0;
    for ( UInt i=0; i < couplings->GetSize(); i++ ) {
      if ( (*couplings)[i]->GetName() == "piezoDirect" ) {
        isPiezo = true;
        idx = i;
      }
    }

    if ( isPiezo ) {

      // found = true;
      ParamNode * nonLinNode = NULL;
      nonLinNode = (*couplings)[idx]->GetParamNode()->Get("nonLinList", false);
      if ( !nonLinNode )
        return false;

      ParamNode* regionNode = NULL;
      regionNode = (*couplings)[idx]->GetParamNode()->
        Get("regionList")->Get("region", "name", regionName);
      if( !regionNode)
        return false;

      // check for nonLin Id
      std::string nonLinId;
      regionNode->Get("nonLinId", nonLinId, false);

      if (nonLinId == "" )
        return false;

      // check for micro-piezo
       isMicroPiezo = nonLinNode->Has("piezoMicroHF", "id", nonLinId);
      return isMicroPiezo;
    }

    return false;

  }

  bool SinglePDE::BelongsPDE2MicroPiezo( ) {
  	
    bool isMicroPiezo = false;
    RegionIdType actRegion;
  	
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      
      // Get current region node
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      if ( IsRegionMicroPiezo( regionName ) ) 
      	isMicroPiezo = true;	
    }
    
    return isMicroPiezo;
  }
  


  // ========================== RegionLoads ==========================
   SinglePDE::RegionLoad::RegionLoad( UInt dim, bool isAxi ) {

     value.Resize( dim );
     value.Init( "0.0");

     isAxi_ = isAxi;
     volume = 0.0;

   }


   VolForceInt * SinglePDE::RegionLoad::GetIntegrator() {
     VolForceInt * forceInt = new VolForceInt( value.GetSize(),
    		 																				phase,
    		 																				isAxi_);

     // Check, if type is "unit"
     bool isUnit;
     if ( type == "total" ) {
       isUnit = false;
     } else {
       isUnit = true;
     }
     forceInt->SetVolForceVector( value,
                                  domain->GetCoordSystem( refCoord),
                                  isUnit, volume );

     return forceInt;


   }

   void SinglePDE::RegionLoad::Print( bool onlyHeader, std::string pdeName ) {
     std::ostringstream out;

     if (onlyHeader) {
       Info->PrintF(pdeName, "The following regions have a region load:\n\n");
       out.clear();
       out << std::setw(15) << "name" << " | "
           << std::setw(15) << "coordSysId" << " | "
           << std::setw(11) << "volume" << " | "
           << std::setw(6) << "type" << " | "
           << std::setw(11) << "load" <<std::endl;
       Info->PrintF(pdeName, out.str().c_str());
       out.str("");
       out << std::setw(90) << std::setfill('-') << ""
           << std::setfill(' ') << std::endl;
       Info->PrintF(pdeName, out.str().c_str());
       out.str("");
     } else {


       // write logging information into info file
       for (UInt k = 0; k < value.GetSize(); k++ ) {
         out.str("");
         if ( k == 0) {
           out << std::setw(15) << name << " | "
               << std::setw(15) << refCoord << " | "
               << std::setw(11) << volume << " | "
               << std::setw(6) << type << "|";
         } else {
           out << std::setw(15) << "" << " | "
               << std::setw(15) << "" << " | "
               << std::setw(15) << "" << " | "
               << std::setw(11) << "" << " | "
               << std::setw(6) << "" << " | ";

         }

         out << std::setw(11) << value[k] << std::endl;

         Info->PrintF(pdeName,out.str().c_str());
       }

     }
   }
} // end of namespace
