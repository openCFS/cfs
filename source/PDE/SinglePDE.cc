#include "PDE/SinglePDE.hh"

#include <fstream>

// use of interpolation
#include <def_use_interpolation.hh>

// for coordinate handling
#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

#include "Utils/EvalIntegrals/BiotSavart.hh"

#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"

// header for logging
#include "DataInOut/Logging/LogConfigurator.hh"

// header for Materialhandling
#include "DataInOut/MaterialHandler.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

// header for Solvestep and assemble
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TransientDriver.hh"
#include "Driver/HarmonicDriver.hh"

// header for iterative coupling
#include "CoupledPDE/PDECoupling.hh"

// header for memento/restart handling
#include "MatVec/vectorSerialization.hh"

// header for resultHandling
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ResultCache.hh"

#include "DataInOut/PostProc.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
//#include "Forms/linearForm.hh"

//feSpaces
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/FeSpaceH1Hi.hh"
using std::string;

//coefFunctions
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

// TEMPORARY
#include "MagEdgePDE.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"


namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(singlepde)
  DEFINE_LOG(singlepde, "singlePde")

  SinglePDE::SinglePDE( Grid *aptgrid, PtrParamNode paramNode ) :
    StdPDE( aptgrid, paramNode ),
    isDirectCoupled_(false),
    isInitialized_(false),
    maxTimeDerivOrder_(0),
    mHandle_(domain->GetMathParser()->GetNewHandle()) // Obtain mathParser handler
  {
    
    // get id for linear system
    std::string systemId = myParam_->Get("systemId")->As<std::string>();
    
    PtrParamNode ls = myParam_->GetParent()
        ->GetParent()->Get("linearSystems",ParamNode::INSERT);
    olasNode_ = ls->GetByVal("system", "id", systemId, ParamNode::INSERT);
    
  }



  // **************
  //   Destructor
  // **************
  SinglePDE::~SinglePDE() {


    // Delete algebraic system only if
    // PDE is not direct coupled
    if ( isDirectCoupled_ == false ) {
      if ( needsAlgsys_ && (algsys_ != NULL)) {
        delete algsys_;
        algsys_ = NULL;
      }
      if (solveStep_) {
        delete solveStep_;
        solveStep_ = NULL;
      }
      if (assemble_) { 
        delete assemble_;
        assemble_ = NULL;
      }
    }


    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      if (it->second != NULL) delete it->second;
    }
    materials_.clear();

  }


  // ********
  //   Init
  // ********
  void SinglePDE::Init( UInt sequenceStep, PtrParamNode base ) {

    sequenceStep_ = sequenceStep;

    infoNode_ = base == NULL ? info->Get("PDE")->Get(pdename_) : base->Get(pdename_);
    infoNode_->Get(ParamNode::HEADER)->Get("sequenceStep")->SetValue(sequenceStep);

    LOG_TRACE(singlepde) << pdename_ << ": Starting Initialization";


    // =====================================================================
    // Get type of analysis
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Obtaining analysis type";
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

    LOG_TRACE(singlepde) << pdename_ << ": Obtaining regions";


    // Obtain regions the pde is defined on
    ParamNodeList regionNodes =
      myParam_->Get("regionList")->GetList("region");

    // output to info-file
    PtrParamNode list = infoNode_->Get(ParamNode::HEADER);

    // output and set subdoms_
    for( UInt i = 0; i < regionNodes.GetSize(); i++ )
    {
      PtrParamNode in_ = list->Get("region");
      in_->Get("name")->SetValue(regionNodes[i]->Get("name")->As<std::string>());
      bool complexMat = false;
      regionNodes[i]->GetValue("complexMaterial",  complexMat, ParamNode::PASS );
      
      RegionIdType actRegionId = ptgrid_->GetRegion().Parse(regionNodes[i]->Get("name")->As<std::string>());
      
      // Check, if region was already defined an issue a warning otherwise
      if( std::find(subdoms_.Begin(), subdoms_.End(), actRegionId ) 
          != subdoms_.End() )  {
        WARN( "The region '" << regionNodes[i]->Get("name")->As<std::string>()
              << "' was already defined for PDE '" << pdename_ 
              << "'. Please remove duplicate entries." );
      }
          
      subdoms_.Push_back( actRegionId );

      complexMatData_[actRegionId] = complexMat;
    }


    // Generate a fitting algebraic system only if PDE is NOT
    // direct coupled
    if( needsAlgsys_ == true ) {
      if ( isDirectCoupled_ == false) {
        olasInfo_ = info->Get("OLAS")->Get(pdename_);
        algsys_ = new AlgebraicSys(olasNode_, olasInfo_);
      }
    }

    // =====================================================================
    // create assemble class
    // =====================================================================

    // Create new assemble class with according analysistype
    if( isDirectCoupled_ == false && needsAlgsys_ == true) {
      assemble_ = new Assemble( algsys_, analysistype_, maxTimeDerivOrder_ );
    }

    //======================================================================
    // trigger the creation of functionDescriptors
    //======================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Define FE-Functions";
    DefineFeFunctions();
    
    // Register all fe functions with the algebraic system
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      shared_ptr<BaseFeFunction> actFct = fncIt->second;
      shared_ptr<FeSpace> actSpace = fncIt->second->GetFeSpace();
      std::string fctName = SolutionTypeEnum.ToString(fncIt->first);
      FeFctIdType fctId = algsys_->ObtainFctId( fctName );
      actFct->SetFctId(fctId);
      fncIt++;
    }


    // =====================================================================
    // trigger definition of available results
    // =====================================================================
    DefinePrimaryResults();
    
    // =====================================================================
    // read in material data
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Reading material information";
    ReadMaterialData();

    // =====================================================================
    // read in damping information
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Reading damping information";
    ReadDampingInformation( );

    // =====================================================================
    // read in NonLinearities
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Initializing non-linearities";
    InitNonLin();

    // Todo: Move this part to the definition of damping
    PtrParamNode in = infoNode_->Get(ParamNode::HEADER);
    for(UInt i = 0; i < subdoms_.GetSize(); i++ )
    {
      PtrParamNode in_ = in->GetByVal("region", "name", domain->GetGrid()->GetRegion().ToString(subdoms_[i]));

      std::map<RegionIdType,DampingType>::const_iterator it = dampingList_.find(subdoms_[i]);
      DampingType dampType = NONE;
      if( it != dampingList_.end() ) {
        dampType = it->second;
      }
      std::string fuck_e2s;
      Enum2String(dampingList_[subdoms_[i]], fuck_e2s);
      in_->Get("damping")->SetValue(fuck_e2s);
    }

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

    LOG_TRACE(singlepde) << pdename_ << ": Reading boundary conditions";
    ReadBCs();
    ReadSpecialBCs();
    //

    // =====================================================================
    // Create time stepping algorithm
    // =====================================================================
    if ( analysistype_ == TRANSIENT &&
         isDirectCoupled_ == false) {
      InitTimeStepping();
    }

    // =====================================================================
    // define the integrators for PDE and initialize eqn object
    // =====================================================================

    // Call initialization of (bi)linear integrators
    LOG_TRACE(singlepde) << pdename_ << ": Defining integrators";
    DefineIntegrators();
    DefineSurfaceIntegrators();
    DefineRhsLoadIntegrators();

    // Print information about defined integrators
    if( needsAlgsys_ == true)
      assemble_->ToInfo(infoNode_->Get(ParamNode::HEADER)->Get("integrators"));

    // now we know about nonlinearities and we can trigger the
    // material objects to perform the approximations of the nonlinear
    // sampled data
    std::map<RegionIdType, BaseMaterial*>::iterator itMat;
    for ( itMat = materials_.begin(); itMat != materials_.end(); itMat++ ) {
      itMat->second->InitApproxCurves();
    }
    
    // =====================================================================
    //  map equations (FeSpaces) and finalize FeFunction (vector creation)
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Mapping Equations";
    
    // Finalize spaces and fefunctions
    fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      shared_ptr<BaseFeFunction> actFct = fncIt->second;
      shared_ptr<FeSpace> actSpace = fncIt->second->GetFeSpace();
      //actFct->SetGrid(shared_ptr<Grid>(ptgrid_));
      actSpace->Finalize();
      actSpace->PreCalcShapeFncs();
      
      // finalize feFunctions
      actFct->Finalize();
      
      // Pass feFctId of primary result also to RHS result
      rhsFeFunctions_[fncIt->first]->SetFctId(actFct->GetFctId());
      rhsFeFunctions_[fncIt->first]->Finalize();
      fncIt++;
    }
    

    if ( analysistype_ == TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain->GetSingleDriver())
        ->GetDeltaT();
      //WARN("Note: The initialization of the timestepping class is currently wrong: "
      //    "The 2nd argument must be the complete SBM-vector of the algebraic system in "
      //    "order to correclty initialize the internal vectors of the timestepping method. "
      //    "In the old implementation it was sufficient to know the number of unknowns. "
      //    "In the current implementation, the SBM-vectors are just defined within the "
      //    "SolveStep classed. Thus maybe the right thing to do is to shift the creation and "
      //    "initialization of the timestepping scheme to the solveStep classes.")
          

      // Call the init function of timescheme of each fefunction
      fncIt= feFunctions_.begin();
      while(fncIt != feFunctions_.end()){
        shared_ptr<BaseFeFunction> actFct = fncIt->second;
        actFct->GetTimeScheme()->Init(actFct->GetSingleVector(),dt);
        fncIt++;
      }

    }

    // =======================================================================
    // Trigger creation of timeDerivative FeFunctions
    // NOTE: Has to be done here after initializaiton of timestepping scheme!
    // =======================================================================
    DefineTimeDerivFeFunctions();

//    // =====================================================================
//    // Set the initial conditions
//    // =====================================================================
//    if ( analysistype_ == TRANSIENT){
//    	SetInitialCondition();
//    }

    // =====================================================================
    // define which solution types have to be saved
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Reading store results";

    DefinePostProcResults();
    ReadStoreResults();
    ReadFieldResults();

    PreparePDE4Computation();

    //! Define step solution driver
    if ( isDirectCoupled_ == false ) {
      LOG_TRACE(singlepde) << pdename_ << ": Defining solveStep class";
      DefineSolveStep();
      
      // check if solve step was defined
      if(! solveStep_) {
        EXCEPTION("No solveStep object defined for PDE '" << pdename_ << "'");
      }
    }

    // Finally set the initialization flag to true
    isInitialized_ = true;
    LOG_TRACE(singlepde) << pdename_ << ": Finished initializaton";
  }


  SubTensorType SinglePDE::GetSubTensorType() const
  {
    if(subType_ == "") return NO_TENSOR;

    SubTensorType stt;
    String2Enum(subType_, stt);
    return stt;
  }

  
  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void SinglePDE::InitNonLin() {

    nonLin_ = false;

    // Check, if "nonLinList" is present
    PtrParamNode nonLinListNode = myParam_->Get("nonLinList", ParamNode::PASS );
    if( nonLinListNode ) { 

      // Get nonlinear types
      ParamNodeList nonLinNodes = nonLinListNode->GetChildren();
      for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

        std::string actTypeString = nonLinNodes[i]->GetName();
        std::string actId = nonLinNodes[i]->Get("id")->As<std::string>();

        NonLinType actType;
        String2Enum( actTypeString, actType );

        //save for each nonlinearity type the id
        nonLinTypes_[actId] = actType;
      }
    

      // Run over all region and set entry in "regionNonLinId"
      ParamNodeList regionNodes = 
        myParam_->Get("regionList")->GetChildren();
      
      RegionIdType actRegionId;
      std::string actRegionName, actNonLinId;
      
      //     if( regionNodes.GetSize() > 0 ) {
      //       Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
      //     }
      
      for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
        //take cae: one region can have more then one nonlinearity!!
        
        // get data
        regionNodes[i]->GetValue( "name", actRegionName );
        regionNodes[i]->GetValue( "nonLinIds", actNonLinId );
        
        if( actNonLinId == "" )
          continue;
        
        typedef boost::tokenizer< boost::char_separator<char> > Tok;
        boost::char_separator<char> sep(";|, ");
        
        Tok tok(actNonLinId, sep);
        
        actRegionId = ptgrid_->GetRegion().Parse( actRegionName );
        
        for(Tok::iterator it=tok.begin(); it!=tok.end(); ++it) {
          std::string nonLinId = (*it);
          
          if(nonLinTypes_.find(nonLinId) == nonLinTypes_.end()) {
            WARN( "NonLinearity with id '" << nonLinId 
                  << "' was not defined in 'nonLinList'");
            continue;
          }
          
          regionNonLinTypes_[actRegionId].Push_back( nonLinTypes_[nonLinId] );
          
          //write info
          std::string nonLinString;
          Enum2String( nonLinTypes_[nonLinId], nonLinString );
          //         Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(), 
          //                       nonLinString.c_str() );
          
          //if one nonlinearity is set, then the whole PDE is set to nonlinear
          nonLin_ = true;
        }
      }

      // Here we need in addition the nonLinMethod_ for the definition
      // of the integrators
      nonLinMethod_ = FIXEDPOINT;
      PtrParamNode nonLinNode = myParam_->Get("nonLinear", ParamNode::PASS );
      if( nonLinNode ) {
        std::string methodString;
        nonLinNode->GetValue(  "method", methodString, ParamNode::PASS );
        nonLinMethod_ = NonLinMethodTypeEnum.Parse(methodString);
      }
    }
  }

   /** can generally be called multiple times. We overwrite old values! Brute force but keeps data size */
   void SinglePDE::WriteGeneralPDEdefines()
   {
    // loads
    PtrParamNode base = infoNode_->Get(ParamNode::HEADER)->Get("loads");

    for(unsigned int i = 0, nn = loads_.GetSize(); i < nn; i++)
    {
      LoadBc const & actBc = *loads_[i];
      EntityList const & actList = *actBc.entities;

      PtrParamNode in = base->GetByVal(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
      in->Get("value")->SetValue(actBc.value);
      in->Get("phase")->SetValue(actBc.phase);
    }

    // Homogeneous Dirichlet BC
    base = infoNode_->Get(ParamNode::HEADER)->Get("homDirichletBC");

    for(unsigned int i = 0, nn = hdBcs_.GetSize(); i < nn; i++)
    {
      HomDirichletBc const & actBc = *hdBcs_[i];
      EntityList const & actList = *actBc.entities;

      PtrParamNode in = base->GetByVal(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
    }

    // Inhomogeneous Dirichlet BC
    base = infoNode_->Get(ParamNode::HEADER)->Get("inhomDirichletBC");

    for(unsigned int i = 0, nn = idBcs_.GetSize(); i < nn; i++)
    {
      InhomDirichletBc const & actBc = *idBcs_[i];
      EntityList const & actList = *actBc.entities;

      PtrParamNode in = base->GetByVal(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
      in->Get("value")->SetValue(actBc.value);
      in->Get("phase")->SetValue(actBc.phase);
    }

    // Inhomogeneous Dirichlet BC read from file
    base = infoNode_->Get(ParamNode::HEADER)->Get("inhomDirichFileBC");

    for(unsigned int i = 0, nn = idFiBcs_.GetSize(); i < nn; i++)
    {
      InhomDirichFileBc const & actBc = *idFiBcs_[i];
      EntityList const & actList = *actBc.entities;

      PtrParamNode in = base->GetByVal(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
      in->Get("inputId")->SetValue(actBc.inputId);
    }

    // Inhomogeneous Neumann BC
    base = infoNode_->Get(ParamNode::HEADER)->Get("inhomNeumannBC");

    for(unsigned int i = 0, nn = inBcs_.GetSize(); i < nn; i++)
    {
      InhomNeumannBc const & actBc = *inBcs_[i];
      EntityList const & actList = *actBc.entities;

      PtrParamNode in = base->GetByVal(actList.listType.ToString(actList.GetType()), "name", actList.GetName());
      in->Get("dof")->SetValue(actBc.result->GetDofName(actBc.dof));
      in->Get("value")->SetValue(actBc.value);
      in->Get("phase")->SetValue(actBc.phase);
    }

    // constraints
    base = infoNode_->Get(ParamNode::HEADER)->Get("constraints");
    // periodic boundary conditions blow this up.
    if(constraints_.GetSize() <= 5 )
    {
      for(unsigned int i = 0, nn = constraints_.GetSize(); i < nn; i++)
      {
        Constraint const & actBc = *constraints_[i];
        EntityList const & masterList = *actBc.masterEntities;
        EntityList const & slaveList = *actBc.slaveEntities;

        PtrParamNode in = base->GetByVal("pair", "master", masterList.GetName());
        in->Get("slave")->SetValue(slaveList.GetName());
        // the names are repeated for the different dofs
        std::string dof = actBc.result->GetDofName(actBc.masterDof);
        if(!in->HasByVal("dof", dof))
          in->Get("dof", ParamNode::APPEND)->SetValue(dof);

        in->Get("periodic")->SetValue(actBc.periodic);
      }
    }
    else
    {
      if(constraints_.GetSize() > 5)
      {
        base->Get("number")->SetValue(constraints_.GetSize());
        base->SetComment("run cfs with -l for list");
      }
    }
  }

  void SinglePDE::ReadStoreResults() {

    ResultSet::iterator it;

    // fetch result node and leave, if none is present
    PtrParamNode resultNode = myParam_->Get("storeResults", ParamNode::PASS);
    if( !resultNode )
      return;

    // Iterate over all availabe results
    for (it = availResults_.begin(); it != availResults_.end(); it++ ) {
       CheckStoreResult(*it);
    }
  }


  bool SinglePDE::CheckStoreResult(shared_ptr<ResultInfo> candidate)
  {

    StdVector<std::string> regionNames, nodeNames, writeResults, actOutDest;
    StdVector<std::string> postProcNames, outDestNames, neighborRegions;
    UInt saveBegin = 0, saveEnd = 0, saveInc = 0;
    std::string quantity, complexFormatString, listElemName, entityName;
    ComplexFormat complexFormat;
    shared_ptr<EntityList> actList;

    EntityList::ListType entityType;
    ResultHandler * resHandler = domain->GetResultHandler();

    // initialize map for relating EntityUnknownType and name of xml-element
    using std::make_pair;
    std::map<ResultInfo::EntityUnknownType, std::string> elemNames;
    std::map<ResultInfo::EntityUnknownType, bool> isHistory;
    elemNames.insert(make_pair(ResultInfo::NODE, "nodeResult"));
    elemNames.insert(make_pair(ResultInfo::ELEMENT, "elemResult"));
    elemNames.insert(make_pair(ResultInfo::SURF_ELEM, "surfElemResult"));
    elemNames.insert(make_pair(ResultInfo::REGION, "regionResult"));
    elemNames.insert(make_pair(ResultInfo::SURF_REGION, "surfRegionResult"));

    isHistory.insert(make_pair(ResultInfo::NODE, false));
    isHistory.insert(make_pair(ResultInfo::ELEMENT, false));
    isHistory.insert(make_pair(ResultInfo::SURF_ELEM, false));
    isHistory.insert(make_pair(ResultInfo::REGION, true));
    isHistory.insert(make_pair(ResultInfo::SURF_REGION, true));
    

    // fetch result node and leave, if none is present
    PtrParamNode resultNode = myParam_->Get("storeResults", ParamNode::PASS);
    if( !resultNode )
      return false;

      // Convert enum
      quantity = SolutionTypeEnum.ToString(candidate->resultType);
      LOG_DBG(singlepde) << "Searching for storeResults of quantity '"
                   << quantity << "'";

      // try to catch possible errors
      try {

        // Get type of result
        std::string xmlElemName = elemNames[candidate->definedOn];
        if( xmlElemName == "" ){
          return false;
        }

        // Remember current result node
        PtrParamNode actResultNode =
          resultNode->GetByVal(xmlElemName, "type", quantity, ParamNode::PASS );

        // Check on which entity type the result is defined on
        switch(candidate->definedOn)
        {
        case ResultInfo::NODE:
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
          break;
        }

        // intialize variables
        neighborRegions.Clear();
        regionNames.Clear();
        //outDestNames.Clear();

        // ========== Look for defineType 'REGION' ==========
        // 1a) Look if result is defined on 'allRegions'

        // if no node was found, continue with next result
        if( !actResultNode) {
          return false;
        }

        // determine complexFormat
        complexFormatString = "amplPhase";
        actResultNode->GetValue("complexFormat", complexFormatString, ParamNode::PASS);
        String2Enum( complexFormatString, complexFormat );

        // otherwise check, if result is to be saved on "allRegions"
        if( actResultNode->Has("allRegions" ) ) {
          ptgrid_->GetRegion().ToString(subdoms_,regionNames);

          PtrParamNode allRegionsNode = actResultNode->Get("allRegions");

          std::string allPostProcName, allOutDestName;

          // fetch postProcNames
          allRegionsNode->GetValue("postProcId", allPostProcName );
          postProcNames.Resize( regionNames.GetSize() );
          postProcNames.Init( allPostProcName );

          //fetch outDestName
          allRegionsNode->GetValue("outputIds", allOutDestName );
          outDestNames.Resize( regionNames.GetSize() );
          outDestNames.Init( allOutDestName );

          // fetch saveBegin, saveEnd and saveInc
          saveBegin = allRegionsNode->Get("saveBegin")->MathParse<UInt>();
          saveEnd = allRegionsNode->Get("saveEnd")->MathParse<UInt>();
          saveInc = allRegionsNode->Get("saveInc")->MathParse<UInt>();

          // fetch writeResult flag
          std::string writeResult;
          allRegionsNode->GetValue("writeResult", writeResult );
          writeResults.Resize( regionNames.GetSize() );
          writeResults.Init( writeResult );

        } else {

          ParamNodeList regionNodes;
          PtrParamNode listNode;
          // 1b) Look for regions the result is defined on
          if(candidate->definedOn == ResultInfo::NODE ||
             candidate->definedOn == ResultInfo::ELEMENT ||
             candidate->definedOn == ResultInfo::REGION ) {
            listNode = actResultNode->Get("regionList", ParamNode::PASS);
            if( listNode )
              regionNodes = listNode->GetList("region");
          } else if(candidate->definedOn == ResultInfo::SURF_ELEM ||
                    candidate->definedOn == ResultInfo::SURF_REGION ) {
            listNode = actResultNode->Get("surfRegionList", ParamNode::PASS);
            if( listNode )
              regionNodes = listNode->GetList("surfRegion");

            // fetch entry with neighboring regions
            for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
              neighborRegions.Push_back( regionNodes[i]->
                                         Get("neighborRegion")->As<std::string>() );
            }
          }

          // only enter, at least one region is present
          if( listNode ) {
            // fetch saveBegin, saveEnd and saveInc
            saveBegin = listNode->Get("saveBegin")->MathParse<UInt>();
            saveEnd = listNode->Get("saveEnd")->MathParse<UInt>();
            saveInc = listNode->Get("saveInc")->MathParse<UInt>();

            // iterate over all regions
            regionNames.Clear();
            postProcNames.Clear();
            outDestNames.Clear();
            writeResults.Clear();
            for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
              regionNames.Push_back( regionNodes[i]->Get("name")->As<std::string>() );
              postProcNames.Push_back( regionNodes[i]->Get("postProcId")->As<std::string>() );
              outDestNames.Push_back( regionNodes[i]->Get("outputIds")->As<std::string>() );
              writeResults.Push_back( regionNodes[i]->Get("writeResult")->As<std::string>() );
            }
          }
        }

        // Check, if any region was found for this result type
        if( regionNames.GetSize() != 0 ) {
          candidate->complexFormat = complexFormat;

          // iterate over all regions
          for( UInt iRegion = 0; iRegion < regionNames.GetSize(); iRegion++ )
          {
            actList = ptgrid_->GetEntityList( entityType, regionNames[iRegion] );
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
            resHandler->RegisterResult( actSol, sequenceStep_, 
                                        saveBegin, saveInc, saveEnd,
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

        PtrParamNode histNode;
        ParamNodeList histEntities;

        if(candidate->definedOn == ResultInfo::NODE ) {
          histNode = actResultNode->Get("nodeList",ParamNode::PASS);
          if( histNode )
            histEntities = histNode->GetList("nodes");
          entityTypeName = "nodes";

        } else if(candidate->definedOn == ResultInfo::ELEMENT ) {
          histNode = actResultNode->Get("elemList",ParamNode::PASS);
          if( histNode )
            histEntities = histNode->GetList("elems");
          entityTypeName = "elements";

        } else if(candidate->definedOn == ResultInfo::SURF_ELEM ) {
          histNode = actResultNode->Get("surfElemList",ParamNode::PASS);
          if( histNode)
            histEntities = histNode->GetList("surfElems");
          entityTypeName = "surfElems";

          // fetch entry with neighboring regions
          // fetch entry with neighboring regions
          for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
            neighborRegions.Push_back( histEntities[i]->
                                       Get("neighborRegion")->As<std::string>() );
          }
        }

        // only proceed, if any history result is defined
        if( histNode && histNode->HasChildren() ) {

          // fetch saveBegin, saveEnd and saveInc
          saveBegin = histNode->Get("saveBegin")->MathParse<UInt>();
          saveEnd = histNode->Get("saveEnd")->MathParse<UInt>();
          saveInc = histNode->Get("saveInc")->MathParse<UInt>();

          // iterate over all regions
          histNames.Clear();
          postProcNames.Clear();
          outDestNames.Clear();
          writeResults.Clear();
          for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
            histNames.Push_back( histEntities[i]->Get("name")->As<std::string>() );
            postProcNames.Push_back( histEntities[i]->Get("postProcId")->As<std::string>() );
            outDestNames.Push_back( histEntities[i]->Get("outputIds")->As<std::string>() );
            writeResults.Push_back( histEntities[i]->Get("writeResult")->As<std::string>() );
          }
        }

        if( histNames.GetSize() > 0 )
        {
          candidate->complexFormat = complexFormat;

          // iterate over all entityNames
          for( UInt i = 0; i < histNames.GetSize(); i++ )
          {
            actList = ptgrid_->GetEntityList( entityType, histNames[i] );
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

            resHandler->RegisterResult( actSol, sequenceStep_, 
                                        saveBegin, saveInc, saveEnd,
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
    LOG_DBG(singlepde) << "WriteResultsInFile() kstep: " <<  kstep
                 << " actTimeFreq: " << actTimeFreq;
    ResultMap::iterator it = resultLists_.begin();
    ResultHandler * resHandler = domain->GetResultHandler();

    // =======================================
    //  Trigger calculation of normal results 
    // =======================================
    
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
    
    
    // ===================================================
    //  Trigger calculation of interpolated field results 
    // ===================================================
    
    // Check for additional field variable
    UInt numFields = fields_.GetSize();
    
    // loop over all fields variables
    for( UInt i = 0; i < numFields; ++i ) {
    
      // call specialized calculation method in sub-class
      FieldAtPoints& fap = fields_[i];
      
      
      // Obtain field resultFunctor object
      SolutionType solType = fap.resultInfo->resultType;
      UInt numDofs = fap.resultInfo->dofNames.GetSize();
      std::map<SolutionType, shared_ptr<BaseFieldFunctor> >::iterator fctIt;
      fctIt = fieldFunctors_.find(solType);
      if( fctIt == fieldFunctors_.end() )  {
        EXCEPTION( "Could not find field functor for result '" 
            << SolutionTypeEnum.ToString(solType) << "'");
      }
       
      // calculate vector entries
      if( isComplex_) {
        shared_ptr<FieldFunctor<Complex> > fct = 
            dynamic_pointer_cast<FieldFunctor<Complex> >(fctIt->second);
        Vector<Complex> temp;
        Vector<Complex>& vec = dynamic_cast<Vector<Complex> &>(*fap.field);
        vec.Resize(fap.elems.GetSize() * numDofs);
        UInt pos = 0;
        LocPointMapped lpm;
        for ( UInt iElem = 0; iElem < fap.elems.GetSize(); ++iElem ) {
          shared_ptr<ElemShapeMap> esm = 
              ptgrid_->GetElemShapeMap( fap.elems[iElem] );
          lpm.Set(fap.locPoints[iElem], esm);
          fct->GetVector(temp, lpm );
          
          for( UInt i = 0; i < numDofs; ++i ) {
            vec[pos++] = temp[i];
          }
        }
      } else {
        shared_ptr<FieldFunctor<Double> > fct = 
            dynamic_pointer_cast<FieldFunctor<Double> >(fctIt->second);
        Vector<Double> temp;
        Vector<Double>& vec = dynamic_cast<Vector<Double> &>(*fap.field);
        vec.Resize(fap.elems.GetSize() * numDofs);
        UInt pos = 0;
        LocPointMapped lpm;
        for ( UInt iElem = 0; iElem < fap.elems.GetSize(); ++iElem ) {
          shared_ptr<ElemShapeMap> esm = 
              ptgrid_->GetElemShapeMap( fap.elems[iElem] );
          lpm.Set(fap.locPoints[iElem], esm);
          fct->GetVector(temp, lpm );

          for( UInt i = 0; i < numDofs; ++i ) {
            vec[pos++] = temp[i];
          }
        }
      }
      
      // print out information
      if( isComplex_ ){
        EXCEPTION("not yet implemented");
      } else {
        // cast solution vector
        Vector<Double>& vec = dynamic_cast<Vector<Double> &>(*(fap.field));
        std::ofstream  out(fap.fileName.c_str(),std::ios::out );

        // obtain result info and print header
        UInt numDofs = fap.resultInfo->dofNames.GetSize();
        
        
        // Loop over all points
        Vector<Double> globPoint;
        for( UInt iPoint = 0; iPoint < fap.locPoints.GetSize(); iPoint++) { 
          
          ElemShapeMap& esm = *(ptgrid_->GetElemShapeMap(fap.elems[iPoint], true));
          esm.Local2Global(globPoint, fap.locPoints[iPoint]);
          // write to file
          out << fap.elems[iPoint]->elemNum << "\t";
          out << globPoint.ToString(0, '\t') << "\t";
          for(UInt j = 0; j < numDofs; ++j ) {
            out << vec[iPoint*numDofs + j] << "\t";
          }
          for(UInt j = 0; j < dim_; ++j ) {
            out << fap.locPoints[iPoint][j] << "\t";
          }
          out << fap.elems[iPoint]->elemNum << "\t";
          out << std::endl;
        }
        out.close();
      
      }
      
    }
  }
  
  
  void SinglePDE::ReadFieldResults() {
    // check, if calculation of field variables is requested at all

    ParamNodeList fieldNodes;
    fieldNodes = myParam_->Get("storeResults")->GetList("field");
    std::string solTypeString;
    
    fields_.Resize(fieldNodes.GetSize());
    // loop over all parts
    for( UInt iPart = 0; iPart <fieldNodes.GetSize(); ++iPart ) {
      PtrParamNode  actNode = fieldNodes[iPart];
      ParamNodeList listNodes = actNode->GetList("list");
      
      FieldAtPoints & actField = fields_[iPart];
      actField.fileName = actNode->Get("fileName")->As<std::string>();
      
      // check for solution type
      solTypeString = actNode->Get("type")->As<std::string>();
      
      SolutionType solType = SolutionTypeEnum.Parse(solTypeString);
      
      // find related result resultinfo
      ResultSet::const_iterator it = availResults_.begin();
      for( ; it != availResults_.end(); ++it ) {
        if( (*it)->resultType == solType ) {
          actField.resultInfo = *it;
          break;
        }
      }
      
      // loop over all components
      StdVector<Double> start(3), stop(3), inc(3);
      StdVector<UInt> numSamples(3);
      start.Init(0);
      stop.Init(0);
      inc.Init(1);
      numSamples.Init(1);
      

      UInt totalPoints = 1;
      std::string comp;
      UInt compIndex;
      for( UInt iComp = 0; iComp < listNodes.GetSize(); iComp++ ) {
        PtrParamNode actCompNode = listNodes[iComp];
        actCompNode->GetValue("comp", comp);
        compIndex = domain->GetCoordSystem("default")->GetVecComponent(comp)-1;
        start[compIndex]=  actCompNode->Get("start")->MathParse<Double>();
        stop[compIndex]=  actCompNode->Get("stop")->MathParse<Double>();
        inc[compIndex] = actCompNode->Get("inc")->MathParse<Double>();
        numSamples[compIndex]  = 
            UInt(floor( (stop[compIndex]-start[compIndex]) / inc[compIndex] ) )+1;
        totalPoints *= numSamples[compIndex];
      }
      // create list

      // generate new vector
      if(isComplex_) {
        actField.field = new Vector<Complex>();
      } else {
        actField.field = new Vector<Double>();
      }
      
      Vector<Double> globPoint(3);
      for( UInt xSample = 0; xSample < numSamples[0]; xSample++ ) {
        Double actX = start[0] + xSample * inc[0];
        for( UInt ySample = 0; ySample < numSamples[1]; ySample++ ) {
          Double actY = start[1] + ySample * inc[1];
          for( UInt zSample = 0; zSample < numSamples[2]; zSample++ ) {
            Double actZ = start[2] + zSample * inc[2];
            globPoint[0] = actX;
            globPoint[1] = actY;
            if( dim_ > 2) {
              globPoint[2] = actZ;
            } 
            else {
             globPoint[2] = 0.0;
            }
//#define USE_INTERPOLATION            
#ifndef USE_INTERPOLATION
             EXCEPTION("Special Results can just be calculated with INTERPOLATION active");
#else
             Vector<Double> locPoint;
             const Elem * ptElem = NULL;
             // now, map global point to localpoint
             ptElem = ptgrid_->GetElemAtGlobalCoord( globPoint, locPoint );
             
             if( !ptElem ) {
               std::string warnStr = "Could not find element at position " 
                   + globPoint.ToString(); 
               WARN( warnStr.c_str());
             } else {
//               std::cerr << "locPoint for globPoint " << globPoint.ToString() 
//                                    << " is " << locPoint.ToString() 
//                                    << " in Elem " << ptElem->elemNum << std::endl;
               
               // check again mapping by performing loc->glob mapping
               shared_ptr<ElemShapeMap>  esm = ptgrid_->GetElemShapeMap(ptElem);
               //esm->Local2Global(globPoint, locPoint);
//               std::cerr << "\tAdditional check loc->glob delivers global point " 
//                   << globPoint.ToString() << std::endl << std::endl;
               
               actField.elems.Push_back(ptElem);
               actField.locPoints.Push_back(LocPoint(locPoint));
             }
#endif
            
          } // z
        } // y
      } // x
    } // loop over <field> entries
  }


  // **********
  //   SetBCs
  // **********
  void SinglePDE::SetBCs() {
    
    // TODO: Is this method really necessary here or can we move it to the SolveStep-class?
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      fncIt->second->ApplyBC();
      fncIt++;
    }
  }

  void SinglePDE::ReadBCs() {


    // fetch "bcsAndLoads" parameter node, if present.
    // otherwise leave
    PtrParamNode bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if( !bcsNode )
      return;

    std::string name, resultName, dof, entType, inputId, value, phase;
    shared_ptr<BaseFeFunction> actFeFunction;


    // =====================================================================
    // homogeneous Dirichlet BC
    // =====================================================================

    // fetch paramnodes for hdbc
    ParamNodeList hdbcNodes = bcsNode->GetList("dirichletHom");

    // iterate over all parameter nodes
    for( UInt i = 0; i < hdbcNodes.GetSize(); i++ ) {

      try {
        // read parameters
        dof.clear();
        hdbcNodes[i]->GetValue( "name", name );
        hdbcNodes[i]->GetValue( "quantity", resultName );
        hdbcNodes[i]->GetValue( "dof", dof, ParamNode::PASS );
        hdbcNodes[i]->GetValue( "entityType", entType );

        shared_ptr<HomDirichletBc> actBc ( new HomDirichletBc );
        EntityList::ListType listType = EntityList::listType.Parse(entType);
        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( listType, name);
        
        // fetch related feFunction object
        actFeFunction = GetFeFunction( SolutionTypeEnum.Parse(resultName));
        
        actBc->entities = actList;
        actBc->result = actFeFunction->GetResultInfo();
        if( dof.empty() ) {
          actBc->dof = 0;
        } else {
          actBc->dof = actFeFunction->GetResultInfo()->GetDofIndex( dof );
        }

        // add definition
        //hdBcs_.Push_back( actBc );
        actFeFunction->AddHomDirichletBc(actBc);
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create homogeneous boundary conditions on '"
                           << name << "'" );
      }
    }

    //=====================================================================
    // inhomogeneous Dirichlet BC
    // =====================================================================

    // fetch paramnodes for idbc
    ParamNodeList idbcNodes = bcsNode->GetList("dirichletInhom");

    // iterate over all parameter nodes
    for( UInt i = 0; i < idbcNodes.GetSize(); i++ ) {

      try {
        // read parameters
        dof.clear();
        idbcNodes[i]->GetValue( "name", name );
        idbcNodes[i]->GetValue( "quantity", resultName );
        idbcNodes[i]->GetValue( "dof", dof, ParamNode::PASS );
        idbcNodes[i]->GetValue( "value", value );
        idbcNodes[i]->GetValue( "phase", phase );
        idbcNodes[i]->GetValue( "entityType", entType );

        // Create inhomogeneous boundary condition
        shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );

        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
                                  name );
        
        // fetch related feFunction object
        actFeFunction = GetFeFunction( SolutionTypeEnum.Parse(resultName));
        
        actBc->entities = actList;
        actBc->result = actFeFunction->GetResultInfo();
        //actBc->eqnMap = eqnMap_;
        if( dof.empty() ) {
          actBc->dof = 0;
        } else {
          actBc->dof = actFeFunction->GetResultInfo()->GetDofIndex( dof );
        }
        actBc->value = value;
        actBc->phase = phase;

        actFeFunction->AddInhomDirichletBc(actBc);
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create inhomogeneous boundary conditions on '"
                           << name << "'" );
      }
     }

    //=====================================================================
    // inhomogeneous Dirichlet BC from File
    // =====================================================================

    // fetch paramnodes for idbc
    ParamNodeList idbcFiNodes = bcsNode->GetList("dirichletFileInhom");

    // iterate over all parameter nodes
    for( UInt i = 0; i < idbcFiNodes.GetSize(); i++ ) {
      EXCEPTION("Dirichlet BC from files not yet adapted");
//
//      try {
//        // read parameters
//        dof.clear();
//        idbcFiNodes[i]->GetValue( "name", name );
//        idbcFiNodes[i]->GetValue( "quantity", resultName );
//        idbcFiNodes[i]->GetValue( "dof", dof, ParamNode::PASS );
//        idbcFiNodes[i]->GetValue( "entityType", entType );
//        idbcFiNodes[i]->GetValue( "inputId", inputId );
//
//        // fetch related resultInfo object
//        actResultInfo = GetResultInfo( SolutionTypeEnum.Parse(resultName) );
//
//        // Create inhomogeneous boundary condition
//        shared_ptr<InhomDirichFileBc> actBc ( new InhomDirichFileBc );
//        if( entType == "nodeList" ) {
//          defineType = EntityList::NAMED_NODES;
//        } else {
//          defineType = EntityList::REGION;
//        }
//
//        shared_ptr<EntityList> actList =
//          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
//                                  name, defineType );
//        actBc->entities = actList;
//        actBc->result = actResultInfo;
//        actBc->eqnMap = eqnMap_;
//        if( dof.empty() ) {
//          actBc->dof = 1;
//        } else {
//          actBc->dof = actResultInfo->GetDofIndex( dof );
//        }
//        actBc->inputId = inputId;
//
//        // add definition
//        idFiBcs_.Push_back( actBc );
//      } catch (Exception & ex ) {
//        RETHROW_EXCEPTION( ex, "Can not create inhomogeneous boundary conditions on '"
//                           << name << "'" );
//      }
     }

    // =====================================================================
    // inhomogeneous Neumann BC
    // =====================================================================

    // fetch paramnodes for inbc
    ParamNodeList inbcNodes = bcsNode->GetList("neumannInhom");

    // iterate over all parameter nodes
    for( UInt i = 0; i < inbcNodes.GetSize(); i++ ) {
      try {
        dof.clear();
        inbcNodes[i]->GetValue( "name", name );
        inbcNodes[i]->GetValue( "quantity", resultName );
        inbcNodes[i]->GetValue( "dof", dof, ParamNode::PASS );
        inbcNodes[i]->GetValue( "value", value );
        inbcNodes[i]->GetValue( "phase", phase );
        inbcNodes[i]->GetValue( "entityType", entType );

        // Create inhomogeneous Neumann boundary condition
        shared_ptr<InhomNeumannBc> actBc ( new InhomNeumannBc );
        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
                                  name );
        
        // fetch related feFunction object
        actFeFunction = GetFeFunction( SolutionTypeEnum.Parse(resultName));
        
        actBc->entities = actList;
        actBc->result =  actFeFunction->GetResultInfo();
        if( dof.empty() ) {
          actBc->dof = 0;
        } else {
          actBc->dof = actFeFunction->GetResultInfo()->GetDofIndex( dof );
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
    ParamNodeList csNodes = bcsNode->GetList("constraint");
    std::string masterDof, slaveDof;

    // iterate over all parameter nodes
    for( UInt i = 0; i < csNodes.GetSize(); i++ ) {
      try {
        csNodes[i]->GetValue( "name", name );
        csNodes[i]->GetValue( "quantity", resultName );
        csNodes[i]->GetValue( "entityType", entType );
        csNodes[i]->GetValue( "masterDof", masterDof );
        csNodes[i]->GetValue( "slaveDof", slaveDof );

        // fetch related resultInfo object
        actFeFunction = GetFeFunction( SolutionTypeEnum.Parse(resultName) );

        // Create constraint condition
        shared_ptr<Constraint> actBc ( new Constraint );
        shared_ptr<EntityList> actList =
          ptgrid_->GetEntityList( EntityList::listType.Parse(entType),
                                  name  );

        actBc->masterEntities = actList;
        actBc->slaveEntities = actList;
        if( masterDof.empty() ) {
          actBc->masterDof = 0;
        } else {
          actBc->masterDof = actFeFunction->GetResultInfo()->GetDofIndex( masterDof );
        }
        if( slaveDof.empty() ) {
          actBc->slaveDof = 0;
        } else {
          actBc->slaveDof = actFeFunction->GetResultInfo()->GetDofIndex( masterDof );
        }

        // add definition
        actFeFunction->AddConstraint(actBc);
      } catch (Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create constraints on '"
                           << name << "'" );
      }
    }

    // =====================================================================
    // Periodic boundary conditions
    // =====================================================================

    // fetch paramnodes for constraint
    ParamNodeList prNodes = bcsNode->GetList("periodic");
    std::string masterName, slaveName;

    // iterate over all parameter nodes
    for( UInt i = 0; i < prNodes.GetSize(); i++ ) {
      try {
        prNodes[i]->GetValue( "master", masterName );
        prNodes[i]->GetValue( "slave", slaveName );
        prNodes[i]->GetValue( "dof", dof );
        prNodes[i]->GetValue( "quantity", resultName );

        // fetch related resultInfo object
        actFeFunction = GetFeFunction( SolutionTypeEnum.Parse(resultName));
        
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
            actBc->masterDof = 0;
          } else {
            actBc->masterDof = actFeFunction->GetResultInfo()->GetDofIndex( dof );
          }
          actBc->slaveDof = actBc->masterDof;
          actBc->result = actFeFunction->GetResultInfo();
          actBc->periodic = true;

          // add definition
          actFeFunction->AddConstraint(actBc);
        }
      } catch (CoupledField::Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create periodic boundary on '"
                           << name << "'" );
      }
    }
  }



  void SinglePDE::ReadRhsExcitation( const std::string& elemName, 
                                     const StdVector<std::string>& compNames,
                                     ResultInfo::EntryType type,
                                     bool isComplex,
                                     StdVector<shared_ptr<EntityList> >& entities, 
                                     StdVector<shared_ptr<CoefFunction> >& coef ) {
    
    // get grip of all elements of that type
    ParamNodeList elems = myParam_->Get("bcsAndLoads")->GetList(elemName);
    
    entities.Resize(elems.GetSize());
    coef.Resize(elems.GetSize());
    for( UInt i = 0; i < elems.GetSize(); ++i ) {

      PtrParamNode xml = elems[i];

      // get entity list, depending on type
      std::string entName = xml->Get("name")->As<std::string>();
      switch( ptgrid_->GetEntityType(entName) ) {
        case EntityList::NAMED_NODES:
          entities[i] = ptgrid_->GetEntityList( EntityList::NODE_LIST, 
                                                entName );
          break;
        case EntityList::REGION:
        case EntityList::NAMED_ELEMS:
          entities[i] = ptgrid_->GetEntityList( EntityList::ELEM_LIST, 
                                               entName );
          break;
        case EntityList::NO_TYPE:
          EXCEPTION("No entities with name '" << entName << "' known");
          break;
      }

     
      
      UInt numComp = compNames.GetSize();
      StdVector<std::string> vals(numComp), phases(numComp);
      vals.Init("0.0");
      phases.Init("0.0");

      // switch type of coef function
      // Note: In case someone request a "vector" valued result and
      // provides no dofNames, we use the scalar parameters.
      if( type == ResultInfo::SCALAR  ||
          (type == ResultInfo::VECTOR && compNames.GetSize() == 0 )  ) {
        // --------------
        //  S C A L A R   
        // --------------

        std::string val, phase;
        xml->GetValue("value", val);
        xml->GetValue("phase", phase);
        std::string real = AmplPhaseToReal(val, phase );

        if( type == ResultInfo::SCALAR) {
          // -- SCALAR case --
          if(!isComplex ) {
            coef[i] = CoefFunction::Generate(Global::REAL, real);  
          } else {
            std::string imag = AmplPhaseToReal(val, phase );
            coef[i] = CoefFunction::Generate(Global::COMPLEX, real, imag);
          } 
        }else {
          // -- VECTOR case --
          // generate coefficient function
          StdVector<std::string> realV, imagV;
          realV = real;
          imagV =  AmplPhaseToReal(val, phase );
          if(!isComplex) {
            StdVector<std::string> valV;
            valV = val;
            coef[i] = CoefFunction::Generate(Global::REAL, valV );
          } else {
            coef[i] = CoefFunction::Generate(Global::COMPLEX,
                                             realV, imagV );
          }
        }
        
      } else if (type == ResultInfo::VECTOR) {
        
        // --------------
        //  V E C T O R  
        // --------------
        // a) all values are given as vector
        if( xml->Has("values") ) {
          std::string valString = xml->Get("values")->As<std::string>();
          SplitStringList(valString, vals, ' ');

          // consistency check
          if( vals.GetSize() != numComp ) {
            EXCEPTION("Boundary condition needs " << numComp << " values!");
          }
          
          // check for phase vector (optional)
          if( xml->Has("phase")) {
            std::string phaseString = xml->Get("phase")->As<std::string>();
            SplitStringList(phaseString, phases, ' ');

            // consistency check
            if( phases.GetSize() != numComp ) {
              EXCEPTION("Boundary condition needs " << numComp 
                         << " phase values!");
            }
          }
        } 
        // b) values are given component-wise
        else if(xml->Has("comp")) {
          ParamNodeList compList = xml->GetList("comp");
          Integer index = 0;
          std::string dof, val, phase;
          for( UInt j = 0; j < compList.GetSize(); ++j  ) {
           compList[j]->GetValue("dof", dof);
           compList[j]->GetValue("value", val);
           compList[j]->GetValue("phase", phase);
           
           // find index
           if( compNames.GetSize() == 0 ) {
             index = 0;
           } else {
             index = compNames.Find(dof);
             if( index == -1 ) {
               EXCEPTION("Could not find component with name '" << dof << "'");
             }
           }
           
           vals[index] = val;
           phases[index] = val;
          } // loop components

        } else {
          EXCEPTION( "No values given for boundary condition '" 
              << elemName << "' on '" << entName << "'" );
        }
        
        // generate coefficient function
        StdVector<std::string> real, imag;
        AmplPhaseToRealImag(vals, phases, real, imag);
        if(!isComplex) {
          coef[i] = CoefFunction::Generate(Global::REAL, vals );
        } else {
          coef[i] = CoefFunction::Generate(Global::COMPLEX, real, imag );
        }
        
      } else if (type == ResultInfo::TENSOR ) {
        // --------------
        //  T E N S O R  
        // --------------
        EXCEPTION("Not yet implemented for tensor-valued boundary conditions");
        //    - all defined: read in both vectors
        //    - components: 

      }
      
      // obtain coordinate system and set it at coefficient function
      std::string coordSysId = "default";
      xml->GetValue("coordSysId", coordSysId, ParamNode::PASS);
      if( coordSysId != "default" ) {
        coef[i]->SetCoordinateSystem( domain->GetCoordSystem(coordSysId) );
      }
    } // for
    // return 
  }

  void SinglePDE::ReadMaterialData() {
    UInt i, numRegions;
    
    // get list of parameter nodes for region definitions
    ParamNodeList regionNodes;

    PtrParamNode regionListNode = param->
      Get("domain")->Get("regionList", ParamNode::PASS );
    if( regionListNode)
      regionNodes = regionListNode->GetList("region");

    numRegions = regionNodes.GetSize();

    // obtain pointer to materialHandler
    MaterialHandler* matLoader = domain->GetMaterialHandler();


    // -------------------
    // NORMAL MATERIALS
    // -------------------
    std::string region, material, composite, refCoordSys;

    // iterate over all regions
    for( i = 0; i < numRegions; ++i ) {

      try{
        // get data from node
        regionNodes[i]->GetValue( "name", region );
        regionNodes[i]->GetValue( "material", material );
        regionNodes[i]->GetValue( "coordSysId", refCoordSys );

        // get regionId
        RegionIdType actRegionId = ptgrid_->GetRegion().Parse( region );

        // if no material is set, continue with next loop run
        if( material.empty() )
          continue;

        // if region is not contained for current pde, simply continue
        // with next loop
        if( subdoms_.Find( actRegionId) < 0 )
          continue;

        // Read data
        materials_[actRegionId] = matLoader->LoadMaterial(material, pdematerialclass_);

        // log the just read material. LoadMaterial() so to say initializes the ToInfo()
        PtrParamNode in = infoNode_->GetByVal("material", "name", material);
        // additional regions are automatically appended
        in->Get("regionList")->GetByVal("region", "name", domain->GetGrid()->GetRegion().ToString(actRegionId));
        materials_[actRegionId]->ToInfo(in);

        // Check for local coordinate system
        if( !refCoordSys.empty() ) {
          CoordSystem * actCoosy =
            domain->GetCoordSystem( refCoordSys);
          materials_[actRegionId]->SetCoordSys( actCoosy );
        }

        // Check for material rotation parameters
        PtrParamNode rotNode = regionNodes[i]->Get("matRotation", ParamNode::PASS );

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
          rotVec[0] = rotNode->Get( "alpha" )->MathParse<Double>();
          rotVec[1] = rotNode->Get( "beta" )->MathParse<Double>(); 
          rotVec[2] = rotNode->Get( "gamma" )->MathParse<Double>();

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
        regionNodes[i]->GetValue( "name", region );
        regionNodes[i]->GetValue( "composite", composite );

        // get regionId
        RegionIdType actRegionId = ptgrid_->GetRegion().Parse( region );

        // if no composite is set, continue with next loop run
        if( composite == "" )
          continue;
        
        if( subdoms_.Find( actRegionId) < 0 )
          continue;

        PtrParamNode in = infoNode_->GetByVal("composite", "region", domain->GetGrid()->GetRegion().ToString(actRegionId));

        // get composite node
        PtrParamNode compNode;
        try {
          compNode = param->Get("domain")->GetByVal("composite", "name", composite);
        } catch( Exception& ex ) {
          RETHROW_EXCEPTION(ex, "No composite material defined with name '" << composite << "'");
        }

        // get laminaNodes
        ParamNodeList laminaNodes = compNode->GetList("lamina");

        // Create new lamina and fill ine materials and thicknesses
        Composite & myMat = compositeMaterials_[actRegionId];
        myMat.name = composite;

        // iterate over all single laminas
        for( UInt j = 0; j < laminaNodes.GetSize(); j++ ) {

          // fetch data for lamina
          std::string lamMaterial;
          Double lamThickness, lamOrientation;

          laminaNodes[j]->GetValue( "material", lamMaterial);
          laminaNodes[j]->GetValue( "thickness", lamThickness);
          laminaNodes[j]->GetValue( "orientation", lamOrientation);

          // Print information
          PtrParamNode lam = in->Get("lamina", ParamNode::APPEND);
          lam->Get("thickness")->SetValue(lamThickness);
          lam->Get("orientation")->SetValue(lamOrientation);

          myMat.thickness.Push_back( lamThickness );
          myMat.orientation.Push_back( lamOrientation );
          myMat.materials.Push_back( matLoader->LoadMaterial( lamMaterial, pdematerialclass_));

          PtrParamNode lm_ = lam->Get("material");
          lm_->Get("name")->SetValue(lamMaterial);
          myMat.materials.Last()->ToInfo(lm_);

        } // over single laminae
      } catch (Exception& ex ) {
        RETHROW_EXCEPTION( ex, "Could not create composite material '"
                           << composite << "'");
      }
    } // over composite

    // once again: loop over all regions and make sure that there is either
    // a normal material or a composite material
    numRegions = subdoms_.GetSize();
    std::map<RegionIdType, BaseMaterial*>::iterator matEnd = materials_.end();
    std::map<RegionIdType, Composite>::iterator
        compEnd = compositeMaterials_.end();
    
    for( i = 0; i < numRegions; ++i ) {
      RegionIdType actRegionId = subdoms_[i];
      if ((materials_.find(actRegionId) == matEnd)
          && (compositeMaterials_.find(actRegionId) == compEnd)) {
        region = ptgrid_->GetRegion().ToString(actRegionId);
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
  // COUPLING SECTION
  // ======================================================

  void SinglePDE::CalcInputCoupling() {


    std::string errMsg;
    StdVector<UInt> * nodes;
    SingleVector * val;
//    Integer eqnNr;
    // UInt couplingDof;

    // Determine maximal allowed equation number for algebraic system
//    Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // Outer loop over all INPUT coupling terms
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {

      //    ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);
      // couplingDof = ptCoupling_->GetInputDof(i);

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
        REFACTOR;
//        ptCoupling_->GetInputNodes(i, nodes);
//
//        for ( UInt dof = 0; dof < couplingDof; dof++ ) {
//          for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
//            eqnNr = eqnMap_->GetNodeEqn( (*nodes)[j], dof+1 );
//            if ( eqnNr != 0 && eqnNr <= maxAllowedEqn ) {
//              algsys_->SetNodeRHS( help[ dof + couplingDof * j ], pdeId_,
//                                   eqnNr );
//            }
//
//#ifndef NDEBUG
//            /*
//            else if ( eqnNr > maxAllowedEqn ) {
//              (*debug) << "SinglePDE::CalcInputCoupling: "
//                       << "(" << pdename_ << ") "
//                       << "Refused to pass "
//                       << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
//                       << "it execeeds numLastFreeDof = " << maxAllowedEqn
//                       << std::endl;
//            }
//            else if ( eqnNr == 0 ) {
//              (*debug) << "SinglePDE::CalcInputCoupling: "
//                       << "(" << pdename_ << ") "
//                       << "Refused to pass node " << (*nodes)[j] << "to SetNodeRHS(), since "
//                       << "it is fixed by hom. Dirichlet BC" << std::endl;
//            }
//            */
//#endif
//
//          }
//        }
        break;

        // -----------------------
        // InhomDirichlet COUPLING
        // -----------------------
      case ID_BC:
        REFACTOR;
//        // How do we want to treat inhomogeneous Dirichlet boundary
//        // conditions?
//        {
//          if ( usePenalty_ == false ) {
//            EXCEPTION( "Cannot use inhom. Dirichlet coupling together with "
//                       << "IDBC elimination, since the equation numbering "
//                       << "object does currently not have the information "
//                       << "required to put those values at the end of the "
//                       << "equation number interval! Please set idbcHandling="
//                       << '"' << "penalty" << '"' << " in your xml-file" );
//          }
//        }
//
//        // Set flag that the boundary conditions have to be incorporated
//        updateCouplingBCs_ = true;
//
//        ptCoupling_->GetInputNodes(i, nodes);
//
//        for ( UInt dof = 0; dof < ptCoupling_->GetInputDof(i); dof++ ) {
//          for ( UInt j = 0; j < nodes->GetSize(); j++) {
//
//            eqnNr = eqnMap_->GetNodeEqn( (*nodes)[j], dof+1 );
//
//            if (eqnNr==0) {
//              // In this case, we can not perform the coupling
//              // The best would be to inform the user ONCE, that some nodes
//              // can not couple. To accomplish this, we would need some way
//              // to remember the warning we have already issued
//              //EXCEPTION( "Coupling node " << (*nodes)[j] << " has a "
//              //           << "zero equation number, so no Dirichlet BC "
//              //           << "can be assigned" );
//        //        std::cerr << "node=" << *nodes)[j]
//            } else {
//            algsys_->SetDirichlet( pdeId_, eqnNr,
//                                   help[dof+j*couplingDof] );
//            
//            }
//          }
//        }
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


  void SinglePDE::WriteRestart()
  {
    EXCEPTION( "Restart-Functionality should be adapted to new"
                "FeSpace / FeFunction structure");
    // prepare output file
//    shared_ptr<SimOutput> restartOutFile;
//    Grid* grid = domain->GetGrid();
//    const std::string simName = progOpts->GetSimName();
//    std::string restartFileName = simName+"_"+pdename_+".restart";
//    PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
//    PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
//    eFiles->SetName("externalFiles");
//    eFiles->SetValue( "false" );
//    h5Node->AddChildNode(eFiles);
//    restartOutFile = shared_ptr<SimOutput>(new SimOutputHDF5(restartFileName, h5Node));
//    restartOutFile->Init(grid, true);
//
//    // time stepping variables
//    const Double& dt = TS_alg_->GetTimeStep();
//    Double timeTmp;
//    UInt lastTimeStep = domain->GetSingleDriver()->GetActStep(pdename_);
//    if (isComplex_)
//    {
//      EXCEPTION("restart file for harmonic results not implemented");
//    }
//
//    shared_ptr<EntityList> entList;
//    ResultMap::iterator it = resultLists_.begin();
//
//    // time stepping vectors in which the results will be stored to write out
//    SBM_Vector solutionTmp;
//    std::map<TIMEStepType, StdVector<shared_ptr<BaseResult> > > outResults;
//    std::map<DERIVType, StdVector<shared_ptr<BaseResult> > > outResults_deriv;
//    std::map<TIMEStepType, SBM_Vector > TsMap = TS_alg_->GetTimeStepMap();
//    /* the number of time steps needed for the algorithm
//     * INFO: TIMESTEP_0 is the new time step */
//    UInt numTimeSteps = TsMap.size();
//    if (TS_alg_->is_SolTimeStep_set(TIMESTEP_0))
//      --numTimeSteps;
//    if (numTimeSteps > lastTimeStep)
//    {
//      EXCEPTION("not enough time steps done to create a feasible restart");
//    }
//
//    // collect all necessary data
//    for (; it != resultLists_.end(); it++)
//    {
//      ResultList& actList = it->second;
//      // iterate over all solutions for each result type
//      for (UInt i = 0; i < actList.GetSize(); ++i)
//      {
//        shared_ptr<ResultInfo> actResultInfo = actList[i]->GetResultInfo();
//        if (TS_alg_->isDeriv(actResultInfo->resultType))
//        {
//          // do not get result if it is a derivative of another result
//          // derivative will be written later
//          continue;
//        }
//        entList = actList[i]->GetEntityList();
//
//        shared_ptr<BaseResult> outResult;
//        outResult = shared_ptr<BaseResult>(new Result<Double>());
//        outResult->SetResultInfo( actResultInfo );
//        outResult->SetEntityList( entList );
//        solutionTmp = TS_alg_->GetOld(TIMESTEP_1);
//        ExtractResult<Double>(outResult, solutionTmp);
//        if (outResults.find(TIMESTEP_1) == outResults.end())
//        {
//          StdVector<shared_ptr<BaseResult> > tmp;
//          outResults[TIMESTEP_1] = tmp;
//        }
//        outResults[TIMESTEP_1].Push_back(outResult);
//        restartOutFile->RegisterResult( outResult, \
//            lastTimeStep - numTimeSteps +1, 1, lastTimeStep, false );
//
//        /* go over all time step except TIMESTEP_1*/
//        std::map<TIMEStepType, SBM_Vector >::iterator itTs = TsMap.begin();
//        for (; itTs != TsMap.end(); ++itTs)
//        {
//          TIMEStepType timeStepType = itTs->first;
//          if (timeStepType == TIMESTEP_0 || timeStepType == TIMESTEP_1)
//          {
//            continue;
//          }
//          shared_ptr<BaseResult> outResult;
//          outResult = shared_ptr<BaseResult>(new Result<Double>());
//          outResult->SetResultInfo( actResultInfo );
//          outResult->SetEntityList( entList );
//
//          solutionTmp = TS_alg_->GetOld(timeStepType);
//          ExtractResult<Double>(outResult, solutionTmp);
//          if (outResults.find(timeStepType) == outResults.end())
//          {
//            StdVector<shared_ptr<BaseResult> > tmp;
//            outResults[timeStepType] = tmp;
//          }
//          outResults[timeStepType].Push_back(outResult);
//        }
//
//        /* go over all derivatives */
//        std::map<DERIVType, Vector<Double> >& DerivMap = TS_alg_->GetDeriveMap();
//        std::map<DERIVType, Vector<Double> >::iterator itDeriv = DerivMap.begin();
//        for (; itDeriv != DerivMap.end(); ++itDeriv)
//        {
//          DERIVType derivType = itDeriv->first;
//          shared_ptr<BaseResult> outResult_solDeriv;
//          outResult_solDeriv = shared_ptr<BaseResult>(new Result<Double>());
//          shared_ptr<ResultInfo> actResultInfo_deriv(new ResultInfo);
//          *actResultInfo_deriv = *actResultInfo;
//
//          outResult_solDeriv->SetResultInfo( actResultInfo_deriv );
//          outResult_solDeriv->SetEntityList( entList );
//          solutionTmp = TS_alg_->GetDeriv(derivType);
//          ExtractResult<Double>(outResult_solDeriv, solutionTmp);
//
//          const SolutionType& tmpSolType = \
//            TS_alg_->mapDerivToSolutionType(actResultInfo_deriv->resultType, derivType);
//          actResultInfo_deriv->resultType = tmpSolType;
//          actResultInfo_deriv->resultName = SolutionTypeEnum.ToString(tmpSolType);
//
//          outResult_solDeriv->SetResultInfo( actResultInfo_deriv );
//          if (outResults_deriv.find(derivType) == outResults_deriv.end())
//          {
//            StdVector<shared_ptr<BaseResult> > tmp;
//            outResults_deriv[derivType] = tmp;
//          }
//          outResults_deriv[derivType].Push_back(outResult_solDeriv);
//          restartOutFile->RegisterResult( outResult_solDeriv, \
//              lastTimeStep - numTimeSteps +1, 1, lastTimeStep, false );
//        }
//      }
//    }
//
//    // start the real storring of the files
//    restartOutFile->BeginMultiSequenceStep(1, analysistype_, lastTimeStep);
//    UInt timeStepTypeAsInt = numTimeSteps;
//    for (UInt i = lastTimeStep - timeStepTypeAsInt +1; i <= lastTimeStep; ++i)
//    {
//      timeTmp = solveStep_->GetActTime() -  (timeStepTypeAsInt -1) * dt;
//      restartOutFile->BeginStep(lastTimeStep - timeStepTypeAsInt +1, timeTmp);
//      writeOutTimeStep(restartOutFile, outResults[(TIMEStepType)timeStepTypeAsInt]);
//      if (timeStepTypeAsInt == TIMESTEP_1)
//      {
//        std::map<DERIVType, StdVector<shared_ptr<BaseResult> > >::iterator itDerivRes = \
//          outResults_deriv.begin();
//        for (; itDerivRes != outResults_deriv.end(); ++itDerivRes)
//        {
//          writeOutTimeStep(restartOutFile, itDerivRes->second);
//        }
//      }
//      restartOutFile->FinishStep();
//      timeStepTypeAsInt -= 1;
//    }
//
//    restartOutFile->FinishMultiSequenceStep();
//    restartOutFile->Finalize();
  }

  void SinglePDE::ReadRestart(UInt &startStep)
  {
    EXCEPTION( "Restart-Functionality should be adapted to new"
                "FeSpace / FeFunction structure");
//    /* initialisation for data reading */
//    std::string simName = progOpts->GetSimName();
//    std::string restartFileName = "results_hdf5/"+simName+"_"+pdename_+".restart.h5";
//    PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
//    PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
//    eFiles->SetName("externalFiles");
//    eFiles->SetValue( "false" );
//    h5Node->AddChildNode(eFiles);
//    shared_ptr<SimInput> input(new SimInputHDF5(restartFileName, h5Node));
//    // read in mesh of input
//    input->InitModule();
//    UInt dim = input->GetDim();
//    Grid* ptGrid = new GridCFS(dim);
//    input->ReadMesh(ptGrid);
//    // FinishInit() does not work if we have named nodes which are added by
//    // coords. This is a not so dirty work around
//    *ptGrid = *ptgrid_;
//    /* ptGrid->FinishInit(); */
//
//    std::map<UInt, BasePDE::AnalysisType> types;
//    std::map<UInt, UInt> numMultiSteps;
//    const bool isHistory = false;
//    input->GetNumMultiSequenceSteps( types, numMultiSteps, isHistory );
//    const UInt lastMultiStep = numMultiSteps.size();
//    if (lastMultiStep != 1)
//    {
//      EXCEPTION("Restart can not handle multiple multistep!")
//    }
//    const UInt lastTimeStep = numMultiSteps[lastMultiStep];
//    startStep = lastTimeStep;
//    std::map<TIMEStepType, Vector<Double> > TsMap = TS_alg_->GetTimeStepMap();
//    // collect all necessary data
//    Vector<Double> tmpVec;
//
//    StdVector<shared_ptr<BaseResult> > inResults;
//
//    ResultMap::iterator it = resultLists_.begin();
//    shared_ptr<EntityList> entList;
//    std::map<shared_ptr<ResultInfo>, std::map<UInt, Double> > resultSteps;
//    for (; it != resultLists_.end(); it++)
//    {
//      ResultList & actList = it->second;
//      // iterate over all solutions for each result type
//      for (UInt i = 0; i < actList.GetSize(); ++i)
//      {
//        shared_ptr<ResultInfo> actResultInfo = actList[i]->GetResultInfo();
//        if (TS_alg_->isDeriv(actResultInfo->resultType))
//        {
//          // do not get result if it is a derivative of another result
//          // derivative will be fetched later
//          continue;
//        }
//        input->GetStepValues( lastMultiStep, actResultInfo,
//            resultSteps[actResultInfo], false);
//
//        // iterate over all regions
//        StdVector<shared_ptr<EntityList> > regions;
//        input->GetResultEntities(lastMultiStep, actResultInfo, regions, isHistory);
//        for (UInt iRegion = 0; iRegion < regions.GetSize(); iRegion++)
//        {
//          // generate new result object and add it to output writer
//          shared_ptr<BaseResult > inResult;
//          if (types[1] != BasePDE::HARMONIC)
//          {
//            inResult  = shared_ptr<BaseResult>( new Result<Double>() );
//          } else {
//            EXCEPTION("restarting over harmonic results ist not implemented");
//          }
//          inResult->SetEntityList( regions[iRegion] );
//          inResult->SetResultInfo( actResultInfo );
//
//          /* read in all time steps */
//          std::map<TIMEStepType, Vector<Double> >::iterator itTs = TsMap.begin();
//          for (; itTs != TsMap.end(); ++itTs)
//          {
//            TIMEStepType timeStepType = itTs->first;
//            if (timeStepType == TIMESTEP_0)
//            {
//              continue;
//            }
//            input->GetResult(lastMultiStep, lastTimeStep - timeStepType +1, \
//                inResult, isHistory);
//            tmpVec = TS_alg_->GetOld(timeStepType);
//            InsertResult<Double>(tmpVec, inResult);
//            if (timeStepType == TIMESTEP_1)
//            {
//              sol_->SetAlgSysVector(tmpVec);
//              algsys_->InitSol(tmpVec);
//            }
//            TS_alg_->SetOld(tmpVec, timeStepType);
//          }
//
//          /* read in all derivatives steps */
//          std::map<DERIVType, Vector<Double> >& DerivMap = TS_alg_->GetDeriveMap();
//          std::map<DERIVType, Vector<Double> >::iterator itDeriv = DerivMap.begin();
//          for (; itDeriv != DerivMap.end(); ++itDeriv)
//          {
//            DERIVType derivType = itDeriv->first;
//            shared_ptr<ResultInfo> actResultInfo_deriv(new ResultInfo);
//            *actResultInfo_deriv = *actList[i]->GetResultInfo();
//            SolutionType tmpSolType = \
//              TS_alg_->mapDerivToSolutionType(actResultInfo_deriv->resultType, derivType);
//
//            actResultInfo_deriv->resultType = tmpSolType;
//            actResultInfo_deriv->resultName = SolutionTypeEnum.ToString(tmpSolType);
//
//            inResult->SetResultInfo( actResultInfo_deriv );
//            input->GetResult(lastMultiStep, lastTimeStep, inResult, isHistory);
//
//            /* revert the changes from actResultInfo_deriv */
//            inResult->SetResultInfo( actResultInfo);
//            tmpVec = TS_alg_->GetDeriv(derivType);
//            InsertResult<Double>(tmpVec, inResult);
//            TS_alg_->SetDeriv(tmpVec, derivType);
//          }
//        }
//      }
//    }
//    if (isIterCoupled_)
//    {
//      CalcOutputCoupling();
//      CalcInputCoupling();
//    }
  }

  void SinglePDE::GetMemento( shared_ptr<PDEMemento>& memento) {

    EXCEPTION("SinglePDE::GetMemento() not wroking yet");
//    // create new memento
//    shared_ptr<PDEMemento> myMemento (new PDEMemento() );
//
//    // first get memento of coupling object
//    if (isIterCoupled_) {
//      ptCoupling_->GetMemento(myMemento->couplingMemento_);
//      myMemento->isIterCoupled_ = true;
//    }
//
//    // then write own data to PDEMemento
//    myMemento->analysisType_ = analysistype_;
//    myMemento->gridFileName_ = progOpts->GetMeshFileStr();
//    myMemento->stepNum_ = domain->GetSingleDriver()->GetActStep(pdename_);
//
//    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {
//
//      // --- Real values --
//      Vector<Double> & solReal =
//        dynamic_cast<NodeStoreSol<Double>&>(*(sol_)).GetAlgSysVector();
//      UInt numDofs = results_[0]->dofNames.GetSize();
//      StdVector<Integer> eqns;
//      for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ){
//
//        // create new entity list
//        shared_ptr<NodeList> actSDList( new NodeList(ptgrid_ ) );
//        actSDList->SetNodesOfRegion( subdoms_[iRegion] );
//
//        // create new vector
//        Vector<Double> * values = new Vector<Double>;
//        Vector<Double> deriv1, deriv2, tn_1;
//
//        values->Resize( actSDList->GetSize() * numDofs );
//        values->Init();
//        if (analysistype_ == TRANSIENT ) {
//          deriv1.Resize( actSDList->GetSize() * numDofs );
//          deriv2.Resize( actSDList->GetSize() * numDofs );
//            tn_1.Resize( actSDList->GetSize() * numDofs );
//        }
//
//        EntityIterator it = actSDList->GetIterator();
//        UInt pos = 0;
//        for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
//          eqnMap_->GetEqns( eqns, *results_[0], it );
//          for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
//            if ( eqns[iDof] != 0 ) {
//              (*values)[numDofs*pos+iDof] = solReal[abs(eqns[iDof])-1];
//              if( analysistype_ == TRANSIENT ) {
//                deriv1[numDofs*pos+iDof] =
//                  TS_alg_->GetDeriv(FIRST_DERIV)[(abs(eqns[iDof])-1)];
//                deriv2[numDofs*pos+iDof] =
//                    TS_alg_->GetDeriv(SECOND_DERIV)[(abs(eqns[iDof]-1))];
//                  tn_1[numDofs*pos+iDof] =
//                    TS_alg_->GetOld(TIMESTEP_1)[(abs(eqns[iDof]-1))];
//              }
//            } else {
//              (*values)[numDofs*pos+iDof] = 0.0;
//              if ( analysistype_ == TRANSIENT ) {
//                deriv1[numDofs*pos+iDof] = 0.0;
//                deriv2[numDofs*pos+iDof] = 0.0;
//                  tn_1[numDofs*pos+iDof] = 0.0;
//              }
//            }
//          }
//        }
//
//        // pass vector to memento object
//        std::string regionName = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
//        myMemento->solution_[regionName] = values;
//        if ( analysistype_ == TRANSIENT ) {
//          myMemento->solDeriv1_[regionName] = deriv1;
//          myMemento->solDeriv2_[regionName] = deriv2;
//          myMemento->sol_tn_1_[regionName] = tn_1;
//        }
//      }
//    } else {
//
//      // --- Complex values --
//
//      // store current frequency
//      myMemento->freq_ = dynamic_cast<HarmonicDriver*>(domain->GetSingleDriver())
//        ->GetActFreq();
//
//      Vector<Complex> & solComp =
//        dynamic_cast<NodeStoreSol<Complex>&>(*(sol_)).GetAlgSysVector();
//      UInt numDofs = results_[0]->dofNames.GetSize();
//      StdVector<Integer> eqns;
//      for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ){
//
//        // create new entity list
//        shared_ptr<NodeList> actSDList( new NodeList(ptgrid_ ) );
//        actSDList->SetNodesOfRegion( subdoms_[iRegion] );
//
//        // create new vector
//        Vector<Complex> * values = new Vector<Complex>;
//        values->Resize( actSDList->GetSize() * numDofs );
//        values->Init();
//        EntityIterator it = actSDList->GetIterator();
//        UInt pos = 0;
//        for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
//          eqnMap_->GetEqns( eqns, *results_[0], it );
//          for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
//            if ( eqns[iDof] != 0 ) {
//              (*values)[numDofs*pos+iDof] = solComp[abs(eqns[iDof])-1];
//            } else {
//              (*values)[numDofs*pos+iDof] = 0.0;
//            }
//          }
//        }
//
//        // pass vector to memento object
//        std::string regionName = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
//        myMemento->solution_[regionName] = values;
//      }
//    }
//
//    // now memento is initialized
//    myMemento->isSet_ = true;
//
//    // pass back memento
//    memento = myMemento;
  }


  void SinglePDE::SetMemento( shared_ptr<PDEMemento>& memento,
                              bool mementoAsDirichlet ) {
    EXCEPTION("SinglePDE::SetMemento not working yet");
//
//    if( isInitialized_ == true ) {
//      EXCEPTION( "SetMemento may only be called, if the method "
//                 << "SinglePDE::Init() was not called yet!" );
//    }
//
//    memento_ = memento;
//    mementoAsDirichlet_ = mementoAsDirichlet;
  }

  void SinglePDE::IncorporateMemento( ) {
    EXCEPTION("SinglePDE::IncorporateMemento not working yet");
//    // if there is no memento present -> leave
//    if( !memento_ ) {
//      return;
//    }
//
//   // if there is no information in the memento just leave
//    if ( memento_->isSet_ == false ) {
//      return;
//    }
//
//    // check that memento is bases on the same grid file
//    if( memento_->gridFileName_ != progOpts->GetMeshFileStr() ) {
//      EXCEPTION( "Error in reading in memento: Memento is based on grid '"
//                 << memento_->gridFileName_
//                 << ", whereas the current simulation is based on grid '"
//                 << progOpts->GetMeshFileStr() );
//    }
//
//    UInt numDofs = results_[0]->dofNames.GetSize();
//    if ( analysistype_ == STATIC 
//        || analysistype_ == TRANSIENT
//        || analysistype_ == EIGENFREQUENCY ) {
//
//      // convert solution to transient StoreSolution type
//      Vector<Double> & solVec =
//        (dynamic_cast<NodeStoreSol<Double> &>(*sol_)).GetAlgSysVector();
//
//
//      Vector<Double> solDeriv1, solDeriv2, solTn_1;
//      
//     // we need derivatives only if we have transient analysis
//      if( analysistype_ == TRANSIENT ) {
//        if( TS_alg_->GetDeriv(FIRST_DERIV).GetSize() != 0 ) {
//          solDeriv1 = TS_alg_->GetDeriv(FIRST_DERIV);
//        } else {
//          solDeriv1.Resize( solVec.GetSize() );
//          solDeriv1.Init();
//        }
//
//        if( TS_alg_->GetDeriv(SECOND_DERIV).GetSize() != 0 ) {
//          solDeriv2 = TS_alg_->GetDeriv(SECOND_DERIV);
//        } else {
//          solDeriv2.Resize( solVec.GetSize() );
//          solDeriv2.Init();
//        }
//
//        if( TS_alg_->is_SolTimeStep_set(TIMESTEP_1) ) {
//          solTn_1 = TS_alg_->GetOld(TIMESTEP_1);
//        } else {
//          solTn_1.Resize( solVec.GetSize() );
//          solTn_1.Init();
//        }
//      }
//
//      
//      if ( memento_->analysisType_ == HARMONIC ) {
//
//        // -> perform complex-to-real adjustment
//
//        // frequency of memento
//        Double freq = memento_->freq_;
//
//        // Iterate over all regions
//        for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ) {
//
//          // Check for related region in memento object
//          std::string name = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
//          if( memento_->solution_.find( name) != memento_->solution_.end() ) {
//
//            // get grip of vector and derivatives of memento
//            Vector<Complex> const & sol =
//              dynamic_cast<const Vector<Complex>& >(*(memento_->solution_[name]) );
//
//            // create entitylist
//            shared_ptr<NodeList> nodes (new NodeList(ptgrid_));
//            nodes->SetNodesOfRegion( subdoms_[iRegion] );
//
//            // iterate over all entries
//            EntityIterator it = nodes->GetIterator();
//            StdVector<Integer> eqns;
//            UInt pos = 0;
//            for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
//              eqnMap_->GetEqns( eqns, *results_[0], it );
//              for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
//                UInt actPos = numDofs*pos+iDof;
//                if ( eqns[iDof] != 0
//                     && (unsigned int) (std::abs(eqns[iDof])-1) < solVec.GetSize() ) {
//                  solVec[eqns[iDof]-1] = sol[actPos].real();
//                  if ( analysistype_ == TRANSIENT ) {
//                    solDeriv1[eqns[iDof]-1] = -2*PI*freq* sol[actPos].imag();
//                    solDeriv2[eqns[iDof]-1] =
//                      -4*PI*PI*freq*freq*sol[actPos].real();
//                  }
//                }
//              }
//            }
//          }
//        }
//
//        // pass derivatives to timestepping algorithm
//        if( analysistype_ == TRANSIENT ) {
//          TS_alg_->SetDeriv( solDeriv1, FIRST_DERIV );
//          TS_alg_->SetDeriv( solDeriv2, SECOND_DERIV );
//        }
//
//      } else {
//
//        // check if derivatives are also needed
//        bool needDeriv = ( analysistype_ == TRANSIENT ) &&
//          (memento_->analysisType_ == TRANSIENT);
//
//        // Iterate over all regions
//        for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ) {
//
//          // Check for related region in memento object
//          std::string name = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
//          if( memento_->solution_.find( name) != memento_->solution_.end() ) {
//
//            // get grip of vector and derivatives of memento
//            Vector<Double> const & sol =
//              dynamic_cast<const Vector<Double>& >(*(memento_->solution_[name]) );
//            Vector<Double> const & deriv1 = memento_->solDeriv1_[name];
//            Vector<Double> const & deriv2 = memento_->solDeriv2_[name];
//            Vector<Double> const & tn_1 = memento_->sol_tn_1_[name];
//
//            // create entitylist
//            shared_ptr<NodeList> nodes (new NodeList(ptgrid_));
//            nodes->SetNodesOfRegion( subdoms_[iRegion] );
//
//            // iterate over all entries
//            EntityIterator it = nodes->GetIterator();
//            StdVector<Integer> eqns;
//            UInt pos = 0;
//            for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
//              eqnMap_->GetEqns( eqns, *results_[0], it );
//              for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
//                if ( eqns[iDof] != 0
//                     && (unsigned int) (std::abs(eqns[iDof])-1) < solVec.GetSize() ) {
//                  solVec[eqns[iDof]-1] = sol[numDofs*pos+iDof];
//                  if ( needDeriv ) {
//                    solDeriv1[eqns[iDof]-1] = deriv1[numDofs*pos+iDof];
//                    solDeriv2[eqns[iDof]-1] = deriv2[numDofs*pos+iDof];
//                      solTn_1[eqns[iDof]-1] = tn_1[numDofs*pos+iDof];
//                  }
//                }
//              }
//            }
//          }
//        }
//
//        // pass derivatives to timestepping algorithm
//        if( needDeriv ) {
//          TS_alg_->SetDeriv( solDeriv1, FIRST_DERIV );
//          TS_alg_->SetDeriv( solDeriv2, SECOND_DERIV );
//          TS_alg_->SetOld( solTn_1, TIMESTEP_1 );
//        }
//      }
//    } else if ( analysistype_ == HARMONIC ) {
//
////      // check value-usage type
////      if( mementoAsDirichlet_ != true) {
////        
////        // nothing to do
////        //EXCEPTION( "For an harmonic simulation only the usage "
////        //           << "of a memento as Drichlet values makes sense!" );
////        return;
////      }
//      
//      if (!mementoAsDirichlet_) {
//        if ( memento_->analysisType_ == STATIC ) {
//          // convert solution to transient StoreSolution type
//          Vector<Complex> & solVec = 
//              (dynamic_cast<NodeStoreSol<Complex> &>(*sol_)).GetAlgSysVector();
//
//          // Iterate over all regions
//          for( UInt iRegion = 0; iRegion < subdoms_.GetSize(); iRegion++ ) {
//
//            // Check for related region in memento object
//            std::string name = ptgrid_->GetRegion().ToString( subdoms_[iRegion] );
//            if( memento_->solution_.find( name) != memento_->solution_.end() ) {
//
//              // get grip of vector and derivatives of memento
//              Vector<Double> const & sol = 
//                  dynamic_cast<const Vector<Double>& >(*(memento_->solution_[name]) );
//              // create entitylist
//           shared_ptr<NodeList> nodes (new NodeList(ptgrid_));
//           nodes->SetNodesOfRegion( subdoms_[iRegion] );
//
//           // iterate over all entries
//           EntityIterator it = nodes->GetIterator();
//           StdVector<Integer> eqns;
//           UInt pos = 0;
//              for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
//                eqnMap_->GetEqns( eqns, *results_[0], it );
//                for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
//                  if ( eqns[iDof] != 0
//                      && (unsigned int) (std::abs(eqns[iDof])-1) < solVec.GetSize() ) {
//                    solVec[eqns[iDof]-1] = Complex(sol[numDofs*pos+iDof],0);
//                  } // if
//                } // loop over dofs
//              } // loop over nodes
//            } // if region found
//          } // loop over regions
//        } // if static
//      }
//      if( mementoAsDirichlet_ == true) {
//        // iterate over all regions of pde
//        for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
//
//          std::string regionName = ptgrid_->GetRegion().ToString( subdoms_[i] );
//
//          // try to find related region in memento object
//          if( memento_->solution_.find( regionName ) !=
//              memento_->solution_.end() ) {
//
//            Vector<Complex> const & regionSol =
//                dynamic_cast<Vector<Complex>&>(*memento_->solution_[regionName]);
//
//            // create entitylist
//            shared_ptr<NodeList> nodes (new NodeList(ptgrid_));
//            nodes->SetNodesOfRegion( subdoms_[i] );
//
//            // iterate over all entries
//            EntityIterator it = nodes->GetIterator();
//            UInt pos = 0;
//            for( it.Begin(); !it.IsEnd(); it++, pos++ ) {
//
//              for( UInt iDof = 0; iDof < numDofs; iDof++ ) {
//                Complex val = regionSol[pos * numDofs + iDof];
//
//                // create idbc-condition and append to class container
//                shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
//                shared_ptr<NodeList> actList (new NodeList(ptgrid_) );
//                StdVector<UInt> nodeList(1);
//                nodeList[0] = it.GetNode();
//                actList->SetNodes(nodeList);
//
//                actBc->entities = actList;
//                actBc->result = results_[0];
//                actBc->eqnMap = eqnMap_;
//                actBc->dof = iDof+1;
//                actBc->value = lexical_cast<std::string>(std::abs( val ) );
//                actBc->phase = lexical_cast<std::string>(std::atan2( val.imag(), 
//                                                                     val.real())
//                *180/PI );
//
//                // append idbc at end of list
//                idBcs_.Push_back( actBc );
//              }
//            }
//          }  
//        } // loop over regions
//      }// memento as Dirichlet
//    } // HARMONIC Analysis

  }

  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
    // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void SinglePDE::GetPMLLayerData(Matrix<Double>& inner,
                                  Matrix<Double>& outer,
                                  RegionIdType actRegion,
                                  std::string& coordSysId )  {

    // fetch coordinate system object
    CoordSystem * coordSys = domain->GetCoordSystem(coordSysId);
    
    // only determine inner region, if not specified by hand using the
    // <propRegion> sub-element
    if ( inner.GetNumCols() != dim_ ) {

      //we have to compute it, since the user has not specified it
      inner.Resize(2,dim_);
      
      // The layout of inner is as follows
      // min( x_min y_min z_min )
      // max( x_max y_max z_max )
      
      // Set maximum values to minus large value and the minimum to large value
      Double largeVal = 1e33;
      for(UInt i = 0; i < dim_; i++ ) {
        inner[0][i] =   largeVal;
        inner[1][i] = - largeVal;
      }

      // Determine list of propagation regions (= all region except PML regions)
      StdVector<RegionIdType> propRegions;
      bool hasRealPropRegion = false;
      for (UInt isd = 0; isd < subdoms_.GetSize(); isd++) {
        RegionIdType aRegion = subdoms_[isd];
        if ( aRegion != actRegion 
         && dampingList_[aRegion] != PML ) {
          hasRealPropRegion = true;
          propRegions.Push_back( subdoms_[isd]);
        }
      }
      // If we have no "real" propagation region, i.e. only one or several
      // PML regions, we can not uniquely determine the damping parameters.
      // In this case the user needs to define the propagation region by hand,
      // so we issue an exception;
      // we can not uniquely 
      if( !hasRealPropRegion ){
        EXCEPTION("The " << pdename_ << "  has only PML damped regions. \nThus the "
            << "automatic calculation of damping parameters does not work!\n"
            << "Use the <propRegion> element instead do define the propagation "
            << "region explicitly!");
      }

      for (UInt isd = 0; isd < propRegions.GetSize(); isd++) {
        StdVector<Elem*> elemssd;
        ptgrid_->GetElems(elemssd, propRegions[isd] );

        for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
          StdVector<UInt> & connecth = elemssd[actEl]->connect;

          Matrix<Double> ptCoord;
          ptgrid_->GetElemNodesCoord(ptCoord, connecth,  false );

          Vector<Double> globCoord(dim_), locCoord(dim_);
          // loop over nodes
          for (UInt i=0; i< ptCoord.GetNumCols(); i++) {

            // convert global coordinate to local coordinate
            for( UInt iDim = 0; iDim < dim_ ; ++iDim )  {
              globCoord[iDim] = ptCoord[iDim][i];
            }
            coordSys->Global2LocalCoord(locCoord, globCoord);

            // determine min / max of propagation region
            for( UInt iDim = 0; iDim < dim_ ; ++iDim )  {
              if ( locCoord[iDim] < inner[0][iDim] )
                inner[0][iDim] = locCoord[iDim];
              if ( locCoord[iDim] > inner[1][iDim] )
                inner[1][iDim] = locCoord[iDim];

            }
          } // loop over nodes
        } // loop over elements
      } // loop over propagation regions
    } // no propagation region specified 

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
      Vector<Double> globCoord(dim_), locCoord(dim_);
      for (UInt i=0; i< ptCoord.GetNumCols(); i++) {
        
        // convert global coordinate to local coordinate
        for( UInt iDim = 0; iDim < dim_ ; ++iDim )  {
          globCoord[iDim] = ptCoord[iDim][i];
        }
        coordSys->Global2LocalCoord(locCoord, globCoord);
        
        // determine extension of PML region
        for( UInt iDim = 0; iDim < dim_ ; ++iDim )  {
          
          if ( locCoord[iDim] < outer[0][iDim] )
            outer[0][iDim] = locCoord[iDim];
         
          if ( locCoord[iDim] > outer[1][iDim] )
            outer[1][iDim] = locCoord[iDim];
        }
      }
    }
    
    // Echo information about PML to info file
    PtrParamNode base = infoNode_->Get(ParamNode::HEADER)
            ->Get("damping")->Get("pml");
    base->Get("coordSysId")->SetValue(coordSysId);
    PtrParamNode prop = base->Get("propRegion");
    PtrParamNode damp = base->Get("dampingRegion");
    std::string comp;
    for( UInt i = 0; i < dim_; ++i ) {
      comp = coordSys->GetDofName(i+1);
      prop->Get(comp)->Get("min")->SetValue(inner[0][i]);
      prop->Get(comp)->Get("max")->SetValue(inner[1][i]);
      damp->Get(comp)->Get("min")->SetValue(outer[0][i]);
      damp->Get(comp)->Get("max")->SetValue(outer[1][i]);
    }
  }


  void SinglePDE::DefineFeFunctions(){
      //This is the default creation of spaces
      //idee: die PDE gibt zum attribute formulation die passenden space zurück
      //DOGMA: PRO UNBEKANNTE EINE FUNCTION UND EIN SPACE
      std::string formulation;
      myParam_->GetValue("feSpaceFormulation",formulation,ParamNode::EX);
      infoNode_->SetComment("List of defined Spaces");
      PtrParamNode feSpaceNode = infoNode_->Get("feSpaces");
      std::map<SolutionType, shared_ptr<FeSpace> > spaces = 
          CreateFeSpaces(formulation, feSpaceNode);

      //loop over all spaces and set an FeFunction
      std::map<SolutionType, shared_ptr<FeSpace> >::iterator spIt = spaces.begin();
      while(spIt != spaces.end()){

        if(feFunctions_.find(spIt->first) != feFunctions_.end()){
          EXCEPTION("It seems that the PDE has created multiple spaces for one result: " << \
                    spIt->first << " This is not how its ought to be!");
        }

        if(analysistype_ == HARMONIC){
          feFunctions_[spIt->first] = shared_ptr<BaseFeFunction >(new FeFunction<Complex>());
          rhsFeFunctions_[spIt->first] = shared_ptr<BaseFeFunction >(new FeFunction<Complex>());
        }else{
          feFunctions_[spIt->first] = shared_ptr<BaseFeFunction >(new FeFunction<Double>());
          rhsFeFunctions_[spIt->first] = shared_ptr<BaseFeFunction >(new FeFunction<Double>());
         
          // Note: in the transient case, we also need fefunctions for the time derivatives
          // ... todo: add initialization
        }
        //let the objects know about each other
        spIt->second->AddFeFunction(feFunctions_[spIt->first]);
        feFunctions_[spIt->first]->SetFeSpace(spIt->second);
        feFunctions_[spIt->first]->SetPDE(this);
        feFunctions_[spIt->first]->SetGrid(ptgrid_);
        
        rhsFeFunctions_[spIt->first]->SetFeSpace(spIt->second);
        rhsFeFunctions_[spIt->first]->SetPDE(this);
        rhsFeFunctions_[spIt->first]->SetGrid(ptgrid_);
        spIt++;
      }
    }


 


  void SinglePDE::writeOutTimeStep(shared_ptr<SimOutput>& outFile, \
      StdVector<shared_ptr<BaseResult> >& outResults)
  {
    for( UInt iRes = 0; iRes < outResults.GetSize(); iRes++)
    {
      try {
        outFile->AddResult(outResults[iRes]);
      } catch (Exception& ex ) {
        std::cerr <<  "\nResult '" << outResults[iRes]->GetResultInfo()->resultName
          << "' in MsStep: " << 1 << ", step: " << 1
          << " could not be converted:\n\n";
        std::cerr << ex.what() << std::endl;
      }
    }
  }
} // end of namespace

