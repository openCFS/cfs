// We put include of CoefFunctionScatteredData.hh and intrin.hh
// in first place to prevent the following error:
// boost/thread/win32/interlocked_read.hpp:61:20: error:
//   '_InterlockedCompareExchange' is not a member of 'boost::detail'
#ifdef WIN32
#include <direct.h>
#endif
#include "Domain/CoefFunction/CoefFunctionScatteredData.hh"
#include "Domain/CoefFunction/CoefFunctionFileData.hh"
#include "Domain/CoefFunction/CoefFunctionFileDataMeas.hh"
#include "Domain/CoefFunction/CoefFunctionPython.hh"

#include "PDE/SinglePDE.hh"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <fstream>
#include <boost/algorithm/string.hpp>

// for coordinate handling
#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

#include "Utils/EvalIntegrals/BiotSavart.hh"

#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"

// header for saving / retrieving the simulation state
#include "DataInOut/SimState.hh"
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

// header for logging
#include "DataInOut/Logging/LogConfigurator.hh"

// header for Materialhandling
#include "DataInOut/ParamHandling/MaterialHandler.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

// header for Solvestep and assemble
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/Assemble.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TransientDriver.hh"
#include "Driver/HarmonicDriver.hh"
#include "Driver/MultiHarmonicDriver.hh"
#include "Driver/Harmonic25DDriver.hh"

// header for resultHandling
#include "DataInOut/ResultHandler.hh"

#include "DataInOut/PostProc.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "CoupledPDE/IterCoupledPDE.hh"

//feSpaces
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/H1/FeSpaceH1Hi.hh"
using std::string;

//coefFunctions
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionGrid.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionDummy.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

// used by Mortar coupling
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"

// used for getting coils from the pde
#include "MagEdgePDE.hh"
#include "MagEdgeSpecialAVPDE.hh"
#include "MagneticPDE.hh"


namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(singlepde, "singlePde")

  SinglePDE::SinglePDE( Grid *aptgrid, PtrParamNode paramNode,
                        PtrParamNode infoNode,
                        shared_ptr<SimState> simState, Domain* domain) :
    StdPDE( aptgrid, paramNode, infoNode, simState, domain ),
    isDirectCoupled_(false),
    isInitialized_(false),
    iterCplPde_(NULL),
    updatedGeo_(false),
    isMaterialComplex_( false ),
    isMultHarm_(false)
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

    // delete materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      if (it->second != NULL) delete it->second;
    }
    materials_.clear();
    
    // Loop overall feFunctions and finalize them
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator feIt;
    for( feIt = feFunctions_.begin(); feIt != feFunctions_.end(); ++feIt)  {
      feIt->second->CleanUp();
    }
    for( feIt = prevFeFunctions_.begin(); feIt != prevFeFunctions_.end(); ++feIt)  {
      feIt->second->CleanUp();
    }
    for( feIt = rhsFeFunctions_.begin(); feIt != rhsFeFunctions_.end(); ++feIt)  {
      feIt->second->CleanUp();
    }

    for( feIt = timeDerivFeFunctions_.begin(); feIt != timeDerivFeFunctions_.end(); ++feIt)  {
      feIt->second->CleanUp();
    }

    feFunctions_.clear();
    prevFeFunctions_.clear();
    rhsFeFunctions_.clear();
    timeDerivFeFunctions_.clear();
    
    // delete external domains
    std::map<shared_ptr<SimState>, Domain* >::iterator inputIt = inputs_.begin();
    for( ; inputIt != inputs_.end(); ++inputIt) {
      inputIt->first->Finalize();
      delete inputIt->second;
    }
  }


  std::string SinglePDE::ToString() const
  {
    std::stringstream ss;
    ss << pdename_ << " s=" << sequenceStep_ << " at=" << BasePDE::analysisType.ToString(analysistype_);
    return ss.str();
  }

  // ********
  //   Init
  // ********
  void SinglePDE::Init_Stage1( UInt sequenceStep, PtrParamNode base ) {

    sequenceStep_ = sequenceStep;

    infoNode_ = base == NULL ? myInfo_->Get("PDE")->Get(pdename_) : base->Get(pdename_);

    LOG_TRACE(singlepde) << pdename_ << ": Starting Initialization";


    // =====================================================================
    // Get type of analysis
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Obtaining analysis type";
    analysistype_ = domain_->GetSingleDriver()->GetAnalysisType();

    if(analysistype_ == MULTIHARMONIC) isMultHarm_ = true;

    isComplex_ = IsComplex();


    // =====================================================================
    // get regions/subdomains for PDE
    // =====================================================================

    LOG_TRACE(singlepde) << pdename_ << ": Obtaining regions";


    // Obtain regions the pde is defined on
    ParamNodeList regionNodes = myParam_->Get("regionList")->GetList("region");

    // output to info-file
    PtrParamNode list = infoNode_->Get(ParamNode::HEADER);

    // output and set regions_
    for( UInt i = 0; i < regionNodes.GetSize(); i++ )
    {
      PtrParamNode in_ = list->Get("region");
      std::string name = regionNodes[i]->Get("name")->As<std::string>();
      in_->Get("name")->SetValue(name);
      bool complexMat = false;
      regionNodes[i]->GetValue("complexMaterial",  complexMat, ParamNode::PASS );
      
      RegionIdType actRegionId = ptGrid_->GetRegion().Parse(name);
      
      // Check, if region was already defined an issue a warning otherwise
      if( std::find(regions_.Begin(), regions_.End(), actRegionId )!= regions_.End() )
        WARN( "The region '" << regionNodes[i]->Get("name")->As<std::string>()
              << "' was already defined for PDE '" << pdename_ 
              << "'. Please remove duplicate entries." );
          
      regions_.Push_back(actRegionId);

      complexMatData_[actRegionId] = complexMat;
    }


    // Generate a fitting algebraic system only if PDE is NOT
    // direct coupled
    if( needsAlgsys_ == true || !simState_->HasInput()) {
      if ( isDirectCoupled_ == false) {
        olasInfo_ = myInfo_->Get("OLAS")->Get(pdename_);
        algsys_ = new AlgebraicSys(olasNode_, olasInfo_, isComplex_, isMultHarm_);
        solStrat_ = algsys_->GetSolStrategy();
      }
    }

    // =====================================================================
    // inform algsys about possible multiharmonics
    // =====================================================================
    if( analysistype_ == MULTIHARMONIC){
      UInt baseFreq = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->baseFreq_;
      UInt numHarm_N = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->numHarmonics_N_;
      UInt numHarm_M = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->numHarmonics_M_;
      UInt numFFT = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->numFFT_;
      bool fullSystem = dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->fullSystem_;
      solStrat_->SetMultHarm(baseFreq, numHarm_N, numHarm_M, numFFT, fullSystem);
    }


    // =====================================================================
    // create assemble class
    // =====================================================================

    // Create new assemble class with according analysistype
    if( isDirectCoupled_ == false && needsAlgsys_ == true) {
      assemble_ = new Assemble( algsys_, analysistype_, this->mp_, myInfo_ );
    }
    
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
    //  read in non-conforming interfaces
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Reading non-conforming interfaces";
    ReadNcInterfaces();

    //======================================================================
    // trigger the creation of functionDescriptors
    //======================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Define FE-Functions";
//    LOG_DBG(singlepde) << "IS1: has MECH_DISPLACEMENT fefunciton? " << (feFunctions_.find(MECH_DISPLACEMENT) != feFunctions_.end());

    DefineFeFunctions();
    
    // Register all fe functions with the algebraic system
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      shared_ptr<BaseFeFunction> actFct = fncIt->second;
      shared_ptr<FeSpace> actSpace = fncIt->second->GetFeSpace();
      std::string fctName = SolutionTypeEnum.ToString(fncIt->first);
      FeFctIdType fctId;
      if( algsys_) {
        fctId = algsys_->ObtainFctId( fctName );
      } else {
        fctId = NO_FCT_ID;
      }
      actFct->SetFctId(fctId);
      fncIt++;
    }
   
    // =====================================================================
    // trigger definition of available results
    // =====================================================================
    DefinePrimaryResults();
    
    // =====================================================================
    // trigger definition of auxiliary results for Mortar ncInterfaces
    // =====================================================================
    DefineNcAuxResults();
    
    // =====================================================================
    // read in NonLinearities
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Initializing non-linearities";
    InitNonLin();

    // =====================================================================
    // read in material dependencies
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Initializing material depenencies";
    InitMaterialDependencies();
   
    // Todo: Move this part to the definition of damping
    PtrParamNode in = infoNode_->Get(ParamNode::HEADER);
    for(UInt i = 0; i < regions_.GetSize(); i++ )
    {
      PtrParamNode in_ = in->GetByVal("region", "name", ptGrid_->GetRegion().ToString(regions_[i]));

      // Not needed at the moment. Commented out due to gcc 4.6.

      std::string fuck_e2s;
      Enum2String(dampingList_[regions_[i]], fuck_e2s);
      in_->Get("damping")->SetValue(fuck_e2s);
    }

    // =====================================================================
    // Create time stepping algorithm
    // =====================================================================
    if ( analysistype_ == TRANSIENT ) {
      InitTimeStepping();
    }

    // =====================================================================
    // trigger definition of available postprocessing results
    // =====================================================================
    DefinePostProcResults();

    // Proceed with initialization stage 2 in the un-coupled case
  }
  
    void SinglePDE::Init_Stage2() {
    
    // =====================================================================
    // read in boundary conditions
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Reading boundary conditions";
    ReadBCs();
    ReadSpecialBCs();

    // =====================================================================
    // define the integrators for PDE and initialize eqn object
    // =====================================================================

    // Call initialization of (bi)linear integrators
    LOG_TRACE(singlepde) << pdename_ << ": Defining integrators";
    DefineIntegrators();
    if ( ncInterfaces_.GetSize() > 0 ) {
      DefineNcIntegrators();
    }
    DefineSurfaceIntegrators();
    DefineRhsLoadIntegrators();

    // Print information about defined integrators
    if( needsAlgsys_ == true && !isDirectCoupled_ )
      assemble_->ToInfo(infoNode_->Get(ParamNode::HEADER)->Get("integrators"));
  }

  void SinglePDE::Init_Stage3() {
    // =====================================================================
    //  map equations (FeSpaces) and finalize FeFunction (vector creation)
    // =====================================================================
    LOG_TRACE(singlepde) << "IS3: " << pdename_ << ": Mapping Equations";
    // Finalize spaces and fefunctions
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      shared_ptr<BaseFeFunction> actFct = fncIt->second;
      shared_ptr<FeSpace> actSpace = fncIt->second->GetFeSpace();
      actFct->SetGrid(ptGrid_);
      actSpace->Finalize();
      actSpace->PreCalcShapeFncs();

      // finalize feFunctions
      actFct->Finalize();
      
      // register FeFunctions with SimState class
      simState_->RegisterFeFct( actFct );
      
      // pass regions of primary function also RHS one
      StdVector< shared_ptr<EntityList> > support =  actFct->GetEntityList();
      LOG_DBG(singlepde) << "IS3: support=" << support.GetSize();
      for( UInt i = 0; i < support.GetSize(); ++i ) {
        LOG_DBG3(singlepde) << "IS3: support[" << i << "]=" << support[i]->GetName() << " size=" << support[i]->GetSize();
        rhsFeFunctions_[fncIt->first]->AddEntityList( support[i] );
      }

      // Pass feFctId of primary result also to RHS result
      rhsFeFunctions_[fncIt->first]->SetFctId(actFct->GetFctId());
      rhsFeFunctions_[fncIt->first]->SetSystem(this->algsys_);
      rhsFeFunctions_[fncIt->first]->Finalize();
      fncIt++;
    }
    
    if ( analysistype_ == TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())->GetDeltaT();

      // Call the init function of timescheme of each fefunction
      fncIt= feFunctions_.begin();
      while(fncIt != feFunctions_.end()){
        shared_ptr<BaseFeFunction> actFct = fncIt->second;
        actFct->GetTimeScheme()->Init(actFct->GetSingleVector(),dt);
        fncIt++;
      }
    }
    //! Define step solution driver
    if ( isDirectCoupled_ == false ) {
      LOG_TRACE(singlepde) << pdename_ << ": Defining solveStep class";
      DefineSolveStep();

      // check if solve step was defined
      if(! solveStep_) {
        EXCEPTION("No solveStep object defined for PDE '" << pdename_ << "'");
      }
    }

    // =====================================================================
    // define which solution types have to be saved
    // =====================================================================
    LOG_TRACE(singlepde) << pdename_ << ": Reading store results";

    FinalizePostProcResults();
    
    // Read result only if a free simulation is performed
    if( !simState_->HasInput()) {
      ReadStoreResults();

      // We might have defined some generic results which depend on the feFunction.
      // Since they get initialized in stage 2, we had to pass the feFunction without
      // knowing the actual coefFunction. We collected these BCs and update their
      // coefFunction know accordingly
      UpdateCoefFuncsForPostProcResults();

      ReadSensorArrayResults();
    }


    // =====================================================================
    // Set the initial conditions
    // =====================================================================
    if( !simState_->HasInput()) {
      ReadInitialConditions();
    }
    
    // Finally set the initialization flag to true
    isInitialized_ = true;
//    LOG_DBG(singlepde) << "IS3: has MECH_DISPLACEMENT fefunciton? " << (feFunctions_.find(MECH_DISPLACEMENT) != feFunctions_.end());
    LOG_TRACE(singlepde) << pdename_ << ": Finished initializaton";

  }


  void SinglePDE::UpdateCoefFuncsForPostProcResults() {
    // loop over a list of dummy coefFunctions and set the true coefFunction
    for (UInt i=0; i<coefsToUpdate_.GetSize(); i++) {
      // init the info to retrieve the coefFunction
      SolutionType solType;
      shared_ptr<EntityList> list;
      std::string pdeName;

      coefsToUpdate_[i]->GetCoefInfo(solType, list, pdeName);

      bool updateGeo;
      PtrCoefFct coef = iterCplPde_->GetCouplingCoefFct(solType, list, pdeName, updateGeo);

      coefsToUpdate_[i]->SetCoef(coef);
    }
  }


  void SinglePDE::SetIterCoupledPDE( IterCoupledPDE* itCplPde ) {
    iterCplPde_ = itCplPde;
    isIterCoupled_ = true;
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

        // thermal radiation BC doesn't have model name attribute
        if (actTypeString != "thermalRadiation") {
          modelName_ = nonLinNodes[i]->Get("model")->As<std::string>();
        }

        //std::cout << "actTypeString " << actTypeString << std::endl;
        //std::cout << "actId " << actId << std::endl;
        
        NonLinType actType;
        String2Enum( actTypeString, actType );

        LOG_DBG(singlepde) << pdename_ << ": INL: set nonLinTypes[" << actId << "]";
        //save for each nonlinearity type the id (Why a map with the string as key?!)
        nonLinTypes_[actId] = actType;
      }
    

      // Run over all region and set entry in "regionNonLinId"
      ParamNodeList regionNodes = 
        myParam_->Get("regionList")->GetChildren();
      
      RegionIdType actRegionId;
      std::string actRegionName, actNonLinId;//, actHysteresis;
      
//           if( regionNodes.GetSize() > 0 ) {
//             Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
//           }
      
      for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
        //take care: one region can have more then one nonlinearity!!
        
        // get data
        regionNodes[i]->GetValue( "name", actRegionName );
        regionNodes[i]->GetValue( "nonLinIds", actNonLinId );
        //regionNodes[i]->GetValue("hysteresis", actHysteresis);

        
        //std::cout << "actRegionName " << actRegionName << std::endl;
        //std::cout << "actNonLinId " << actNonLinId << std::endl;
        if( actNonLinId == "" )
          continue;
        
        typedef boost::tokenizer< boost::char_separator<char> > Tok;
        boost::char_separator<char> sep(";|, ");
        
        Tok tok(actNonLinId, sep);
        
        actRegionId = ptGrid_->GetRegion().Parse( actRegionName );
        //std::cout << "actRegionId " << actRegionId << std::endl;
        
        
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

          //Here i have to do something. O
          //Original Code:
           if ( nonLinTypes_[nonLinId] == HYSTERESIS ){
             isHysteresis_ = true;
            //or nonLinTypes_[nonLinId] == HYSTERESIS_FIXPOINT )// enum removed
            //fixpoint iteration can be selected via evaluationParameter flag
            // > see stdSolveStep for more details
          } else {
            // new flag used by std solvestep hyst
            // it indicates, if there are additional nonlinearities, which are not due to
            // hysteresis (e.g. nonlinear bhcurve in a different region)
            // this information is needed for the residual computation
            // pure hysteretic case uses the linear system matrix during residual computation
            // nonlinear bh-curve used the nonlinear system matrix
            nonLinNonHyst_ = true;
          }

        }
      }

      if(isHysteresis_ == false){
        // Here we need in addition the nonLinMethod_ for the definition
        // of the integrators
        nonLinMethod_ = FIXEDPOINT;
        PtrParamNode nonLinNode = solStrat_->GetNonLinNode();
        // NEW: additionally check if nonLinearity is used at all for some region
        // otherwise we do not have to search for nonlinear methods
        if(( nonLinNode ) && (nonLin_ == true)) {
          std::string methodString;
          nonLinNode->GetValue(  "method", methodString, ParamNode::PASS );
          nonLinMethod_ = NonLinMethodTypeEnum.Parse(methodString);
        }
      } else {
        // > read in during solveStep hyst
      }
    }
  }


  // *********************************
  //  Initialize material depenencies
  // *********************************
  void SinglePDE::InitMaterialDependencies() {

    matDepend_ = false;

    // Check, if "matDependencyList" is present
    PtrParamNode matDepListNode = myParam_->Get("matDependencyList", ParamNode::PASS );
    if( matDepListNode ) {

      // Get nonlinear types
      ParamNodeList matDepNodes = matDepListNode->GetChildren();
      for( UInt i = 0; i < matDepNodes.GetSize(); i++ ) {

        std::string actTypeString = matDepNodes[i]->GetName();
        std::string actId = matDepNodes[i]->Get("id")->As<std::string>();

        NonLinType actType;
        String2Enum( actTypeString, actType );

        //save for each nonlinearity type the id
        matDepTypes_[actId] = actType;
      }


      // Run over all region and set entry in "regionNonLinId"
      ParamNodeList regionNodes =
        myParam_->Get("regionList")->GetChildren();

      RegionIdType actRegionId;
      std::string actRegionName, actMatDepId;

      for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
        //take care: one region can have more then one material dependency!!

        // get data
        regionNodes[i]->GetValue( "name", actRegionName );
        regionNodes[i]->GetValue( "matDependIds", actMatDepId );

        if( actMatDepId == "" )
          continue;

        typedef boost::tokenizer< boost::char_separator<char> > Tok;
        boost::char_separator<char> sep(";|, ");

        Tok tok(actMatDepId, sep);

        actRegionId = ptGrid_->GetRegion().Parse( actRegionName );

        for(Tok::iterator it=tok.begin(); it!=tok.end(); ++it) {
          std::string matDepId = (*it);

          if(matDepTypes_.find(matDepId) == matDepTypes_.end()) {
            WARN( "Material depenency with id '" << matDepId
                  << "' was not defined in 'matDependencyList'");
            continue;
          }

          regionMatDepTypes_[actRegionId].Push_back( matDepTypes_[matDepId] );

          //write info
          std::string matDepString;
          Enum2String( matDepTypes_[matDepId], matDepString );

          //if one nonlinearity is set, then the whole PDE is set to nonlinear
          matDepend_ = true;
        }
      }
    }
  }

   /** can generally be called multiple times. We overwrite old values! Brute force but keeps data size */
   void SinglePDE::WriteGeneralPDEdefines() {
     
     PtrParamNode feFctNode = infoNode_->Get("feFunctions");
     
     // loop over all feFunctions
     std::map<SolutionType, shared_ptr<BaseFeFunction> >::const_iterator it = feFunctions_.begin();
     for( ; it != feFunctions_.end(); ++it ) {
       
       SolutionType solType = it->first;
       shared_ptr<BaseFeFunction> feFct = it->second;
       
       std::string solName = SolutionTypeEnum.ToString(solType);
       PtrParamNode actNode = feFctNode->Get(solName);
       
       // === Homogeneous Dirichlet BC ===
       PtrParamNode hdbcNode = feFctNode->Get("homDirichletBC");
       HdBcList hdbcs = feFct->GetHomDirichletBCs();
       
       HdBcList::iterator hdbcIt = hdbcs.Begin();
       for( ; hdbcIt != hdbcs.End(); ++hdbcIt ) {
         HomDirichletBc & actBc = *(*hdbcIt);
         EntityList const & actList = *actBc.entities;
        
         PtrParamNode in = hdbcNode->GetByVal(actList.listType.ToString(actList.GetType()), 
                                              "name", actList.GetName());
         std::string dofString;
         std::set<UInt>::const_iterator dofIt = actBc.dofs.begin();
         for(; dofIt != actBc.dofs.end(); ++dofIt) {
           dofString += feFct->GetResultInfo()->dofNames[*dofIt] + " "; 
         }
         boost::trim(dofString);
         in->Get("dofs")->SetValue(dofString);
       }

       // === Inhomogeneous Dirichlet BC ===
       PtrParamNode idbcNode = feFctNode->Get("inhomDirichletBC");
       IdBcList idbcs = feFct->GetInHomDirichletBCs();

       IdBcList::iterator idbcIt = idbcs.Begin();
       for( ; idbcIt != idbcs.End(); ++idbcIt ) {
         InhomDirichletBc & actBc = *(*idbcIt);
         EntityList const & actList = *actBc.entities;

         PtrParamNode in = idbcNode->GetByVal(actList.listType.ToString(actList.GetType()), 
                                              "name", actList.GetName());
         std::string dofString;
         std::set<UInt>::const_iterator dofIt = actBc.dofs.begin();
         for(; dofIt != actBc.dofs.end(); ++dofIt) {
           dofString += feFct->GetResultInfo()->dofNames[*dofIt] + " "; 
         }
         boost::trim(dofString);
         in->Get("dofs")->SetValue(dofString);
         in->Get("value")->SetValue(actBc.value->ToString());
       }
     }
  }

  void SinglePDE::ReadStoreResults() {

    ResultSet::iterator it;

    // fetch result node and leave, if none is present
    PtrParamNode resultNode = myParam_->Get("storeResults", ParamNode::PASS);
    if( !resultNode )
      return;

    // Iterate over all available results
    for (it = availResults_.begin(); it != availResults_.end(); it++ ) {
       CheckStoreResult(*it);
    }
  }


  bool SinglePDE::CheckStoreResult(shared_ptr<ResultInfo> candidate) {

    StdVector<std::string> regionNames, nodeNames, writeResults, actOutDest;
    StdVector<std::string> postProcNames, outDestNames, neighborRegions, writeAsHistResult;
    UInt saveBegin = 0, saveEnd = 0, saveInc = 0;
    std::string quantity, complexFormatString, listElemName, entityName;
    ComplexFormat complexFormat;
    shared_ptr<EntityList> actList;

    EntityList::ListType entityType;
    ResultHandler * resHandler = domain_->GetResultHandler();

    // initialize map for relating EntityUnknownType and name of xml-element
    using std::make_pair;
    std::map<ResultInfo::EntityUnknownType, std::string> elemNames;
    std::map<ResultInfo::EntityUnknownType, bool> isHistory;
    elemNames.insert(make_pair(ResultInfo::NODE, "nodeResult"));
    elemNames.insert(make_pair(ResultInfo::ELEMENT, "elemResult"));
    elemNames.insert(make_pair(ResultInfo::SURF_ELEM, "surfElemResult"));
    elemNames.insert(make_pair(ResultInfo::REGION, "regionResult"));
    elemNames.insert(make_pair(ResultInfo::REGION_AVERAGE, "regionAverageResult"));
    elemNames.insert(make_pair(ResultInfo::SURF_REGION, "surfRegionResult"));
    elemNames.insert(make_pair(ResultInfo::COIL, "coilResult"));

    isHistory.insert(make_pair(ResultInfo::NODE, false));
    isHistory.insert(make_pair(ResultInfo::ELEMENT, false));
    isHistory.insert(make_pair(ResultInfo::SURF_ELEM, false));
    isHistory.insert(make_pair(ResultInfo::REGION, true));
    isHistory.insert(make_pair(ResultInfo::REGION_AVERAGE, true));
    isHistory.insert(make_pair(ResultInfo::SURF_REGION, true));
    isHistory.insert(make_pair(ResultInfo::COIL, true));


    // fetch result node and leave, if none is present
    PtrParamNode resultNode = myParam_->Get("storeResults", ParamNode::PASS);
    if( !resultNode )
      return false;

    // Convert enum
    quantity = SolutionTypeEnum.ToString(candidate->resultType);
    LOG_DBG(singlepde) << pdename_ << ": Searching for storeResults of quantity '" << quantity << "'";
    // try to catch possible errors
    try {

      // Get type of result
      std::string xmlElemName = elemNames[candidate->definedOn];
      if( xmlElemName == "" ){
        return false;
      }
      LOG_DBG(singlepde) << pdename_ << ":   xmlElemName = " << xmlElemName;

      // Remember current result node
      PtrParamNode actResultNode =
          resultNode->GetByVal(xmlElemName, "type", quantity, ParamNode::PASS );
      LOG_DBG(singlepde) << pdename_ << ":   quantity = " << quantity;

      // Check on which entity type the result is defined on
      switch(candidate->definedOn) {
      case ResultInfo::NODE:
        entityType = EntityList::NODE_LIST;
        LOG_DBG(singlepde) << pdename_ << ":   -> defined on nodes";
        break;
      case ResultInfo::REGION_AVERAGE:
        LOG_DBG(singlepde) << pdename_ << ":   -> defined on REGION_AVERAGE";
        entityType = EntityList::NAME_LIST;
        break;
      case ResultInfo::REGION:
        LOG_DBG(singlepde) << pdename_ << ":   -> defined on REGION";
        entityType = EntityList::NAME_LIST;
        break;
      case ResultInfo::SURF_REGION:
        entityType = EntityList::NAME_LIST;
        break;
      case ResultInfo::SURF_ELEM:
        entityType = EntityList::SURF_ELEM_LIST;
        break;
      case ResultInfo::ELEMENT:
        entityType = EntityList::ELEM_LIST;
        break;
      case ResultInfo::COIL:
        entityType = EntityList::COIL_LIST;
        break;
      default:
        EXCEPTION("Type of 'definedOn' was not found");
        break;
      }

      // intialize variables
      neighborRegions.Clear();
      regionNames.Clear();

      // ========== Look for defineType 'REGION' ==========
      // if no node was found, continue with next result
      if( !actResultNode) {
        LOG_DBG(singlepde) << pdename_ << ":   not found";
        return false;
      } else {
        LOG_DBG(singlepde) << pdename_ << ": ------- FOUND -------";
      }

      // determine complexFormat
      complexFormatString = "amplPhase";
      actResultNode->GetValue("complexFormat", complexFormatString, ParamNode::PASS);
      String2Enum( complexFormatString, complexFormat );

      // otherwise check, if result is to be saved on "allRegions"
      if( actResultNode->Has("allRegions" ) ) {
        ptGrid_->GetRegion().ToString(regions_,regionNames);

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
            candidate->definedOn == ResultInfo::REGION_AVERAGE ||
            candidate->definedOn == ResultInfo::REGION ) {
          listNode = actResultNode->Get("regionList", ParamNode::PASS);
          if( listNode )
            regionNodes = listNode->GetList("region");
        }
        else if(candidate->definedOn == ResultInfo::SURF_ELEM ||
            candidate->definedOn == ResultInfo::SURF_REGION ) {
          listNode = actResultNode->Get("surfRegionList", ParamNode::PASS);
          if( listNode )
            regionNodes = listNode->GetList("surfRegion");

//          // fetch entry with neighboring regions
//          for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
//            std::string str = regionNodes[i]->Get("neighborRegion")->As<std::string>();
//            neighborRegions.Push_back( str );
//          }
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
            // if ( candidate->resultType == MAG_FORCE_MAXWELL_DENSITY ||
            //      candidate->resultType == MAG_FORCE_MAXWELL ||
            //    candidate->resultType == MAG_FORCE_VWP) {
            //   neighborRegions.Push_back( regionNodes[i]->Get("neighborRegion")->As<std::string>());
            // }
          }
        }
      }

      // Check, if any region was found for this result type
      if( regionNames.GetSize() != 0 ) {
        candidate->complexFormat = complexFormat;

        // iterate over all regions
        for( UInt iRegion = 0; iRegion < regionNames.GetSize(); iRegion++ )
        {
          actList = ptGrid_->GetEntityList( entityType, regionNames[iRegion] );
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


          // try to get result functor
          shared_ptr<ResultFunctor> fnc;
          if( candidate->definedOn == ResultInfo::REGION_AVERAGE ) {
            LOG_DBG(singlepde) << pdename_ << "  -> regionAverage result we need to get functor from fieldAverageFunctors_";
            if(fieldAverageFunctors_.find(candidate->resultType) == fieldAverageFunctors_.end()){
              LOG_DBG(singlepde) << pdename_ << "    not found in fieldAverageFunctors_, trying resultFunctors_";
              EXCEPTION( "No result functor defined for results of type '" << quantity << "' (field average) - this should not happen - check DefineFieldResult");
            } else {
              fnc = fieldAverageFunctors_[candidate->resultType];
            }
          } else {
            if(resultFunctors_.find(candidate->resultType) == resultFunctors_.end()){
              LOG_DBG(singlepde) << pdename_ << "  -> no functors found";
              EXCEPTION( "No result functor defined for results of type '" << quantity << "'");
            } else {
              fnc = resultFunctors_[candidate->resultType];
            }
          }

          // if ( candidate->resultType == MAG_FORCE_MAXWELL_DENSITY ||
          //      candidate->resultType == MAG_FORCE_MAXWELL ||
          //    candidate->resultType == MAG_FORCE_VWP) {
          //   std::string neighborReg =  neighborRegions[iRegion];
          //   RegionIdType surfRegionId = ptGrid_->GetRegion().Parse( regionNames[iRegion] );
          //   RegionIdType volNeighborRegionId = ptGrid_->GetRegion().Parse( neighborReg );
          //   fnc->GetCoefFct()->SetVolNeighborRegionId(surfRegionId,volNeighborRegionId);
          //   }

          // update sequence step for result handler (also done in ResultHandler::BeginMultiSequenceStep, but AFTER reading the xml for the current step, hence we would be one step behind regarding reading e.g. postproc results from the xml)
          // nevertheless, it appears that we have to keep it in ResultHandler::BeginMultiSequenceStep as well since not setting the sequence step in this routine breaks stuff
          resHandler->SetSequenceStep(sequenceStep_);
          // pass result to result  
          resHandler->RegisterResult( actSol, fnc, sequenceStep_,
              saveBegin, saveInc, saveEnd,
              actOutDest,
              postProcNames[iRegion], writeResult,
              isHistory[candidate->definedOn] );
        }
      }


      // ========== Look for defineType node/elemList/coilList (history) ==========

      std::string entityTypeName;
      StdVector<std::string> histNames;
      neighborRegions.Clear();
      writeAsHistResult.Clear();

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
        for( UInt i = 0; i < histEntities.GetSize(); i++ ) {
          std::string str = histEntities[i]->Get("neighborRegion")->As<std::string>();
          if ( str != "" )
            neighborRegions.Push_back( str );
        }
      } else if(candidate->definedOn == ResultInfo::COIL ) {
        histNode = actResultNode->Get("coilList", ParamNode::PASS);
        if( histNode )
          histEntities = histNode->GetList("coil");
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
          std::string nameType = "name";
          if( candidate->definedOn == ResultInfo::COIL )
            nameType = "id";
          histNames.Push_back( histEntities[i]->Get(nameType)->As<std::string>() );
          postProcNames.Push_back( histEntities[i]->Get("postProcId")->As<std::string>() );
          outDestNames.Push_back( histEntities[i]->Get("outputIds")->As<std::string>() );
          writeResults.Push_back( histEntities[i]->Get("writeResult")->As<std::string>() );
          writeAsHistResult.Push_back( histEntities[i]->Get("writeAsHistResult")->As<std::string>() );
        }
      }

      if ( histNames.GetSize() > 0 ) {
        candidate->complexFormat = complexFormat;

        // iterate over all entityNames
        for( UInt i = 0; i < histNames.GetSize(); i++ )
        {
          if( candidate->definedOn != ResultInfo::COIL ){
            actList = ptGrid_->GetEntityList( entityType, histNames[i] );
          } else {
            // The grid does not know about coils beause depending on the space used
            // we don't know if we need approximation in space, e.g. with the FeSpaceConst.
            // But we know that we want only one result per coil, not for each element in the coil.
            shared_ptr<Coil> actCoil;
            if( pdename_ == "magneticEdge" ){
              MagEdgePDE* askThePDE = dynamic_cast<MagEdgePDE*>(this);
              actCoil = askThePDE->GetCoilById( histNames[i] );
            } else if( pdename_ == "magneticEdgeSpecialAV" ){
              MagEdgeSpecialAVPDE* askThePDE = dynamic_cast<MagEdgeSpecialAVPDE*>(this);
              actCoil = askThePDE->GetCoilById( histNames[i] );
            } else {
              MagneticPDE* askThePDE = dynamic_cast<MagneticPDE*>(this);
              actCoil = askThePDE->GetCoilById( histNames[i] );
            }
            shared_ptr<CoilList> singleCoilList( new CoilList( ptGrid_ ) );
            singleCoilList->AddCoil( actCoil );
            actList = singleCoilList;
          }
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
          bool writeAsHistoryResult = ( writeAsHistResult[i] == "yes"  ? true : false );

          // try to get result functor
          shared_ptr<ResultFunctor> fnc;
          if( resultFunctors_.find(candidate->resultType) ==
              resultFunctors_.end() ) {
            EXCEPTION( "No result functor defined for results of type '"
                << quantity << "'");
          }

          fnc = resultFunctors_[candidate->resultType];
//          if ( neighborRegions.GetSize() != 0 ) {
//            std::string neighborReg =  neighborRegions[i];
//            RegionIdType actRegionId = ptGrid_->GetRegion().Parse( neighborReg );
//            fnc->GetCoefFct()->SetNeighborRegionId(actRegionId);
//          }

          resHandler->RegisterResult( actSol, fnc, sequenceStep_,
              saveBegin, saveInc, saveEnd,
              actOutDest, postProcNames[i],
              writeResult, writeAsHistoryResult);

        }
      }
    } catch( Exception &ex ) {
      RETHROW_EXCEPTION(ex, "Could not determine storeResults for quantity '"
          << quantity << "' within pde '" << pdename_ << "'" );
    }
    return true;
  }
  
  void SinglePDE::FinalizePostProcResults() {
    {
      // 1) Associate all stiffness related coeffunctions and result functors
      //    with the bilinearform
      std::multimap<RegionIdType, BaseBDBInt*>::iterator stiffIt = bdbInts_.begin();
      for(; stiffIt != bdbInts_.end(); ++stiffIt ) {
        RegionIdType region = stiffIt->first;
        BaseBDBInt* bdb = stiffIt->second;
        if( !bdb)
          continue;

        // 1) pass it to all coefficient functions related to stiffness
        // Although multiple integrators can be stored in bdbInts_, only one can be passed to the coefFunction
        // This can be controlled by assigning the name of the requested integrator
        std::set<shared_ptr<CoefFunctionFormBased> >::iterator stiffCoefIt;
        for( stiffCoefIt = stiffFormCoefs_.begin();
            stiffCoefIt != stiffFormCoefs_.end(); ++stiffCoefIt) {
          (*stiffCoefIt)->AddIntegrator(bdb, region);
        }
        // 2) pass it to all result functors related to stiffness
        std::set<shared_ptr<ResultFunctor> >::iterator stiffFuncIt;
        for( stiffFuncIt = stiffFormFunctors_.begin();
            stiffFuncIt != stiffFormFunctors_.end(); ++stiffFuncIt) {
          (*stiffFuncIt)->AddIntegrator(bdb, region);
        }
        // 3) set region to to all surfCoefFcts
        std::map<shared_ptr<CoefFunctionSurf>, PtrCoefFct >::iterator surfCoefIt;
        for( surfCoefIt = surfCoefFcts_.begin(); 
            surfCoefIt != surfCoefFcts_.end(); ++surfCoefIt ) {
          surfCoefIt->first->AddVolumeCoef(region, surfCoefIt->second);
        }
      }


      // 1.1) Auxiliary BDB Integrator, associate all stiffness related
      // coeffunctions and result functors with the bilinearform
      std::map<RegionIdType, BaseBDBInt*>::iterator stiffItAux1 = bdbIntsAux1_.begin();
      for(; stiffItAux1 != bdbIntsAux1_.end(); ++stiffItAux1 ) {
        RegionIdType region = stiffItAux1->first;
        BaseBDBInt* bdb = stiffItAux1->second;
        if( !bdb)
          continue;

        // 1) pass it to all coefficient functions related to stiffness
        std::set<shared_ptr<CoefFunctionFormBased> >::iterator stiffCoefIt;
        for( stiffCoefIt = stiffFormCoefsAux1_.begin();
            stiffCoefIt != stiffFormCoefsAux1_.end(); ++stiffCoefIt) {
          (*stiffCoefIt)->AddIntegrator(bdb, region);
        }
        // 2) pass it to all result functors related to stiffness
        std::set<shared_ptr<ResultFunctor> >::iterator stiffFuncIt;
        for( stiffFuncIt = stiffFormFunctorsAux1_.begin();
            stiffFuncIt != stiffFormFunctorsAux1_.end(); ++stiffFuncIt) {
          (*stiffFuncIt)->AddIntegrator(bdb, region);
        }
        // 3) set region to to all surfCoefFcts
        std::map<shared_ptr<CoefFunctionSurf>, PtrCoefFct >::iterator surfCoefIt;
        for( surfCoefIt = surfCoefFctsAux1_.begin();
            surfCoefIt != surfCoefFctsAux1_.end(); ++surfCoefIt ) {
          surfCoefIt->first->AddVolumeCoef(region, surfCoefIt->second);
        }
      }
    }

    // 2) Associate all mass related coeffunctions and result functors
    //    with the bilinearform
    std::map<RegionIdType, BaseBDBInt*>::iterator massIt = massInts_.begin();
    for( ; massIt != massInts_.end(); ++massIt ) {

      RegionIdType region = massIt->first;
      BaseBDBInt* mass = massIt->second;

      // check, that mass integrator is defined at all
      if( !mass)
        continue;

      // 1) pass it to all coefficient functions related to mass
      std::set<shared_ptr<CoefFunctionFormBased> >::iterator massCoefIt;
      for( massCoefIt = massFormCoefs_.begin();
          massCoefIt != massFormCoefs_.end(); ++massCoefIt) {
        (*massCoefIt)->AddIntegrator(mass, region);
      }
      // 2) pass it to all result functors related to mass
      std::set<shared_ptr<ResultFunctor> >::iterator massFuncIt;
      for( massFuncIt = massFormFunctors_.begin();
          massFuncIt != massFormFunctors_.end(); ++massFuncIt) {
        (*massFuncIt)->AddIntegrator(mass, region);
      }
    }
    
    // 3) Pass regions of primary FeFunction to all time derivatives
    if( analysistype_ != STATIC && analysistype_ != BUCKLING ) {
      // Workaround for transient simulation:
      // We have to pass directly the time-derivative vector of the time integration
      // scheme to the feFunctions related to time derivative results
      std::map<SolutionType, shared_ptr<BaseFeFunction> >::const_iterator it;
      for( it = timeDerivFeFunctions_.begin(); it != timeDerivFeFunctions_.end(); ++it ) {
        shared_ptr<BaseFeFunction> primFeFct = feFunctions_[timeDerivPrimaryResults_[it->first]];
        shared_ptr<BaseFeFunction> derivFeFct = it->second;
        if( !primFeFct || !derivFeFct ) {
          EXCEPTION( "The time derivative information for PDE '" << pdename_  << "' is not set correctly!" );
        }
        LOG_DBG(singlepde) << "FPPR: " << pdename_ << " derivFeFct =" << derivFeFct->ToString() << " complex:" << derivFeFct->IsComplex();

        // Now loop over all regions of the primary FeFunction and pass the
        // all regions of the primary one
        StdVector< shared_ptr<EntityList> > support =  primFeFct->GetEntityList();
        for( UInt i = 0; i < support.GetSize(); ++i ) {
          derivFeFct->AddEntityList( support[i] );
        }
        
        derivFeFct->Finalize();
        derivFeFct->SetPDE(this);
        UInt timeDerivOrder = timeDerivOrder_[it->first];
        if( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ||  analysistype_ == INVERSESOURCE || analysistype_ == EIGENFREQUENCY || analysistype_ == EIGENVALUE || analysistype_ == HARMONIC25D) {
          FeFunction<Complex> & cDerivFct = dynamic_cast<FeFunction<Complex>& >(*derivFeFct);
          shared_ptr<FeFunction<Complex> > cPrimFct = dynamic_pointer_cast<FeFunction<Complex> >(primFeFct);
          cDerivFct.SetTimeDerivOrder( timeDerivOrder, cPrimFct );
        } else {
          primFeFct->GetTimeScheme()->SetTimeDerivVector(timeDerivOrder, derivFeFct->GetSingleVector());
          simState_->RegisterFeFct(derivFeFct);
        }
      }
    }
  }


  PtrCoefFct SinglePDE::GetCoefFct( SolutionType type ) {
    LOG_DBG(singlepde) <<  pdename_ << ": Requesting coefficient function for solution type " << SolutionTypeEnum.ToString(type);
    
    PtrCoefFct ret;
    // HACK!!!!!!!!!!!!!!!!!
    if ( type == MAG_MAGNETIZATION ) {
      if ( domain_->GetRegion4Hyst() != NO_REGION_ID )
        ret = matModelCoefm_[domain_->GetRegion4Hyst()];
    }
    else {
      // 1) look in fieldCoefs
      if ( fieldCoefs_.find(type) == fieldCoefs_.end() ) {
        if( matCoefs_.find(type) == matCoefs_.end() ) {
//        EXCEPTION( "No coefficient function for result type '"
//            << SolutionTypeEnum.ToString( type ) << "' found");
        } else {
          ret = matCoefs_[type];
       }
      } else {
        ret = fieldCoefs_[type];
      }
    }
    if( !ret ) {
      LOG_DBG(singlepde) << pdename_ << ": \t=> NOT FOUND";
    } else {
      LOG_DBG(singlepde) << pdename_ << ": \t=> SUCCESS!";
    }    
    return ret;
  }
  
  SubTensorType SinglePDE::GetSubTensorType() const
  {
    if(subType_ == "") return NO_TENSOR;

    SubTensorType stt;
    String2Enum(subType_, stt);
    return stt;
  }

  void SinglePDE::WriteResultsInFile(const UInt kstep, const Double actTimeFreq)
  {
    LOG_DBG(singlepde) << pdename_  << ": WriteResultsInFile() kstep: " <<  kstep << " actTimeFreq: " << actTimeFreq;
    
    // shared amongst e.g. WriteResults, UpdateResults, FinishStep
    shared_ptr<Timer> timer = domain->GetInfoRoot()->Get(ParamNode::HEADER)->Get("results/timer")->AsTimer();
    timer->Start();

    // ===================================================
    //  Trigger calculation of interpolated field results 
    // ===================================================
    
    // Check for additional field variable
    UInt numFields = sensors_.GetSize();
    
    // loop over all fields variables
    for( UInt i = 0; i < numFields; ++i ) {
    
      // call specialized calculation method in sub-class
      FieldAtPoints& fap = sensors_[i];
      
      
      // Obtain field resultFunctor object
      SolutionType solType = fap.resultInfo->resultType;
      StdVector<std::string> dofNames = fap.resultInfo->dofNames;
      UInt numDofs = dofNames.GetSize();
      std::string solTypeString;
      solTypeString = SolutionTypeEnum.ToString(solType);
      std::map<SolutionType, PtrCoefFct >::iterator fctIt;
      fctIt = fieldCoefs_.find(solType);
      if( fctIt == fieldCoefs_.end() )  {
        EXCEPTION( "Could not find field functor for result '" 
            << SolutionTypeEnum.ToString(solType) << "'");
      }
      PtrCoefFct fct =  fctIt->second;
      
      // calculate vector entries
      if( isComplex_) {
        Vector<Complex> temp;
        Vector<Complex>& vec = dynamic_cast<Vector<Complex> &>(*fap.field);
        vec.Resize(fap.elems.GetSize() * numDofs);
        UInt pos = 0;
        LocPointMapped lpm;
        for ( UInt iElem = 0; iElem < fap.elems.GetSize(); ++iElem ) {
          shared_ptr<ElemShapeMap> esm =
              ptGrid_->GetElemShapeMap( fap.elems[iElem], true );
          lpm.Set(fap.locPoints[iElem], esm, 0.0);
          fct->GetVector(temp, lpm );
          
          for( UInt i = 0; i < numDofs; ++i ) {
            vec[pos++] = temp[i];
          }
        }
      } else {
        Vector<Double> temp;
        Vector<Double>& vec = dynamic_cast<Vector<Double> &>(*fap.field);
        vec.Resize(fap.elems.GetSize() * numDofs);
        UInt pos = 0;
        LocPointMapped lpm;
        for ( UInt iElem = 0; iElem < fap.elems.GetSize(); ++iElem ) {
          shared_ptr<ElemShapeMap> esm =
              ptGrid_->GetElemShapeMap( fap.elems[iElem], true );
          lpm.Set(fap.locPoints[iElem], esm, 0.0);
          fct->GetVector(temp, lpm );

          for( UInt i = 0; i < numDofs; ++i ) {
            vec[pos++] = temp[i];
          }
        }
      }


      std::ofstream  out((fap.fileName+"-"+lexical_cast<std::string>(kstep)).c_str(),
                          std::ios::out );

      // Ensure that no precision is lost
      out.precision(15);

      Vector<Double> globPoint, globPointcSys;
      
      StdVector<std::string> globCoordNames;
      StdVector<std::string> locCoordNames;
      for(UInt i = 0; i < dim_; ++i ) {
        globCoordNames.Push_back(fap.coordSys->GetDofName(i+1));
      }
      locCoordNames.Push_back("xi");
      locCoordNames.Push_back("eta");
      locCoordNames.Push_back("zeta");      
      std::string delim = "\t";
      if(fap.csv) 
      {
        delim = std::string(1,fap.delim);
      }
      
      // print out information
      if( isComplex_ ){
        // cast solution vector
        Vector<Complex>& vec = dynamic_cast<Vector<Complex> &>(*(fap.field));

        // Write header line with descriptions of columns
        if(fap.csv) 
        {
          out << "origElemNum" << delim;        
          for(UInt j = 0; j < dim_; ++j ) {
            out << "globCoord-" << globCoordNames[j] << delim;
          }
          for(UInt j = 0; j < numDofs; ++j ) {
            out << solTypeString << "-real" << "-" << dofNames[j] << delim;
          }
          for(UInt j = 0; j < numDofs; ++j ) {
            out << solTypeString << "-imag" << "-" << dofNames[j] << delim;
          }
          for(UInt j = 0; j < dim_-1; ++j ) {
            out << "locCoord-" << locCoordNames[j] << delim;
          }
          out << "locCoord-" << locCoordNames[dim_-1] << std::endl;
        }
      
        // Loop over all points
        for( UInt iPoint = 0; iPoint < fap.locPoints.GetSize(); iPoint++) { 
          
          shared_ptr<ElemShapeMap> esm =
              ptGrid_->GetElemShapeMap(fap.elems[iPoint], true);
          esm->Local2Global(globPoint, fap.locPoints[iPoint]);
          
          fap.coordSys->Global2LocalCoord(globPointcSys, globPoint);
          
          // write to file
          out << fap.elems[iPoint]->elemNum << delim;
          out << globPointcSys.ToString(TS_PLAIN, delim) << delim;
          for(UInt j = 0; j < numDofs; ++j ) {
            out << vec[iPoint*numDofs + j].real() << delim;
          }
          for(UInt j = 0; j < numDofs; ++j ) {
            out << vec[iPoint*numDofs + j].imag() << delim;
          }
          for(UInt j = 0; j < dim_-1; ++j ) {
            out << fap.locPoints[iPoint][j] << delim;
          }
          out << fap.locPoints[iPoint][dim_-1] << std::endl;
        }

      } else {
        // cast solution vector
        Vector<Double>& vec = dynamic_cast<Vector<Double> &>(*(fap.field));

        // Write header line with descriptions of columns
        if(fap.csv) 
        {
          out << "origElemNum" << delim;        
          for(UInt j = 0; j < dim_; ++j ) {
            out << "globCoord-" << globCoordNames[j] << delim;
          }
          for(UInt j = 0; j < numDofs; ++j ) {
            out << solTypeString << "-" << dofNames[j] << delim;
          }
          for(UInt j = 0; j < dim_-1; ++j ) {
            out << "locCoord-" << locCoordNames[j] << delim;
          }
          out << "locCoord-" << locCoordNames[dim_-1] << std::endl;
        }
        
        // Loop over all points
        for( UInt iPoint = 0; iPoint < fap.locPoints.GetSize(); iPoint++) { 
          
          shared_ptr<ElemShapeMap> esm =
              ptGrid_->GetElemShapeMap(fap.elems[iPoint], true);
          esm->Local2Global(globPoint, fap.locPoints[iPoint]);
          
          fap.coordSys->Global2LocalCoord(globPointcSys, globPoint);
          // write to file
          out << fap.elems[iPoint]->elemNum << delim;
          out << globPointcSys.ToString(TS_PLAIN, delim) << delim;
          for(UInt j = 0; j < numDofs; ++j ) {
            out << vec[iPoint*numDofs + j] << delim;
          }
          for(UInt j = 0; j < dim_-1; ++j ) {
            out << fap.locPoints[iPoint][j] << delim;
          }
          out << fap.locPoints[iPoint][dim_-1] << std::endl;
        }
      }

      out.close();
    }
    timer->Stop();
  }
  
  
  void SinglePDE::ReadSensorArrayResults() {

    ParamNodeList sensorNodes;
    sensorNodes = myParam_->Get("storeResults")->GetList("sensorArray");
    std::string solTypeString;
    static std::map< SolutionType, bool> warningPrinted;

    if (sensorNodes.GetSize()==0)
      return;

    sensors_.Resize(sensorNodes.GetSize());
    // loop over all parts
    for( UInt iPart = 0; iPart <sensorNodes.GetSize(); ++iPart ) {
      PtrParamNode  actNode = sensorNodes[iPart];
      
      FieldAtPoints & actField = sensors_[iPart];
      actField.fileName = actNode->Get("fileName")->As<std::string>();

      /* check if directory-path for sensor array file exists */
      // search for last Slash in fileName
      int idx_lastSlash = actField.fileName.find_last_of("/");
      // if idx_lastSlash = -1 -> "/" not found, else position of the last slash
      // if there is a "/" in the filename -> save directory is not "." -> check if it exists
      if ( idx_lastSlash != -1){
        // get directory name
        std::string directoryName;
        directoryName = actField.fileName.substr(0,idx_lastSlash);
        // ensure errno is cleared and call mkdir with the directory name
        errno = 0;
        int mkdir_call;
        #ifndef WIN32
          mkdir_call = mkdir( directoryName.c_str(), S_IRWXU | S_IRGRP | S_IWGRP | S_IROTH );
        #else
          mkdir_call = _mkdir( directoryName.c_str() );
        #endif

        if ( mkdir_call == -1 && errno == EEXIST ){
          // directory exists, do nothing
          errno = 0;
        } else if ( mkdir_call == 0 ){
          // directory didn't exist but was created, do nothing
        } else{
          // directory didn't exist, and couldn't be created -> raise exception
          EXCEPTION("The directory: '" << directoryName << "' to save the sensor arrays doesn't exist and couldn't be created! Please create it by hand!" );
        }

        /* working alternative version
        // check existence
        struct stat directory_attr;
        int stat_call;
        errno = 0;
        stat_call = stat( directoryName.c_str(), &directory_attr);

        if( stat_call < 0 ){
          // stat call failed
          if( errno == ENOENT ){
            EXCEPTION("The directory: '" << directoryName << "' to save the sensor arrays does not exist!" );
          } else {
            // other errno
            EXCEPTION("The 'stat'-call of directory: '" << directoryName << "' to save the sensor arrays failed with errno=" << errno << " ! Check the directory path!" );
            }
        } else {
          // stat call successful -> do nothing
        }

        if(  !(directory_attr.st_mode & S_IFDIR) ){
          EXCEPTION("The specified path: '" << directoryName << "' to save the sensor arrays is not a directory!");
          }
        else{
          // is directory do nothing //
          }
          */
      } else {
        // no slash in filename -> do nothing
      }
      actField.csv = actNode->Get("csv")->As<bool>();
      std::string coordSysId = actNode->Get("coordSysId")->As<std::string>();
      actField.coordSys = domain_->GetCoordSystem(coordSysId);
      
      std::string delim = actNode->Get("delimiter")->As<std::string>();
      if(actField.csv && delim.length() == 0) 
      {
        actField.delim = ',';
      }
      else 
      {        
        actField.delim = delim[0];
      }
      
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
      
      //array of global sensor coordinates
      StdVector< Vector<Double> > globPoints;

      //get entity list of current pde
      StdVector<shared_ptr<EntityList> > lists;
      StdVector<RegionIdType>::iterator regIt = regions_.Begin();
      for(; regIt != regions_.End(); regIt++ ) {
        shared_ptr<ElemList> newList(new ElemList(ptGrid_));
        newList->SetRegion(*regIt);
        lists.Push_back(newList);
      }

      // create list
      // generate new vector
      if(isComplex_) {
        actField.field = new Vector<Complex>();
      } else {
        actField.field = new Vector<Double>();
      }

      if(actNode->Has("parametric")){
        // define sensors according to parametric line definitions
        ParamNodeList listNodes = actNode->Get("parametric")->GetList("list");
        // loop over all components
        StdVector<Double> start(3), stop(3), inc(3);
        StdVector<UInt> numSamples(3);
        start.Init(0);
        stop.Init(0);
        inc.Init(1);
        numSamples.Init(1);
        std::string comp;
        UInt compIndex;
        for( UInt iComp = 0; iComp < listNodes.GetSize(); iComp++ ) {
          PtrParamNode actCompNode = listNodes[iComp];
          actCompNode->GetValue("comp", comp);
          compIndex = actField.coordSys->GetVecComponent(comp)-1;
          start[compIndex]=  actCompNode->Get("start")->MathParse<Double>();
          stop[compIndex]=  actCompNode->Get("stop")->MathParse<Double>();
          inc[compIndex] = actCompNode->Get("inc")->MathParse<Double>();
          numSamples[compIndex] = UInt(floor( (stop[compIndex]-start[compIndex]) / inc[compIndex] ) )+1;
        }

        globPoints.Resize( numSamples[0] *
                           numSamples[1] *
                           numSamples[2] );
        UInt pIdx = 0;

        for( UInt xSample = 0; xSample < numSamples[0]; xSample++ ) {
          Double actX = start[0] + xSample * inc[0];
          for( UInt ySample = 0; ySample < numSamples[1]; ySample++ ) {
            Double actY = start[1] + ySample * inc[1];
            for( UInt zSample = 0; zSample < numSamples[2]; zSample++ ) {
              Double actZ = start[2] + zSample * inc[2];

              // transform global point w.r.t. to coordinate system
              // to global point w.r.t. to global cartesian system
              Vector<Double> globPointcSys;
              globPointcSys.Resize(dim_);

              globPointcSys[0] = actX;
              globPointcSys[1] = actY;
              if( dim_ > 2) {
                globPointcSys[2] = actZ;
              }
              actField.coordSys->Local2GlobalCoord(globPoints[pIdx++],
                                                   globPointcSys);
            } // z
          } // y
        } // x
      }else if(actNode->Has("coordinateFile")){
        globPoints.Reserve(200);

        PtrParamNode coordFileNode = actNode->Get("coordinateFile");
        std::string inFileName = coordFileNode->Get("fileName")->As<std::string>();
        std::string delim = coordFileNode->Get("delimiter")->As<std::string>();
        std::string comment = coordFileNode->Get("commentCharacter")->As<std::string>();
        UInt xCol = coordFileNode->Get("xCoordColumn",ParamNode::PASS)->As<UInt>();
        UInt yCol = coordFileNode->Get("yCoordColumn",ParamNode::PASS)->As<UInt>();
        UInt zCol = coordFileNode->Get("zCoordColumn",ParamNode::PASS)->As<UInt>();

        if(xCol == 0 || yCol ==0 || zCol == 0){
          EXCEPTION("Read coordinate file for sensor array: column indices need to be one based.");
        }

        if(comment.size()>1 || delim.size() > 1){
          WARN("Read coordinate file for sensor array: Comment and delimiter strings need to be single characters!");
        }


        if(!boost::filesystem::exists(inFileName)){
          EXCEPTION("Read coordinate file for sensor array: Could not find coordinate file \"" + inFileName + "\" to read sensor positions!");
          continue;
        }

        std::fstream coordFile(inFileName.c_str(),std::ios::in);

        std::string curLine;
        coordFile >> std::ws;
        UInt lineCounter = 0;
        while(std::getline (coordFile,curLine)){
          lineCounter++;
          //ignore leading whitespace
          string::size_type pos = 0;
          while (pos < curLine.size() && std::isspace(curLine[pos], std::locale()))
            pos++;

          curLine.erase(0, pos);

          //check for comment character
          if(curLine.at(0) == comment.at(0)){
            continue;
          }

          //tokenize line with tokenizer
          typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
          boost::char_separator<char> sep(delim.c_str());
          tokenizer tokens(curLine, sep);

          UInt numNumbers = std::distance(tokens.begin(),tokens.end());

          //ignore empty lines
          if(numNumbers==0){
            continue;
          }

          //ignore invalid lines and print a warning
          //here we check for dimension
          if( (dim_ == 2 && numNumbers < 2) ||
              (numNumbers < xCol) ||
              (numNumbers < yCol) ||
              (numNumbers < zCol) ){
            WARN("Read coordinate file for sensor array: Invalid coordinate definition at line: " << lineCounter << " in file : " << inFileName );
          }

          //finally read in the tokens
          Vector<Double> curCoord(dim_);
          tokenizer::iterator tokIter=tokens.begin();
          for(UInt i=0;i<dim_ && tokIter!=tokens.end();i++,tokIter++){
            try{
              curCoord[i] = boost::lexical_cast<Double>(*tokIter);
            }catch(const boost::bad_lexical_cast &e){
              EXCEPTION("Read coordinate file for sensor array: Error reading coordinates in line: " << lineCounter << ". " << *tokIter << " The line was:\n" << curLine << e.what());
            }
          }
          Vector<Double> globPointcSys;
          actField.coordSys->Local2GlobalCoord(globPointcSys,
                                               curCoord);
          globPoints.Push_back(globPointcSys);
        }
        globPoints.Trim();
      }else{
        EXCEPTION("Could not find valid sensor coordinate definition in xml file although tag was given.");
      }
      StdVector< LocPoint > locPoints;
      StdVector< const Elem* > elems;

      // now, map global points to local points restricted to regions of this PDE.
      ptGrid_->GetElemsAtGlobalCoords( globPoints,
                                       locPoints,
                                       elems,
                                       lists,
                                       0.0, 1.0e-2,
                                       false );

      for(UInt i=0, n=globPoints.GetSize(); i<n; i++) {
        const Elem* ptElem = elems[i];
        
        if( !ptElem ) {
          bool wP = !warningPrinted[actField.resultInfo->resultType];
          if( wP ) {
            std::stringstream sstr;
            sstr << "Could not find element at position " 
                 << globPoints[i].ToString()
                 << " for evaluation of field values for "
                 << solTypeString << ".";
            WARN( sstr.str() );
            warningPrinted[actField.resultInfo->resultType] = true;
          }
        } else {
          //               std::cerr << "locPoint for globPoint " << globPoint.ToString() 
          //                                    << " is " << locPoint.ToString() 
          //                                    << " in Elem " << ptElem->elemNum << std::endl;
               
               // check again mapping by performing loc->glob mapping
          shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(ptElem);
          //esm->Local2Global(globPoint, locPoint);
          //               std::cerr << "\tAdditional check loc->glob delivers global point " 
          //                   << globPoint.ToString() << std::endl << std::endl;
          
          actField.elems.Push_back(ptElem);
          actField.locPoints.Push_back(locPoints[i]);
        }
      }

      if(warningPrinted[actField.resultInfo->resultType]) {
        std::stringstream sstr;
        sstr << "Could not find " << (globPoints.GetSize()-actField.locPoints.GetSize())
             << " locations for evaluation of field values for "
             << solTypeString << ".";
        WARN( sstr.str() );
      }
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
      fncIt->second->ApplyLoads();
      fncIt++;
    }
  }


  void SinglePDE::ReadInitialConditions() {
    
    LOG_TRACE(singlepde) << pdename_ << ": Reading initial conditions";
    PtrParamNode icNode = myParam_->Get("initialValues", ParamNode::PASS );
    if( !icNode )
      return;
    
    // create info node for initial conditions
    PtrParamNode icInfo = infoNode_->Get("initialConditions");
    
    // ===========================
    //  1) Initial State
    // ===========================
    // Initial state denotes, that all FeFunctions and their time derivatives are
    // initialized from the same PDE of either a previous sequenceStep 
    // (multiSequence analysis) or loaded from an external HDF file (e.g. displacement,
    // velocity and acceleration).
    
    PtrParamNode isNode = icNode->Get("initialState", ParamNode::PASS );
    if( isNode ) {
      LOG_TRACE(singlepde) << pdename_ << ": Reading initial state";
      
      PtrParamNode isInfo = icInfo->Get("initialState");
      
      // Ensure, that we have a static or transient analysis
      if( !( analysistype_ == STATIC || analysistype_ == TRANSIENT ) ) {
        WARN( "Initial conditions are only meaningful in a transient analysis and "
            << "will be omitted for this type of analysis" );
      }
      
      PtrParamNode srcNode = isNode->GetChild();
      if( srcNode->GetName() == "sequenceStep" ) {
        UInt sequenceStep = srcNode->Get("index")->As<UInt>();
        Integer stepNum = srcNode->Get("step")->As<Integer>();
        bool extrapolateStatic = srcNode->Get("extrapolateStatic")->As<bool>();
        
        Domain * inDomain = NULL;
        // create SimState (for input)
        boost::shared_ptr<SimState> inState(new SimState(true, domain_));
        
        try {
          LOG_DBG(singlepde) << pdename_ << ": Use initial condition from sequenceStep " << sequenceStep;
        
          // create new simState from current hdf file
          std::string fileName = simState_->GetOutputWriter()->GetFileName().string();

          // create new param and info node (without logging to console) for the
          // newly created Domain object
          PtrParamNode node(new ParamNode());
          PtrParamNode infoNode = ParamNode::GenerateWriteNode("", "",ParamNode::APPEND); // empty filename means we don't write and ignore ParamNode::ToFile()
          boost::shared_ptr<SimInputHDF5> in;
          in.reset(new SimInputHDF5(fileName, node, infoNode));
          inState->SetInputHdf5Reader(in);

          // Get grid map of own domain, as the grids can be re-used
          SimState::GridMap gridMap = domain_->GetGridMap();

          // Obtain temporary Domain object, from which the initial state is read 
          // in. As this generates inferior logging output, we make a visual
          // break.
          LOG_DBG(singlepde) << pdename_ 
              << ": Obtaining Domain from simState object";
          LOG_DBG(singlepde) << pdename_ << ": =================================="; 
          LOG_DBG(singlepde) << pdename_ << ":  BEGIN OUTPUT OF TEMPORARY DOMAIN ";
          LOG_DBG(singlepde) << pdename_ << ": ==================================";            

          inDomain = inState->GetDomain(sequenceStep, gridMap);

          LOG_DBG(singlepde) << pdename_ << ": ==================================="; 
          LOG_DBG(singlepde) << pdename_ << ":  END OF OUTPUT OF TEMPORARY DOMAIN ";
          LOG_DBG(singlepde) << pdename_ << ": ===================================";

          // try to get last step number
          if( stepNum == -1 ) {
            Double stepVal = 0.0;
            UInt lastStepNum = 0;
            inState->GetLastStepNum(sequenceStep, lastStepNum, stepVal);
            stepNum = lastStepNum;
          }
          LOG_DBG(singlepde) << pdename_ << ": Step number to be read in: " << stepNum;

          // log to info node
          PtrParamNode seqInfo = isInfo->Get("sequenceStep");
          seqInfo->Get("index")->SetValue(sequenceStep);
          seqInfo->Get("stepNum")->SetValue(stepNum);

          // update to last step number
          inState->SetInterpolation(SimState::CONSTANT, mp_, analysistype_, 0);
          inState->UpdateToStep(sequenceStep, stepNum);

          // Obtain same PDE from new domain
          SinglePDE * inPDE = inDomain->GetSinglePDE(pdename_);

          // -------------------
          //  Primary Unknowns
          // -------------------
          // Loop over all feFunctions of this pde
          LOG_DBG(singlepde) << pdename_  << ": Transferring primary FeFunctions";
          std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
          for(; fncIt != feFunctions_.end(); fncIt++ ) {
            SolutionType solType = fncIt->first;
            shared_ptr<BaseFeFunction> & myFct = fncIt->second;
            shared_ptr<BaseFeFunction>  inFct = inPDE->GetFeFunction( solType );
            LOG_DBG3(singlepde) << pdename_ << ": \t" << SolutionTypeEnum.ToString(solType);

            // Now transfer results from inFct to myFct
            myFct->InitFromFeFunction( inFct );
          }

          // --------------------
          //  Time Derivative(s)
          // --------------------
          // 2) Loop over all time derivative functions (only if analysis type of both PDEs is
          //    transient )
          if( inPDE->GetAnalysisType() == TRANSIENT && analysistype_ == TRANSIENT ) {
            LOG_DBG(singlepde) << pdename_  << ": Transferring time derivative FeFunctions";

            fncIt= timeDerivFeFunctions_.begin();
            for(; fncIt != timeDerivFeFunctions_.end(); fncIt++ ) {
              SolutionType solType = fncIt->first;
              shared_ptr<BaseFeFunction> & myFct = fncIt->second;
              shared_ptr<BaseFeFunction>  inFct = inPDE->GetFeFunction( solType );
              LOG_DBG3(singlepde) << pdename_ << ": \t" << SolutionTypeEnum.ToString(solType);

              // Now transfer results from inFct to myFct
              myFct->InitFromFeFunction( inFct );
            }
          } // if TRANSIENT

          if( inPDE->GetAnalysisType() == STATIC && analysistype_ == TRANSIENT ) {
            // Set the extrapolateStatic bool for each feFunction
            LOG_DBG(singlepde) << pdename_  << ": Extrapolating old solutions based on the static solution";
            fncIt= feFunctions_.begin();
            while(fncIt != feFunctions_.end()){
              shared_ptr<BaseFeFunction> actFct = fncIt->second;
              actFct->GetTimeScheme()->ModifyInit(extrapolateStatic);
              fncIt++;
            }
          }

          // Cleanup everything, so that temporary memory needed for domain gets freed
          in.reset();
          // important: This deletes the internal references to the
          inState->Finalize(); 
          inState.reset();
          delete inDomain;
          
        } catch (Exception& e) {
          if( inState ) { 
            inState->Finalize(); 
            inState.reset();
          }
          delete inDomain;
          RETHROW_EXCEPTION(e, "Could not transfer initial state for Physic '" << pdename_ 
                            << "' from sequenceStep " << sequenceStep );
        }

      } else if( srcNode->GetName() == "externalFile" ) {
        EXCEPTION( "No implemented yet")
            
      } else {
        EXCEPTION( "Unknown type of source for initial state.")
      }
      
    } // if initialState
    
    // ===========================
    //  2) Initial condition
    // ===========================
    // Here only specific FeFunction(s) gets initialized, e.g. only the displacement,
    PtrParamNode ifNode = icNode->Get("initialField", ParamNode::PASS );
    if( ifNode ) {
      LOG_TRACE(singlepde) << pdename_ << ": Reading initial condition";
      //get scalar or vector element

      //read which quantity to initialize
      std::string quantityStr = ifNode->Get("quantity",ParamNode::EX)->As<std::string>();
      SolutionType solType = SolutionTypeEnum.Parse(quantityStr);
      shared_ptr<ResultInfo> aResult = this->feFunctions_[solType]->GetResultInfo();

      //get every region which has this ID
      std::string idStr = ifNode->Get("id",ParamNode::PASS)->As<std::string>();
      ParamNodeList regionList = myParam_->Get("regionList")->GetListByVal("region","initialFieldId",idStr);

      for(UInt aNode = 0; aNode < regionList.GetSize(); aNode++){
        // create new entity list
        RegionIdType actRegion = ptGrid_->GetRegion().Parse(regionList[aNode]->Get("name")->As<std::string>());

        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( actRegion );

        //Create a CoefFunction to process the initial field
        PtrCoefFct regionCoef;
        std::set<UInt> definedDofs;
        bool coefUpdateGeo;
        if(aResult->dofNames.GetSize()>1){
          ReadUserFieldValues( actSDList, ifNode->Get("vector"), aResult->dofNames, aResult->entryType,
                               isComplex_, regionCoef, definedDofs, coefUpdateGeo );
        }else{
          ReadUserFieldValues( actSDList, ifNode->Get("scalar"), aResult->dofNames, aResult->entryType,
                                         isComplex_, regionCoef, definedDofs, coefUpdateGeo );
        }
        this->feFunctions_[solType]->AddExternalDataSource(regionCoef,actSDList);
      }

      //Trigger the feFunction to fill itself from the field
      this->feFunctions_[solType]->ApplyExternalData();
    }
    

    ifNode = icNode->Get("initialFieldD1", ParamNode::PASS );
    if( ifNode ) {
      LOG_TRACE(singlepde) << pdename_ << ": Reading initial condition";
      //get scalar or vector element

      //read which quantity to initialize
      std::string quantityStr = ifNode->Get("quantity",ParamNode::EX)->As<std::string>();
      SolutionType solType = SolutionTypeEnum.Parse(quantityStr);
      shared_ptr<ResultInfo> aResult = this->timeDerivFeFunctions_[solType]->GetResultInfo();

      //get every region which has this ID
      std::string idStr = ifNode->Get("id",ParamNode::PASS)->As<std::string>();
      ParamNodeList regionList = myParam_->Get("regionList")->GetListByVal("region","initialFieldD1Id",idStr);

      for(UInt aNode = 0; aNode < regionList.GetSize(); aNode++){
        // create new entity list
        RegionIdType actRegion = ptGrid_->GetRegion().Parse(regionList[aNode]->Get("name")->As<std::string>());

        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( actRegion );

        //Create a CoefFunction to process the initial field
        PtrCoefFct regionCoef;
        std::set<UInt> definedDofs;
        bool coefUpdateGeo;
        if(aResult->dofNames.GetSize()>1){
          ReadUserFieldValues( actSDList, ifNode->Get("vector"), aResult->dofNames, aResult->entryType,
                               isComplex_, regionCoef, definedDofs, coefUpdateGeo );
        }else{
          ReadUserFieldValues( actSDList, ifNode->Get("scalar"), aResult->dofNames, aResult->entryType,
                                         isComplex_, regionCoef, definedDofs, coefUpdateGeo );
        }
        this->timeDerivFeFunctions_[solType]->AddExternalDataSource(regionCoef,actSDList);

      }

      //Trigger the feFunction to fill itself from the field
      this->timeDerivFeFunctions_[solType]->ApplyExternalData();
    }

    LOG_TRACE(singlepde) << pdename_ << ": Finished reading initial conditions";
  }

  void SinglePDE::ReadBCs() {
    LOG_TRACE(singlepde) << pdename_ << ": Reading boundary conditions";


    // fetch "bcsAndLoads" parameter node, if present.
    // otherwise leave
    PtrParamNode bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if( !bcsNode )
      return;

    std::string name, resultName, entType, inputId, value, phase;
    shared_ptr<BaseFeFunction> actFeFunction;

    // =====================================================================
    // homogeneous Dirichlet BC
    // =====================================================================
    // iterate over all available result types
    std::map<SolutionType,std::string>::const_iterator hdbcIt;
    hdbcIt = hdbcSolNameMap_.begin();
    for( ; hdbcIt != hdbcSolNameMap_.end(); ++hdbcIt ) {
      
      // get for each solutiontype the corresponding element name for the
      // homogeneous Dirichlet Bc
      actFeFunction = GetFeFunction( hdbcIt->first );
      std::string elemName = hdbcIt->second;
      shared_ptr<ResultInfo> info = actFeFunction->GetResultInfo();
      StdVector<std::string>  dofNames = info->dofNames;
      // additional check: if we have a vector-valued function approximation,
      // we might have scalar unknowns
      if( actFeFunction->GetFeSpace()->GetNumDofs() == 1 ) {
        dofNames.Clear();
        dofNames.Push_back("_");
      }
      ParamNodeList hdbcNodes = bcsNode->GetList(elemName);

      // iterate over all parameter nodes
      for( UInt i = 0; i < hdbcNodes.GetSize(); i++ ) {
          hdbcNodes[i]->GetValue( "name", name );

          shared_ptr<HomDirichletBc> actBc ( new HomDirichletBc );
          shared_ptr<EntityList> actList;
          switch( ptGrid_->GetEntityType(name) ) {
            case EntityList::NAMED_NODES:
              actList = ptGrid_->GetEntityList( EntityList::NODE_LIST, name); 
              break;
            case EntityList::REGION:
            case EntityList::NAMED_ELEMS:
              actList = ptGrid_->GetEntityList( EntityList::ELEM_LIST, name );
              break;
            case EntityList::NO_TYPE:
              EXCEPTION("No entities with name '" << name << "' known");
              break;
          }
          
          // Read defined excitation
          std::set<UInt> definedDofs;
          PtrCoefFct coef;
          bool coefUpdatedGeo = false;
          ResultInfo::EntryType entryType = ResultInfo::VECTOR;
          if(dofNames.GetSize() == 0 || dofNames.GetSize() == 1 ) {
            entryType = ResultInfo::SCALAR;
          }
          ReadUserFieldValues( actList, hdbcNodes[i], dofNames, entryType,
                               isComplex_, coef, definedDofs,  coefUpdatedGeo );
          
          // ensure, that only the default coordinate system is used
          if( coef->GetCoordinateSystem() ) {
            if( coef->GetCoordinateSystem()->GetName() != "default" ) {
              EXCEPTION( "Dirichlet boundary conditions can only be defined on "
                        << "the default Cartesian system" );
            }
          }
          actBc->entities = actList;
          actBc->result = actFeFunction->GetResultInfo();
          actBc->dofs = definedDofs;
          
          // add definition to feFunction
          actFeFunction->AddHomDirichletBc(actBc);
      } // loop: hdbcs
    } // loop: solutiontypes
    
    
    //=====================================================================
    // inhomogeneous Dirichlet BC
    // =====================================================================
    // iterate over all available result types
    
    // loop over timederivative
    std::map<SolutionType,std::string>::const_iterator idbcIt;
    std::map<SolutionType,std::string>::const_iterator idbcIt_end;
    
    for(UInt timeDeriv = 0; timeDeriv < 3; timeDeriv++){
      if(timeDeriv == 0){
      idbcIt = idbcSolNameMap_.begin();
      idbcIt_end = idbcSolNameMap_.end();
        }else if(timeDeriv == 1){
      idbcIt = idbcSolNameMapD1_.begin();
      idbcIt_end = idbcSolNameMapD1_.end();
        }else if(timeDeriv == 2){
      idbcIt = idbcSolNameMapD2_.begin();
      idbcIt_end = idbcSolNameMapD2_.end();    
        } else {
      EXCEPTION("Max timederiv = 2");
        }
      for( ; idbcIt != idbcIt_end; ++idbcIt ) {

        // get for each solutiontype the corresponding element name for the
        // homogeneous Dirichlet Bc
        actFeFunction = GetFeFunction( idbcIt->first );
        std::string elemName = idbcIt->second;
        shared_ptr<ResultInfo> info = actFeFunction->GetResultInfo();
        StdVector<std::string>  dofNames = info->dofNames;
        // additional check: if we have a vector-valued function approximation,
        // we might have scalar unknowns
  //      if( actFeFunction->GetFeSpace()->GetNumDofs() == 1 ) {
  //        dofNames.Clear();
  //        dofNames.Push_back("_");
  //      }
        ParamNodeList idbcNodes = bcsNode->GetList(elemName);

        // iterate over all parameter nodes
        for( UInt i = 0; i < idbcNodes.GetSize(); i++ ) {
          idbcNodes[i]->GetValue( "name", name );

          shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
          shared_ptr<EntityList> actList;
          switch( ptGrid_->GetEntityType(name) ) {
            case EntityList::NAMED_NODES:
              actList = ptGrid_->GetEntityList( EntityList::NODE_LIST, name); 
              break;
            case EntityList::REGION:
            case EntityList::NAMED_ELEMS:
              actList = ptGrid_->GetEntityList( EntityList::ELEM_LIST, name );
              break;
            case EntityList::NO_TYPE:
              EXCEPTION("No entities with name '" << name << "' known");
              break;
          }

          // Read defined excitation
          std::set<UInt> definedDofs;
          PtrCoefFct coef;
          bool updatedGeo;
          //only needed for MultiHarm IDBC
          PtrCoefFct harm;
          ResultInfo::EntryType entryType = ResultInfo::VECTOR;
          if(dofNames.GetSize() == 0 || dofNames.GetSize() == 1 ) {
            entryType = ResultInfo::SCALAR;
          }

          ReadUserFieldValues( actList, idbcNodes[i], dofNames, entryType,
                                isComplex_, coef, definedDofs, updatedGeo, harm );
          
          // ensure, that only the default coordinate system is used
          if( coef->GetCoordinateSystem() ) {
            if( coef->GetCoordinateSystem()->GetName() != "default" ) {
              EXCEPTION( "Dirichlet boundary conditions can only be defined on "
                  << "the default Cartesian system" );
            }
          }

          actBc->entities = actList;
          actBc->result = actFeFunction->GetResultInfo();
          if( actFeFunction->GetFeSpace()->GetNumDofs() == 1 ) {
            actBc->dofs.insert(0);
          } else {
            actBc->dofs = definedDofs;        
          }

          actBc->value = coef;
          actBc->updatedGeo = updatedGeo;
          actBc->timeDerivOrder = timeDeriv;
          actBc->harm = harm;

          // add definition to feFunction
          actFeFunction->AddInhomDirichletBc(actBc);
        } // loop: idbcs
      } // loop: solutiontypes
    } // loop: timederiv
    
    
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
        csNodes[i]->GetValue( "masterDof", masterDof );
        csNodes[i]->GetValue( "slaveDof", slaveDof );

        // fetch related resultInfo object
        actFeFunction = GetFeFunction( SolutionTypeEnum.Parse(resultName) );

        // Create constraint condition
        shared_ptr<Constraint> actBc ( new Constraint );
        shared_ptr<EntityList> actList;
        switch( ptGrid_->GetEntityType(name) ) {
          case EntityList::NAMED_NODES:
            actList = ptGrid_->GetEntityList( EntityList::NODE_LIST, name); 
            break;
          case EntityList::REGION:
          case EntityList::NAMED_ELEMS:
            actList = ptGrid_->GetEntityList( EntityList::ELEM_LIST, name );
            break;
          case EntityList::NO_TYPE:
            EXCEPTION("No entities with name '" << name << "' known");
            break;
        }

        actBc->masterEntities = actList;
        actBc->slaveEntities = actList;
        actBc->name = name;
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

    // iterate over all parameter nodes
    for( UInt i=0, numNodes=prNodes.GetSize(); i < numNodes; ++i ) {
      try {
        ReadPeriodicBC(prNodes[i]);
      } catch (CoupledField::Exception & ex ) {
        RETHROW_EXCEPTION( ex, "Can not create periodic boundary on '"
                           << prNodes[i]->Get("primary")->As<std::string>() << "/"
                           << prNodes[i]->Get("secondary")->As<std::string>() << "'" );
      }
    }
    LOG_TRACE(singlepde) << pdename_ << ": Finished reading boundary conditions";
  }

  void SinglePDE::ReadPeriodicBC(PtrParamNode prNode) {
    bool allCoordsFree = false;
    Double tol = 1e-6;
    UInt i, nDims, nFixed;
    std::string masterName, slaveName, resultName, dof, coordStr, coordSysId;
    StdVector<UInt> fixedCoords;
    shared_ptr<BaseFeFunction> actFeFunction;
    CoordSystem *coordSys = NULL;

    // Read attribute values
    prNode->GetValue( "primary", masterName );
    prNode->GetValue( "secondary", slaveName );
    prNode->GetValue( "dof", dof, ParamNode::PASS );
    prNode->GetValue( "quantity", resultName );
    prNode->GetValue( "fixedCoords", coordStr, ParamNode::PASS );
    prNode->GetValue( "tolerance", tol, ParamNode::PASS );
    prNode->GetValue( "coordSysId", coordSysId );

    // fetch related resultInfo object
    actFeFunction = GetFeFunction(SolutionTypeEnum.Parse(resultName));
    
    // get entitylists
    NodeList masterList( ptGrid_ ), slaveList( ptGrid_ );
    masterList.SetNamedNodes( masterName );
    slaveList.SetNamedNodes( slaveName );

    // ensure, that both lists have the same length
    if( masterList.GetSize() != slaveList.GetSize() ) {
      EXCEPTION( "Node lists '" << masterName << "' and '"
                 << slaveName << "' have different size" );
    }

    // Get coordinate system
    coordSys = domain_->GetCoordSystem( coordSysId );
    assert( coordSys );
    nDims = coordSys->GetDim();
    
    // Find out which coordinates are fixed
    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";|, ");
    Tok tok(coordStr, sep);
    if ( tok.begin() == tok.end() ) {
      allCoordsFree = true;
      nFixed = 0;
    }
    else {
      for ( Tok::iterator tokIt = tok.begin(); tokIt != tok.end(); ++tokIt ) {
        for ( i = 0; i < nDims; ++i ) {
          if ( coordSys->GetDofName(i+1) == *tokIt ) {
            fixedCoords.Push_back( i );
            break;
          }
        }
        if ( i == nDims ) {
          EXCEPTION("'" << *tokIt
                    << "' is not a valid component name in coordinate system '"
                    << coordSys->GetName() << "'");
        }
      }
      allCoordsFree = fixedCoords.IsEmpty();
      nFixed = fixedCoords.GetSize();
    }
    

    // Brute force algorithm:
    // iterate over all master nodes and try to find "nearest"
    // node in slave list with respect to centers of gravity
    // get center of gravity of bounding box for each list
    Vector<Double> mLoc, sLoc, diff, tmp;
    Vector<Double> mMin (dim_), sMin (dim_), mMax (dim_), sMax (dim_), mCOG, sCOG;
    EntityIterator masterIt = masterList.GetIterator();
    EntityIterator slaveIt = slaveList.GetIterator();
    for( UInt i=0; i<dim_; i++) {
      mMin[i] = 1e42;
      sMin[i] = 1e42;
      mMax[i] = -1e42;
      sMax[i] = -1e42;
    }
    // get bounding box of master nodes
    for( masterIt.Begin(); !masterIt.IsEnd(); masterIt++ ) {
      // obtain nodal coordinate
      ptGrid_->GetNodeCoordinate( mLoc, masterIt.GetNode() );
      for( UInt i=0; i<dim_; i++) {
        if( mLoc[i] < mMin[i]) {
          mMin[i] = mLoc[i];
        }
        if( mLoc[i] > mMax[i] ) {
          mMax[i] = mLoc[i];
        }
      }
    }
    // get bounding box of slave nodes
    for( slaveIt.Begin(); !slaveIt.IsEnd(); slaveIt++ ) {
      // obtain nodal coordinate
      ptGrid_->GetNodeCoordinate( sLoc, slaveIt.GetNode() );
      for( UInt i=0; i<dim_; i++) {
        if( sLoc[i] < sMin[i]) {
          sMin[i] = sLoc[i];
        }
        if( sLoc[i] > sMax[i] ) {
          sMax[i] = sLoc[i];
        }
      }
    }
    mCOG = (mMax + mMin);
    sCOG = (sMax + sMin);
    for( UInt i=0; i<dim_; i++) {
      mCOG[i] = mCOG[i]/2;
      sCOG[i] = sCOG[i]/2;
    }

    Double minDist, dist, minFixed, fixedDiff;
    StdVector<UInt> nodes(2);
    for( masterIt.Begin(); !masterIt.IsEnd(); masterIt++ ) {

      minDist = 1e42;
      minFixed = minDist;

      // obtain nodal coordinate
      ptGrid_->GetNodeCoordinate( mLoc, masterIt.GetNode() );
      if ( !allCoordsFree && coordSysId != "default" ) {
        // if all coordinates are free, the coordinate system is irrelevant
        coordSys->Global2LocalCoord(tmp, mLoc);
        mLoc = tmp;
      }
      
      // initialize node pair
      nodes.Init();
      nodes[0] = masterIt.GetNode();

      // iterate over all slave nodes and find the one with minimum
      // distance
      for( slaveIt.Begin(); !slaveIt.IsEnd(); slaveIt++ ) {
        ptGrid_->GetNodeCoordinate( sLoc, slaveIt.GetNode() );
        if ( !allCoordsFree ) {
          if ( coordSysId != "default" ) {
            // obtain coordinates in given coordinate system
            coordSys->Global2LocalCoord(tmp, sLoc);
            sLoc = tmp;
          }
          // first make sure that fixed coordinates match
          for ( i = 0; i < nFixed; ++i ) {
            // absolute difference
            fixedDiff = abs(mLoc[fixedCoords[i]] - sLoc[fixedCoords[i]]);
            // compute relative difference, if possible (i.e. coordinate != 0)
            if ( mLoc[fixedCoords[i]] > 1e-14 ) fixedDiff /= mLoc[fixedCoords[i]];
            // reject node if it exceeds the tolerance
            if ( fixedDiff > tol )  {
              if ( fixedDiff < minFixed ) {
                // store best guess, so we can display a hint for the
                // tolerance, that would be needed to find a match
                minFixed = fixedDiff;
              }
              break;
            }
          }
          // skip to next slave node, because a fixed coordinate did not match
          if ( i < nFixed ) continue;
        }

        // calculate distance between master and slave node with respect to centers of gravity
        diff = mLoc - sLoc - (mCOG - sCOG);
        dist = diff.NormL2();
        if( dist < minDist ) {
          // store slave node with least distance
          minDist = dist;
          nodes[1] = slaveIt.GetNode();
        }
      }
      
      if ( nodes[1] == 0 ) {
        EXCEPTION("Could not find a matching slave node for master node "
                  << nodes[0] << ".\nRelative deviation of nearest miss: "
                  << minFixed);
      }
      else if ( nodes[0] == nodes[1] ) {
        // Ignore nodes that are contained both in master and slave regions.
        // That usually happens in rotationally periodic setups with nodes, that
        // are located on the axis of rotation.
        continue;
      }
      
      shared_ptr<NodeList> nodePair(new NodeList( ptGrid_ ) );
      nodePair->SetNodes( nodes );

      // create new constraint condition
      shared_ptr<Constraint> actBc ( new Constraint );
      actBc->masterEntities = nodePair;
      actBc->slaveEntities = nodePair;
      if( dof.empty() ) {
        if (actFeFunction->GetResultInfo()->entryType != ResultInfo::SCALAR) {
          EXCEPTION("The 'dof' attribute of the periodic boundary is missing. "
                    << "It is mandatory for vectorial unknowns.");
        }
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
    
  }

  void SinglePDE::ReadMaterialDependency( const std::string& elemName,
                                     const StdVector<std::string>& compNames,
                                     ResultInfo::EntryType type,
                                     bool isComplex,
                                     shared_ptr<EntityList>& entity,
                                     PtrCoefFct& coef,
                                     bool& updateGeo ) {

    // get grip of all elements of that type
    if( !myParam_->Has("matDependencyList") )
      return;

    PtrParamNode xml = myParam_->Get("matDependencyList")->Get(elemName);
    std::set<UInt> definedDofs;
    ReadUserFieldValues(entity,xml,compNames,type,isComplex,coef,
                            definedDofs, updateGeo );
  }

  void SinglePDE::ReadEntities( const std::string& elemName,
                                     const StdVector<std::string>& compNames,
                                     ResultInfo::EntryType type,
                                     StdVector<shared_ptr<EntityList> >& entities,
                                     StdVector<PtrParamNode>& xmls,
                                     StdVector<PtrCoefFct>& coef,
                                     bool& updateGeo,
                                     PtrParamNode input) {

    // get grip of all elements of that type
    if(!input && !myParam_->Has("bcsAndLoads") )
      return;

    ParamNodeList elems = !input ? myParam_->Get("bcsAndLoads")->GetList(elemName) : input->GetList(elemName);

    // necessary for constraints on displacements
    UInt end = 0;
    if (elemName == "displacement_constraint") {
      assert(elems.GetSize() == 1);
      PtrParamNode xml = elems[0];
      // read number of nodes where displacement constraint is applied
      end = xml->Get("multiple_nodes")->As<int>();
    } else {
      end = elems.GetSize();
    }

    // allocate
    entities.Resize(elems.GetSize());
    xmls.Resize(elems.GetSize());
    coef.Resize(elems.GetSize());

    for( UInt i = 0; i < end; ++i ) {
      PtrParamNode xml = elems[i];
      bool hasName = xml->Has("name");
      bool hasRegionList = xml->Has("regionList");

      if (hasName && hasRegionList) {
        EXCEPTION(elemName << " element contains name attribute and regionList element, both are not allowed together");
      } else if (!hasName && !hasRegionList) {
        EXCEPTION(elemName << " element contains neither name attribute nor regionList element, exactly one these is required");
      }

      std::string entName;
      if (hasRegionList) {
        StdVector<PtrParamNode> regs = xml->Get("regionList")->GetList("region");
        if (regs.GetSize() > 1) {
          StdVector<RegionIdType> regionTypes;
          for(UInt r=0;r<regs.GetSize();++r) {
            std::string regName = regs[r]->Get("name")->As<std::string>();
            regionTypes.Push_back(ptGrid_->GetRegion().Parse(regName));
          }
          RegionList* regionList = new RegionList(ptGrid_);
          regionList->SetRegions(regionTypes);
          shared_ptr<EntityList> entList(regionList);
          entities[i] = entList;
        } else if (regs.GetSize() == 1) {
          entName = regs[0]->Get("name")->As<std::string>();
          ptGrid_->GetRegion().Parse(entName);
          hasName = true;
        } else {
          EXCEPTION(elemName << " element contains regionList without any regions");
        }
      }
      if (hasName) {
        // get entity list, depending on type
        entName = xml->Get("name")->As<std::string>();
        try {
          // determine list type: In case we have have surface elements, generate explicitly
          // a surface element list
          EntityList::ListType listType = EntityList::ELEM_LIST;
          if( ptGrid_->GetEntityDim( entName ) == ptGrid_->GetDim() - 1) {
            listType = EntityList::SURF_ELEM_LIST;
          }

          switch( ptGrid_->GetEntityType(entName) ) {
            case EntityList::NAMED_NODES:
              entities[i] = ptGrid_->GetEntityList( EntityList::NODE_LIST, entName);
              break;
            case EntityList::REGION:
            case EntityList::NAMED_ELEMS:
              entities[i] = ptGrid_->GetEntityList( listType, entName );
              break;
            case EntityList::NO_TYPE:
              EXCEPTION("No entities with name '" << entName << "' known");
              break;
          }
        } catch (Exception& e) {
          RETHROW_EXCEPTION(e, pdename_ << ": Could not read definition for '" << elemName
                            << "' on entities '" << entName <<"'");
        }
      }
      std::set<UInt> definedDofs;
      xmls[i] = xml;
    } // loop: elements
  }

  void SinglePDE::ReadRhsExcitation( const std::string& elemName, 
                                     const StdVector<std::string>& compNames,
                                     ResultInfo::EntryType type,
                                     bool isComplex,
                                     StdVector<shared_ptr<EntityList> >& entities, 
                                     StdVector<PtrCoefFct >& coef,
                                     bool& updateGeo,
                                     PtrParamNode input) {
    // read entities and allocate xmls, coef
    StdVector<PtrParamNode> xmls;
    ReadEntities( elemName, compNames, type, entities, xmls, coef, updateGeo, input);
    // read field values
    for( UInt i = 0; i < xmls.GetSize(); ++i ) {
      std::set<UInt> definedDofs;
      ReadUserFieldValues(entities[i],xmls[i],compNames,type,isComplex,coef[i],definedDofs, updateGeo );
    } // loop: elements
  }

  void SinglePDE::ReadRhsExcitation( const std::string& elemName,
                                       const StdVector<std::string>& compNames,
                                       ResultInfo::EntryType type,
                                       StdVector<shared_ptr<EntityList> >& entities,
                                       StdVector<PtrCoefFct >& coef,
                                       bool& updateGeo,
                                       PtrParamNode input) {
    StdVector<PtrParamNode> xmls;
    ReadEntities( elemName, compNames, type, entities, xmls, coef, updateGeo, input);
    for( UInt i = 0; i < xmls.GetSize(); ++i ) {
      std::set<UInt> definedDofs;
      ReadUserFieldValues(entities[i],xmls[i],compNames,type,coef[i],definedDofs, updateGeo );
    } // loop: elements
  }

  void SinglePDE::ReadVolumeRegions( const std::string& elemName, StdVector<std::string>& volumeRegions){
    ParamNodeList elems = myParam_->Get("bcsAndLoads")->GetList(elemName);
    // read the Volume Region from each node
    volumeRegions.Resize(elems.GetSize());
    for( UInt i = 0; i < elems.GetSize(); ++i ) {
      PtrParamNode xml = elems[i];
      std::string volRegName;
      xml->GetValue( "volumeRegion", volRegName, ParamNode::PASS );
      volumeRegions[i] = volRegName;
    }
  }

  void SinglePDE::ReadRhsExcitation( const std::string& elemName,
                                  const StdVector<std::string>& compNames,
                                  ResultInfo::EntryType type,
                                  bool isComplex,
                                  StdVector<shared_ptr<EntityList> >& entities,
                                  StdVector<PtrCoefFct>& coef,
                                  bool& updateGeo,
                                  StdVector<std::string>& volumeRegions){
      // read volume regions
      ReadVolumeRegions(elemName, volumeRegions);
      // read the rest
      ReadRhsExcitation(elemName,compNames,type,isComplex,entities,coef,updateGeo);
  }

  void SinglePDE::ReadRhsExcitation( const std::string& elemName,
                                  const StdVector<std::string>& compNames,
                                  ResultInfo::EntryType type,
                                  StdVector<shared_ptr<EntityList> >& entities,
                                  StdVector<PtrCoefFct>& coef,
                                  bool& updateGeo,
                                  StdVector<std::string>& volumeRegions){
      // read volume regions
      ReadVolumeRegions(elemName, volumeRegions);
      // read the rest
      ReadRhsExcitation(elemName,compNames,type,entities,coef,updateGeo);
  }

  template<typename T>
  void SinglePDE::ReadUserHistValues( PtrParamNode valueNode,
                           ResultInfo::EntryType type,
                           Vector<T>& resV,
                           std::string regionName){

    // some checks
    if( !valueNode->Has("sequenceStep") ) EXCEPTION("History results can only be read for sequence steps!")
    if( type != ResultInfo::EntryType::SCALAR ) EXCEPTION("History result currently needs to be a scalar!")

    PtrParamNode esNode = valueNode->Get("sequenceStep");
    PtrParamNode qNode = esNode->Get("quantity");
    PtrParamNode tfm = esNode->Get("timeFreqMapping");
    std::string tfmString =  tfm->GetChild()->GetName();
    if( tfmString != "constant") EXCEPTION("Only 'constant' interpolation allowed");

    // obtain fileId and SequenceStep
    std::string fileId;
    UInt sequenceStep = 0;
    std::string quantityName = qNode->Get("name")->As<std::string>();
    std::string pdeName = qNode->Get("pdeName")->As<std::string>();
    SolutionType solType = SolutionTypeEnum.Parse(quantityName);

    Domain * inDomain = NULL;
    // create SimState (for input)
    boost::shared_ptr<SimState> inState(new SimState(true, domain_));
    shared_ptr<SimInput> reader;
    shared_ptr<SimInputHDF5> in;
    sequenceStep = esNode->Get("index")->As<UInt>();
    // create new simState from current hdf file
    if( !simState_->GetOutputWriter() ){
      // Sometimes the writer is not yet set if using initial values and external data.
      // Therefore the SimState is instructed to create it now.
      shared_ptr<SimOutputHDF5> writer;
      simState_->SetOutputHdf5Writer( writer );
    }
    std::string fileName = simState_->GetOutputWriter()->GetFileName().string();
    // create new param and info node (without logging to console) for the
    // newly created Domain object
    PtrParamNode node(new ParamNode());
    PtrParamNode infoNode = ParamNode::GenerateWriteNode("", "", ParamNode::APPEND); // empty filename means we don't write and ignore ParamNode::ToFile()
    in.reset(new SimInputHDF5(fileName, node, infoNode));
    inState->SetInputHdf5Reader(in);

    // Get grid map of own domain, as the grids can be re-used
    SimState::GridMap gridMap = domain_->GetGridMap();

    inDomain = inState->GetDomain(sequenceStep, gridMap);

    // Obtain same PDE from new domain
    SinglePDE * inPDE = inDomain->GetSinglePDE(pdeName);

    // remeber input simState and domain
    inputs_[inState] = inDomain;
    // set domain to one specific step
    inState->SetInterpolation( SimState::CONSTANT, mp_, analysistype_ , 0 );
    UInt stepNum = tfm->Get("constant")->Get("step")->As<UInt>();
    inState->UpdateToStep(sequenceStep, stepNum);
    shared_ptr<ResultInfo> rInfo = inPDE->GetResultInfo(solType);
    inPDE->ReadStoreResults();
    ResultMap rMap = inPDE->GetResults();
    ResultList rList = rMap[rInfo];
    shared_ptr<BaseResult> res = rList[0];
    in->SetTempRegionName(regionName);
    in->GetResult(sequenceStep, stepNum, res, true);
    resV = dynamic_cast<Result<T>&>(*res).GetVector();
    inDomain->GetResultHandler()->FinishMultiSequenceStep();
    in->ResetTempRegionName();
  }

  void SinglePDE::ReadUserFieldValues( shared_ptr<EntityList> list,
                                         PtrParamNode valueNode,
                                         const StdVector<std::string>& compNames,
                                         ResultInfo::EntryType type,
                                         PtrCoefFct & coef,
                                         std::set<UInt>& definedDofs,
                                         bool& updateGeo){
    bool isComplex;
    valueNode->GetValue("isComplex", isComplex, ParamNode::PASS );
    ReadUserFieldValues(list, valueNode, compNames, type, isComplex, coef, definedDofs, updateGeo);
  }

  void SinglePDE::
  ReadUserFieldValues( shared_ptr<EntityList> list,
                            PtrParamNode valueNode,
                            const StdVector<std::string>& compNames,
                            ResultInfo::EntryType type,
                            bool isComplex,
                            PtrCoefFct & coef,
                            std::set<UInt>& definedDofs,
                            bool& updateGeo){
    PtrCoefFct harm;
    ReadUserFieldValues(list, valueNode, compNames, type, isComplex, coef, definedDofs, updateGeo, harm);
  }

  void SinglePDE::ReadUserFieldValues( shared_ptr<EntityList> list,
                                       PtrParamNode valueNode,
                                       const StdVector<std::string>& compNames,
                                       ResultInfo::EntryType type,
                                       bool isComplex,
                                       PtrCoefFct & coef,
                                       std::set<UInt>& definedDofs,
                                       bool& updateGeo,
                                       PtrCoefFct & harm){

    UInt numComp = compNames.GetSize();
    StdVector<std::string> vals(numComp), phases(numComp);
    std::string dofString = "";

    vals.Init("0.0");
    phases.Init("0.0");
    definedDofs.clear();

    // switch type of coef function
    if( valueNode->Has("grid") ) {
      // ====================
      //  EXTERNAL GRID DATA 
      // ====================
      shared_ptr<RegionList> regions;
      if (list->GetType() == EntityList::REGION_LIST) {
        regions = boost::static_pointer_cast<RegionList>(list);
      } else {
        RegionList* regionList = new RegionList(ptGrid_);
        regionList->SetRegion(list->GetRegion());
        regions = shared_ptr<RegionList>(regionList);
      }
      if(!isComplex) {
        coef = CoefFunctionGrid::Generate(domain_, Global::REAL, infoNode_ , valueNode->Get("grid"),
                                          regions,type);
        //this is hardcoded so far. should be changed or generated depending on the type
        //of grid (nodal or higher order)
        //coef.reset(new CoefFunctionNodalGrid<Double>(valueNode->Get("grid")));
      } else {
        coef = CoefFunctionGrid::Generate(domain_, Global::COMPLEX, infoNode_ , valueNode->Get("grid"),
                                          regions,type);
        //coef.reset(new CoefFunctionNodalGrid<Complex>(valueNode->Get("grid")));
      }
      //read in the defined dofs
      std::string dofString = valueNode->Get("grid")->Get("dofs")->As<std::string>();
      std::istringstream iss(dofString);
      do{
          string sub;
          iss >> sub;
          if(sub=="all"){
            // add all dofs to the definedDofs
            for( UInt i = 0; i < numComp; ++i ) {
              definedDofs.insert(i);
            }
            break;
          }else{
            UInt index = compNames.Find(sub);
            definedDofs.insert(index);
          }

      } while (iss);
    
      // here we assume no updated geometry
      updateGeo = false;
      
    } else if( valueNode->Has("externalSimulation") ||
               valueNode->Has("sequenceStep") ) {
      // ===========================================
      //  EXTERNAL SIMULATION / MULTISEQUENCE STEP
      // ============================================
      PtrParamNode esNode;
      if( valueNode->Has("externalSimulation") ) {
        esNode = valueNode->Get("externalSimulation");
      } else {
        esNode = valueNode->Get("sequenceStep");
      }
      PtrParamNode qNode = esNode->Get("quantity");
      PtrParamNode tfm = esNode->Get("timeFreqMapping");

      // obtain fileId and SequenceStep
      std::string fileId;
      UInt sequenceStep = 0; 
      std::string quantityName = qNode->Get("name")->As<std::string>();
      std::string pdeName = qNode->Get("pdeName")->As<std::string>();
      SolutionType solType = SolutionTypeEnum.Parse(quantityName);
      
      //read in the defined dofs
      dofString = qNode->Get("dofs")->As<std::string>();
      
      try {
        Domain * inDomain = NULL;

        // create SimState (for input)
        boost::shared_ptr<SimState> inState(new SimState(true, domain_));
        shared_ptr<SimInput> reader;
        shared_ptr<SimInputHDF5> in;

        if( valueNode->Has("externalSimulation") ) {
          sequenceStep = esNode->Get("sequenceStep")->As<UInt>();
          fileId = esNode->Get("inputId")->As<std::string>();
          reader = domain_->GetResultHandler()->GetInputReader(fileId);
          try {
            in = dynamic_pointer_cast<SimInputHDF5>(reader);
          } catch (...) {
            EXCEPTION( "Reader with id'" << fileId << "' has not HDF5 format." );
          }
        } else { 
          
          sequenceStep = esNode->Get("index")->As<UInt>();
          // create new simState from current hdf file
          if( !simState_->GetOutputWriter() ){
            // Sometimes the writer is not yet set if using initial values and external data.
            // Therefore the SimState is instructed to create it now.
            shared_ptr<SimOutputHDF5> writer;
            simState_->SetOutputHdf5Writer( writer );
          }
          std::string fileName = simState_->GetOutputWriter()->GetFileName().string();
          // create new param and info node (without logging to console) for the
          // newly created Domain object
          PtrParamNode node(new ParamNode());
          PtrParamNode infoNode = ParamNode::GenerateWriteNode("", "", ParamNode::APPEND); // empty filename means we don't write and ignore ParamNode::ToFile()
          in.reset(new SimInputHDF5(fileName, node, infoNode));
        }

        inState->SetInputHdf5Reader(in);

        // Get grid map of own domain, as the grids can be re-used
        SimState::GridMap gridMap = domain_->GetGridMap();

        // Obtain temporary Domain object, from which the initial state is read 
        // in. As this generates inferior logging output, we make a visual
        // break.
        LOG_DBG(singlepde) << pdename_ 
            << ": Obtaining Domain from simState object";
        LOG_DBG(singlepde) << pdename_ << ": =================================="; 
        LOG_DBG(singlepde) << pdename_ << ":  BEGIN OUTPUT OF TEMPORARY DOMAIN ";
        LOG_DBG(singlepde) << pdename_ << ": ==================================";            

        inDomain = inState->GetDomain(sequenceStep, gridMap);

        LOG_DBG(singlepde) << pdename_ << ": ==================================="; 
        LOG_DBG(singlepde) << pdename_ << ":  END OF OUTPUT OF TEMPORARY DOMAIN ";
        LOG_DBG(singlepde) << pdename_ << ": ===================================";

        // Obtain same PDE from new domain
        SinglePDE * inPDE = inDomain->GetSinglePDE(pdeName);

        // Check type of interpolation
        std::string tfmString =  tfm->GetChild()->GetName();
        if( tfmString == "constant") {

          // set domain to one specific step
          inState->SetInterpolation( SimState::CONSTANT, mp_, analysistype_ , 0 );
          UInt stepNum = tfm->Get("constant")->Get("step")->As<UInt>();
          inState->UpdateToStep(sequenceStep, stepNum);

        } else if( tfmString == "continuous" ) {

          // read parameters for continuous interpolation
          PtrParamNode contNode = tfm->Get("continuous");
          std::string interpolString = contNode->Get("interpolation")->As<std::string>();
          Double offset = contNode->Get("offset")->As<Double>();
          SimState::InterpolType interpolType = 
              SimState::InterpolTypeEnum.Parse(interpolString);
          inState->SetInterpolation( interpolType, mp_, analysistype_ , offset );

        } else {
          EXCEPTION( "Time / frequency mapping of type '" << tfmString 
                     << "' not known.");
        }


        // Return coefficient function
        coef = inPDE->GetCoefFct(solType);
        if( !coef ) {
          EXCEPTION( "Quantity '" << quantityName << "' is not computable by physic '"
                     << pdeName << "'.")
        }
        
        // remember input simState and domain
        inputs_[inState] = inDomain;
        
        // Check dimensionality of coefficient function
        CoefFunction::CoefDimType dimType = coef->GetDimType();
        if( !( 
            (dimType == CoefFunction::SCALAR && type == ResultInfo::SCALAR ) ||
            (dimType == CoefFunction::VECTOR && type == ResultInfo::VECTOR ) ||
            (dimType == CoefFunction::TENSOR && type == ResultInfo::TENSOR ) ) ) {
          EXCEPTION( "Quantity '" << quantityName << "' is of type '" 
                     << CoefFunction::coefDimType.ToString(coef->GetDimType())
                     << "' but type '" 
                     << ResultInfo::EntryTypeEnum_.ToString(type) 
                     << "' is required." );
        }
      } catch (Exception& e) {
        RETHROW_EXCEPTION(
            e, "Could not obtain quantity '" << quantityName << "' from physic '"
            << pdeName << "' from external simulation with id '" <<  fileId << "'.\n"
            << "Please check, if desired quantity and physic are defined for the "
            << "requested sequence step " << sequenceStep << ".");
      }
      // add the defined components
      for( UInt i = 0; i < numComp; ++i ) {
        definedDofs.insert(i);
      }
      
      
    } else if( valueNode->Has("coupling") ) {
      // ====================
      //  ITERATIVE COUPLING 
      // ====================
      // read pdeName and quantity
      std::string pdeName, quantity;
      PtrParamNode cplNode =valueNode->Get("coupling"); 
      cplNode->GetValue("pdeName", pdeName);
      cplNode->Get("quantity")->GetValue("name", quantity);
      
      // try to access SinglePDE and acquire result from there
      if( !iterCplPde_ ) {
        EXCEPTION( "Can not get quantity '" << quantity << "' from physic '"
                   << pdeName << "', as no coupling object is defined.");
      }

      // check if we can parse the quantity (by passing the false flag, we just get an
      // invalid solution type but not an error)
      // we do this before checking for generic postProcessing results (although they might not be defined already)
      // since they COULD be defined already be a PDE read before this one
      SolutionType solType = SolutionTypeEnum.TryParse(quantity, INVALID_SOLUTION_TYPE);
      if( solType == INVALID_SOLUTION_TYPE ) {
        // direct conversion was not successful, check if we can add the result based on postProcessing results

        // get the postProcList and access all names of type function
        PtrParamNode ppListNode = myParam_->GetParent()->GetParent()->Get("postProcList",ParamNode::PASS);
        std::string funcName;

        if( ppListNode ) {
          ParamNodeList ppListNodeChildren = ppListNode->GetChildren();

          // loop over all postProc definitions and check, if a function is defined
          for( UInt i = 0; i < ppListNodeChildren.GetSize(); i++ ) {

            // get all children of one postProc definition
            ParamNodeList ppNodeChildren = ppListNodeChildren[i]->GetChildren();
            
            // we only consider the function type postProcResults
            ParamNodeList postProcFuncs = ParamNode::GetListByChild(ppNodeChildren, "function");

            // loop over all function-type children and check their name
            for( UInt u = 0; u < postProcFuncs.GetSize(); u++ ) {
              postProcFuncs[u]->GetValue("resultName",funcName);

              if( funcName==quantity ) {
                // the quantity we need is defined as a postProcResult but not yet parseable
                // add it to the environment and get a usable solution type
                AddGenericSolution(quantity, domain_);
              }
            }
          }
        }
        // finally, try to parse the quantity again
        // if it does not work, the quantity is not defined (neither pre-defined nor via the postProcList)
        solType = SolutionTypeEnum.Parse(quantity);  
      } 

      // now we know that we have a quantity for sure (not an invalid one)
      // in the next step, we have to check if they are generic, because then we
      // still need to define a dummy coefFunction since the real one will not be
      // available yet
      if( solType == GENERIC_RESULT_0 || solType == GENERIC_RESULT_1 ||
          solType == GENERIC_RESULT_2 || solType == GENERIC_RESULT_3 ||
          solType == GENERIC_RESULT_4 || solType == GENERIC_RESULT_5 ||
          solType == GENERIC_RESULT_6 || solType == GENERIC_RESULT_7 ||
          solType == GENERIC_RESULT_8 || solType == GENERIC_RESULT_9 ) {

        // For generic results, we temporarely define a coefFunctionDummy
        // where the feFunction is temporarly set during initialization.
        // Afterwards, we replace it with the actual coefFunction when 
        // the postProc result is actually defined.
        // We have to do this since the initialization procedure calls 
        // this in stage 2, while we get the postproc stuff in stage 3

        coef.reset(new CoefFunctionDummy(solType, list, pdeName));
        // we have to set the coordsys to pass the check below
        std::string coordSysId = "default";
        coef->SetCoordinateSystem( domain_->GetCoordSystem(coordSysId) );

        // set the geometry update acoording to the feFunction
        iterCplPde_->GetUpdateGeoForPDE(solType, list, pdeName, updateGeo);

        // append to the vector of coefs which will get updated later on
        coefsToUpdate_.push_back(coef);
      
      } else {
        // the result is already known (non-generic)
        // --> define the coefFunction
        coef = iterCplPde_->GetCouplingCoefFct(solType, list, pdeName, updateGeo);
      }
      
      // add the defined components
      for( UInt i = 0; i < numComp; ++i ) {
        definedDofs.insert(i);
      }

      //read in the defined dofs
      dofString = cplNode->Get("quantity")->Get("dofs")->As<std::string>();

    }
    else if(valueNode->Has("scatteredData"))
    {
      PtrParamNode scatteredDataNode = valueNode->Get("scatteredData");
      if(type == ResultInfo::SCALAR){
        if(isComplex)
          coef.reset(new CoefFunctionScatteredData<Complex, 1>(scatteredDataNode));
        else
          coef.reset(new CoefFunctionScatteredData<Double, 1>(scatteredDataNode));
      }
      else if(type == ResultInfo::VECTOR){
        if( dim_ == 2 )
        {
          if(isComplex)
            coef.reset(new CoefFunctionScatteredData<Complex, 2>(scatteredDataNode));
          else
            coef.reset(new CoefFunctionScatteredData<Double, 2>(scatteredDataNode));
        } else  {
          if(isComplex)
            coef.reset(new CoefFunctionScatteredData<Complex, 3>(scatteredDataNode));
          else
            coef.reset(new CoefFunctionScatteredData<Double, 3>(scatteredDataNode));
        }
      }
      else
        EXCEPTION("TENSOR not implemented yet!");
    }
    else if(valueNode->Has("fileData"))
    {
      if(type == ResultInfo::TENSOR)
        EXCEPTION("TENSOR not implemented yet!");

      PtrParamNode pnfd = valueNode->Get("fileData");
      int mydim = type == ResultInfo::SCALAR ? 1 : dim_;
      if (compNames[0] == "Hx" || compNames[0] == "Bx") {
        coef.reset(new CoefFunctionFileDataMeas(domain_, pdename_, pnfd, mydim));  
      } else {
        coef.reset(new CoefFunctionFileData(pnfd, mydim));
      }
    }
    else if(valueNode->Has("python"))
    {
      unsigned int mydim = type == ResultInfo::SCALAR ? 1 : dim_;
      coef.reset(new CoefFunctionPython(valueNode->Get("python"), mydim));
    }

    else
    {
      // ======================================
      //  STANDARD EXPLICIT BOUNDARY CONDITION 
      // ======================================
      // Note: In case someone request a "vector" valued result and
      // provides no dofNames, we use the scalar parameters.
      if( type == ResultInfo::SCALAR  ||
          (type == ResultInfo::VECTOR && compNames.GetSize() <= 1 )  ) {
        // --------------
        //  S C A L A R
        // --------------
        // in the scalar case, the component with index 0 is always defined
        definedDofs.insert(0);
        std::string val = "0.0";
        std::string phase = "0.0";
        std::string harmonic = "0";
        // allow for optional value / phase attributes 
        valueNode->GetValue("value", val, ParamNode::PASS );
        valueNode->GetValue("phase", phase, ParamNode::PASS );
        //for multiharmonic case
        valueNode->GetValue("harmonic", harmonic, ParamNode::PASS);
        std::string real = AmplPhaseToReal(val, phase );
      
        if( type == ResultInfo::SCALAR) {
          // -- SCALAR case --
          if(!isComplex ) {
            coef = CoefFunction::Generate(mp_, Global::REAL, real);
          } else {
            std::string imag = AmplPhaseToImag(val, phase );
            coef = CoefFunction::Generate(mp_, Global::COMPLEX, real, imag);
            harm = CoefFunction::Generate(mp_, Global::INTEGER, harmonic);
          }
        }else {
          // -- VECTOR case --
          // generate coefficient function
          StdVector<std::string> realV, imagV;
          realV = real;
          imagV =  AmplPhaseToImag(val, phase );
          if(!isComplex) {
            StdVector<std::string> valV;
            valV = val;
            coef = CoefFunction::Generate(mp_, Global::REAL, valV );
          } else {
            coef = CoefFunction::Generate(mp_, Global::COMPLEX, realV, imagV);
          }
        }
      }
      else if (type == ResultInfo::VECTOR)
      {
        CoordSystem * coordSys = NULL;
        std::string coordSysId = "default";
        valueNode->GetValue("coordSysId", coordSysId, ParamNode::PASS);
        if( coordSysId != "default" ) {
          coordSys = domain_->GetCoordSystem(coordSysId);
        }
        
        // --------------
        //  V E C T O R
        // --------------
        // a) all values are given as vector
        if( valueNode->Has("values") ) {
          std::string valString = valueNode->Get("values")->As<std::string>();
          SplitStringList(valString, vals, ' ');

          // consistency check
          if( vals.GetSize() != numComp ) {
            EXCEPTION("Boundary condition needs " << numComp << " values!");
          }

          // check for phase vector (optional)
          if( valueNode->Has("phase")) {
            std::string phaseString = valueNode->Get("phase")->As<std::string>();
            SplitStringList(phaseString, phases, ' ');

            // consistency check
            if(phases.GetSize() != numComp )
              EXCEPTION("Boundary condition needs " << numComp << " phase values!");
          }
          
          // add all dofs to the definedDofs
          for(UInt i = 0; i < compNames.GetSize(); ++i )
            definedDofs.insert(i);
        }
        // b) values are given component-wise
        else if(valueNode->Has("comp")) {
          ParamNodeList compList = valueNode->GetList("comp");
          Integer index = 0;
          std::string dof;
          std::string val = "0.0";
          std::string phase = "0.0";
          for( UInt j = 0; j < compList.GetSize(); ++j  ) {
           compList[j]->GetValue("dof", dof);
           compList[j]->GetValue("value", val, ParamNode::PASS );
           compList[j]->GetValue("phase", phase, ParamNode::PASS);

           // find index
           if( compNames.GetSize() == 0 ) {
             index = 0;
             definedDofs.insert(0);
           } else {
             // try to map found component name to coordinate-system local one
             if(coordSys) {
               index = coordSys->GetVecComponent(dof)-1; 
             } else {
               index = compNames.Find(dof);
             }
             
             if( index == -1 ) {
               EXCEPTION("Could not find component with name '" << dof << "'");
             }
             definedDofs.insert(UInt(index));
           }

           vals[index] = val;
           phases[index] = phase;
          } // loop components

        } else {
          // This branch will only be handled for homogeneous values, i.e. it
          // is okay to have no value given. In this case
          EXCEPTION( "No values given for boundary condition '"
              << valueNode->GetParent()->GetName() << "' on '" 
              << list->GetName() << "'" );
        }

        // generate coefficient function
        StdVector<std::string> real, imag;
        AmplPhaseToRealImag(vals, phases, real, imag);
        if(!isComplex) {
          coef = CoefFunction::Generate(mp_, Global::REAL, vals );
        } else {
          coef = CoefFunction::Generate(mp_, Global::COMPLEX, real, imag );
        }

      } else if (type == ResultInfo::TENSOR ) {
        // --------------
        //  T E N S O R
        // --------------
        EXCEPTION("Not yet implemented for tensor-valued boundary conditions");
        //    - all defined: read in both vectors
        //    - components:

      }
      
      // explicitly defined boundary conditions are assumed to comply with
      // the formulation of the PDE, i.e. if the PDE has total Lagrangian 
      // formulation, so will be all the BCs and source terms.
      updateGeo = updatedGeo_;
    }


    // obtain coordinate system and set it at coefficient function
    std::string coordSysId = "default";
    valueNode->GetValue("coordSysId", coordSysId, ParamNode::PASS);
    if( coordSysId != "default" ) {
      coef->SetCoordinateSystem( domain_->GetCoordSystem(coordSysId) );
    }

    //defined in xml, set dof names
    coef->SetDofNames(dofString);

    // return 
  }

  void SinglePDE::ReadMaterialData() {
    UInt i, numRegions;
    
    // get list of parameter nodes for region definitions
    ParamNodeList regionNodes;

    PtrParamNode regionListNode = domain_->GetParamRoot()->
      Get("domain")->Get("regionList", ParamNode::PASS );
    if( regionListNode)
      regionNodes = regionListNode->GetList("region");

    numRegions = regionNodes.GetSize();

    // obtain pointer to materialHandler
    MaterialHandler* matLoader = domain_->GetMaterialHandler();


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
        RegionIdType actRegionId = ptGrid_->GetRegion().Parse( region );

        // if no material is set, continue with next loop run
        if( material.empty() )
          continue;

        // if region is not contained for current pde, simply continue
        // with next loop
        if( regions_.Find( actRegionId) < 0 )
          continue;

        // Read data
        materials_[actRegionId] = matLoader->LoadMaterial(material, pdematerialclass_);

        // log the just read material. LoadMaterial() so to say initializes the ToInfo()
        PtrParamNode in = infoNode_->GetByVal("material", "name", material);
        // additional regions are automatically appended
        in->Get("regionList")->GetByVal("region", "name", ptGrid_->GetRegion().ToString(actRegionId));
        materials_[actRegionId]->ToInfo(in);

        // Check for local coordinate system
        if( !refCoordSys.empty() ) {
          CoordSystem * actCoosy =
            domain_->GetCoordSystem( refCoordSys);
          materials_[actRegionId]->SetCoordSys( actCoosy );
        }

        // Check for material rotation parameters
        PtrParamNode rotNode = regionNodes[i]->Get("matRotation", ParamNode::INSERT );

        Vector<Double> rotVec (3);
        rotVec.Init();

        // NOTE: If no rotation is specified and the dimension is
        // 2D, -> material is rotated by
        // alpha = -90 and gamma = -90 degree,
        // so that we pick by default the yz-plane
        if( !rotNode->HasChildren() ) {
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
        RegionIdType actRegionId = ptGrid_->GetRegion().Parse( region );

        // if no composite is set, continue with next loop run
        if( composite == "" )
          continue;
        
        if( regions_.Find( actRegionId) < 0 )
          continue;

        PtrParamNode in = infoNode_->GetByVal("composite", "region", 
                                              ptGrid_->GetRegion().ToString(actRegionId));

        // get composite node
        PtrParamNode compNode;
        try {
          compNode = domain_->GetParamRoot()->Get("domain")->GetByVal("composite", "name", composite);
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
    numRegions = regions_.GetSize();
    std::map<RegionIdType, BaseMaterial*>::iterator matEnd = materials_.end();
    std::map<RegionIdType, Composite>::iterator
        compEnd = compositeMaterials_.end();
    
    for( i = 0; i < numRegions; ++i ) {
      RegionIdType actRegionId = regions_[i];
      if ((materials_.find(actRegionId) == matEnd)
          && (compositeMaterials_.find(actRegionId) == compEnd)) {
        region = ptGrid_->GetRegion().ToString(actRegionId);
        EXCEPTION("Region '" << region << "' has no material assigned.");
      }
    }
  }

  void SinglePDE::ReadNcInterfaces() {
    PtrParamNode nciListNode = myParam_->Get("ncInterfaceList", ParamNode::PASS);
    if ( !nciListNode )
      return;
    ParamNodeList nciNodes = nciListNode->GetList("ncInterface");
    ParamNodeList::iterator nciIt = nciNodes.Begin(),
                            endIt = nciNodes.End();
    // iterate over the defined interfaces...
    for ( ; nciIt != endIt; ++nciIt ) {
      PtrParamNode nciNode = (*nciIt);
      NcInterfaceInfo newIface;
      // get the interface ids and types
      newIface.interfaceId = ptGrid_->GetNcInterfaceId( nciNode->Get("name")->As<std::string>() );
      newIface.type = ncCouplingType_.Parse( nciNode->Get("formulation", ParamNode::INSERT)->As<std::string>() );
      // check for lagrange multiplyer
      if (newIface.type == NC_MORTAR) {
        // actually the lagrangeMultType us unused!!
        newIface.lagrangeMultType = LM_STANDARD;
      } else {
        // get nitsche factor
        nciNode->GetValue( "nitscheFactor", newIface.nitscheFactor, ParamNode::INSERT );
      }
      if (nciNode->Has( "thinLayer")){ // check if thin layer element is within xml file
        newIface.thinLayer = true;
        nciNode->GetValue( "thinLayer/layerThickness", newIface.layerThickness, ParamNode::INSERT ); // parameter for thin layer formulation non conforming interface condition
        nciNode->GetValue( "thinLayer/layerMaterial", newIface.layerMaterial, ParamNode::INSERT ); // parameter for thin layer formulation non conforming interface condition
        if(newIface.type != NC_NITSCHE){
          EXCEPTION("Thin layer formulation is only available with Nitsche type coupling");
        }
      }
      else {
        newIface.thinLayer = false;
      }
      ncInterfaces_.Push_back(newIface);
    }
  }

  // ======================================================
  // GET /SET  METHODS
  // ======================================================

  //! Activate the direct coupling
  void SinglePDE::SetDirectCoupling () {
    isDirectCoupled_ = true;
  }


  void SinglePDE::UpdateToSolStrategy() {
    
    // this is hopefully a general way to update all information related
    // to a step update in a multistep solution process
    
    // update all feFunctions 
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt = 
        feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      FeSpace* actSpace = fncIt->second->GetFeSpace().get();
      actSpace->UpdateToSolStrategy();
      fncIt++;
    }
  }

  
  void SinglePDE::DefineFieldResult( PtrCoefFct coef, shared_ptr<ResultInfo> res ) {

    LOG_DBG(singlepde) << pdename_ << ": Defining field result " << SolutionTypeEnum.ToString(res->resultType);

    // create new result functor based on coefficient
    shared_ptr<ResultFunctor> func;
    if( isComplex_ ) {
      func.reset(new FieldCoefFunctor<Complex>(coef, res));
    } else {
      func.reset(new FieldCoefFunctor<Double>(coef, res));
    }
    
    // insert result to list of available results and field functors
    resultFunctors_[res->resultType] = func;
    fieldCoefs_[res->resultType] = coef;
    availResults_.insert(res);

    // define the averaged result of the field result
    shared_ptr<ResultInfo> avgResInfo(new ResultInfo);
    avgResInfo->definedOn = ResultInfo::REGION_AVERAGE; // on the region, because we integrate it over the volume
    // copy the rest
    avgResInfo->resultType = res->resultType ;
    avgResInfo->dofNames = res->dofNames;
    avgResInfo->unit = res->unit;
    avgResInfo->entryType = res->entryType;
    //avgResInfo->name = ;
    availResults_.insert( avgResInfo );
    if ( res->GetFeFunction().expired() && (feFunctions_.size() <= 1) ) { // feFunction was not set to result info in PDE, but there is only one anyway
      res->SetFeFunction( feFunctions_.rbegin()->second );//we set it
    }
    if ( res->GetFeFunction().expired() ) {
      LOG_DBG(singlepde) << "Developer Info: use SetFeFct() to make averaged results work for '" << SolutionTypeEnum.ToString(res->resultType) <<"'";
    } else {
      shared_ptr<BaseFeFunction> feFct = shared_ptr<BaseFeFunction>(res->GetFeFunction());
      shared_ptr<ResultFunctor> avgResFunctor;
      if(isComplex_) {
        avgResFunctor.reset(new ResultFunctorIntegrate<Complex>(coef, feFct, avgResInfo));
        dynamic_pointer_cast< ResultFunctorIntegrate<Complex> >(avgResFunctor)->SetAveraged(true);
      }
      else {
        avgResFunctor.reset(new ResultFunctorIntegrate<Double>(coef, feFct, avgResInfo));
        dynamic_pointer_cast< ResultFunctorIntegrate<Double> >(avgResFunctor)->SetAveraged(true);
      }
      fieldAverageFunctors_[avgResInfo->resultType] = avgResFunctor;
    }
  }
  
  void SinglePDE::DefineFieldResult(SolutionType solType, ResultInfo::EntryType entryType, ResultInfo::EntityUnknownType definedOn, const std::string& dofNames, bool fromOptimization)
  {
    shared_ptr<ResultInfo> ri(new ResultInfo);
    ri->resultType = solType;
    ri->entryType = entryType;
    ri->definedOn = definedOn;
    ri->dofNames = dofNames;
    ri->fromOptimization = fromOptimization;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), ri);
  }

  void SinglePDE::SetFieldCoef( PtrCoefFct coef, SolutionType resultType)
  {
    fieldCoefs_[resultType] = coef;
  }

  void SinglePDE::DefineTimeDerivResult( SolutionType derivSolType,
                                         UInt timeDerivOrder,
                                         SolutionType primSolType ) {
    
    // only define time derivatives in transient or harmonic case 
    if ( analysistype_ == STATIC || analysistype_ == BUCKLING){
      return;
    }
      
    // get primary feFunction
    shared_ptr<BaseFeFunction> primFeFct = feFunctions_[primSolType];
    
    // obtain resultinfos
    shared_ptr<ResultInfo> primRes = GetResultInfo( primSolType );
    shared_ptr<ResultInfo> derivRes = GetResultInfo( derivSolType);
    
    // ensure that both result types could be found
    if( !primRes ) {
      EXCEPTION( "Could not find result information for result type '"
          << SolutionTypeEnum.ToString(primSolType) << "'" );
    }
    if( !derivRes ) {
      EXCEPTION( "Could not find result information for result type '"
          << SolutionTypeEnum.ToString(derivSolType) << "'" );
    }
    
    // create new fe function
    shared_ptr<BaseFeFunction> derivFeFct;
    if( isComplex_) {
      derivFeFct.reset(new FeFunction<Complex>(domain_->GetMathParser()));
    }  else {
      derivFeFct.reset(new FeFunction<Double>(domain_->GetMathParser()));
    }
    // copy information (entitylists, grid, space, fctId)
    StdVector< shared_ptr<EntityList> > entList = primFeFct->GetEntityList();
    for(UInt i =0;i<entList.GetSize();i++){
      derivFeFct->AddEntityList(entList[i]);
    }
    derivFeFct->SetFeSpace( primFeFct->GetFeSpace());
    derivFeFct->SetGrid(ptGrid_);
    derivFeFct->SetResultInfo(derivRes);
    derivFeFct->SetFctId( primFeFct->GetFctId());
    //derivFeFct->Finalize();
    
    // Note:
    // The initialization of the time derivative function has to be performed
    // at a later step, as the primary FeFct and the TimeStepping algorithm
    // are not initialized yet (see SinglePDE::FinalizeInit()); 
    
    
    // insert the fefunction to the time derivative section
    timeDerivFeFunctions_[derivSolType] = derivFeFct;
    timeDerivPrimaryResults_[derivSolType] = primRes->resultType;
    timeDerivOrder_[derivSolType] = timeDerivOrder;
    
    // define field result
    DefineFieldResult( derivFeFct, derivRes );
    
  }

  void SinglePDE::DefineFeFunctions(){
    //This is the default creation of spaces
    //idee: die PDE gibt zum attribute formulation die passenden space zurueck
    //DOGMA: PRO UNBEKANNTE EINE FUNCTION UND EIN SPACE
    std::string formulation;
    myParam_->GetValue("feSpaceFormulation",formulation,ParamNode::EX);
    PtrParamNode feSpaceNode = infoNode_->Get("feSpaces");
    std::map<SolutionType, shared_ptr<FeSpace> > spaces = CreateFeSpaces(formulation, feSpaceNode);
    CreateNcFeSpaces(spaces, formulation, feSpaceNode);
    
    //loop over all spaces and set an FeFunction
    std::map<SolutionType, shared_ptr<FeSpace> >::iterator spIt = spaces.begin();
    while(spIt != spaces.end()){

      if(feFunctions_.find(spIt->first) != feFunctions_.end()){
        EXCEPTION("It seems that the PDE has created multiple spaces for one result: " << \
            spIt->first << " This is not how its ought to be!");
      }

      if( isComplex_ ) {
        feFunctions_[spIt->first].reset(new FeFunction<Complex>(domain_->GetMathParser()));
        rhsFeFunctions_[spIt->first].reset(new FeFunction<Complex>(domain_->GetMathParser()));
      }else{
        feFunctions_[spIt->first].reset(new FeFunction<Double>(domain_->GetMathParser()));
        rhsFeFunctions_[spIt->first].reset(new FeFunction<Double>(domain_->GetMathParser()));
        // Note: in the transient case, we also need fefunctions for the time derivatives
      }
      //let the objects know about each other
      spIt->second->AddFeFunction(feFunctions_[spIt->first]);
      feFunctions_[spIt->first]->SetFeSpace(spIt->second);
      feFunctions_[spIt->first]->SetPDE(this);
      feFunctions_[spIt->first]->SetGrid(ptGrid_);

      rhsFeFunctions_[spIt->first]->SetFeSpace(spIt->second);
      rhsFeFunctions_[spIt->first]->SetPDE(this);
      rhsFeFunctions_[spIt->first]->SetGrid(ptGrid_);
      spIt++;
    }
  }
  
  void SinglePDE::CreateNcFeSpaces(std::map<SolutionType, shared_ptr<FeSpace> > &spaces,
                                   const std::string &formulation, PtrParamNode infoNode) {
    // create FeSpaces for Lagrange multipliers of ncInterfaces (if any)
    StdVector<NcInterfaceInfo>::iterator ncIt = ncInterfaces_.Begin(),
                                         endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      if ( ncIt->type == NC_MORTAR ) {
        FeSpace::SpaceType type;
        if ( formulation == "default" ) {
          type = FeSpace::H1;
        } else {
          type = FeSpace::SpaceTypeEnum.Parse(formulation);
          if ( type != FeSpace::H1 ) {
            EXCEPTION("Mortar ncInterfaces currently support only H1 elements");
          }
        }
        spaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, infoNode, type, ptGrid_);
        spaces[LAGRANGE_MULT]->SetLagrSurfSpace();
        spaces[LAGRANGE_MULT]->Init(solStrat_);
        break; // One FeSpace for the Lagrange multiplier is enough, so exit loop
        // TODO: This might change, when we allow different FeSpaces for
        // different regions of the same PDE.
      }
    }
  }

  void SinglePDE::DefineNcAuxResults() {
    // create results for Lagrange multipliers of ncInterfaces (if any)
    StdVector<NcInterfaceInfo>::iterator ncIt = ncInterfaces_.Begin(),
                                                endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      if ( ncIt->type == NC_MORTAR ) {
        shared_ptr<ResultInfo> lm( new ResultInfo() );
        lm->resultType = LAGRANGE_MULT;
        lm->complexFormat = results_[0]->complexFormat;
        lm->definedOn = results_[0]->definedOn;
        lm->dofNames = results_[0]->dofNames;
        lm->entryType = results_[0]->entryType;

        feFunctions_[LAGRANGE_MULT]->SetResultInfo(lm);
        lm->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
        results_.Push_back(lm);
        DefineFieldResult(feFunctions_[LAGRANGE_MULT], lm);
        
        shared_ptr<ResultInfo> lm_d1( new ResultInfo() );
        lm_d1->resultType = LAGRANGE_MULT_DERIV_1;
        lm_d1->complexFormat = lm->complexFormat;
        lm_d1->definedOn = lm->definedOn;
        lm_d1->dofNames = lm->dofNames;
        lm_d1->entryType = lm->entryType;
        availResults_.insert(lm_d1);
        DefineTimeDerivResult(lm_d1->resultType, 1, lm->resultType);
        
        shared_ptr<ResultInfo> lm_d2( new ResultInfo() );
        lm_d2->resultType = LAGRANGE_MULT_DERIV_2;
        lm_d2->complexFormat = lm->complexFormat;
        lm_d2->definedOn = lm->definedOn;
        lm_d2->dofNames = lm->dofNames;
        lm_d2->entryType = lm->entryType;
        availResults_.insert(lm_d2);
        DefineTimeDerivResult(lm_d2->resultType, 2, lm->resultType);
        
        break; // one result is enough for all ncInterfaces
      }
    }
  }
  
  template<UInt DIM, UInt D_DOF>
  void SinglePDE::DefineMortarCoupling(SolutionType solType,
                                       NcInterfaceInfo &iface)
  {
    // Get interface from grid and cast to MortarInterface class
    shared_ptr<BaseNcInterface> ncIf =
        ptGrid_->GetNcInterface(iface.interfaceId);
    MortarInterface *mortarIf = dynamic_cast<MortarInterface*>(&(*ncIf));
    assert(mortarIf);
    
    // currently we have a moving formulation only for acoustics
    updatedGeo_ = updatedGeo_ || ncIf->IsMoving();

    // create ElemLists for slave surface and intersection
    shared_ptr<SurfElemList> elMaster(new SurfElemList(ptGrid_)),
                             elSlave(new SurfElemList(ptGrid_));
    elMaster->SetRegion(mortarIf->GetPrimarySurfRegion());
    elSlave->SetRegion(mortarIf->GetSecondarySurfRegion());
    shared_ptr<NcSurfElemList> elMortar = mortarIf->GetElemList();
    
    // make sure there are appropriate FeFunctions for both unknowns
    if ( feFunctions_.find(solType) == feFunctions_.end() ) {
      EXCEPTION("FeFunction of " << SolutionTypeEnum.ToString(solType)
                << " not found");
    }
    if ( feFunctions_.find(LAGRANGE_MULT) == feFunctions_.end() ) {
      EXCEPTION("FeFunction of Lagrange multiplier not found");
    }
    
    // add element lists to FeFunctions
    feFunctions_[solType]->AddEntityList(elMaster);
    feFunctions_[solType]->AddEntityList(elSlave);
    feFunctions_[LAGRANGE_MULT]->AddEntityList(elSlave);
    
    // For transient case we need to initialize the TimeStepping
    // of the Lagrange multiplier
    if ( analysistype_ == TRANSIENT ) {
      TimeSchemeGLM* tsSol = dynamic_cast<TimeSchemeGLM*>(
          feFunctions_[solType]->GetTimeScheme().get());
      assert( tsSol );
      shared_ptr<TimeSchemeGLM> tsCopy( new TimeSchemeGLM( *tsSol ));
      feFunctions_[LAGRANGE_MULT]->SetTimeScheme( tsCopy );
    }
    
    // Set the same approximation for the Lagrange multiplier as for the
    // primary unknown in the slave region
    shared_ptr<FeSpace> mortarSpace = feFunctions_[LAGRANGE_MULT]->GetFeSpace();
    std::string slaveVolName = ptGrid_->GetRegion()
        .ToString(mortarIf->GetSecondaryVolRegion());
    PtrParamNode slaveRegNode = myParam_->Get("regionList", ParamNode::EX)
        ->GetByVal("region", "name", slaveVolName, ParamNode::EX);
    std::string polyId = slaveRegNode->Get("polyId")->As<std::string>();
    std::string integId = slaveRegNode->Get("integId")->As<std::string>();
    
    mortarSpace->SetRegionApproximation( mortarIf->GetSecondarySurfRegion(),
                                         polyId, integId);

    // create a mass integrator on the slave surface (conforming grid)
    PtrCoefFct unity = CoefFunction::Generate( mp_, Global::REAL, "1" );
    BiLinearForm *massInt = nullptr;
    massInt = new BBInt<>( new IdentityOperator<FeH1,DIM,D_DOF>, unity, 1.0, updatedGeo_);
    massInt->SetName("MortarMassInt");

    // create a context to put the mass integrator into the stiffness matrix
    BiLinFormContext *massContext = new BiLinFormContext(massInt, STIFFNESS);
    massContext->SetEntities(elSlave, elSlave);
    massContext->SetFeFunctions( feFunctions_[solType], feFunctions_[LAGRANGE_MULT] );
    massContext->SetCounterPart(true);
    assemble_->AddBiLinearForm(massContext);
    
    // create a non-conforming mass integrator
    BiLinearForm *ncInt = nullptr;
    ncInt = new SurfaceMortarABInt<>( new IdentityOperator<FeH1,DIM,D_DOF>(),
        new IdentityOperator<FeH1,DIM,D_DOF>(),
        unity, -1.0,
        mortarIf->GetPrimaryVolRegion(),
        mortarIf->GetSecondaryVolRegion(),
        mortarIf->IsCoplanar(),
        updatedGeo_ );
    ncInt->SetName("MortarNcInt");
    
    NcBiLinFormContext *ncContext = new NcBiLinFormContext(ncInt, STIFFNESS);
    ncContext->SetEntities(elMortar, elMortar);
    ncContext->SetFeFunctions( feFunctions_[solType], feFunctions_[LAGRANGE_MULT] );
    ncContext->SetCounterPart(true);
    assemble_->AddBiLinearForm(ncContext);
    ncIf->RegisterIntegrator(ncContext);
    ncContext->SetMotion(ncIf->IsMoving());

    // check for eulerian formulation of moving grid
    DefineEulerianSystem<DIM, D_DOF>(solType, iface);
  }
  
  template<UInt DIM, UInt D_DOF>
  void SinglePDE::DefineNitscheCoupling( SolutionType solType,
                                         NcInterfaceInfo &iface,
                                         shared_ptr<CoefFunctionMulti> additionalCoef )
  {
    shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(iface.interfaceId);
    MortarInterface * nitscheIf = dynamic_cast<MortarInterface*>(ncIf.get());
    assert(nitscheIf);
    
    //in case of Nitsche coupling edge/face information is required
    this->ptGrid_->MapEdges();
    this->ptGrid_->MapFaces();

    // currently we have a moving formulation only for acoustics
    updatedGeo_ = updatedGeo_ || ncIf->IsMoving();

    // create new entity list
    shared_ptr<ElemList> actSDList = ncIf->GetElemList();

    //we set here the penalty factor
    Double beta = iface.nitscheFactor;

    Double assignedFactor = 0.0; // used to set beta / factors within integrator definition below
    Double alphaHeat = 0.0;
    bool isPenalty = true;

    // check if thin layer formulation is set, read and calculate parameters
    bool isThinLayer = false;
    if (iface.thinLayer){
      isThinLayer = true;
      Double layerThickness = iface.layerThickness; // not yet used
      string layerMaterial = iface.layerMaterial; // not yet used
      if (solType == HEAT_TEMPERATURE){
        // Read material data
        MaterialHandler* matLoader = domain->GetMaterialHandler();
        BaseMaterial* mat = matLoader->LoadMaterial(layerMaterial, THERMIC);
        std::map<MaterialClass, std::map<MaterialType, PtrCoefFct> > tl_materials;
        tl_materials[THERMIC][HEAT_CONDUCTIVITY_SCALAR] = mat->GetScalCoefFnc( HEAT_CONDUCTIVITY_SCALAR, Global::REAL );
        double tl_heatConductivity = std::stod(tl_materials[THERMIC][HEAT_CONDUCTIVITY_SCALAR]->ToString());
        alphaHeat = tl_heatConductivity / layerThickness; // calculated value for heat transfer coefficient
        //WARN("Thin layer formulation used with alphaHeat = " << alphaHeat);
      }
      else
        EXCEPTION("Thin layer is currently implemented only for HeatPDE"); // Handle not implemented formulation
    }


    //possible material parameter and adaption of penalty term
    PtrCoefFct factor;
    if ( solType == HEAT_TEMPERATURE && !isThinLayer) {
      factor = materials_[nitscheIf->GetPrimaryVolRegion()]->GetScalCoefFnc( HEAT_CONDUCTIVITY_SCALAR, Global::REAL );
    }
    else if ( solType == ELEC_POTENTIAL && pdename_  != "elecQuasistatic") {
      if(additionalCoef){
        factor = materials_[nitscheIf->GetPrimaryVolRegion()]->GetScalCoefFnc( MAG_CONDUCTIVITY_SCALAR, Global::REAL );
      }else{
        factor = materials_[nitscheIf->GetPrimaryVolRegion()]->GetScalCoefFnc( ELEC_CONDUCTIVITY_SCALAR, Global::REAL );
      }
    }
    else if ( solType == MAG_POTENTIAL) {
      //TODO Clean this up
      PtrCoefFct permeability, reluctivity, permeabilityM, permeabilityS, factorM, factorS, factorAdd;
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      PtrCoefFct constOneC = CoefFunction::Generate( mp_, Global::COMPLEX, "1.0");
      PtrCoefFct constTwo = CoefFunction::Generate( mp_, Global::REAL, "2.0");

      if(additionalCoef){
      // Per convention, master is the linear and slave the nonlinear region
      // Get master and slave permeabilities
      //permeabilityM = materials_[nitscheIf->GetPrimaryVolRegion()]->GetScalCoefFnc( MAG_PERMEABILITY, Global::REAL );
      // Compute the reluctivities
      //factorM = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, constOne, permeabilityM, CoefXpr::OP_DIV));
      //factorS = additionalCoef;
      // And compute the average
      //factorAdd = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp( mp_, factorM, factorS, CoefXpr::OP_ADD ) );
      //factor = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, additionalCoef, constTwo, CoefXpr::OP_DIV));

      reluctivity = additionalCoef->GetRegionCoef(nitscheIf->GetSecondaryVolRegion());
      factor = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOneC, reluctivity, CoefXpr::OP_MULT));


      //permeability = materials_[nitscheIf->GetSecondaryVolRegion()]->GetScalCoefFnc( MAG_PERMEABILITY, Global::REAL );
      //factor = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOneC, permeability, CoefXpr::OP_DIV));

      }else{
        permeability = materials_[nitscheIf->GetSecondaryVolRegion()]->GetScalCoefFnc( MAG_PERMEABILITY_SCALAR, Global::REAL );
        factor = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, constOne, permeability, CoefXpr::OP_DIV));
      }

    }
    else if ( solType == ACOU_PRESSURE || solType == ACOU_POTENTIAL ) {
      EXCEPTION("Nitsche coupling for acoustics is now implemented in the AcousticPDE class!!!");
    }
    else if ( solType == ELEC_POTENTIAL && pdename_  == "elecQuasistatic") {
      factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      if ( isMaterialComplex_ ) {
        // ----- alpha = electricConductivity + i*omega*permittivity
        PtrCoefFct omega = CoefFunction::Generate( mp_,  Global::COMPLEX, "0", "2*pi*f");
        
        //get characteristic conductivity
        shared_ptr<CoefFunction > condCoef;
        condCoef = materials_[nitscheIf->GetPrimaryVolRegion()]->GetTensorCoefFnc(ELEC_CONDUCTIVITY_TENSOR, FULL,
                                                                                 Global::REAL);
        StdVector< Matrix<Double> > mat;
        Double matVal = 0.0;
        StdVector<Vector<Double> > points(1);
        Vector<Double> p1(DIM);
        p1.Init();
        points[0]= p1;
        condCoef->GetTensorValuesAtCoords(points, mat, this->ptGrid_);
        for (UInt i = 0, numRows = mat[0].GetNumRows(); i < numRows; ++i) {
          matVal += mat[0][i][i];
        }
        matVal /= (Double) mat[0].GetNumRows();
        std::string strCond = std::to_string(matVal);
        PtrCoefFct constCond = CoefFunction::Generate( mp_, Global::REAL, strCond);

        //get characteristic permittivity
        shared_ptr<CoefFunction > permCoef;
        permCoef = materials_[nitscheIf->GetPrimaryVolRegion()]->GetTensorCoefFnc(ELEC_PERMITTIVITY_TENSOR, FULL, 
                                                                                 Global::REAL);
        p1.Init();
        points[0]= p1;
        permCoef->GetTensorValuesAtCoords(points, mat, this->ptGrid_);
        matVal = 0.0;
        for (UInt i = 0, numRows = mat[0].GetNumRows(); i < numRows; ++i) {
          matVal += mat[0][i][i];
        }
        matVal /= (Double) mat[0].GetNumRows();
        std::string strPerm  = std::to_string(matVal);
        PtrCoefFct constPerm = CoefFunction::Generate( mp_, Global::REAL, strPerm);
        
        //shared_ptr<CoefFunction > matCoef;
        shared_ptr<CoefFunction > omegaPermCoef;
        omegaPermCoef = CoefFunction::Generate( mp_, Global::COMPLEX,
                                                CoefXprBinOp(mp_, constPerm, omega, CoefXpr::OP_MULT));
        factor = CoefFunction::Generate( mp_, Global::COMPLEX,
                                         CoefXprBinOp(mp_, constCond, omegaPermCoef, CoefXpr::OP_ADD));
        }
    }
    else
      factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");


    //notation> assume the test function is called v
    BiLinearForm *penalty_v1_u1 = nullptr;
    BiLinearForm *penalty_v1_u2 = nullptr;
    BiLinearForm *penalty_v2_u2 = nullptr;
    //now bilinear forms related to the normal derivatives
    //du1 refers to the normal derivative directing from 1 to 2
    BiLinearForm *flux_dv1_u1 = nullptr;
    BiLinearForm *flux_dv1_u2 = nullptr;
    BiLinearForm *flux_v1_du1 = nullptr;
    BiLinearForm::CouplingDirection curcpl;

     // in case of mechanical PDE, we need the material tensor
    shared_ptr<CoefFunction > coefMech;
    SubTensorType tensorType = NO_TENSOR;
    if ( solType == MECH_DISPLACEMENT ) {
      if ( subType_ == "3d" )
        tensorType = FULL;
      else if ( subType_ == "axi" )
        tensorType = AXI;
      else if ( subType_ == "planeStrain" )
        tensorType = PLANE_STRAIN;
      else if ( subType_ == "planeStress" )
        tensorType = PLANE_STRESS;

      //get correct scaling of penalty term
      Double matVal = 0.0;
      StdVector<Vector<Double> > points(1);
      Vector<Double> p1(DIM);
      p1.Init();
      points[0]= p1;

      if ( isMaterialComplex_ ) {
        coefMech = materials_[nitscheIf->GetPrimaryVolRegion()]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, Global::COMPLEX);

        StdVector< Matrix<Complex> > mat;
        coefMech->GetTensorValuesAtCoords(points, mat, this->ptGrid_);

        for (UInt i = 0, numRows = mat[0].GetNumRows(); i < numRows; ++i) {
          matVal += mat[0][i][i].real();
        }

        matVal /= (Double) mat[0].GetNumRows();
      }
      else {
        coefMech = materials_[nitscheIf->GetPrimaryVolRegion()]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, Global::REAL);

        StdVector< Matrix<Double> > mat;
        coefMech->GetTensorValuesAtCoords(points, mat, this->ptGrid_);

        for (UInt i = 0, numRows = mat[0].GetNumRows(); i < numRows; ++i) {
          matVal += mat[0][i][i];
        }

        matVal /= (Double) mat[0].GetNumRows();
      }

      //scale the penalty value
      beta *= matVal;
    }

    // NOTE: the algebraic system sets the system matrix to
    // nonSym  if any bilinear form with the same fctID1 and fctID2 is nonSym.
    // We set here the symmetric flag to true in the constructor
    // of the SurfaceNitscheABInt even though the bilinear form itself is
    // not symmetric. Nitsche formulation is basically sym due to the
    // set counterpart directive for the context.

    // NOTE: when using edge elements, we overwrite the D_DOF, when handing it
    // over to the surface-operator because edge elements have vectorial shape
    // functions BUT scalar DOF's!! Therefore we treat it here as a vector valued
    // DOF and perform the scalar multiplication later on in the CalcElemMatrix of
    // the integrator


    // first, assign BOperators that include solution u1 and and test funciton v1 of the primary side (upper-left part of the coupling matrix)
    curcpl = BiLinearForm::PRIM_PRIM;
    // primary-primary part of the penalty BOperator (v1*u1)
    if ( isMaterialComplex_) {
      penalty_v1_u1 = new SurfaceNitscheABInt<Complex,Complex>
                     (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      factor, beta, curcpl, updatedGeo_, true, true);
    } else {
      if(pdename_ == "magneticEdge" || pdename_ == "magneticEdgeMixedAV" || pdename_ == "magneticEdgeSpecialAV") {
        if(additionalCoef) {
          // multiharmonic case
          penalty_v1_u1 = new SurfaceNitscheABInt<Complex,Complex>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, beta, curcpl, updatedGeo_, true, true);
        } else {
          penalty_v1_u1 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, beta, curcpl, updatedGeo_, true, true);
        }
      } else {
        if(isThinLayer) { // Nitsche term with thinLayer formulaiton
          assignedFactor = alphaHeat; 
          isPenalty = false;
        } else {
          assignedFactor = beta;
          isPenalty = true;
        }
        penalty_v1_u1 = new SurfaceNitscheABInt<Double,Double>
                       (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                        new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                        factor, assignedFactor, curcpl, updatedGeo_, true, isPenalty);
      }
    }
    // primary-primary of the symmetrization term, (dv1/dn*u1)
    if ( solType == MECH_DISPLACEMENT ) {
      flux_dv1_u1 = new SurfaceNitscheABInt<Double,Double>
                   (new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
                    new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                    factor, -1.0, curcpl, updatedGeo_, true);
                    flux_dv1_u1->SetBCoefFunctionOpA(coefMech);
    } else {
      if (isMaterialComplex_) {
        flux_dv1_u1 = new SurfaceNitscheABInt<Complex,Complex>
                     (new SurfaceNormalDerivOperator<FeH1,DIM,D_DOF>(),
                      new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      factor, -1.0, curcpl, updatedGeo_, true);
      } else {
        if(pdename_ == "magneticEdge" || pdename_ == "magneticEdgeMixedAV" || pdename_ == "magneticEdgeSpecialAV") {
          if(additionalCoef) {
            // multiharmonic case
            flux_dv1_u1 = new SurfaceNitscheABInt<Complex,Complex>
                         (new SurfaceCurlNormalOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, -1.0, curcpl, updatedGeo_, true);
          } else {
            flux_dv1_u1 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceCurlNormalOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, -1.0, curcpl, updatedGeo_, true);
          }
        } else {
          if(!isThinLayer) {
            flux_dv1_u1 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceNormalDerivOperator<FeH1,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                          factor, -1.0, curcpl, updatedGeo_, true);
          }
        }
      }
    }

    // primary-primary part of the flux term, (v1*du1/dn)
    if ( solType == MECH_DISPLACEMENT ) {
      flux_v1_du1 = new SurfaceNitscheABInt<Double,Double>
                   (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                    new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
                    factor, -1.0, curcpl, updatedGeo_, true);
      flux_v1_du1->SetBCoefFunctionOpB(coefMech);
    } else {
      if ( isMaterialComplex_) {
        flux_v1_du1 = new SurfaceNitscheABInt<Complex,Complex>
                     (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      new SurfaceNormalDerivOperator<FeH1,DIM,D_DOF>(),
                      factor, -1.0, curcpl, updatedGeo_, true);
      } else {
        if(pdename_ == "magneticEdge" || pdename_ == "magneticEdgeMixedAV" || pdename_ == "magneticEdgeSpecialAV") {
          if(additionalCoef) {
            // multiharmonic case
            flux_v1_du1 = new SurfaceNitscheABInt<Complex,Complex>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceCurlNormalOperator<FeHCurl,DIM,D_DOF>(),
                          factor, -1.0, curcpl, updatedGeo_, true);
          } else {
            flux_v1_du1 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceCurlNormalOperator<FeHCurl,DIM,D_DOF>(),
                          factor, -1.0, curcpl, updatedGeo_, true);
          }
        } else {
          if(!isThinLayer) {
            flux_v1_du1 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                          new SurfaceNormalDerivOperator<FeH1,DIM,D_DOF>(),
                          factor, -1.0, curcpl, updatedGeo_, true);
          }
        }
      }
    }

    // (upper-right part and lower-left part of the coupling matrix)
    curcpl = BiLinearForm::PRIM_SEC;
    // mixed part of penalty term (-v1*u2)
    if (isMaterialComplex_) {
      penalty_v1_u2 = new SurfaceNitscheABInt<Complex,Complex>
                     (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      factor, beta * -1.0, curcpl, updatedGeo_, true, true);
    } else {
      if(pdename_ == "magneticEdge" || pdename_ == "magneticEdgeMixedAV" || pdename_ == "magneticEdgeSpecialAV") {
        if(additionalCoef){
          // multiharmonic case
          penalty_v1_u2 = new SurfaceNitscheABInt<Complex,Complex>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                factor, beta * -1.0, curcpl, updatedGeo_, true, true);
        } else {
          penalty_v1_u2 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                factor, beta * -1.0, curcpl, updatedGeo_, true, true);
        }
      } else {
        if(isThinLayer) {
          assignedFactor = alphaHeat * -1.0; // thin layer formulation
          isPenalty = false;
        } else {
          assignedFactor = beta * -1.0;
          isPenalty = true;
        }
        penalty_v1_u2 = new SurfaceNitscheABInt<Double,Double>
                       (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                        new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                        factor, assignedFactor, curcpl, updatedGeo_, true, isPenalty);
      }
    }
    // mixed part of symmetrization term (dv1/dn*u2)
    if ( solType == MECH_DISPLACEMENT ) {
      flux_dv1_u2 = new SurfaceNitscheABInt<Double,Double>
                   (new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
                    new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                    factor, 1.0, curcpl, updatedGeo_, true);
      flux_dv1_u2->SetBCoefFunctionOpA(coefMech);
    } else {
      if ( isMaterialComplex_) {
        flux_dv1_u2 = new SurfaceNitscheABInt<Complex,Complex>
                     (new SurfaceNormalDerivOperator<FeH1,DIM,D_DOF>(),
                      new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      factor, 1.0, curcpl, updatedGeo_, true);
      } else {
        if(pdename_ == "magneticEdge" || pdename_ == "magneticEdgeMixedAV" || pdename_ == "magneticEdgeSpecialAV") {
          if(additionalCoef) {
            // multiharmonic case
            flux_dv1_u2 = new SurfaceNitscheABInt<Complex,Complex>
                         (new SurfaceCurlNormalOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, 1.0, curcpl, updatedGeo_, true);
          } else {
            flux_dv1_u2 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceCurlNormalOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, 1.0, curcpl, updatedGeo_, true);
          }

        } else {
          if(!isThinLayer) {
            flux_dv1_u2 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceNormalDerivOperator<FeH1,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                          factor, 1.0, curcpl, updatedGeo_, true);
          }
        }
      }
    }

    // finally, we need to define the penalty term living purely on the secondary domain (lower-right part of the coupling matrix)
    // (no fluxes here as the secondary fluxes were expressed via the primary ones)
    curcpl = BiLinearForm::SEC_SEC;
    if (isMaterialComplex_) {
      penalty_v2_u2 = new SurfaceNitscheABInt<Complex,Complex>
                     (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                      factor, beta, curcpl, updatedGeo_, true, true);
    } else {
      if(pdename_ == "magneticEdge" || pdename_ == "magneticEdgeMixedAV" || pdename_ == "magneticEdgeSpecialAV"){
        if(additionalCoef) {
          // multiharmonic case
          penalty_v2_u2 = new SurfaceNitscheABInt<Complex,Complex>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, beta, curcpl, updatedGeo_, true, true);
        } else {
          penalty_v2_u2 = new SurfaceNitscheABInt<Double,Double>
                         (new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          new SurfaceIdentityOperator<FeHCurl,DIM,D_DOF>(),
                          factor, beta, curcpl, updatedGeo_, true, true);
        }
      } else {
        if (isThinLayer) {
          assignedFactor = alphaHeat; //thin layer formulation
          isPenalty = false;
        } else {
          assignedFactor = beta;
          isPenalty = true;
        }
        penalty_v2_u2 = new SurfaceNitscheABInt<Double,Double>
                       (new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                        new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
                        factor, assignedFactor, curcpl, updatedGeo_, true, isPenalty);
      }
    }

    // Nitsche coupling matrix is a type of stiffness matrix
    FEMatrixType targetMatrix = STIFFNESS;
    if(ncIf->IsMoving()){
      targetMatrix = STIFFNESS_UPDATE;
    }

    // now the BOperators are set, so define the contexts...
    SurfaceBiLinFormContext *penalty_v1_u1_Context = nullptr;
    SurfaceBiLinFormContext *flux_dv1_u1_Context   = nullptr;
    SurfaceBiLinFormContext *flux_v1_du1_Context   = nullptr;
    SurfaceBiLinFormContext *penalty_v2_u2_Context = nullptr;
    SurfaceBiLinFormContext *penalty_v1_u2_Context = nullptr;
    SurfaceBiLinFormContext *flux_dv1_u2_Context   = nullptr;
    curcpl = BiLinearForm::PRIM_PRIM;
    penalty_v1_u1_Context = new SurfaceBiLinFormContext(penalty_v1_u1, targetMatrix, curcpl);
    if (flux_dv1_u1){ flux_dv1_u1_Context = new SurfaceBiLinFormContext(flux_dv1_u1  , targetMatrix, curcpl);}
    if (flux_v1_du1){ flux_v1_du1_Context = new SurfaceBiLinFormContext(flux_v1_du1  , targetMatrix, curcpl);}
    curcpl = BiLinearForm::SEC_SEC;
    penalty_v2_u2_Context = new SurfaceBiLinFormContext(penalty_v2_u2, targetMatrix, curcpl);
    curcpl = BiLinearForm::PRIM_SEC;
    penalty_v1_u2_Context = new SurfaceBiLinFormContext(penalty_v1_u2, targetMatrix, curcpl);
    if (flux_dv1_u2){ flux_dv1_u2_Context = new SurfaceBiLinFormContext(flux_dv1_u2  , targetMatrix, curcpl);}
    curcpl = BiLinearForm::SEC_PRIM;
    // assign motion to the contexts
    penalty_v1_u1_Context->SetMotion(ncIf->IsMoving());
    penalty_v2_u2_Context->SetMotion(ncIf->IsMoving());
    penalty_v1_u2_Context->SetMotion(ncIf->IsMoving());
    if (flux_dv1_u1_Context){ flux_dv1_u1_Context->SetMotion(ncIf->IsMoving());}
    if (flux_v1_du1_Context){ flux_v1_du1_Context->SetMotion(ncIf->IsMoving());}
    if (flux_dv1_u2_Context){ flux_dv1_u2_Context->SetMotion(ncIf->IsMoving());}
    // assign names to the operators
    penalty_v2_u2->SetName("penalty_v2_u2");
    penalty_v1_u2->SetName("penalty_v1_u2");
    penalty_v1_u1->SetName("penalty_v1_u1");
    if (flux_dv1_u1){ flux_dv1_u1->SetName("flux_dv1_u1");}
    if (flux_v1_du1){ flux_v1_du1->SetName("flux_v1_du1");}
    if (flux_dv1_u2){ flux_dv1_u2->SetName("flux_dv1_u2");}

    penalty_v1_u1_Context->SetEntities(actSDList,actSDList);
    penalty_v2_u2_Context->SetEntities(actSDList,actSDList);
    penalty_v1_u2_Context->SetEntities(actSDList,actSDList);
    if (flux_dv1_u1_Context){ flux_dv1_u1_Context->SetEntities(actSDList,actSDList);}
    if (flux_v1_du1_Context){ flux_v1_du1_Context->SetEntities(actSDList,actSDList);}
    if (flux_dv1_u2_Context){ flux_dv1_u2_Context->SetEntities(actSDList,actSDList);}

    penalty_v1_u1_Context->SetFeFunctions( feFunctions_[solType], feFunctions_[solType] );
    penalty_v2_u2_Context->SetFeFunctions( feFunctions_[solType], feFunctions_[solType] );
    penalty_v1_u2_Context->SetFeFunctions( feFunctions_[solType], feFunctions_[solType] );
    if (flux_dv1_u1_Context){ flux_dv1_u1_Context->SetFeFunctions( feFunctions_[solType], feFunctions_[solType] );}
    if (flux_v1_du1_Context){ flux_v1_du1_Context->SetFeFunctions( feFunctions_[solType], feFunctions_[solType] );}
    if (flux_dv1_u2_Context){ flux_dv1_u2_Context->SetFeFunctions( feFunctions_[solType], feFunctions_[solType] );}

    penalty_v1_u2_Context->SetCounterPart(true);
    if (flux_dv1_u2_Context){ flux_dv1_u2_Context->SetCounterPart(true);}

    assemble_->AddBiLinearForm( penalty_v1_u1_Context );
    assemble_->AddBiLinearForm( penalty_v2_u2_Context );
    assemble_->AddBiLinearForm( penalty_v1_u2_Context );
    if (flux_dv1_u1_Context){ assemble_->AddBiLinearForm( flux_dv1_u1_Context );}
    if (flux_v1_du1_Context){ assemble_->AddBiLinearForm( flux_v1_du1_Context );}
    if (flux_dv1_u2_Context){ assemble_->AddBiLinearForm( flux_dv1_u2_Context );}

    ncIf->RegisterIntegrator( penalty_v1_u1_Context );
    ncIf->RegisterIntegrator( penalty_v2_u2_Context );
    ncIf->RegisterIntegrator( penalty_v1_u2_Context );
    if (flux_dv1_u1_Context){ ncIf->RegisterIntegrator( flux_dv1_u1_Context );}
    if (flux_v1_du1_Context){ ncIf->RegisterIntegrator( flux_v1_du1_Context );}
    if (flux_dv1_u2_Context){ ncIf->RegisterIntegrator( flux_dv1_u2_Context );}

    // check for eulerian formulation of moving grid
    DefineEulerianSystem<DIM, D_DOF>(solType, iface);
  }


  template<UInt DIM, UInt D_DOF>
  void SinglePDE::DefineEulerianSystem(SolutionType solType, NcInterfaceInfo &iface) {
    shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(iface.interfaceId);
    MortarInterface* ncInterface = dynamic_cast<MortarInterface*>(ncIf.get());
    assert(ncInterface);

    // check for vailidity
    if ((ncIf->IsMoving() && ncInterface->IsEulerian()) != true)
      return;
    if (IsComplex() == true)
      EXCEPTION("Eulerian formulation is only available for real-valued problems.")
    if (solType != ACOU_PRESSURE && solType != ACOU_POTENTIAL)
      EXCEPTION("Eulerian formulation is only tested for acoustics. Please consider adding a testcase to enable this feature for your PDE!")

    std::cout << "Setting up convective terms for Eulerian formulation of moving grid." << std::endl;
    // get the moving region from the mortar interface
    RegionIdType movingRegion = ncInterface->GetMovingVolRegion();
    // get the grid velocity
    PtrCoefFct gridVelocity = ncInterface->GetGridVelocity();
    // create the coefficient of the mass integrator
    PtrCoefFct massCoef = this->massInts_[movingRegion]->GetCoef();
    shared_ptr<ElemList> movingELems(new ElemList(this->ptGrid_));
    movingELems->SetRegion(movingRegion);
    // create bilinear forms
    BiLinearForm *eulerDampInt = new ABInt<Double>(new IdentityOperator<FeH1,DIM,D_DOF>(),
                                                    new ConvectiveOperator<FeH1,DIM,D_DOF>(),
                                                    massCoef, 1.0, true);
    BiLinearForm *eulerDampIntTrans = new ABInt<Double>(new ConvectiveOperator<FeH1,DIM,D_DOF>(),
                                                        new IdentityOperator<FeH1,DIM,D_DOF>(),
                                                        massCoef, -1.0, true);
    BiLinearForm *eulerStiffInt = new BBInt<Double>(new ConvectiveOperator<FeH1,DIM,D_DOF>(),
                                                    massCoef, -1.0, true);
    // assign coefFunction
    eulerDampInt->SetBCoefFunctionOpB(gridVelocity);
    eulerDampIntTrans->SetBCoefFunctionOpA(gridVelocity);
    eulerStiffInt->SetBCoefFunctionOpB(gridVelocity);
    // set names
    eulerDampInt->SetName("eulerDampingInt");
    eulerDampIntTrans->SetName("eulerDampingIntTransposed");
    eulerStiffInt->SetName("eulerStiffnessInt");
    // create bilinear contexts
    BiLinFormContext *eulerDampContext = new BiLinFormContext(eulerDampInt, DAMPING);
    BiLinFormContext *eulerDampContextTrans = new BiLinFormContext(eulerDampIntTrans, DAMPING);
    BiLinFormContext *eulerStiffContext = new BiLinFormContext(eulerStiffInt, STIFFNESS);
    // add the element list to the contexts
    eulerDampContext->SetEntities(movingELems, movingELems);
    eulerDampContextTrans->SetEntities(movingELems, movingELems);
    eulerStiffContext->SetEntities(movingELems, movingELems);
    // add the fe functions to the contexts
    eulerDampContext->SetFeFunctions(feFunctions_[solType], feFunctions_[solType]);
    eulerDampContextTrans->SetFeFunctions(feFunctions_[solType], feFunctions_[solType]);
    eulerStiffContext->SetFeFunctions(feFunctions_[solType], feFunctions_[solType]);
    // add the contexts to the assembler
    assemble_->AddBiLinearForm(eulerDampContext);
    assemble_->AddBiLinearForm(eulerDampContextTrans);
    assemble_->AddBiLinearForm(eulerStiffContext);
  }


  template void SinglePDE::ReadUserHistValues(PtrParamNode, ResultInfo::EntryType, Vector<Double>&, std::string);
  template void SinglePDE::ReadUserHistValues(PtrParamNode, ResultInfo::EntryType, Vector<Complex>&, std::string);
  template void SinglePDE::DefineMortarCoupling<2,1>(SolutionType,NcInterfaceInfo&);
  template void SinglePDE::DefineMortarCoupling<2,2>(SolutionType,NcInterfaceInfo&);
  template void SinglePDE::DefineMortarCoupling<3,1>(SolutionType,NcInterfaceInfo&);
  template void SinglePDE::DefineMortarCoupling<3,3>(SolutionType,NcInterfaceInfo&);
  template void SinglePDE::DefineNitscheCoupling<2,1>(SolutionType,NcInterfaceInfo&, shared_ptr<CoefFunctionMulti> additionalCoef);
  template void SinglePDE::DefineNitscheCoupling<2,2>(SolutionType,NcInterfaceInfo&, shared_ptr<CoefFunctionMulti> additionalCoef);
  template void SinglePDE::DefineNitscheCoupling<3,1>(SolutionType,NcInterfaceInfo&, shared_ptr<CoefFunctionMulti> additionalCoef);
  template void SinglePDE::DefineNitscheCoupling<3,3>(SolutionType,NcInterfaceInfo&, shared_ptr<CoefFunctionMulti> additionalCoef);
  template void SinglePDE::DefineEulerianSystem<2,1>(SolutionType,NcInterfaceInfo&);
  template void SinglePDE::DefineEulerianSystem<2,2>(SolutionType,NcInterfaceInfo&);
  template void SinglePDE::DefineEulerianSystem<3,1>(SolutionType,NcInterfaceInfo&);
  template void SinglePDE::DefineEulerianSystem<3,3>(SolutionType,NcInterfaceInfo&);
} // end of namespace
