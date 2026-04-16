// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "StdPDE.hh"

#include "MatVec/Vector.hh"
#include "MatVec/generatematvec.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TransientDriver.hh"
#include "Driver/MultiSequenceDriver.hh"

#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "DataInOut/SimState.hh"

// headers for Paramhandling
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/Assemble.hh"

// header for logging
#include "DataInOut/Logging/LogConfigurator.hh"

#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"

#include "Utils/Timer.hh"

namespace CoupledField {

// declare logging stream
  DEFINE_LOG(stdPde, "stdPde")

  StdPDE::StdPDE(Grid *aptgrid, PtrParamNode paramNode, PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain ) :
    BasePDE(paramNode, infoNode, simState, domain),
    ptGrid_(aptgrid),
    subType_(),
    nonLin_(false),
    nonLinMaterial_(false),
    nonLinTotalFormulation_(false),
    isHysteresis_(false),
    modelName_("nonlinearCurve"),
    matDepend_(false),
    nonLinMethod_(FIXEDPOINT),
    pdematerialclass_(NO_CLASS),
    isIterCoupled_(false),
    diagMass_(false),
    needsAlgsys_(true),
    analysistype_(BasePDE::NO_ANALYSIS),
    dim_(ptGrid_->GetDim()), 
    isaxi_(ptGrid_->IsAxi()),
    isComplex_(false),    
    assemble_(NULL),
    solveStep_(NULL),
    algsys_(NULL)
  {
  }
  
  StdPDE::~StdPDE()
  {

    // delete all Fe functions explicitly
    feFunctions_.clear();
    prevFeFunctions_.clear();
    timeDerivFeFunctions_.clear();
    rhsFeFunctions_.clear();
    
  }

  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  void StdPDE::CreateMatrices_Solver() {
    
    
    // ==============================
    // create and initialize matrices 
    // ==============================

    // create algebraic system and intialize matrices
    algsys_->CreateLinSys();
    algsys_->InitMatrix();
    
    // Check for analysistype
    if ( analysistype_ != EIGENFREQUENCY && analysistype_ != BUCKLING && analysistype_ != EIGENVALUE) {
      
      // create solver and preconditioner
      algsys_->CreateSolver();
      algsys_->CreatePrecond();
      
    } else {
      // create eigenvalue solver
      algsys_->CreateEigenSolver();
    }
        
    // now reset AlgebraicSystem 
    algsys_->InitRHS();
  	algsys_->InitSol();
  }
  
  // ======================================================
  // ALGSYS SECTION (SOLVER, ...)
  // ======================================================
  void StdPDE::DefineAlgSys() 
  {
    LOG_TRACE(stdPde) << pdename_ << ": Defining Algsys";

    // Immediately leave if external simState is provided
    if (simState_->HasInput())
      return;
    
    // Trigger writing of info file
    myInfo_->GetRoot()->ToFile("", true );
    
    // First check if the PDE needs an algebraic system at all
    if( needsAlgsys_ == false ) {
      return;
    }

    // defining the graph can be expensive. Most expensive are (ordered by cost):
    // - AlgebraicSys::GraphSetupDone()
    // - Assemble::SetupMatrixGraph()
    // when we initialize the BaseGraph instances below, they are implicitly started 
   
    // ==============================================
    //   DEFINE GRAPH AND SBM BLOCKS
    // ==============================================
    // vector containing SBM-block definitions (length: numSbmBlocks)
    StdVector<AlgebraicSys::SBMBlockDef > sbmBlocks;

    // Introduce new mapping:
    //typedef StdVector<std::pair<FeFctIdType, Integer> > MinorBlockDef;
    //std::map<UInt,StdVector <MinorBlockDef > > minorBlocks;
    //          ^        ^
    //          |        |
    //        sbm   minor blockNum

    // Structure for mapping of minor blocks 
    std::map<UInt,StdVector<std::set<Integer> > > minorBlocks;


    StdVector<UInt> sbmInd(0);
    if(algsys_->GetSolStrategy()->IsMultHarm()){
      UInt M = algsys_->GetSolStrategy()->GetNumHarmM();
      UInt a = (domain->GetDriver()->IsFullSystem())? (M+1) : ((M-1)/2 + 1);
      // same as ComputeIndex method in GraphManager, here with a lambda function
      auto ComputeIndex = [](UInt a, UInt b ) { return (domain->GetDriver()->GetNumFreq()) * a + b;};

      // store the sbm-indices of the nnz sbm-blocks
      for( UInt iRow = 0; iRow < domain->GetDriver()->GetNumFreq(); ++iRow ) {
        sbmInd.Push_back( ComputeIndex(iRow, iRow) );
        for( UInt iCol = iRow + 1; iCol < iRow + a ; ++iCol ) {
          if( iCol < domain->GetDriver()->GetNumFreq()){
            sbmInd.Push_back( ComputeIndex(iRow, iCol) );
          }
        }
      }

      // register it at algsys
      algsys_->SetNnzSBMInd(sbmInd);
    }// endif is multiharmonic

    // -----------------------------------------------------------
    //  1) Register FeFunctions with Algebraic System
    // -----------------------------------------------------------
    // Note: currently we use the same graph for all matrix
    // types (STIFFNESS, MASS, SYSTEM...)
    UInt numFcts = feFunctions_.size();
    bool useDistinctGraphs = false;
    algsys_->GraphSetupInit( numFcts,  useDistinctGraphs );

    LOG_DBG(stdPde) << pdename_ << ": Registering functions";
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fctIt;
    for( fctIt = feFunctions_.begin(); fctIt != feFunctions_.end(); fctIt++ ) {
      FeSpace & feSpace = *(fctIt->second->GetFeSpace());
      FeFctIdType fctId = fctIt->second->GetFctId();
      shared_ptr<ResultInfo> res = fctIt->second->GetResultInfo();
      std::string resultName = SolutionTypeEnum.ToString(res->resultType);

      feSpace.GetOlasMappings( sbmBlocks, minorBlocks);
      LOG_DBG(stdPde) << pdename_ << ":\tfctId #" << fctId
          << ", Type: " << resultName << ", #Equations: " 
          << feSpace.GetNumEquations();
      // Only register function, if there are any equations at all
      if(feSpace.GetNumEquations() > 0 ) {
        algsys_->RegisterFct( fctId, feSpace.GetNumEquations(),
                              feSpace.GetNumFreeEquations() );
      } else {
        LOG_DBG(stdPde) << pdename_ << ":\tfctId #" << fctId
                  << ", Omitting registration, as no equations are defined!";
      }
    }


    // ---------------------------------------
    //  2) Define SBM-Blocks and minor blocks
    // ---------------------------------------
    LOG_DBG(stdPde) << pdename_ << ": Defining SBM-blocks";

    UInt numBlocks = sbmBlocks.GetSize();

    if( solStrat_->IsMultHarm() && sbmBlocks.GetSize() > 1 )EXCEPTION("No submatrices allowed for multiharmonic analysis");

    // security check: ensure that at least one block is defined
    if (numBlocks == 0 ) {
      EXCEPTION( "There are no SBM blocks defined!" );
    }

    // Loop over blocks and register them at OLAS
    // This also holds for multiharmonic case, because there we only have
    // one set of equations, which are spread over the different blocks later on
    Integer sbmIndex = -1;
    for( UInt i = 0; i < numBlocks; ++i ) {

      // register block. In addition we check, if this is the inner block
      // and static condensation is activated
      bool isInnerBlock = solStrat_->UseStaticCondensation() && (i == numBlocks-1);

      if( solStrat_->IsMultHarm() && i != 0 ) EXCEPTION("Only one block allowed in multiharmonic algsys!");

      sbmIndex = algsys_->DefineSBMMatrixBlock( sbmBlocks[i], isInnerBlock );


      if( minorBlocks.size() != 0 && sbmIndex != -1) {
        StdVector<std::set<Integer> >& sbmSubBlocks = minorBlocks[i];

        // check if minor blocks are defined at all
        if(sbmSubBlocks.GetSize() == 0 ) continue;

        // Warning: Sub-matrix block is currently restricted to
        // just one FeFunction
        if( feFunctions_.size() > 1){
          WARN( "Sub-matrix block is currently just working for "
              << "1 FeFunction!"
              << "However unless you don't use any sub-block specific stuff like "
              << "special preconditioners it could work." );
          continue;
        }
        FeFctIdType fctId = feFunctions_.begin()->second->GetFctId();

        algsys_->RegisterSubMatrixBlocks(i, sbmSubBlocks.GetSize());
        // loop over all minor blocks
        for(UInt j = 0; j < sbmSubBlocks.GetSize(); ++j ) {
          UInt blockSize = sbmSubBlocks[j].size();
          StdVector<FeFctIdType> fctIds(blockSize);
          fctIds.Init(fctId);
          StdVector<Integer> eqns(blockSize);
          std::set<Integer>::iterator it = sbmSubBlocks[j].begin();
          UInt pos = 0;

          // loop over all eqns
          for( UInt k = 0; k < blockSize; ++k ) {
            eqns[pos++] = *it++;
          }
          algsys_->DefineSubMatrixBlocks(i,j, fctIds, eqns);
        }
      } // if block is defined at all
    } // loop over blocks

    // Finalize registration of blocks, which includes the generation of the graph manager
    algsys_->FinishRegistration();


    // Trigger writing of info file
    myInfo_->GetRoot()->ToFile("", true );


    // -----------------------------------
    //  3) Setup Sparsity Patterns
    // -----------------------------------
    LOG_DBG(stdPde) << pdename_ << ": Setting sparsity pattern";
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator it1;
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator it2;
    for( it1 = feFunctions_.begin(); it1 != feFunctions_.end(); ++it1 ) {
      for( it2 = feFunctions_.begin(); it2 != feFunctions_.end(); ++it2 ) {
        FeFctIdType fctId1 = it1->second->GetFctId();
        FeFctIdType fctId2 = it2->second->GetFctId();
        // assemble upper diagonal blocks including diagonal
        LOG_DBG(stdPde)<<pdename_<<":\tset graph for fctIds #"<< fctId1<<" and #"<<fctId2<<std::endl;
        assemble_->SetupMatrixGraph(fctId1, fctId2);
      } // it2
    } // it1

    // finish the assembly of the matrix graph
    if( solStrat_->IsMultHarm() ){
      algsys_->GraphSetupDoneMH();
    }else{
      algsys_->GraphSetupDone();
    }

    algsys_->GetGraphManager()->StopGraphTimers(); 

    // create matrices and solver object, if PDE is not direct coupled
    CreateMatrices_Solver();

    //Pass the system to every feFunction
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      fncIt->second->SetSystem(algsys_);
      fncIt++;
    }

    if(progOpts->DoEquationMapping())
      CreateEquationMapFile();

    // Trigger writing of info file
    myInfo_->GetRoot()->ToFile(true);

  }

   void StdPDE::CreateEquationMapFile()
   {
     // to be called after DefineAlgSys()

     std::string name =progOpts->GetSimName();

     if(domain->GetMultiSequenceDriver() != NULL)
       name += "_sequence_" + boost::lexical_cast<std::string>(domain->GetMultiSequenceDriver()->GetActSequenceStep());

     if(domain->GetSinglePDEs().GetSize() > 1)
       name += "_" + GetName(); // in case of more than one pde name it

     name += ".map";

     std::ofstream* out = new std::ofstream(name.c_str());
     if(out == NULL)
       throw Exception("cannot open equation mapping file " + name + " for writing");

     std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
     fncIt= feFunctions_.begin();
     while(fncIt != feFunctions_.end()){
       // Print equation information
       fncIt->second->GetFeSpace()->PrintEqnMap(out);
       fncIt++;
     }
     out->close();
     delete out;
   }

  BaseSolveStep * StdPDE::GetSolveStep()
  {
    assert(solveStep_!=nullptr);
    return solveStep_;
  }

  DampingType StdPDE::GetDamping(RegionIdType reg_id) const
  {
    std::map<RegionIdType,DampingType>::const_iterator it = dampingList_.find(reg_id);

    return it != dampingList_.end() ? it->second : NONE;
  }



  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================


  shared_ptr<ResultInfo> StdPDE::GetResultInfo( SolutionType solType ) {
    
    shared_ptr<ResultInfo> res;
    ResultSet::const_iterator it = availResults_.begin();
    for( ; it != availResults_.end(); ++it )  {
      if( (*it)->resultType == solType) {
        res = *it;
        break;
      }
    }

    // check, if resultinfo was found
    if( !res.get() ) {
      EXCEPTION( "A result with solutionType '" << SolutionTypeEnum.ToString(solType)
                 << "' was not found for " << pdename_ );
    }

    return res;
  }


  //============================================================================================
  //FeFunction Methods
  //============================================================================================

  shared_ptr<BaseFeFunction> StdPDE::GetFeFunction(SolutionType solType) {

    shared_ptr<BaseFeFunction> feFct;    
    if( feFunctions_.find(solType) != feFunctions_.end() ){
      feFct = feFunctions_[solType];
    }
    
    if( timeDerivFeFunctions_.find(solType) != timeDerivFeFunctions_.end() ){
      feFct = timeDerivFeFunctions_[solType];
    }
    
    if( !feFct )  {
      EXCEPTION( "A FeFunction descriptor with solutionType '" << SolutionTypeEnum.ToString(solType)
                 << "' was not found for " << pdename_ );
    }
    return feFct;

  }

  // **********
  // Geometry Information
  // **********
  void StdPDE::SetGeomInfo() {
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator FncIt= this->feFunctions_.begin();

    while(FncIt != this->feFunctions_.end()){
      FncIt->second->ApplyGeomInfo();
      FncIt++;
    }
  }

  // **********
  // SetRhsLoads
  // **********
  void StdPDE::SetRhsValues() {

    //do the same for RHS
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator rFncIt= this->rhsFeFunctions_.begin();
    while(rFncIt != this->rhsFeFunctions_.end()) {
    	rFncIt->second->ApplyLoads();
    	rFncIt++;
    }
  }

  // **********
  // Hysteresis
  // **********
  void StdPDE::SetPreviousHystVals(bool lastTS, bool forceMemoryLock){
    if ( isHysteresis_ ){
        //set current values to previous values for hysteresis operator
        //needed for the next time step
        std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
        std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
        for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
          /*
           * Note: If locked = true, overwrite = false
           */
          it->second->SetPreviousHystVals(lastTS, forceMemoryLock);
        }
     }
  }

  void StdPDE::CheckSaturationOfHystOperators(Double& lastTSSatAvg, Double& lastItSatAvg, Double& curItSatAvg,
      Double& oppositeDirAsTSAvg, Double& oppositeDirAsItAvg){

    Double lastTSSatAvgTMP,lastItSatAvgTMP,curItSatAvgTMP,oppositeDirAsTSAvgTMP,oppositeDirAsItAvgTMP;
    lastTSSatAvgTMP = 0.0;
    lastItSatAvgTMP = 0.0;
    curItSatAvgTMP  = 0.0;
    oppositeDirAsTSAvgTMP = 0.0;
    oppositeDirAsItAvgTMP = 0.0;
    UInt cnt = 0;
    if ( isHysteresis_ ){
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;

      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
        cnt++;
        it->second->checkSaturationStateAllElements(lastTSSatAvgTMP,lastItSatAvgTMP,curItSatAvgTMP,
          oppositeDirAsTSAvgTMP,oppositeDirAsItAvgTMP);

        LOG_DBG(stdPde) << "Saturation information for region " << it->first;
        LOG_DBG(stdPde) << lastTSSatAvgTMP*100.0 << " percent of elements were saturated during last timestep";
        LOG_DBG(stdPde) << lastItSatAvgTMP*100.0 << " percent of elements were saturated during last iteration";
        LOG_DBG(stdPde) << curItSatAvgTMP*100.0 << " percent of elements are currently saturated";
        LOG_DBG(stdPde) << oppositeDirAsTSAvgTMP*100.0 << " percent of elements changed direction of saturation by more than 90 degree compared to last timestep";
        LOG_DBG(stdPde) << oppositeDirAsItAvgTMP*100.0 << " percent of elements changed direction of saturation by more than 90 degree compared to last iteration";
      }
    } else {
      cnt = 1;
    }
    lastTSSatAvg = lastTSSatAvgTMP/( (Double) cnt);
    lastItSatAvg = lastItSatAvgTMP/( (Double) cnt);
    curItSatAvg  = curItSatAvgTMP/( (Double) cnt);
    oppositeDirAsTSAvg = oppositeDirAsTSAvgTMP/( (Double) cnt);
    oppositeDirAsItAvg = oppositeDirAsItAvgTMP/( (Double) cnt);

  }

  void StdPDE::EstimateCurrentSlopeForHysteresis(Double steppingLength, Double scaling){

    if ( isHysteresis_ ){
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;

      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
        it->second->ActiveOneShotSlopeEstimation(steppingLength, scaling);
      }
    }
  }

  bool StdPDE::MaterialTensorsHystDependent(){
    /*
     * aim: find out if material tensors are dependent on hysteresis
     * this is true, in case of piezos and magnetostriction as the coupling
     * tensors depend on hysteresis;
     * in electrostatics and magnetics, eps/nu do not depend on hysteresis, so
     * this function returns false here
     * NOTE: we just check the material tensors for dependency; the system matrix
     * might still depend on hysteresis in the case of deltaMatrix formulations
     */
    if ( isHysteresis_ ){
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;

      /*
       * Idea: iterate over all hysteresis regions (for each region we may have a
       * different coefFctHyst) and check if the coupling tensors are set or not;
       * if they are not set, the material tensors are indepenedent on hysteresis;
       * if they are set, we have to assume a dependency
       */
      bool dependency = false;
      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
        dependency = it->second->couplingTensorSet();
        if(dependency == true){
          break;
        }
      }
      return dependency;
    } else {
      // no hysteresis > no dependency
      return false;
    }
  }

  void StdPDE::TestInversionOfHystOperator(PtrParamNode testNode){
    if ( isHysteresis_ ){
        //set current values to previous values for hysteresis operator
        //needed for the next time step
        std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
        std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;

          // set flag with with the corresponding flagname
          for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
            /*
             * Note: If locked = true, overwrite = false
             */
            it->second->TestInversion(testNode,this->GetInfoNode());
          }
     }
  }

  void StdPDE::SetDoubleFlagInCoefFncHyst(std::string flagName, Double newState){
    if ( isHysteresis_ ){
      //set current values to previous values for hysteresis operator
      //needed for the next time step
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;

      // set flag with with the corresponding flagname
      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
        /*
         * Note: If locked = true, overwrite = false
         */
        it->second->SetDoubleFlag(flagName,newState);
      }
    }
  }

  void StdPDE::SetFlagInCoefFncHyst(std::string flagName, Integer newState){

    if ( isHysteresis_ ){
        //set current values to previous values for hysteresis operator
        //needed for the next time step
        std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
        std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;

        if(flagName == "outputDebugInfos"){
          for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
            /*
             * Note: If locked = true, overwrite = false
             */
            std::cout << it->second->ToString() << std::endl;
          }
        } else {
          // set flag with with the corresponding flagname
          for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
            /*
             * Note: If locked = true, overwrite = false
             */
            it->second->SetFlag(flagName,newState);
          }
        }
     } else {
			// check for coupled case if one of the pdes is hysteretic
			//if(){}
		 }
  }
} // end of namespace


EnumTuple StdPDE::ncTypeTuples_[] =
{
    EnumTuple(StdPDE::NC_NITSCHE, "Nitsche"),
    EnumTuple(StdPDE::NC_MORTAR, "Mortar"),
    EnumTuple(StdPDE::NC_NONE, "none")
};

Enum<StdPDE::NcCouplingType> StdPDE::ncCouplingType_ = 
    Enum<StdPDE::NcCouplingType>("Type of non-conforming formulation used",
                                  sizeof(ncTypeTuples_) / sizeof(EnumTuple),
                                  ncTypeTuples_);

EnumTuple StdPDE::lmTypeTuples_[] =
{
    EnumTuple(StdPDE::LM_STANDARD, "standard"),
    EnumTuple(StdPDE::LM_DUAL_DISCONTINUOUS, "dualDiscont"),
    EnumTuple(StdPDE::LM_DUAL_CUBIC, "dualCubic")
};

Enum<StdPDE::LagrangeMultType> StdPDE::lmType_ =
    Enum<StdPDE::LagrangeMultType>("Type of ansatz functions for Lagrange Multiplier",
                                   sizeof(lmTypeTuples_) / sizeof(EnumTuple),
                                   lmTypeTuples_);
