// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "StdPDE.hh"

#include "MatVec/Vector.hh"
#include "MatVec/generatematvec.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TransientDriver.hh"

#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

// headers for Paramhandling
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/Assemble.hh"

// header for logging
#include "DataInOut/Logging/LogConfigurator.hh"

#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"

namespace CoupledField {

// declare logging stream
  DECLARE_LOG(stdPde)
  DEFINE_LOG(stdPde, "stdPde")

  StdPDE::StdPDE(Grid *aptgrid, PtrParamNode paramNode ) :
    BasePDE(paramNode),
    ptGrid_(aptgrid),
    subType_(),
    numCouplingBcs_(0),
    nonLin_(false),
    nonLinMaterial_(false),
    isHysteresis_(false),
    isIterCoupled_(false),
    updateCouplingBCs_(false),
    diagMass_(false),
    needsAlgsys_(true),
    isAlwaysStatic_(false),
    dim_(ptGrid_->GetDim()), 
    isaxi_(param->Get("domain")->Get("geometryType")->As<std::string>() == "axi"),
    isComplex_(false),    
    needSolPrev_(false),
    isIncrFormulation_(false),
    updatedLagrangeForm_(false),
    assemble_(NULL),
    solveStep_(NULL),
    algsys_(NULL)
  {
  }
  
  StdPDE::~StdPDE() 
  {
    
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
    if ( analysistype_ != EIGENFREQUENCY ) {
      
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
   
    // Trigger writing of info file
    info->ToFile("", true );
    
    // First check if the PDE needs an algebraic system at all
    if( needsAlgsys_ == false ) {
      return;
    }


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
    //        sbm   min blockNum

    // Structure for mapping of minor blocks 
    std::map<UInt,StdVector<std::set<Integer> > > minorBlocks;

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
    for( fctIt = feFunctions_.begin(); 
        fctIt != feFunctions_.end(); 
        fctIt++ ) {
      FeSpace & feSpace = *(fctIt->second->GetFeSpace());
      FeFctIdType fctId = fctIt->second->GetFctId();
      shared_ptr<ResultInfo> res = fctIt->second->GetResultInfo();
      std::string resultName = SolutionTypeEnum.ToString(res->resultType);

      feSpace.GetOlasMappings( sbmBlocks, minorBlocks);
      LOG_DBG(stdPde) << pdename_ << ":\tfctId #" << fctId
          << ", Type: " << resultName << ", #Equations: " 
          << feSpace.GetNumEquations();
      algsys_->RegisterFct( fctId, feSpace.GetNumEquations(),
                            feSpace.GetNumFreeEquations() );
    }

    // ---------------------------------------
    //  2) Define SBM-Blocks and minor blocks
    // ---------------------------------------

    LOG_DBG(stdPde) << pdename_ << ": Defining SBM-blocks";

    UInt numBlocks = sbmBlocks.GetSize();
    // security check: ensure that at least one block is defined
    if (numBlocks == 0 ) {
      EXCEPTION( "There are no SBM blocks defined!" );
    }

    // Loop over blocks and register them at OLAS
    Integer sbmIndex = -1;
    for( UInt i = 0; i < numBlocks; ++i ) {

      // register block. In addition we check, if this is the inner block
      // and static condensation is activated
      bool isInnerBlock = solStrat_->UseStaticCondensation() &&
          (i == numBlocks-1);
      sbmIndex = algsys_->DefineSBMMatrixBlock( sbmBlocks[i], isInnerBlock );
      if( minorBlocks.size() != 0 && sbmIndex != -1) {
        StdVector<std::set<Integer> >& sbmSubBlocks = minorBlocks[i];

        // check if minor blocks are defined at all
        if(sbmSubBlocks.GetSize() == 0 ) continue;

        // Warning: Sub-matrix block is currently restricted to 
        // just one FeFunction
        if( feFunctions_.size() > 1){
          EXCEPTION( "Sub-matrix block is currently just working for "
              << "1 FeFunction!" );
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

    // Finalize registration of blocks
    algsys_->FinishRegistration();


    // Trigger writing of info file
    info->ToFile("", true );

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
        LOG_DBG(stdPde) << pdename_ << ":\tset graph for fctIds #"
            << fctId1 << " and # " << fctId2 << std::endl;
        assemble_->SetupMatrixGraph(fctId1, fctId2);
      } // it2
    } // it1

    // finish the assembly of the matrix graph
    algsys_->GraphSetupDone();

    // create matrices and solver object, if PDE is not direct coupled
    CreateMatrices_Solver();

    //Pass the system to every feFunction
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt= feFunctions_.begin();
    fncIt= feFunctions_.begin();
    while(fncIt != feFunctions_.end()){
      fncIt->second->SetSystem(algsys_);

      // Print equation information
      //fncIt->second->GetFeSpace()->PrintEqnMap();
      fncIt++;
    }
    //exit(0);
    // Trigger writing of info file
    info->ToFile("", true );

  }

  BaseSolveStep * StdPDE::GetSolveStep() {
    
    return solveStep_;
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

  shared_ptr<BaseFeFunction> StdPDE::GetFeFunction( SolutionType solType ) {
    
    //TODO> We need to find a more failsafe way to store the entity names associated with a 
    //FeFunction
    shared_ptr<BaseFeFunction> feFct;
    SolutionType mySolType;
    if( feFunctions_.find(solType) == feFunctions_.end()){
      EXCEPTION( "A FeFunction descriptor with solutionType '" << SolutionTypeEnum.ToString(solType)
                 << "' was not found for " << pdename_ );
    }else{
      mySolType = solType;
    }
    if(feFunctions_.find(mySolType) != feFunctions_.end()){
      feFct = feFunctions_[mySolType];
    }else{
      EXCEPTION("StdPDE::GetFeFunction: Could not find the corresponding FeFunction for the given entity name\n \
                         Did you specify all Regions, Surfregions and NamedNodes in the xml?");
    }
    return feFct;

  }

} // end of namespace


static EnumTuple ncTypeTuples[] =
{
 EnumTuple(StdPDE::NITSCHE, "Nitsche"),
 EnumTuple(StdPDE::MORTAR, "Mortar"),
 EnumTuple(StdPDE::NONE, "None")
};

Enum<StdPDE::NcCouplingType>StdPDE::ncCouplingType_ = \
    Enum<StdPDE::NcCouplingType>("Type of non-conforming formulation used",
                                  sizeof(ncTypeTuples) / sizeof(EnumTuple),
                                  ncTypeTuples);
