// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <fstream>

#include "algsys/olascomm.hh"
#include "algsys/standardsys.hh"

#include "utils/utils.hh"
#include "graph/graphmanagersimple.hh"
#include "graph/graphmanagerstdmat.hh"
#include "matvec/matvec.hh"
#include "solver/solver.hh"
#include "precond/precond.hh"

#include "algsys/generateentrymanipulator.hh"
#include "algsys/baseentrymanipulator.hh"
#include "algsys/generateidbchandler.hh"
#include "algsys/baseidbchandler.hh"

#include "General/exception.hh"

#include "DataInOut/Logging/cfslog.hh"

using CoupledField::Exception;

DECLARE_LOG(stdSys)
DEFINE_LOG(stdSys, "stdSys")

namespace OLAS {
 

  // ***********************
  //   Default Constructor
  // ***********************
  StandardSystem::StandardSystem(ParamNode* pn) : BaseSystem(pn) {


    totalSize_           = 0;
    precond_             = NULL;
    rhs_                 = NULL;
    sol_                 = NULL;
    solComm_             = NULL;
    patternPool_         = NULL;
    eigenValues_         = NULL;
    algSysType_          = STANDARD_SYSTEM;

    // Initialize pointer to finite-element matrices
    NewArray( sysmat_, StdMatrix*, MAX_NUM_FE_MATRICES );
    for ( UInt i = 1; i <= MAX_NUM_FE_MATRICES; i++ ) {
      sysmat_[i] = NULL;
    }

    // Set flag for insertion of penalty terms into matrix
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      assembleDirichletToSysMat_ = true;
    }
    else {
      assembleDirichletToSysMat_ = false;
    }
  }


  // **************
  //   Destructor
  // **************
  StandardSystem::~StandardSystem() {


    delete precond_;
    delete sol_;
    delete solComm_;
    delete rhs_;
    delete eigenValues_;

    // Delete finite-element matrices
    for ( UInt i = 1; i <= MAX_NUM_FE_MATRICES; i++ ) {
      delete sysmat_[i];
    }
    DeleteArray( sysmat_ );

    // After matrices we may delete the pattern pool
    delete patternPool_;
  }


  // *****************************************
  //   Trigger setup phase of preconditioner
  // *****************************************
  void StandardSystem::SetupPrecond( ) {
    precond_->Setup( *sysmat_[SYSTEM] );
  }


  // *********************************
  //   Trigger setup phase of solver
  // *********************************
  void StandardSystem::SetupSolver( ) {
    solver_->Setup( *sysmat_[SYSTEM] );
  }
  
  // **************************************
  //   Trigger setup phase of EigenSolver
  // **************************************
  void StandardSystem::SetupEigenSolver( UInt numFreq, Double shift, 
                                         bool isQuadratic ) {

    // Determine if a generalized or a quadratic eigenvalue 
    // problem has to be solved
    // In the latter case, a damping matrix has to be present, 
    // otherwise we issue a warning
    bool dampPresent = (matrixTypes_.find( DAMPING) != matrixTypes_.end());
    
    if( isQuadratic == true ) {
      if( dampPresent == false ) {
        (*error) << "To solve a quadratic eigenvalue problem, a damping "
                 << "matrix has to be present!";
        Error( __FILE__, __LINE__ );
      }
     
      // Setup the quadratic eigenvalue solver
      eigenSolver_->Setup( *sysmat_[STIFFNESS], *sysmat_[MASS],
                           *sysmat_[DAMPING], numFreq, shift );
    } else {
      if( dampPresent == true ) {
        (*warning) << "Although a damping matrix is present, only a generalized "
                   << "eigenvalue will be solved (as you have specified)!";
        Warning( __FILE__, __LINE__ );
      }
      
      // Setup the eigenvalue solver
      eigenSolver_->Setup( *sysmat_[STIFFNESS], *sysmat_[MASS],
                           numFreq, shift );
    }
    
    // Determine some basic properties and create vectors
    // for eigenvalues and related error bounds
    MatrixEntryType eType;
    if( isQuadratic == true ) { 
      eType = COMPLEX;
    } else { 
      eType = DOUBLE;
    }

    MatrixStorageType sType = sysmat_[SYSTEM]->GetStorageType();
    UInt totalSize = totalSize_;
    totalSize *= blockSize_;

    BaseVector *bVec = GenerateSparseVectorObject( sType, eType, 1, totalSize );
    BaseVector *errVec = GenerateSparseVectorObject( sType, eType, 1, totalSize );
    eigenValues_ = dynamic_cast<SparseVector*>( bVec );
    eigenValError_ = dynamic_cast<SparseVector*>( errVec );
    
  }

  // ***********************
  //   Solve linear system
  // ***********************
  void StandardSystem::Solve(const std::string& comment) {


    // If the penalty formulation is used and we have inhomogeneous
    // Dirichlet boundary conditions, then the righ-hand side is
    // "contaminated" with penalty terms
    if ( numDirichletValues_ > 0 &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      myParams_.SetValue( "RHSwithPenalty", true );
      LOG_DBG(stdSys) << "Solve: rhs with penalty";
    }

#ifdef DEBUG_STANDARDSYSTEM
    (*debug) << "Systemmatrix:" << std::endl;
    sysmat_[SYSTEM]->Print(*debug);
    (*debug) << "RHS" << std::endl;
    rhs_->Print(*debug);
#endif

#ifdef PROFILING
    Double t1 = Profiler::GetRealTime();
#endif

    // Iterative solvers require an initial guess and in the penalty case
    // we should insert the Dirichlet values into it
    if ( dynamic_cast<BaseIterativeSolver*>(solver_) != NULL &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      idbcHandler_->SetDofsToIDBC( sol_ );
      (*cla) << " Inserted Dirichlet values into initial guess"
             << std::endl;
      LOG_DBG(stdSys) << "Solve: Inserted Idbc values into initial guess for iterative solver"; 
    }

    // Perform a simple sanity check
    if ( sysmat_[SYSTEM] == NULL || precond_ == NULL || rhs_ == NULL ||
         sol_ == NULL ) {
      Error( "Detected NULL pointer where there should be none!", __FILE__,
             __LINE__ );
    }

    // Assume that everything will go well
    myReport_.SetValue( "solutionIsOkay", true );

    // Now modifiy the right-hand side vector
    LOG_DBG(stdSys) << "Solve: add idbc to rhs";
    idbcHandler_->AddIDBCToRHS( rhs_ );
    
    // check if we do export stuff
    ParamNode* els = xml  != NULL && xml->Has("exportLinSys") ? xml->Get("exportLinSys") : NULL;
    std::string file;
    std::string base;

    // need it common even when exclusive solution
    if(els) {
      std::ostringstream os;
      os << els->Get("baseName")->AsString();
      if(comment != "") os << "_" << comment;
      base = os.str();  
    }

    // check if we do not only want the solution
    if(els && els->Get("solution")->AsString() != "exclusive")
    {
      // two formats. The harwell-boing format includes the rhs!
      if(els->Get("format")->AsString()  == "harwell-boeing") 
      {
        sysmat_[SYSTEM]->ExportHarwellBoeing(base+".hb", *rhs_);

        if(els->Has("damping", true) && sysmat_[DAMPING] != NULL)
          sysmat_[DAMPING]->ExportHarwellBoeing(base+"_damping.hb", *rhs_);
        
        if(els->Has("auxiliary", true) && sysmat_[AUXILIARY] != NULL)
          sysmat_[AUXILIARY]->ExportHarwellBoeing(base+"_aux.hb", *rhs_);
      } 
      else // classical (default) matrix-market 
      {
        sysmat_[SYSTEM]->Export((base+".mtx").c_str() );

        if(els->Has("damping", true) && sysmat_[DAMPING] != NULL)
          sysmat_[DAMPING]->Export((base+"_damping.mtx").c_str() );
        
        if(els->Has("auxiliary", true) && sysmat_[AUXILIARY] != NULL)
          sysmat_[AUXILIARY]->Export((base+"_aux.mtx").c_str() );

        // rhs is only in harwell-boing included    
        rhs_->Export((base+".vec").c_str() );
      }    
    }

    if(els && els->Has("initialGuess", true))
      sol_->Export((base+"_intial_guess.vec").c_str());
    
    LOG_TRACE2(stdSys) << "before solve: euclidian norm of initial guess = " << sol_->NormEuclid(); 
    LOG_DBG(stdSys) << "euclidian norm of rhs = " << rhs_->NormEuclid();
    
    // ---------------------------
    // This is the expensize part! solve the system
    solver_->Solve( *sysmat_[SYSTEM], *precond_, *rhs_, *sol_ );
    
    LOG_TRACE2(stdSys) << "after solve: euclidian norm of solution = " << sol_->NormEuclid();

    if(els && els->Get("solution")->AsString() != "no")
      sol_->Export((base+".sol.vec").c_str());

    LOG_DBG(stdSys) << "Solve: remove idbc from rhs";
    // Now de-modifiy the right-hand side vector
    idbcHandler_->RemoveIDBCFromRHS( rhs_ );
    

    // Check that solution went fine, if not issue a warning
    if ( myReport_.GetBoolValue( "solutionIsOkay" ) == false ) {
      // std::cerr << "Life's a piece of shit, when you look at it ...\n";
      (*warning) << "Solver reports a problem! Consult .las file for "
                 << "further diagnostics!";
      Warning( __FILE__, __LINE__ );
    }

#ifdef PROFILING    
    Double t2 = Profiler::GetRealTime();
    (*cla)  << "solution time: " << t2-t1 << " seconds " << std::endl;
    Profiler::WriteReport();
#endif

    LOG_DBG(stdSys) << "Solve invalidate flag for solution buffer, why not also fro rhs?";
    // Set invaldiation flag for assembling routine
    assemble_->InvalidateSolBuffer();
  }

  // ******************************
  //   Calculate Eigenfrequencies
  // ******************************
  UInt StandardSystem::
  CalcEigenFrequencies( const Double * &frequencies, const Double* &err  ) {

    
    // Trigger calculation of eigenvalues
    UInt numConverged = 
      eigenSolver_->CalcEigenFrequencies( *eigenValues_, *eigenValError_ );
    
    // Hard coded cast for double values
    RefCast( *eigenValues_, Vector<Double>, valVec );
    frequencies = valVec.GetPointer();
    frequencies++;

    RefCast( *eigenValError_, Vector<Double>, errVec );
    err = errVec.GetPointer();
    err++;

    // Return number of converged eigenvalues
    return numConverged;
    
  }

  UInt StandardSystem::
  CalcEigenFrequencies( const Complex * &frequencies, const Double* &err  ) {

    
    // Check, if eigenvalue solver is quadratic, as only in this case
    // this method is well-defined
    
    if( eigenSolver_->IsQuadratic() == false ) {
      (*error) << "When solving a generalized eigenvalue problem, only "
               << "real-valued results are obtained! Use the second "
               << "CalcEigenFrequencies()-method!";
      Error( __FILE__, __LINE__ );
    }

    // Trigger calculation of eigenvalues
    UInt numConverged = 
      eigenSolver_->CalcEigenFrequencies( *eigenValues_, *eigenValError_ );
    
    // Hard coded cast for double values
    RefCast( *eigenValues_, Vector<Complex>, valVec );
    frequencies = valVec.GetPointer();
    frequencies++;

    RefCast( *eigenValError_, Vector<Double>, errVec );
    err = errVec.GetPointer();
    err++;

    // Return number of converged eigenvalues
    return numConverged;
    
  }

  
  // ******************************
  //   Calculate Eigenmodes
  // ******************************
  void StandardSystem::CalcEigenMode( UInt numMode ) {
    
    Vector<Double> & solHelp =
      dynamic_cast<Vector<Double> &> (*sol_);
    eigenSolver_->CalcEigenMode( numMode, solHelp );

    // Set invaldiation flag for assembling routine
    assemble_->InvalidateSolBuffer();
  }


  // ************************************
  //   Creation of matrices and vectors
  // ************************************
  void StandardSystem::CreateLinSys() {


    MatrixEntryType   entryType;
    MatrixStorageType storageType;

    myParams_.GetEnumValue( "MatrixEntryType"  , entryType   );
    myParams_.GetEnumValue( "MatrixStorageType", storageType );

    // get the matrix graph of the system matrix
    BaseGraph *graph = graphManager_->GetGraph();
    nne_  = graph->GetNNE();
    size_ = graph->GetSize();

    // Determine total number of equation numbers for this system
    // including those for fixed dofs
    totalSize_ = 0;
    for ( UInt i = 1; i <= numPDEs_; i++ ) {
      totalSize_ += sizePerPDE_[i];
    }

    LOG_DBG(stdSys)  << "CreateLinSys() size=" << size_ << " nne=" << nne_ << " totalSize=" << totalSize_; 
    LOG_DBG3(stdSys) << " Matrix Graph=" << graph->ToString();

    // ------------------------
    //  Re-set numLastFreeDof
    // ------------------------
    DeleteArray( numLastFreeDof_ );
    NewArray( numLastFreeDof_, UInt, 1 );
    numLastFreeDof_[1] = graph->GetSize();

    LOG_DBG(stdSys) << " CreateLinSys(): Re-set numLastFreeDof to " << numLastFreeDof_[1];

    // --------------------------------------------
    //  Treatment of Dirichlet Boundary Conditions
    // --------------------------------------------

    // Check the case we have (elimination vs. penalty) and generate
    // an appropriate object for handling inhomogeneous Dirichlet
    // boundary conditions
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      idbcHandler_ = GenerateIDBC_HandlerObject( numDirichletValues_,
                                                 blockSize_, entryType );
    }
    else {
      idbcHandler_ = GenerateIDBC_HandlerObject( matrixTypes_,
                                                 graphManager_,
                                                 1,
                                                 numDirichletValues_,
                                                 entryType,
                                                 false );
    }


    // ------------------------
    //  Generation of matrices
    // ------------------------
    std::set<FEMatrixType>::iterator it;
    
    // Print information about matrix types to .las file
    (*cla) << "\n StandardSystem: The following matrices are created:\n ";
    for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ ) {
      (*cla) << " " << Enum2String(*it) << '\n';
    }
    (*cla) << std::endl;

    //check for solver
    SolverType solver;
    myParams_.GetEnumValue( "Solver", solver );

    // In the case of SRCS_Matrices we can save memory by sharing
    // sparsity patterns
    bool sharePatternPossible = false;
    bool sharePattern = false;
    PatternIdType patternID = NO_PATTERN_ID;
    if ( storageType == SPARSE_SYM && solver != DIAGSOLVER ) {
      sharePatternPossible = true;
      patternPool_ = New PatternPool;
    }

    // actual matrix generation
    for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ ) 
    {
      MatrixStorageType actStorageType;

      LOG_DBG(stdSys) << "Generating matrix[" << *it << "]" << " MatrixEntryType=" << entryType
      << " MatrixStorageType= " << storageType << " dofs= " << blockSize_
      << " rows/cols= " << size_ << " nnz=" << nne_;

      actStorageType = storageType;
      if ( solver == DIAGSOLVER ) {
        if ( *it == SYSTEM || *it == MASS ) {
          // just allocate a diagonal matrix
          actStorageType = DIAG;
        }
      }

      LOG_DBG(stdSys) << "calling GenerateStdMatrixObject" << std::endl;

      sysmat_[*it] = GenerateStdMatrixObject( entryType, actStorageType,
          blockSize_, size_, size_,
          nne_ );
      
      assert(sysmat_[*it] != NULL);
      sysmat_[*it]->name = Enum2String(*it);

      LOG_DBG(stdSys) << "created matrix " << sysmat_[*it]->ToString(); 
      
      // Generate sparsity pattern once, but share it afterwards
      // if this is possible
      if ( sharePattern == false ) 
      {
        sysmat_[*it]->SetSparsityPattern( (*graph) );
        if ( sharePatternPossible == true ) 
        {
          sharePattern = true;
          patternID = sysmat_[*it]->TransferPatternToPool( patternPool_ );
        }
      }
      else 
        sysmat_[*it]->SetSparsityPattern( patternPool_, patternID );
    }


    // ---------------------------------
    //  Generation of remaining objects
    // ---------------------------------

    // Generate vectors
    rhs_ = dynamic_cast<SparseVector *>
      (GenerateVectorObject( *(sysmat_[SYSTEM]) ));
    sol_ = dynamic_cast<SparseVector *>
      (GenerateVectorObject( *(sysmat_[SYSTEM]) ));

    // Generate communication buffers
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) 
      GenerateCFSTransferBuffers();

    if ( rhs_ == NULL || sol_ == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

    #ifdef DEBUG_STANDARDSYSTEM
    // Seens to be somehow buggy, as I got an Segemtation fault
    // when running it with block-matrices
    //(*debug) << "System Matrix connectivity: " << std::endl;
    //sysmat_[SYSTEM]->Print(*debug);
    (*debug) << "rhs: " << std::endl;
    rhs_->Print(*debug);
    (*debug) << "sol_: " << std::endl;
    sol_->Print(*debug);
#endif

    // Create the assemble object
    assemble_ = GenerateEntryManipulatorObject( entryType, blockSize_ );


    // -----------------
    //  Memory clean-up
    // -----------------

    // At this point, hopefully the graph object is no longer
    // required by anyone, so release pointer and delete manager
    // to free memory
    graph = NULL;
    delete graphManager_;
    graphManager_ = NULL;
  }


  // ************************
  //   Create solver object
  // ************************
  void StandardSystem::CreateSolver()
  {
    SolverType solver;
    myParams_.GetEnumValue( "Solver", solver );
    solver_ = GenerateSolverObject( *(sysmat_[SYSTEM]), solver, xml, &myParams_, &myReport_ );
  }


  // ********************************
  //   Create preconditioner object
  // ********************************

  // Das Interface muessen wir uns noch ueberlegen
  void StandardSystem::CreatePrecond() {

    
    PrecondType precond;
    myParams_.GetEnumValue("Precond", precond);
    precond_ = GenerateStdPrecondObject( *(sysmat_[SYSTEM]), precond,
                                         &myParams_, &myReport_ );

  }

  // ********************************
  //   Create eigensolver object
  // ********************************

  void StandardSystem::CreateEigenSolver() {

    
    EigenSolverType egSolver;
    myParams_.GetEnumValue("EigenSolver", egSolver);
    eigenSolver_ = GenerateEigenSolverObject( *(sysmat_[SYSTEM]), egSolver,
                                              xml, 
                                              &myParams_, &myReport_ );

  }


  // =======================================================================
  //   Methods for passing data
  // =======================================================================



  // *******************
  //   SetFEMatrixType
  // *******************
  void StandardSystem::SetFEMatrixType(const FEMatrixType matType,
                                       const PdeIdType identifier1, 
                                       const PdeIdType identifier2) 
  {
    LOG_DBG2(stdSys) << "SetFEMatrixType matType=" << matType << " pde1=" 
                     << identifier1 << " pde2=" << identifier2;
    if(matType != NOTYPE) matrixTypes_.insert(matType);
  }


  // ****************
  //   SetBlockSize
  // ****************
  void StandardSystem::SetBlockSize( const PdeIdType identifier,
                                     const UInt bs ) {
    blockSize_ = bs;
  }



  // =======================================================================
  //   Methods for assembling, setting BCs and rhs
  // =======================================================================
  


  // **************
  //   InitMatrix
  // **************
  void StandardSystem::InitMatrix( FEMatrixType matrixID,
                                   const PdeIdType identifierPDE ) 
  {

    LOG_TRACE2(stdSys) << "InitMatrix matrixID=" << matrixID << " identifierPDE=" << identifierPDE;
    
    // If matrix specified init this one
    if ( matrixID != NOTYPE ) 
    {
      assert(sysmat_[matrixID] != NULL);
      sysmat_[matrixID]->Init();
      idbcHandler_->InitMatrix(matrixID);
    }

    // Otherwise init all matrices
    else {
      std::set<FEMatrixType>::iterator it;
      for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ ) {
        sysmat_[*it]->Init();
        idbcHandler_->InitMatrix(*it);
      }
    }

    // Do we need to re-insert penalty values into matrix?
    if ( ( matrixID == NOTYPE || matrixID == SYSTEM ) &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      assembleDirichletToSysMat_ = true;
    }
  }


  // ***********************
  //   InitRHS (Version 1)
  // ***********************
  void StandardSystem::InitRHS( const PdeIdType identifierPDE ) {
    rhs_->Init();

    // Invalidate rhs buffer of assemble
    assemble_->InvalidateRhsBuffer();
  }


  // ***********************
  //   InitRHS (Version 2)
  // ***********************
  void StandardSystem::InitRHS( const Double *newRHS ) {
    assemble_->InitRHS( rhs_, newRHS );
  }

  void StandardSystem::InitRHS( const Vector<Complex>* newRHS ) {
    assemble_->InitRHS( rhs_, newRHS );
  }
  

  // **********************
  //   InitSol (Version1)
  // **********************
  void StandardSystem::InitSol( const PdeIdType identifierPDE ) {
    sol_->Init();

    // Invalidate sol buffer of assemble
    assemble_->InvalidateSolBuffer();
  }

  // ***********************
  //   InitSol (Version 2)
  //   for doubles
  // ***********************
  void StandardSystem::InitSol( const Double *newSol, const UInt size ) {
    
//     try{
// 	    std::cerr << "\n Initializing solution vector to initial conditions " << std::endl;
	
	    if( size == sol_->GetSize() ){
//	    	std::cout << "\n InitSol: " ;
//	    	for(UInt i=0; i<size; i++){
//	    		std::cout <<  newSol[i] << "," ;
//	    	}
//	    	std::cout << std::endl;
	    	assemble_->InitSol( sol_, newSol );
	    }
// 	    else{
// #ifdef DEBUG_STANDARDSYSTEM
//           std::cout << "\n aux:" << size <<" != sol_:"<< sol_->GetSize() << std::endl;
//           (*error) << "size of initial solution vector doesn't match "
//                    << "with the size of the initialized solution vector!";
//           Error( __FILE__, __LINE__ );
// #endif
// 	    }
//     }
//     catch(Exception& e) {
//       throw e.ToString();
//     }
        
  }
  


  // ********************
  //   SetElementMatrix 
  // ********************
  void StandardSystem::SetElementMatrix( FEMatrixType matrixID,
                                         Double *elemMat, 
                                         PdeIdType idPDE1,
                                         Integer *connect1,
                                         Integer elemSize1,
                                         PdeIdType idPDE2,
                                         Integer *connect2,
                                         Integer elemSize2,
                                         bool setCounterPart ) {


    // then normal assembly:
    // element matrix and its transposed on the counter  
    // part of the global matrix when setCounterPart is 
    // true (upper and lower global assembly).
    assemble_->SetElementMatrix( matrixID, idPDE1, idPDE2,
                                 sysmat_[matrixID],
                                 sysmat_[matrixID],
                                 idbcHandler_, elemMat,
                                 connect1, (UInt)elemSize1,
                                 connect2, (UInt)elemSize2,
                                 numLastFreeDof_[1],
                                 numLastFreeDof_[1],
                                 setCounterPart );



  }


  // *****************
  //   SetElementRHS
  // *****************
  void StandardSystem::SetElementRHS( Double *elemRHS, 
                                      const PdeIdType identifierPDE,
                                      Integer *connect, UInt length ) {


    // Delegate work to EntryManipulator
    assemble_->SetElementRHS( rhs_, elemRHS, connect, length, size_ );
  }


  // ***********************
  //   SetNodeRHS (Double)
  // ***********************
  void StandardSystem::SetNodeRHS( Double val, const PdeIdType identifierPDE,
                                   Integer eqnNum ) {


    // Perform consistency check
#ifdef DEBUG_STANDARDSYSTEM
    if ( eqnNum > size_ || eqnNum < 1 ) {
      (*error) << "StandardSystem::SetNodeRHS: Inconsistency detected:"
               << "\n size      = " << size_
               << "\n totalSize = " << totalSize_
               << "\n eqnNum    = " << eqnNum
               << "\n val       = " << val
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to entry manipulation object
    assemble_->SetNodeRHS( rhs_, val, eqnNum  );
  }

  void StandardSystem::SetNodeRHS( Complex val, const PdeIdType identifierPDE,
                                   Integer eqnNum  ) {

    
    // Perform consistency check
#ifdef DEBUG_STANDARDSYSTEM
    if ( eqnNum > size_ || eqnNum < 1 ) {
      (*error) << "StandardSystem::SetNodeRHS: Inconsistency detected:"
               << "\n size      = " << size_
               << "\n totalSize = " << totalSize_
               << "\n eqnNum    = " << eqnNum
               << "\n val       = " << val
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'";
      Error( __FILE__, __LINE__ );
    }
#endif
    
    // Delegate work to entry manipulation object
    assemble_->SetNodeRHS( rhs_, val, eqnNum );
  }


  // ************************
  //   SetNodeRHS (Complex)
  // ************************

  // **************
  //   Update RHS
  // **************  
  void StandardSystem::UpdateRHS( FEMatrixType matrix_id, Double *fup ) {
    if ( numDirichletValues_ > 0 &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {
      (*warning) << "UpdateRHS may fail in combo with "
                 << "idbcHandling = elimination";
      Warning( __FILE__, __LINE__ );
    }
    assemble_->UpdateRHS( rhs_, sysmat_[matrix_id], fup );
  }

  // *************************************
  //   Remove IDBC Information from matrix
  // *************************************

  void StandardSystem::RemoveIDBCInfoFromMatrix() const {

    std::set<FEMatrixType>::iterator it;
    BaseVector *rhs;
    for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ ) {
      //sysmat_[*it]->Init();
      idbcHandler_->InitMatrix(*it);
      idbcHandler_->RemoveIDBCFromRHS( rhs );
    }
    
  }


  // *************************
  //   SetDirichlet (Double)
  // *************************
  void StandardSystem::SetDirichlet(const PdeIdType pdeID,
                                     Integer eqnNum, const Double &val  ) {


    // Perform some consistency checks
#ifdef DEBUG_STANDARDSYSTEM

    // Some situation dependend modifiers
    bool usingPenalty = myParams_.GetBoolValue( "UsingPenaltyFormulation" );
    UInt maxVal = usingPenalty == true ? size_ : size_ + numDirichletValues_;
    UInt minVal = usingPenalty == true ? 1 : size_ + 1;

    if ( eqnNum > maxVal || eqnNum < minVal ) {
      (*error) << "StandardSystem::SetDirichlet: Inconsistency detected:"
               << "\n pdeID  = " << pdeID
               << "\n eqnNum = " << eqnNum
               << "\n val    = " << val
               << "\n numBC  = " << numDirichletValues_
               << "\n size   = " << size_
               << "\n minVal = " << minVal
               << "\n maxVal = " << maxVal
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( pdeID, eqnNum, val );
  }


  // **************************
  //   SetDirichlet (Complex)
  // **************************
  void StandardSystem::SetDirichlet( const PdeIdType pdeID,
                                     Integer eqnNum, const Complex &val ) {


    // Perform some consistency checks
#ifdef DEBUG_STANDARDSYSTEM

    // Some situation dependend modifiers
    bool usingPenalty = myParams_.GetBoolValue( "UsingPenaltyFormulation" );
    UInt maxVal = usingPenalty == true ? size_ : size_ + numDirichletValues_;
    UInt minVal = usingPenalty == true ? 1 : size_ + 1;

    if ( eqnNum > maxVal || eqnNum < minVal || ) {
      (*error) << "StandardSystem::SetDirichlet: Inconsistency detected:"
               << "\n pdeID  = " << pdeID
               << "\n eqnNum = " << eqnNum
               << "\n val    = " << val
               << "\n numBC  = " << numDirichletValues_
               << "\n size   = " << size_
               << "\n minVal = " << minVal
               << "\n maxVal = " << maxVal;
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( pdeID, eqnNum, val );
  }


  // *******************
  //   BuildInDirichlet
  // *******************
  void StandardSystem::BuildInDirichlet() {


    // If necessary modify matrix diagonal for penalty approach
    if ( assembleDirichletToSysMat_ == true ) {
      idbcHandler_->AdaptSystemMatrix( *(sysmat_[SYSTEM]) );
      assembleDirichletToSysMat_ = false;
    }
  }


  // ****************************
  //   ConstructEffectiveMatrix
  // ****************************
  void StandardSystem::ConstructEffectiveMatrix( const factorMap &matFactors ) 
  {
    LOG_DBG(stdSys) << "ConstructEffectiveMatrix() matFactors = " << matFactors.size()
                    << " matrixTypes: " << matrixTypes_.size();

    factorMap::const_iterator it;
    StdMatrix *sys = sysmat_[SYSTEM];

    // It's okay, if there are no factors, if there is only a system
    // matrix and no other non auxiliary ones
    if(matFactors.empty() == true) 
    {
      // one element is auxiliary
      bool has_aux = std::find(matrixTypes_.begin(), matrixTypes_.end(), AUXILIARY) != matrixTypes_.end();

      if(sys != NULL && (matrixTypes_.size() == 1 || (matrixTypes_.size() == 2 && has_aux)))
      {

        LOG_TRACE(stdSys) << "ConstructEffectiveMatrix: Nothing to be done, "
                          << "since there is only one (non auxiliary) FE matrix";

        // Also assemble the effective auxilliary system matrix for moving
        // IDBCs to the right-hand side
        idbcHandler_->BuiltSystemMatrix( matFactors );

        // Now we are done
        return;
      }
      else {
        (*warning) << "StandardSystem::ConstructEffectiveMatrix: "
                   << "Map with factors is empty, but there are "
                   << matrixTypes_.size() << " FE matrices in the game!";
        Warning( __FILE__, __LINE__ );
      }
    }

    // Intialize system matrix before adding the weighted
    // matrices to it
    this->InitMatrix( SYSTEM );

    for ( it = matFactors.begin(); it != matFactors.end(); it++ )
    {
      LOG_DBG(stdSys) << " matFactor: " << Enum2String((*it).first) << ":" << (*it).second;      
      if (sysmat_[(*it).first] != NULL  && (*it).second != 0.0 ) 
        sys->Add((*it).second, *sysmat_[(*it).first] );
    }
    LOG_DBG3(stdSys) << "effective Matrix: " << sys->ToString();

    // Also assemble the effective auxilliary system matrix for moving
    // IDBCs to the right-hand side
    idbcHandler_->BuiltSystemMatrix( matFactors );
  }


  // ******************************
  //   GetSolutionVal (Version 1)
  // ******************************
  Integer StandardSystem::
  GetSolutionVal( Double* &ptSol,
                  const PdeIdType identifierPDE ) {


    // Elimination case
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {

      // Will not work
      if ( myParams_.GetBoolValue( "Parallel" ) ) {
        (*error) << "Parallel only supported in combo with penalty!";
        Error( __FILE__, __LINE__ );
      }

      // Sequential setting: We copy the solution and add the Dirichlet values
      Vector<Double> *dVec1 = dynamic_cast<Vector<Double>*>( sol_ );
      Vector<Double> *dVec2 = dynamic_cast<Vector<Double>*>( solComm_ );
      Double *data1 = dVec1->GetPointer();
      Double *data2 = dVec2->GetPointer();

      for ( int i = 1; i <= size_; i++ ) {
        data2[i] = data1[i];
      }
      idbcHandler_->SetDofsToIDBC( solComm_ );

      ptSol = data2 + 1;
    }

    // Penalty case
    else {
      assemble_->GetSolutionVal( sol_, ptSol, identifierPDE );
    }

    // Compute total size
    return totalSize_ * blockSize_;
  }


  // ******************************
  //   GetSolutionVal (Version 2)
  // ******************************
  Integer StandardSystem::
  GetSolutionVal( Complex* &ptSol,
                  const PdeIdType identifierPDE ) {


    // Elimination case
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {

      // Will not work
      if ( myParams_.GetBoolValue( "Parallel" ) ) {
        (*error) << "Parallel only supported in combo with penalty!";
        Error( __FILE__, __LINE__ );
      }

      // Sequential setting: We copy the solution and add the Dirichlet values
      Vector<Complex> *cVec1 = dynamic_cast<Vector<Complex>*>( sol_ );
      Vector<Complex> *cVec2 = dynamic_cast<Vector<Complex>*>( solComm_ );
      Complex *data1 = cVec1->GetPointer();
      Complex *data2 = cVec2->GetPointer();

      for ( int i = 1; i <= size_; i++ ) {
        data2[i] = data1[i];
      }
      idbcHandler_->SetDofsToIDBC( solComm_ );

      ptSol = data2 + 1;
    }

    // Penalty case
    else {
      assemble_->GetSolutionVal( sol_, ptSol, identifierPDE );
    }

    // Compute total size
    return totalSize_ * blockSize_;
  }
  

  // *************************
  //   GetRHSVal (Version 1)
  // *************************
  Integer StandardSystem::GetRHSVal( Double* &ptRhs,
                                     const PdeIdType identifierPDE ) {
   //  if ( numDirichletValues_ > 0 &&
//          myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {
//       (*error) << "GetRHSVal does not work for idbcHandling=elimination";
//       Error( __FILE__, __LINE__ );
//     }
    assemble_->GetRHSVal( rhs_, ptRhs, identifierPDE );
    return size_ * blockSize_;
  }


  // *************************
  //   GetRHSVal (Version 2)
  // *************************
  Integer StandardSystem::GetRHSVal( Complex* &ptRhs,
                                     const PdeIdType identifierPDE ) {
//     if ( numDirichletValues_ > 0 &&
//          myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {
//       (*error) << "GetRHSVal does not work for idbcHandling=elimination";
//       Error( __FILE__, __LINE__ );
//     }
    assemble_->GetRHSVal( rhs_, ptRhs, identifierPDE );
    return size_ * blockSize_;
  }
  

  // *********
  //   Print
  // *********
  void StandardSystem::Print( FEMatrixType matrix_id) const {
    sysmat_[matrix_id]->Print(*cla);
  }

  // **********
  //   Export
  // **********
  void StandardSystem::Export( FEMatrixType type, Char *filename,
                               Char *comment ) const {
    sysmat_[type]->Export( filename, comment );
  }


  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  void StandardSystem::AddToDiagMatrixEntry( FEMatrixType matrixID,
                                             const PdeIdType pdeID,
                                             Integer eqnNum,
                                             Double *val ) {


    // Delegate work to implementation in assemble class
    assemble_->AddToDiagMatrixEntry( sysmat_[matrixID], eqnNum, val );
  }
  
  // ******************
  //   GetMatrixEntry
  // ******************
  void StandardSystem::GetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, 
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum, 
                                       Double & val ) {


    // Delegate work to implementation in assemble class
    assemble_->GetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum, 
                               colEqnNum, numLastFreeDof_[1],
                               numLastFreeDof_[1] );
    
  }
  
  void StandardSystem::GetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, 
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum, 
                                       Complex & val ) {
    

    // Delegate work to implementation in assemble class
    assemble_->GetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum, 
                               colEqnNum,  numLastFreeDof_[1],
                               numLastFreeDof_[1] );
  }

  // ******************
  //   GetMatrixEntry
  // ******************

  void StandardSystem::SetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, 
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum, 
                                       Double val, bool setCounterPart ) {

     
    // Delegate work to implementation in assemble class
    assemble_->SetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum,
                               colEqnNum,  numLastFreeDof_[1],
                               numLastFreeDof_[1], setCounterPart );
  }
  void StandardSystem::SetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, 
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum,
                                       Complex val, bool setCounterPart ) {

     // Delegate work to implementation in assemble class
    assemble_->SetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum, 
                               colEqnNum, numLastFreeDof_[1],
                               numLastFreeDof_[1], setCounterPart );
  }

  // ******************************
  //   GenerateCFSTransferBuffers
  // ******************************
  void StandardSystem::GenerateCFSTransferBuffers() {


    BaseVector *bVec = NULL;

    // Determine total number of equation numbers for this system
    // including those for fixed dofs
    UInt totalSize = totalSize_;
    totalSize *= blockSize_;

    // Determine some basic properties
    MatrixStorageType sType = sysmat_[SYSTEM]->GetStorageType();
    MatrixEntryType eType = sysmat_[SYSTEM]->GetEntryType();

    // Generate vector for passing solution back to CFS++
    bVec = GenerateSparseVectorObject( sType, eType, 1, totalSize );
    solComm_ = dynamic_cast<SparseVector*>( bVec );

  }

}
