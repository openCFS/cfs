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

#include "utils/Exception.hh"

namespace OLAS {
 

  // ***********************
  //   Default Constructor
  // ***********************
  StandardSystem::StandardSystem() : BaseSystem() {

    ENTER_FCN( "StandardSystem::StandardSystem" );

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

    ENTER_FCN( "StandardSystem::~StandardSystem" );

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
    ENTER_FCN("StandardSystem::SetupPrecond");
    precond_->Setup( *sysmat_[SYSTEM] );
  }


  // *********************************
  //   Trigger setup phase of solver
  // *********************************
  void StandardSystem::SetupSolver( ) {
    ENTER_FCN( "StandardSystem::SetupSolver" );
    solver_->Setup( *sysmat_[SYSTEM] );
  }
  
  // **************************************
  //   Trigger setup phase of EigenSolver
  // **************************************
  void StandardSystem::SetupEigenSolver( UInt numFreq, Double shift, 
                                         bool isQuadratic ) {
    ENTER_FCN( "StandardSystem::SetupEigenSolver" );

    // Determine if a generalized or a quadratic eigenvalue 
    // problem has to be solved.
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

    BaseVector *bVec = GenerateStdVectorObject( sType, eType, 1, totalSize );
    BaseVector *errVec = GenerateStdVectorObject( sType, eType, 1, totalSize );
    eigenValues_ = dynamic_cast<StdVector*>( bVec );
    eigenValError_ = dynamic_cast<StdVector*>( errVec );
    
  }

  // ***********************
  //   Solve linear system
  // ***********************
  void StandardSystem::Solve(int step) {

    ENTER_FCN( "StandardSystem::Solve" );

    // If the penalty formulation is used and we have inhomogeneous
    // Dirichlet boundary conditions, then the righ-hand side is
    // "contaminated" with penalty terms
    if ( numDirichletValues_ > 0 &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      myParams_.SetValue( "RHSwithPenalty", false );
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
    idbcHandler_->AddIDBCToRHS( rhs_ );

    // Export linear system
    bool doExport = myParams_.GetBoolValue( "exportLinSys" );
    std::string file;
    std::string base;

    // need it common even when exclusive solution
    if(doExport) {
      std::ostringstream os;
      os << myParams_.GetStringValue("exportLinSysBaseName");
      if(step != -1) os << "_" << step;
      base = os.str();  
    }

    if (doExport && myParams_.GetStringValue("exportLinSysSolution") != "exclusive") 
    {
      // two formats
      if(myParams_.GetStringValue("exportLinSysFormat") == "harwell-boeing") 
      {
          file = base +  ".hb";
          sysmat_[SYSTEM]->ExportHarwellBoeing(file, *rhs_);
          (*cla) << " Export system matrix and RHS to " << file << std::endl;
      } 
      else // classical (default) matrix-market 
      {
          file = base + ".mtx";
          sysmat_[SYSTEM]->Export( file.c_str() );
          (*cla) << " Export system matrix to " << file << std::endl;
    
          file = base + ".vec";
          rhs_->Export( file.c_str() );
          (*cla) << " Export system right hand side to " << file << " euclidian norm=" << rhs_->NormEuclid() << std::endl;
      }    
    }


    // Trigger solution
    // as long as olas is not merges with cfs and I'm not able to let cfs know the olas
    // exception catch it manually and transfer
    try
    {
       solver_->Solve( *sysmat_[SYSTEM], *precond_, *rhs_, *sol_ );
    }
    catch(Exception& e)
    {
        throw e.ToString(); // cfs knows this type
    }   
    

    if (doExport && myParams_.GetStringValue("exportLinSysSolution") != "no") {
      file = base + ".sol.vec";
      sol_->Export(file.c_str());
      (*cla) << " Export system solution to " << file << " euclidian norm=" << sol_->NormEuclid() << std::endl;
    }

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

    // Set invaldiation flag for assembling routine
    assemble_->InvalidateSolBuffer();
  }

  // ******************************
  //   Calculate Eigenfrequencies
  // ******************************
  UInt StandardSystem::
  CalcEigenFrequencies( const Double * &frequencies, const Double* &err  ) {

    ENTER_FCN( "StandardSystem::CalcEigenFrequencies" );
    
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

    ENTER_FCN( "StandardSystem::CalcEigenFrequencies" );
    
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
    
    ENTER_FCN( "StandardSystem::CalcEigenMode" );
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

    ENTER_FCN( "StandardSystem::CreateLinSys" );

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

#ifdef DEBUG_STANDARDSYSTEM  
    (*debug) << "size " << size_ << std::endl;
    (*debug) << "nne  " << nne_ << std::endl;
    (*debug) << "Matrix Graph: " << std::endl;
    graph->Print(*debug);
#endif


    // ------------------------
    //  Re-set numLastFreeDof
    // ------------------------
    DeleteArray( numLastFreeDof_ );
    NewArray( numLastFreeDof_, UInt, 1 );
    numLastFreeDof_[1] = graph->GetSize();

#ifdef DEBUG_STANDARDSYSTEM
    (*debug) << " StandardSystem::CreateLinSys:"
             << " Re-set numLastFreeDof to " << numLastFreeDof_[1]
             << std::endl;
#endif

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

#ifdef DEBUG_STANDARDSYSTEM
    (*debug) << "+++ CreateLinSys +++" << std::endl;
#endif

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
    for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ ) {

      MatrixStorageType actStorageType;

#ifdef DEBUG_STANDARDSYSTEM
      (*debug) << "Generating matrix[" << *it << "]"
               << "\nMatrixEntryType     = " << entryType
               << "\nMatrixStorageType   = " << storageType
               << "\nDegrees of freedom  = " << blockSize_
               << "\nNumber of rows      = " << size_
               << "\nNumber of columns   = " << size_
               << "\nNumber of non-zeros = " << nne_
               << std::endl;
#endif

      actStorageType = storageType;
      if ( solver == DIAGSOLVER ) {
	if ( *it == SYSTEM || *it == MASS ) {
	  // just allocate a diagonal matrix
	  actStorageType = DIAG;
	}
      }

#ifdef DEBUG_STANDARDSYSTEM
        (*debug) << "calling GenerateStdMatrixObject" << std::endl;
#endif

        sysmat_[*it] = GenerateStdMatrixObject( entryType, actStorageType,
                                                blockSize_, size_, size_,
                                                nne_ );
	// }

#ifdef DEBUG_STANDARDSYSTEM
      if ( sysmat_[*it] == NULL ) {
        Error( "Matrix generation failed!", __FILE__, __LINE__ );
      }
#endif

      // Generate sparsity pattern once, but share it afterwards
      // if this is possible
      if ( sharePattern == false ) {
        sysmat_[*it]->SetSparsityPattern( (*graph) );
        if ( sharePatternPossible == true ) {
          sharePattern = true;
          patternID = sysmat_[*it]->TransferPatternToPool( patternPool_ );
        }
      }
      else {
        sysmat_[*it]->SetSparsityPattern( patternPool_, patternID );
      }
    }

   
    // ---------------------------------
    //  Generation of remaining objects
    // ---------------------------------

    // Generate vectors
    rhs_ = dynamic_cast<StdVector *>
      (GenerateVectorObject( *(sysmat_[SYSTEM]) ));
    sol_ = dynamic_cast<StdVector *>
      (GenerateVectorObject( *(sysmat_[SYSTEM]) ));

    // Generate communication buffers
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {
      GenerateCFSTransferBuffers();
    }

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
     // remove this when merging olas with cfs!!
     try
     {
        ENTER_FCN("StandardSystem::CreateSolver");
        
        SolverType solver;
        myParams_.GetEnumValue( "Solver", solver );
        solver_ = GenerateSolverObject( *(sysmat_[SYSTEM]), solver, &myParams_,
                                        &myReport_ );
     }
     catch(Exception& e)
     {
        // cfs doesn't know our olas exception
        throw e.ToString();
     }
  }


  // ********************************
  //   Create preconditioner object
  // ********************************

  // Das Interface muessen wir uns noch ueberlegen
  void StandardSystem::CreatePrecond() {

    ENTER_FCN("StandardSystem::CreatePrecond");
    
    PrecondType precond;
    myParams_.GetEnumValue("Precond", precond);
    precond_ = GenerateStdPrecondObject( *(sysmat_[SYSTEM]), precond,
                                         &myParams_, &myReport_ );

  }

  // ********************************
  //   Create eigensolver object
  // ********************************

  void StandardSystem::CreateEigenSolver() {

    ENTER_FCN("StandardSystem::CreateEigenSolver");
    
    EigenSolverType egSolver;
    myParams_.GetEnumValue("EigenSolver", egSolver);
    eigenSolver_ = GenerateEigenSolverObject( *(sysmat_[SYSTEM]), egSolver,
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
                                       const PdeIdType identifier2) {
    ENTER_FCN( "StandardSystem::SetFEMatrixType" );

    if( matType != NOTYPE ) {
      matrixTypes_.insert( matType );
    }
  }


  // ****************
  //   SetBlockSize
  // ****************
  void StandardSystem::SetBlockSize( const PdeIdType identifier,
                                     const UInt bs ) {
    ENTER_FCN( "StandardSystem::SetNumDof" );
    blockSize_ = bs;
  }



  // =======================================================================
  //   Methods for assembling, setting BCs and rhs
  // =======================================================================
  


  // **************
  //   InitMatrix
  // **************
  void StandardSystem::InitMatrix( FEMatrixType matrixID,
                                   const PdeIdType identifierPDE ) {

    ENTER_FCN( "StandardSystem::InitMatrix" );

    // If matrix specified init this one
    if ( matrixID != NOTYPE ) {
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
    ENTER_FCN("StandardSystem::InitRHS");
    rhs_->Init();

    // Invalidate rhs buffer of assemble
    assemble_->InvalidateRhsBuffer();
  }


  // ***********************
  //   InitRHS (Version 2)
  // ***********************
  void StandardSystem::InitRHS( const Double *newRHS ) {
    ENTER_FCN( "StandardSystem::InitRHS" );
    assemble_->InitRHS( rhs_, newRHS );
  }


  // ***********
  //   InitSol
  // ***********
  void StandardSystem::InitSol( const PdeIdType identifierPDE ) {
    ENTER_FCN("StandardSystem::InitSol");
    sol_->Init();

    // Invalidate sol buffer of assemble
    assemble_->InvalidateSolBuffer();
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

    ENTER_IFCN( "StandardSystem::SetElementMatrix" );


    // The following section is obsolete alltogether, as
    // in CFS++ the flag 'setCounterPart' is really only
    // set for non-diagonal matrices, where the transposed
    // has to be set as well. Previously the default for
    // 'setCounterPart' was by default 'true', even for
    // on-diagonal matrices.

    // // Set flag for setting the symmetric counter-part of
//     // the element matrix
//     if ( setCounterPart == true ) {
//       if ( idPDE1 == idPDE2 ) {

//         // NOTE: Changed the default to true, since we have no way to handle
//         //       it in CFS++ now. Thus, warning was commented out
//         // (*warning) << "StandardSystem::SetElementMatrix: Refusing to set "
//         //            << "symmetric counterpart for PDE with ID = "
//         //            << identifierPDE1;
//         // Warning( __FILE__, __LINE__ );
//         setCounterPart = false;
//       }
//     }

    // Delegate remaining work to assemble
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

    ENTER_IFCN( "StandardSystem::SetElementRHS" );

    // Delegate work to EntryManipulator
    assemble_->SetElementRHS( rhs_, elemRHS, connect, length, size_ );
  }


  // ***********************
  //   SetNodeRHS (Double)
  // ***********************
  void StandardSystem::SetNodeRHS( Double val, const PdeIdType identifierPDE,
                                   Integer eqnNum, Integer comp ) {

    ENTER_IFCN( "StandardSystem::SetNodeRHS" );

    // Perform consistency check
#ifdef DEBUG_STANDARDSYSTEM
    if ( eqnNum > size_ || eqnNum < 1 ) {
      (*error) << "StandardSystem::SetNodeRHS: Inconsistency detected:"
               << "\n size      = " << size_
               << "\n totalSize = " << totalSize_
               << "\n eqnNum    = " << eqnNum
               << "\n comp      = " << comp
               << "\n val       = " << val
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to entry manipulation object
    assemble_->SetNodeRHS( rhs_, val, eqnNum, comp );
  }

  void StandardSystem::SetNodeRHS( Complex val, const PdeIdType identifierPDE,
                                   Integer eqnNum, Integer comp ) {

    ENTER_IFCN( "StandardSystem::SetNodeRHS" );
    
    // Perform consistency check
#ifdef DEBUG_STANDARDSYSTEM
    if ( eqnNum > size_ || eqnNum < 1 ) {
      (*error) << "StandardSystem::SetNodeRHS: Inconsistency detected:"
               << "\n size      = " << size_
               << "\n totalSize = " << totalSize_
               << "\n eqnNum    = " << eqnNum
               << "\n comp      = " << comp
               << "\n val       = " << val
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'";
      Error( __FILE__, __LINE__ );
    }
#endif
    
    // Delegate work to entry manipulation object
    assemble_->SetNodeRHS( rhs_, val, eqnNum, comp );
  }


  // ************************
  //   SetNodeRHS (Complex)
  // ************************

  // **************
  //   Update RHS
  // **************  
  void StandardSystem::UpdateRHS( FEMatrixType matrix_id, Double *fup ) {
    ENTER_FCN( "StandardSystem::UpdateRHS" );
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

    ENTER_FCN("StandardSystem::RemoveIDBCFromMatrix" );
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
  void StandardSystem::SetDirichlet( UInt bcNum, const PdeIdType pdeID,
                                     Integer eqnNum, const Double &val,
                                     UInt comp ) {

    ENTER_IFCN( "StandardSystem::SetDirichlet" );

    // add offset to the boundary value number
    bcNum += bcOffsets_[ pdeID ];

    // Perform some consistency checks
#ifdef DEBUG_STANDARDSYSTEM

    // Some situation dependend modifiers
    bool usingPenalty = myParams_.GetBoolValue( "UsingPenaltyFormulation" );
    UInt maxVal = usingPenalty == true ? size_ : size_ + numDirichletValues_;
    UInt minVal = usingPenalty == true ? 1 : size_ + 1;

    if ( eqnNum > maxVal || eqnNum < minVal || comp < 1 || comp > 4 ||
         bcNum > numDirichletValues_ ) {
      (*error) << "StandardSystem::SetDirichlet: Inconsistency detected:"
               << "\n bcNum  = " << bcNum
               << "\n pdeID  = " << pdeID
               << "\n eqnNum = " << eqnNum
               << "\n val    = " << val
               << "\n comp   = " << comp
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
    idbcHandler_->SetIDBC( pdeID, eqnNum, comp, bcNum, val );
  }


  // **************************
  //   SetDirichlet (Complex)
  // **************************
  void StandardSystem::SetDirichlet( UInt bcNum, const PdeIdType pdeID,
                                     Integer eqnNum, const Complex &val,
                                     UInt comp ) {

    ENTER_IFCN( "StandardSystem::SetDirichlet" );

    // add offset to the boundary value number
    bcNum += bcOffsets_[ pdeID ];

    // Perform some consistency checks
#ifdef DEBUG_STANDARDSYSTEM

    // Some situation dependend modifiers
    bool usingPenalty = myParams_.GetBoolValue( "UsingPenaltyFormulation" );
    UInt maxVal = usingPenalty == true ? size_ : size_ + numDirichletValues_;
    UInt minVal = usingPenalty == true ? 1 : size_ + 1;

    if ( eqnNum > maxVal || eqnNum < minVal || comp < 1 || comp > 4 ||
         bcNum > numDirichletValues_ ) {
      (*error) << "StandardSystem::SetDirichlet: Inconsistency detected:"
               << "\n bcNum  = " << bcNum
               << "\n pdeID  = " << pdeID
               << "\n eqnNum = " << eqnNum
               << "\n val    = " << val
               << "\n comp   = " << comp
               << "\n numBC  = " << numDirichletValues_
               << "\n size   = " << size_
               << "\n minVal = " << minVal
               << "\n maxVal = " << maxVal;
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( pdeID, eqnNum, comp, bcNum, val );
  }


  // *******************
  //   BuildInDirichlet
  // *******************
  void StandardSystem::BuildInDirichlet() {

    ENTER_IFCN( "StandardSystem::BuildInDirichlet" );

    // If necessary modify matrix diagonal for penalty approach
    if ( assembleDirichletToSysMat_ == true ) {
      idbcHandler_->AdaptSystemMatrix( *(sysmat_[SYSTEM]) );
      assembleDirichletToSysMat_ = false;
    }
  }


  // ****************************
  //   ConstructEffectiveMatrix
  // ****************************
  void StandardSystem::
  ConstructEffectiveMatrix( const factorMap &matFactors ) {

    ENTER_FCN( "StandardSystem::ConstructEffectiveMatrix" );

    factorMap::const_iterator it;
    StdMatrix *sys = sysmat_[SYSTEM];

    // It's okay, if there are no factors, if there is only a system
    // matrix and no other ones
    if ( matFactors.empty() == true ) {
      if ( matrixTypes_.size() == 1 && sys != NULL ) {
#ifdef DEBUG_STANDARDSYSTEM
        (*debug) << " ConstructEffectiveMatrix: Nothing to be done, since"
                 << " there is only one FE = system matrix\n\n";
#endif

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

#ifdef DEBUG_STANDARDSYSTEM
    (*debug) << "ConstructEffectiveMatrix with factors: " << std::endl;
    
    for ( it = matFactors.begin(); it != matFactors.end(); it++ ) {
      (*debug) << Enum2String((*it).first) << ":" 
               << (*it).second << std::endl;
    }
#endif

    // Intialize system matrix before adding the weighted
    // matrices to it
    this->InitMatrix( SYSTEM );

    for ( it = matFactors.begin(); it != matFactors.end(); it++ ) {
      if ( sysmat_[(*it).first] != NULL  && (*it).second != 0.0 ) {
        sys->Add( (*it).second, *sysmat_[(*it).first] );
      }
    }
    
#ifdef DEBUG_STANDARDSYSTEM
    (*debug) << "effective Matrix: " << std::endl;
    sys->Print(*debug);
#endif

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

    ENTER_FCN("BaseSystem::GetSolutionVal");

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

    ENTER_FCN( "BaseSystem::GetSolutionVal" );

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
    ENTER_FCN( "StandardSystem::GetRHSVal" );
    if ( numDirichletValues_ > 0 &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {
      (*error) << "GetRHSVal does not work for idbcHandling=elimination";
      Error( __FILE__, __LINE__ );
    }
    assemble_->GetRHSVal( rhs_, ptRhs, identifierPDE );
    return size_ * blockSize_;
  }


  // *************************
  //   GetRHSVal (Version 2)
  // *************************
  Integer StandardSystem::GetRHSVal( Complex* &ptRhs,
                                     const PdeIdType identifierPDE ) {
    ENTER_FCN("StandardSystem::GetRHSVal");
    if ( numDirichletValues_ > 0 &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {
      (*error) << "GetRHSVal does not work for idbcHandling=elimination";
      Error( __FILE__, __LINE__ );
    }
    assemble_->GetRHSVal( rhs_, ptRhs, identifierPDE );
    return size_ * blockSize_;
  }
  

  // *********
  //   Print
  // *********
  void StandardSystem::Print( FEMatrixType matrix_id) const {
    ENTER_FCN( "StandardSystem::Print" );
    sysmat_[matrix_id]->Print(*cla);
  }

  // ********
  //  returns system matrix (if stored in scrs!) as a vector
  //  containing all upper triangle entries 
  // ********
  void StandardSystem::GetFullSystemMatrixAsVec(Complex*  &vec ) const {
    ENTER_FCN("StandardSystem::GetFullSystemMatrix");
     
    sysmat_[1]->CopySCRSMatrix2Vec(vec);

    // needs further programming for crs, ... matrices                         
  }



  // **********
  //   Export
  // **********
  void StandardSystem::Export( FEMatrixType type, Char *filename,
                               Char *comment ) const {
    ENTER_FCN( "StandardSystem::Export" );
    sysmat_[type]->Export( filename, comment );
  }


  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  void StandardSystem::AddToDiagMatrixEntry( FEMatrixType matrixID,
                                             const PdeIdType pdeID,
                                             Integer eqnNum, UInt dof,
                                             Double *val ) {

    ENTER_IFCN( "StandardSystem::AddToDiagMatrixEntry" );

    // Perform consistency checks
#ifdef DEBUG_STANDARDSYSTEM
    if ( dof > (UInt)blockSize_ || dof < 1 ) {
      (*error) << "StandardSystem::AddToDiagMatrixEntry: Parameter dof = "
               << dof << ", but must be in [1,blockSize_] with "
               << "blockSize_ = " << blockSize_;
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to implementation in assemble class
    assemble_->AddToDiagMatrixEntry( sysmat_[matrixID], eqnNum, dof, val );
  }
  
  // ******************
  //   GetMatrixEntry
  // ******************
  void StandardSystem::GetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, UInt rowDof,
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum, UInt colDof,
                                       Double & val ) {

    ENTER_IFCN( "StandardSystem::GetMatrixEntry");

    // Perform consistency checks
#ifdef DEBUG_STANDARDSYSTEM
    if ( dof > (UInt)blockSize_ || rowDof < 1 || colDof < 1) {
      (*error) << "StandardSystem::GetMatrixEntry: Parameter rowDof = "
               << rowDof << " and colDof = " << colDof 
               << ", but must be in [1,blockSize_] with "
               << "blockSize_ = " << blockSize_;
      Error( __FILE__, __LINE__ );
    }
#endif
    
    // Delegate work to implementation in assemble class
    assemble_->GetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum, rowDof,
                               colEqnNum, colDof,  numLastFreeDof_[1],
                               numLastFreeDof_[1] );
    
  }
  
  void StandardSystem::GetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, UInt rowDof,
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum, UInt colDof,
                                       Complex & val ) {
    ENTER_IFCN( "StandardSystem::GetMatrixEntry");
    

    // Perform consistency checks
#ifdef DEBUG_STANDARDSYSTEM
    if ( dof > (UInt)blockSize_ || rowDof < 1 || colDof < 1) {
      (*error) << "StandardSystem::GetMatrixEntry: Parameter rowDof = "
               << rowDof << " and colDof = " << colDof 
               << ", but must be in [1,blockSize_] with "
               << "blockSize_ = " << blockSize_;
      Error( __FILE__, __LINE__ );
    }
#endif
    
    // Delegate work to implementation in assemble class
    assemble_->GetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum, rowDof,
                               colEqnNum, colDof,  numLastFreeDof_[1],
                               numLastFreeDof_[1] );
  }

  // ******************
  //   GetMatrixEntry
  // ******************

  void StandardSystem::SetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, UInt rowDof,
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum, UInt colDof,
                                       Double val, bool setCounterPart ) {
    ENTER_IFCN( "StandardSystem::SetMatrixEntry");

    // Perform consistency checks
#ifdef DEBUG_STANDARDSYSTEM
    if ( dof > (UInt)blockSize_ || rowDof < 1 || colDof < 1) {
      (*error) << "StandardSystem::SetMatrixEntry: Parameter rowDof = "
               << rowDof << " and colDof = " << colDof 
               << ", but must be in [1,blockSize_] with "
               << "blockSize_ = " << blockSize_;
      Error( __FILE__, __LINE__ );
    }
#endif
    
    // Delegate work to implementation in assemble class
    assemble_->SetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum, rowDof,
                               colEqnNum, colDof,  numLastFreeDof_[1],
                               numLastFreeDof_[1], setCounterPart );
  }
  void StandardSystem::SetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum, UInt rowDof,
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum, UInt colDof,
                                       Complex val, bool setCounterPart ) {
    ENTER_IFCN( "StandardSystem::SetMatrixEntry");

        // Perform consistency checks
#ifdef DEBUG_STANDARDSYSTEM
    if ( dof > (UInt)blockSize_ || rowDof < 1 || colDof < 1) {
      (*error) << "StandardSystem::SetMatrixEntry: Parameter rowDof = "
               << rowDof << " and colDof = " << colDof 
               << ", but must be in [1,blockSize_] with "
               << "blockSize_ = " << blockSize_;
      Error( __FILE__, __LINE__ );
    }
#endif
    
    // Delegate work to implementation in assemble class
    assemble_->SetMatrixEntry( matrixID, rowPdeID, colPdeID, 
                               sysmat_[matrixID], idbcHandler_,
                               val, rowEqnNum, rowDof,
                               colEqnNum, colDof,  numLastFreeDof_[1],
                               numLastFreeDof_[1], setCounterPart );
  }

  // ******************************
  //   GenerateCFSTransferBuffers
  // ******************************
  void StandardSystem::GenerateCFSTransferBuffers() {

    ENTER_FCN( "StandardSystem::GenerateCFSTransferBuffers" );

    BaseVector *bVec = NULL;

    // Determine total number of equation numbers for this system
    // including those for fixed dofs
    UInt totalSize = totalSize_;
    totalSize *= blockSize_;

    // Determine some basic properties
    MatrixStorageType sType = sysmat_[SYSTEM]->GetStorageType();
    MatrixEntryType eType = sysmat_[SYSTEM]->GetEntryType();

    // Generate vector for passing solution back to CFS++
    bVec = GenerateStdVectorObject( sType, eType, 1, totalSize );
    solComm_ = dynamic_cast<StdVector*>( bVec );

  }

}
