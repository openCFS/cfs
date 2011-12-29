// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <algorithm>
#include <complex>
#include <fstream>
#include <set>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/basevector.hh"
#include "MatVec/generatematvec.hh"
#include "MatVec/patternpool.hh"
#include "MatVec/stdmatrix.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"
#include "OLAS/algsys/baseidbchandler.hh"
#include "OLAS/algsys/generateidbchandler.hh"
#include "OLAS/algsys/standardsys.hh"
#include "OLAS/graph/basegraph.hh"
#include "OLAS/graph/basegraphmanager.hh"
#include "OLAS/precond/baseprecond.hh"
#include "OLAS/precond/generateprecond.hh"
#include "OLAS/solver/baseEigensolver.hh"
#include "OLAS/solver/basesolver.hh"
#include "OLAS/solver/generateEigensolver.hh"
#include "OLAS/solver/generatesolver.hh"
#include "Utils/Timer.hh"
#include "boost/algorithm/string/replace.hpp"

namespace CoupledField {
template <class TYPE> class Matrix;
}  // namespace CoupledField

using std::string;

DECLARE_LOG(stdSys)
DEFINE_LOG(stdSys, "stdSys")

namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  StandardSystem::StandardSystem(PtrParamNode pn) : BaseSystem(pn) {


    totalSize_           = 0;
    precond_             = NULL;
    rhs_                 = NULL;
    sol_                 = NULL;
    solComm_             = NULL;
    patternPool_         = NULL;
    eigenValues_         = NULL;
    algSysType_          = STANDARD_SYSTEM;

    // Initialize pointer to finite-element matrices
    sysmat_[SYSTEM]     = NULL;
    sysmat_[STIFFNESS]  = NULL;
    sysmat_[DAMPING]    = NULL;
    sysmat_[CONVECTION] = NULL;
    sysmat_[MASS]       = NULL;
    sysmat_[AUXILIARY]  = NULL;

    if ( usingPenalty_ ) {
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
    std::map<FEMatrixType, StdMatrix*>::iterator iter;
    for (iter = sysmat_.begin(); iter != sysmat_.end(); ++iter)
    {
      delete iter->second;
    }
    sysmat_.clear();

    // After matrices we may delete the pattern pool
    delete patternPool_;
  }


  // *****************************************
  //   Trigger setup phase of preconditioner
  // *****************************************
  void StandardSystem::SetupPrecond() 
  {
    info->ToFile();
    precond_->Setup( *sysmat_[SYSTEM]);
  }


  // *********************************
  //   Trigger setup phase of solver
  // *********************************
  void StandardSystem::SetupSolver(PtrParamNode analysis_id) {
    solver_->GetSetupTimer()->Start();
    solver_->Setup(*sysmat_[SYSTEM], analysis_id);
    solver_->GetSetupTimer()->Stop();
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
    bool massPresent = (matrixTypes_.find( MASS) != matrixTypes_.end());

    if( isQuadratic == true ) {
      if( dampPresent == false ) {
        EXCEPTION("To solve a quadratic eigenvalue problem, a damping " \
                  << "matrix has to be present!");
      }

      // Setup the quadratic eigenvalue solver
      eigenSolver_->Setup( *sysmat_[STIFFNESS], *sysmat_[MASS],
                           *sysmat_[DAMPING], numFreq, shift );
    } else {
      if( dampPresent == true ) {
        WARN("Although a damping matrix is present, only a generalized "
             << "eigenvalue will be solved (as you have specified)!");
      }
      
      if( massPresent == true ) {
        // Setup the eigenvalue solver for generalized EV problem
        eigenSolver_->Setup( *sysmat_[STIFFNESS], *sysmat_[MASS],
                             numFreq, shift );
      } else {
        // Setup the eigenvalue solver for standard EV problem
        eigenSolver_->Setup( *sysmat_[STIFFNESS], 
                             numFreq, shift );
      }
    }

    // Determine some basic properties and create vectors
    // for eigenvalues and related error bounds
    BaseMatrix::EntryType eType;
    if( isQuadratic == true ) {
      eType = BaseMatrix::COMPLEX;
    } else {
      eType = BaseMatrix::DOUBLE;
    }

    BaseMatrix::StorageType sType = sysmat_[SYSTEM]->GetStorageType();
    UInt totalSize = totalSize_;
    totalSize *= blockSize_;

    BaseVector *bVec = GenerateSingleVectorObject( sType, eType, totalSize );
    BaseVector *errVec = GenerateSingleVectorObject( sType, eType, totalSize );
    eigenValues_ = dynamic_cast<SingleVector*>( bVec );
    eigenValError_ = dynamic_cast<SingleVector*>( errVec );

  }

  // ***********************
  //   Solve linear system
  // ***********************
  void StandardSystem::Solve(PtrParamNode analysis_id) 
  {
    info->ToFile(); // write current info state
    solver_->GetSolveTimer()->Start();
    
    // Check, if condition number is to be calculated
    bool calcCapa = false;
    xml_->Get("matrix", ParamNode::INSERT)
    ->GetValue("calcConditionNumber", calcCapa, ParamNode::INSERT );
    
    if( calcCapa ) {
      Double condNumber = 0.0;
      Vector<Double> ev, err;
      BaseEigenSolver * evs = 
        GenerateEigenSolverObject( *(sysmat_[SYSTEM]), xml_, olasInfo_->Get("solve_eigen") );
      PtrParamNode in = systemInfo_->Get(ParamNode::PROCESS)->Get("conditionNumber", ParamNode::APPEND);
      in->Get("analysis_id")->SetValue(analysis_id);
      try {
        evs->CalcConditionNumber( *(sysmat_[SYSTEM]), condNumber,
                                  ev, err );
        in->Get("value")->SetValue(condNumber);
        in = in->Get("extremalEigenValues");

        for( UInt i = 0; i < ev.GetSize(); i++ )  {
          PtrParamNode t = in->Get("eigenvalue", ParamNode::APPEND);
          t->Get("value")->SetValue(ev[i]);
          t->Get("tolerance")->SetValue(err[i]);
        }

      } catch (Exception& ex ) {
        WARN("Calculation of condition number for system matrix did not converge");
      
        in->Get("value")->SetValue(0.0);
      }
      
      delete evs;
      // in the end prevent-relcaulation of evs by re-setting the value for CalcConditionNumber
      xml_->Get("matrix")->Get("calcConditionNumber")->SetValue( boost::any(false) );
    }
   
    // If the penalty formulation is used and we have inhomogeneous
    // Dirichlet boundary conditions, then the righ-hand side is
    // "contaminated" with penalty terms
    if ( numDirichletValues_ > 0 && usingPenalty_ ) {
      solver_->SetUsingPenalty(true);
      LOG_DBG(stdSys) << "Solve: rhs with penalty";
    }

    // Iterative solvers require an initial guess and in the penalty case
    // we should insert the Dirichlet values into it
    if ( dynamic_cast<BaseIterativeSolver*>(solver_) != NULL &&
         usingPenalty_ ) {
      idbcHandler_->SetDofsToIDBC( sol_ );
      PtrParamNode in = systemInfo_->Get(ParamNode::HEADER)->Get("idbc");
      in->Get("penalty")->SetValue(true);
      in->SetComment("Inserted Dirichlet values into initial guess");
      LOG_DBG(stdSys) << "Solve: Inserted Idbc values into initial guess for iterative solver";
    }

    // Perform a simple sanity check
    if ( sysmat_[SYSTEM] == NULL || precond_ == NULL || rhs_ == NULL ||sol_ == NULL ) {
      EXCEPTION( "Detected NULL pointer where there should be none!" );
    }

    // Assume that everything will go well
    PtrParamNode out = olasInfo_->Get(ParamNode::PROCESS)->Get("solver");
    out->Get("solutionIsOkay")->SetValue(true);

    // Now modifiy the right-hand side vector
    LOG_DBG(stdSys) << "Solve: add idbc to rhs";
    idbcHandler_->AddIDBCToRHS( rhs_ );

    // check if we do export stuff
    PtrParamNode els =  xml_->Get("exportLinSys", ParamNode::PASS);
    string file;
    string base;
    
    // need it common even when exclusive solution
    if(els) {
      string id = analysis_id->Get("analysis_id")->As<std::string>();
      boost::replace_all(id, ":", "_");
      string name = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
      base = name + "_" + id;
      systemInfo_->Get(ParamNode::PROCESS)->Get("exportLinearSystem", ParamNode::APPEND)->Get("name")->SetValue(base);
    }

    // check if we do not only want the solution
    if(els && els->Get("solution")->As<std::string>() != "exclusive")
    {
      // two formats. The harwell-boing format includes the rhs!
      if(els->Get("format")->As<std::string>()  == "harwell-boeing")
      {
        sysmat_[SYSTEM]->ExportHarwellBoeing(base+".hb", *rhs_);

        if(els->HasByVal("damping", true) && sysmat_[DAMPING] != NULL)
          sysmat_[DAMPING]->ExportHarwellBoeing(base+"_damping.hb", *rhs_);

        if(els->HasByVal("auxiliary", true) && sysmat_[AUXILIARY] != NULL)
          sysmat_[AUXILIARY]->ExportHarwellBoeing(base+"_aux.hb", *rhs_);
      }
      else // classical (default) matrix-market
      {
        sysmat_[SYSTEM]->Export((base+".mtx").c_str() );

        if(els->HasByVal("stiffness", true) && sysmat_[STIFFNESS] != NULL)
          sysmat_[STIFFNESS]->Export((base+"_stiffness.mtx").c_str() );

        if(els->HasByVal("damping", true) && sysmat_[DAMPING] != NULL)
          sysmat_[DAMPING]->Export((base+"_damping.mtx").c_str() );

        if(els->HasByVal("mass", true) && sysmat_[MASS] != NULL)
          sysmat_[MASS]->Export((base+"_mass.mtx").c_str() );

        if(els->HasByVal("auxiliary", true) && sysmat_[AUXILIARY] != NULL)
          sysmat_[AUXILIARY]->Export((base+"_aux.mtx").c_str() );

        // rhs is only in harwell-boing included
        rhs_->Export((base+".vec").c_str() );
      }
    }

    if(els && els->HasByVal("initialGuess", true))
      sol_->Export((base+"_intial_guess.vec").c_str());

    LOG_TRACE2(stdSys) << "before solve: euclidian norm of initial guess = " << sol_->NormL2();
    LOG_DBG3(stdSys) << "initial guess = " << sol_->ToString();
    LOG_DBG(stdSys) << "euclidian norm of rhs = " << rhs_->NormL2();
    LOG_DBG3(stdSys) << "rhs = " << rhs_->ToString();

    // ---------------------------
    // This is the expensize part! solve the system
    solver_->Solve( *sysmat_[SYSTEM], *precond_, *rhs_, *sol_, analysis_id);

    LOG_TRACE2(stdSys) << "after solve: euclidian norm of solution = " << sol_->NormL2();
    LOG_DBG3(stdSys) << "after solve: solution = " << sol_->ToString();

    if(els && els->Get("solution")->As<std::string>() != "no")
      sol_->Export((base+".sol.vec").c_str());

    LOG_DBG(stdSys) << "Solve: remove idbc from rhs";
    // Now de-modifiy the right-hand side vector
    idbcHandler_->RemoveIDBCFromRHS( rhs_ );


    // Check that solution went fine, if not issue a warning
    if ( out->Get("solutionIsOkay")->As<bool>() == false ) {
      // std::cerr << "Life's a piece of shit, when you look at it ...\n";
      WARN("Solver reports a problem! Consult .info.xml file for "
           << "further diagnostics!");
    }

    solver_->GetSolveTimer()->Stop();
    LOG_DBG(stdSys) << "Solve invalidate flag for solution buffer, why not also fro rhs?";
  }

  // ******************************
  //   Calculate Eigenfrequencies
  // ******************************
  void StandardSystem::
  CalcEigenFrequencies(  Vector<Double>& frequencies,
                         Vector<Double>& err ) {

    // Trigger calculation of eigenvalues
    eigenSolver_->CalcEigenFrequencies( *eigenValues_, *eigenValError_ );

    // Hard coded cast for double values
    Vector<Double>& valVec = dynamic_cast<Vector<Double>&>(*eigenValues_);
    frequencies = valVec;

    Vector<Double>& errVec = dynamic_cast<Vector<Double>&>(*eigenValError_);
    err = errVec;

  }

  void StandardSystem::
  CalcEigenFrequencies(  Vector<Complex>& frequencies,
                         Vector<Double>& err ) {

    
    // Check, if eigenvalue solver is quadratic, as only in this case
    // this method is well-defined

    if( eigenSolver_->IsQuadratic() == false ) {
      EXCEPTION("When solving a generalized eigenvalue problem, only " \
               << "real-valued results are obtained! Use the second " \
               << "CalcEigenFrequencies()-method!");
    }

    // Trigger calculation of eigenvalues
    eigenSolver_->CalcEigenFrequencies( *eigenValues_, *eigenValError_ );

    // Hard coded cast for double values
    Vector<Complex>& valVec = dynamic_cast<Vector<Complex>&>(*eigenValues_);
    frequencies = valVec;

    Vector<Double>& errVec = dynamic_cast<Vector<Double>&>(*eigenValError_);
    err = errVec;
  }


  // ******************************
  //   Calculate Eigenmodes
  // ******************************
  void StandardSystem::CalcEigenMode( UInt numMode ) {

    Vector<Double> & solHelp =
      dynamic_cast<Vector<Double> &> (*sol_);
    eigenSolver_->CalcEigenMode( numMode, solHelp );
    
    // DEBUG: Make check, if mode is really an eigenmode
    // by calculating (A - lambda * I) * u = 0
    // Vector<Double>& valVec = dynamic_cast<Vector<Double>&>(*eigenValues_);
//    Double lambda = (2.0*PI*valVec[numMode])*(2.0*PI*valVec[numMode]);
    
//    StdMatrix * mat = 
//    CopyStdMatrixObject( *sysmat_[STIFFNESS] );
//    // calculate (A - lambda * I)
//    std::cerr << "lambda = " << lambda << std::endl;
//    std::cerr << "mat before is " << mat->ToString() << std::endl;
//    Double diag = 0.0;
//    for( UInt i = 0; i < mat->GetNumRows(); i++ ) {
//      mat->GetDiagEntry(i, diag);
//      diag = diag - lambda;
//      mat->SetDiagEntry(i,diag);
//    }
//    std::cerr << "mat after is " << mat->ToString() << std::endl;
//    
//    Vector<Double> rhs(mat->GetNumRows());
//    mat->Mult( solHelp, rhs);
//    std::cerr << "Check for " << numMode << "th eigenmode gives\n";
//    std::cerr << rhs << std::endl;
  }


  // ************************************
  //   Creation of matrices and vectors
  // ************************************
  void StandardSystem::CreateLinSys() {


    BaseMatrix::EntryType   entryType;
    BaseMatrix::StorageType storageType;

    std::string entryTypeStr = "double";
    std::string storageTypeStr = "";
    
    PtrParamNode matrixNode = xml_->Get("matrix", ParamNode::INSERT);
    matrixNode->GetValue("entry", entryTypeStr, ParamNode::INSERT);
    matrixNode->GetValue("storage", storageTypeStr, ParamNode::INSERT);
    
    entryType = BaseMatrix::entryType.Parse( entryTypeStr );
    storageType = BaseMatrix::storageType.Parse( storageTypeStr );

    // get the matrix graph of the system matrix
    BaseGraph *graph = graphManager_->GetGraph();
    nne_  = graph->GetNNE();
    size_ = graph->GetSize();

    // Determine total number of equation numbers for this system
    // including those for fixed dofs
    totalSize_ = 0;
    for ( UInt i = 0; i < numPDEs_; i++ ) {
      totalSize_ += sizePerPDE_[i];
    }

    LOG_DBG(stdSys)  << "CreateLinSys() size=" << size_ << " nne=" << nne_ << " totalSize=" << totalSize_;
    LOG_DBG3(stdSys) << " Matrix Graph=" << graph->ToString();

    // ------------------------
    //  Re-set numLastFreeDof
    // ------------------------
    delete [] ( numLastFreeDof_ );  numLastFreeDof_  = NULL;
    NEWARRAY( numLastFreeDof_, UInt, 1 );
    numLastFreeDof_[0] = graph->GetSize();

    LOG_DBG(stdSys) << " CreateLinSys(): Re-set numLastFreeDof to " << numLastFreeDof_[0];

    // --------------------------------------------
    //  Treatment of Dirichlet Boundary Conditions
    // --------------------------------------------

    // Check the case we have (elimination vs. penalty) and generate
    // an appropriate object for handling inhomogeneous Dirichlet
    // boundary conditions
    if ( usingPenalty_ ) {
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

    //check for solver
    BaseSolver::SolverType solver;
    std::string solverStr;
    PtrParamNode solverNode = xml_->Get("solver", ParamNode::INSERT);
    solverNode->GetValue("type", solverStr, ParamNode::INSERT);
    solver = BaseSolver::solverType.Parse(solverStr);

    // In the case of SRCS_Matrices we can save memory by sharing
    // sparsity patterns
    bool sharePatternPossible = false;
    bool sharePattern = false;
    PatternIdType patternID = NO_PATTERN_ID;
    if ( storageType == BaseMatrix::SPARSE_SYM &&
         solver != BaseSolver::DIAGSOLVER ) {
      sharePatternPossible = true;
      patternPool_ = new PatternPool;
    }

    // actual matrix generation
    for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ )
    {
      BaseMatrix::StorageType actStorageType;

      LOG_DBG(stdSys) << "Generating matrix[" << *it << "]" << " MatrixEntryType=" << entryType
      << " MatrixStorageType= " << storageType << " dofs= " << blockSize_
      << " rows/cols= " << size_ << " nnz=" << nne_;

      actStorageType = storageType;
      if ( solver == BaseSolver::DIAGSOLVER ) {
        if ( *it == SYSTEM || *it == MASS ) {
          // just allocate a diagonal matrix
          actStorageType = BaseMatrix::DIAG;
        }
      }

      LOG_DBG(stdSys) << "calling GenerateStdMatrixObject" << std::endl;

      sysmat_[*it] = GenerateStdMatrixObject( entryType, actStorageType,
                                              size_, size_, nne_ );

      assert(sysmat_[*it] != NULL);

      string tmp;
      Enum2String( *it, tmp );

      LOG_DBG(stdSys) << "created matrix " << tmp;

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
    rhs_ = dynamic_cast<SingleVector *>
      (GenerateVectorObject( *(sysmat_[SYSTEM]) ));
    sol_ = dynamic_cast<SingleVector *>
      (GenerateVectorObject( *(sysmat_[SYSTEM]) ));

    // Generate communication buffers
    if ( !usingPenalty_ )
      GenerateCFSTransferBuffers();

    if ( rhs_ == NULL || sol_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    // Create the assemble object
    assemble_ =  new BaseEntryManipulator();
    //GenerateEntryManipulatorObject( entryType, blockSize_ );

    // Print information about matrix types
    PtrParamNode inm = info->Get("OLAS")->Get(ParamNode::HEADER)->Get("matrices");
    for ( it = matrixTypes_.begin(); it != matrixTypes_.end(); it++ )
    {
      string tmp;
      Enum2String( *it, tmp ); // such fucking non-sense :(
      PtrParamNode in = inm->Get(tmp);
      StdMatrix* mat = sysmat_[*it];
      in->Get("rows")->SetValue(mat->GetNumRows());
      in->Get("cols")->SetValue(mat->GetNumCols());
      in->Get("nnz")->SetValue(mat->GetNnz());
      in->Get("storageType")->SetValue(BaseMatrix::storageType.ToString(mat->GetStorageType()));
      in->Get("structueType")->SetValue(BaseMatrix::structureType.ToString(mat->GetStructureType()));
    }

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
    solver_ = GenerateSolverObject( *(sysmat_[SYSTEM]), xml_, olasInfo_ );
  }


  // ********************************
  //   Create preconditioner object
  // ********************************

  // Das Interface muessen wir uns noch ueberlegen
  void StandardSystem::CreatePrecond()
  {
    precond_ = GenerateStdPrecondObject( *(sysmat_[SYSTEM]), xml_, olasInfo_ );
  }

  // ********************************
  //   Create eigensolver object
  // ********************************

  void StandardSystem::CreateEigenSolver() 
  {
    eigenSolver_ = GenerateEigenSolverObject( *(sysmat_[SYSTEM]), xml_, olasInfo_->Get("eigenSolver"));
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
         usingPenalty_ ) {
      assembleDirichletToSysMat_ = true;
    }
  }


  // ***********************
  //   InitRHS (Version 1)
  // ***********************
  void StandardSystem::InitRHS( const PdeIdType identifierPDE ) {
    rhs_->Init();
  }


  // ***********************
  //   InitRHS (Version 2)
  // ***********************
  void StandardSystem::InitRHS( const BaseVector& newRHS ) {
    if( newRHS.GetEntryType() == BaseMatrix::COMPLEX ) {
      assemble_->InitRHS( rhs_, dynamic_cast<const Vector<Complex>&>(newRHS) );
    } else {
      assemble_->InitRHS( rhs_, dynamic_cast<const Vector<Double>&>(newRHS) );
    }
  }



  // **********************
  //   InitSol (Version1)
  // **********************
  void StandardSystem::InitSol( const PdeIdType identifierPDE ) {
    sol_->Init();
  }

  // ***********************
  //   InitSol (Version 2)
  //   for doubles
  // ***********************
  void StandardSystem::InitSol( const BaseVector& newSol ) {

//     try{
    // 	    std::cerr << "\n Initializing solution vector to initial conditions " << std::endl;

//	    	std::cout << "\n InitSol: " ;
//	    	for(UInt i=0; i<size; i++){
//	    		std::cout <<  newSol[i] << "," ;
//	    	}
    //	    	std::cout << std::endl;
    if( newSol.GetEntryType() == BaseMatrix::COMPLEX ) {
      assemble_->InitSol( sol_, dynamic_cast<const Vector<Complex>&>(newSol) );
    } else {
      assemble_->InitSol( sol_, dynamic_cast<const Vector<Double>&>(newSol) );
    }
// 	    else{
// #ifdef DEBUG_STANDARDSYSTEM
//           EXCEPTION(size <<" != sol_:"<< sol_->GetSize() << std::endl;
//                     << "size of initial solution vector doesn't match "
//                     << "with the size of the initialized solution vector!");
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
                                         const Matrix<Double>& elemMat,
                                         PdeIdType idPDE1,
                                         const StdVector<Integer>& eqnNrs1,
                                         PdeIdType idPDE2,
                                         const StdVector<Integer>& eqnNrs2,
                                         bool setCounterPart) {


    // then normal assembly:
    // element matrix and its transposed on the counter
    // part of the global matrix when setCounterPart is
    // true (upper and lower global assembly).
    assemble_->SetElementMatrix( matrixID, idPDE1, idPDE2,
                                 sysmat_[matrixID],
                                 sysmat_[matrixID],
                                 idbcHandler_, elemMat,
                                 eqnNrs1,
                                 eqnNrs2,
                                 numLastFreeDof_[0],
                                 numLastFreeDof_[0],
                                 setCounterPart );

  }
  
  void StandardSystem::SetElementMatrix( FEMatrixType matrixID, 
                                         const Matrix<Complex>& elemMat,
                                         PdeIdType idPDE1,
                                         const StdVector<Integer>& eqnNrs1,
                                         PdeIdType idPDE2,
                                         const StdVector<Integer>& eqnNrs2,
                                         bool setCounterPart) {
    // then normal assembly:
       // element matrix and its transposed on the counter
       // part of the global matrix when setCounterPart is
       // true (upper and lower global assembly).
       assemble_->SetElementMatrix( matrixID, idPDE1, idPDE2,
                                    sysmat_[matrixID],
                                    sysmat_[matrixID],
                                    idbcHandler_, elemMat,
                                    eqnNrs1,
                                    eqnNrs2,
                                    numLastFreeDof_[0],
                                    numLastFreeDof_[0],
                                    setCounterPart );
  }


  // *****************
  //   SetElementRHS
  // *****************
  void StandardSystem::SetElementRHS( const Vector<Double>& elemRHS, 
                                      const PdeIdType idPDE,
                                      StdVector<Integer>& eqnNrs ) {
    assemble_->SetElementRHS( rhs_, elemRHS, eqnNrs, size_ );
  }
  void StandardSystem::SetElementRHS( const Vector<Complex>& elemRHS, 
                                      const PdeIdType idPDE,
                                      StdVector<Integer>& eqnNrs ) {
    assemble_->SetElementRHS( rhs_, elemRHS, eqnNrs,  size_ );
   }


  // ***********************
  //   SetNodeRHS (Double)
  // ***********************
  void StandardSystem::SetNodeRHS( Double val, const PdeIdType identifierPDE,
                                   Integer eqnNum ) {


    // Perform consistency check
#ifdef DEBUG_STANDARDSYSTEM
    if ( eqnNum > size_ || eqnNum < 1 ) {
      EXCEPTION("StandardSystem::SetNodeRHS: Inconsistency detected:"
               << "\n size      = " << size_
               << "\n totalSize = " << totalSize_
               << "\n eqnNum    = " << eqnNum
               << "\n val       = " << val
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'");
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
      EXCEPTION("StandardSystem::SetNodeRHS: Inconsistency detected:"
               << "\n size      = " << size_
               << "\n totalSize = " << totalSize_
               << "\n eqnNum    = " << eqnNum
               << "\n val       = " << val
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'");
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
  void StandardSystem::UpdateRHS( FEMatrixType matrix_id, const BaseVector& fup ) {
    if ( numDirichletValues_ > 0 &&
         !usingPenalty_ ) {
      WARN("UpdateRHS may fail in combo with "
           << "idbcHandling = elimination");
    }
    if( fup.GetEntryType() == BaseMatrix::COMPLEX ) {
      assemble_->UpdateRHS( rhs_, sysmat_[matrix_id],
                            dynamic_cast<const Vector<Complex>&>(fup) );
    } else {
    const Vector<Double> & tmp = dynamic_cast<const Vector<Double>&>(fup);
    assemble_->UpdateRHS( rhs_, sysmat_[matrix_id], tmp );
    }

  }

  // *************************************
  //   Remove IDBC Information from matrix
  // *************************************

  void StandardSystem::RemoveIDBCInfoFromMatrix() const {

    std::set<FEMatrixType>::iterator it;
    BaseVector *rhs = NULL;
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
    UInt maxVal = usingPenalty_ == true ? size_ : size_ + numDirichletValues_;
    UInt minVal = usingPenalty_ == true ? 1 : size_ + 1;

    if ( eqnNum > maxVal || eqnNum < minVal ) {
      EXCEPTION("StandardSystem::SetDirichlet: Inconsistency detected:"
                << "\n pdeID  = " << pdeID
                << "\n eqnNum = " << eqnNum
                << "\n val    = " << val
                << "\n numBC  = " << numDirichletValues_
                << "\n size   = " << size_
                << "\n minVal = " << minVal
                << "\n maxVal = " << maxVal
                << "\n SystemName is '"
                << myParams_.GetStringValue( "SystemName" ) << "'");
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
    UInt maxVal = usingPenalty_ == true ? size_ : size_ + numDirichletValues_;
    UInt minVal = usingPenalty_ == true ? 1 : size_ + 1;

    if ( eqnNum > maxVal || eqnNum < minVal || ) {
      EXCEPTION("StandardSystem::SetDirichlet: Inconsistency detected:"
                << "\n pdeID  = " << pdeID
                << "\n eqnNum = " << eqnNum
                << "\n val    = " << val
                << "\n numBC  = " << numDirichletValues_
                << "\n size   = " << size_
                << "\n minVal = " << minVal
                << "\n maxVal = " << maxVal);
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
        WARN("StandardSystem::ConstructEffectiveMatrix: "
             << "Map with factors is empty, but there are "
             << matrixTypes_.size() << " FE matrices in the game!");
      }
    }

    // Intialize system matrix before adding the weighted
    // matrices to it
    this->InitMatrix( SYSTEM );

    for ( it = matFactors.begin(); it != matFactors.end(); it++ )
    {
#ifndef NDEBUG
      string tmp;
      Enum2String( (*it).first, tmp );
      LOG_DBG(stdSys) << " matFactor: " << tmp << ":" << (*it).second;
#endif
      if (sysmat_[(*it).first] != NULL  && (*it).second != 0.0 )
        sys->Add((*it).second, *sysmat_[(*it).first] );
    }
    LOG_DBG3(stdSys) << "effective Matrix: " << sys->ToString();

    // Also assemble the effective auxilliary system matrix for moving
    // IDBCs to the right-hand side
    idbcHandler_->BuiltSystemMatrix( matFactors );
  }


  // ******************************
  //   GetSolutionVal
  // ******************************
  void StandardSystem::
  GetSolutionVal( SingleVector& ptSol,
                  const PdeIdType identifierPDE ) {
    

    // Elimination case
    if ( !usingPenalty_ ) {

      if( ptSol.GetEntryType() == BaseMatrix::COMPLEX ) {
        // Sequential setting: We copy the solution and add the Dirichlet values
        const Vector<Complex> * dVec1 = dynamic_cast<Vector<Complex>*>( sol_ );
        Vector<Complex> * retVec = dynamic_cast<Vector<Complex>*>( &ptSol );
        (*retVec) = (*dVec1);
        idbcHandler_->SetDofsToIDBC( retVec );
      } else {
        // Sequential setting: We copy the solution and add the Dirichlet values
        const Vector<Double> & dVec1 = dynamic_cast<Vector<Double>&>( *sol_ );
        Vector<Double> & retVec = dynamic_cast<Vector<Double>& >( ptSol );
        retVec.Resize( totalSize_ );
        for ( int i = 0; i < size_; i++ ) {
          retVec[i] = dVec1[i];
        }
        idbcHandler_->SetDofsToIDBC( &retVec );
      }
    }

    // Penalty case
    else {
      if( ptSol.GetEntryType() == BaseMatrix::COMPLEX ) {
       Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptSol );
        retVec = dynamic_cast<Vector<Complex>&>(*sol_);
      } else {
        Vector<Double> & retVec = dynamic_cast<Vector<Double>&>( ptSol );
        retVec = dynamic_cast<Vector<Double>&>(*sol_);
      }
    }

  }


  // ************
  //   GetRHSVal 
  // ************
  void StandardSystem::GetRHSVal( SingleVector &ptRhs,
                                  const PdeIdType identifierPDE ) {
    if ( numDirichletValues_ > 0 && ! usingPenalty_ ) {
      WARN("GetRHSVal does not work for idbcHandling=elimination");
    }
    if( ptRhs.GetEntryType() == BaseMatrix::COMPLEX ) {
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>&>( ptRhs );
      retVec = dynamic_cast<Vector<Complex>&>(*rhs_);
    } else {
      Vector<Double> & retVec = dynamic_cast<Vector<Double>&>( ptRhs );
      retVec = dynamic_cast<Vector<Double>&>(*rhs_);
    }
  }

  // *********
  //   Print
  // *********
  void StandardSystem::Print( FEMatrixType matrix_id) const {
    *cla << sysmat_.find(matrix_id)->second->ToString();
  }

  // **********
  //   Export
  // **********
  void StandardSystem::Export( FEMatrixType type, char *filename,
                               char *comment ) const {
    sysmat_.find(type)->second->Export( filename, comment );
  }


  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  void StandardSystem::AddToDiagMatrixEntry( FEMatrixType matrixID,
                                             const PdeIdType pdeID,
                                             Integer eqnNum,
                                             Double val ) {

    WARN( "Adapt me");
    // Delegate work to implementation in assemble class
    //assemble_->AddToDiagMatrixEntry( sysmat_[matrixID], eqnNum, val );
  }
  
  void StandardSystem::AddToDiagMatrixEntry( FEMatrixType matrixID,
                                             const PdeIdType pdeID,
                                             Integer eqnNum,
                                             Complex val )  {
                                     
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
                               colEqnNum, numLastFreeDof_[0],
                               numLastFreeDof_[0] );

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
                               colEqnNum,  numLastFreeDof_[0],
                               numLastFreeDof_[0] );
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
                               colEqnNum,  numLastFreeDof_[0],
                               numLastFreeDof_[0], setCounterPart );
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
                               colEqnNum, numLastFreeDof_[0],
                               numLastFreeDof_[0], setCounterPart );
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
    BaseMatrix::StorageType sType = sysmat_[SYSTEM]->GetStorageType();
    BaseMatrix::EntryType eType = sysmat_[SYSTEM]->GetEntryType();

    // Generate vector for passing solution back to CFS++
    bVec = GenerateSingleVectorObject( sType, eType, totalSize );
    solComm_ = dynamic_cast<SingleVector*>( bVec );

  }

}
