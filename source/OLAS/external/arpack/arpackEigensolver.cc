// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <limits>
#include <string.h>

#include "MatVec/stdmatrix.hh"
#include "MatVec/generatematvec.hh"
#include "OLAS/algsys/olasparams.hh"

#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/baseprecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/basesolver.hh"

#include "arpackEigensolver.hh"

namespace CoupledField {

  ArpackEigenSolver::ArpackEigenSolver( ParamNode * xml,
                                        OLAS_Params *myParams,
                                        OLAS_Report *myReport )
    : BaseEigenSolver( myParams, myReport ){

    interface_    = NULL;
    matrixA_      = NULL;
    matrixB_      = NULL;
    matrixD_      = NULL;
    arpackSolver_ = NULL;
    solver_       = NULL;
    precond_      = NULL;
    which_        = NULL;
    xml_          = xml;
    
    isGeneralized_ = false;
    shiftAndInvert_ = false;
  }

  ArpackEigenSolver::~ArpackEigenSolver() {

    delete solver_;
    delete precond_;
    delete interface_;
    delete arpackSolver_;
    delete[] which_;

    matrixA_ = NULL;
    matrixB_ = NULL;
    matrixD_ = NULL;
  }
  void ArpackEigenSolver::Setup(const  BaseMatrix & mat,
                                UInt numFreq, Double freqShift ) {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = false;
    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(mat);
    if ( matrixA_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    UInt size = matrixA_->GetNumRows();

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

   

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver();

    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = (myParams_->GetStringValue( "ARPACK_which" ));
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
    Setup( interface_, size, numFreq_, freqShift_, which_ , "I", shiftAndInvert_);

    // Set additional parameters for tolerance, number of Arnoldi vectors and
    // number of iterations
    Double tol = myParams_->GetDoubleValue("ARPACK_tolerance");
    Integer maxIt = myParams_->GetIntValue("ARPACK_maxIt");
    Integer numVec = myParams_->GetIntValue("ARPACK_numVec");

    if (tol > 0.0)
      arpackSolver_->SetTolerance(tol);
    if (maxIt > 0)
      arpackSolver_->SetIterations(maxIt);
    if (numVec > 0)
      arpackSolver_->SetNumVectors(numVec);

    // Check trace settings
    if ( myParams_->GetBoolValue( "ARPACK_logging" ) == true ) {
      arpackSolver_->DebugOn();
    }

    // Print log-info about EigenSolver
    PrintInfo();

    // Create solver
    SolverType solver;
    myParams_->GetEnumValue( "Solver", solver );
    // killme! check the ParamNode parameter
    solver_ = GenerateSolverObject( *matrixA_, solver, xml_, myParams_,
                                    myReport_ );

    // Create preconditioner
    PrecondType precond;
    myParams_->GetEnumValue("Precond", precond);

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixA_ );
      precond_ = GenerateStdPrecondObject( mat, precond,
                                           myParams_, myReport_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->Setup( solver_, precond_ );

  }
  
  
  void ArpackEigenSolver::Setup(const  BaseMatrix & stiffMat,
                                const  BaseMatrix & massMat,
                                UInt numFreq, Double freqShift ) {

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = true;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);
    if ( matrixA_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);
    if ( matrixB_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    UInt size = matrixA_->GetNumRows();

    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;

    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, matrixB_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver();

    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = (myParams_->GetStringValue( "ARPACK_which" ));
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
        Setup( interface_, size, numFreq_, freqShift_, which_ , "G", shiftAndInvert_);

    // Set additional parameters for tolerance, number of Arnoldi vectors and
    // number of iterations
    Double tol = myParams_->GetDoubleValue("ARPACK_tolerance");
    Integer maxIt = myParams_->GetIntValue("ARPACK_maxIt");
    Integer numVec = myParams_->GetIntValue("ARPACK_numVec");

    if (tol > 0.0)
        arpackSolver_->SetTolerance(tol);
    if (maxIt > 0)
        arpackSolver_->SetIterations(maxIt);
    if (numVec > 0)
        arpackSolver_->SetNumVectors(numVec);

    // Check trace settings
    if ( myParams_->GetBoolValue( "ARPACK_logging" ) == true ) {
       arpackSolver_->DebugOn();
    }

    // Print log-info about EigenSolver
    PrintInfo();

    // Create solver
    SolverType solver;
    myParams_->GetEnumValue( "Solver", solver );
    // killme! check the ParamNode parameter
    solver_ = GenerateSolverObject( *matrixB_, solver, xml_, myParams_,
                                    myReport_ );

    // Create preconditioner
    PrecondType precond;
    myParams_->GetEnumValue("Precond", precond);

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixB_ );
      precond_ = GenerateStdPrecondObject( mat, precond,
                                           myParams_, myReport_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->Setup( solver_, precond_ );

  }

  void ArpackEigenSolver::Setup(const  BaseMatrix & stiffMat,
                                const  BaseMatrix & massMat,
                                const  BaseMatrix & dampMat,
                                UInt numFreq, Double freqShift ) {

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = true;
    isGeneralized_ = true;

    // Temporary hard coded debug information
    // std::cerr << "--- Setup of QUADRATIC EV problem! ---\n";

    // Copy matrix references and convert them to StdMatrices
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);
    if ( matrixA_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);
    if ( matrixB_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    matrixD_ = & dynamic_cast<const StdMatrix&>(dampMat);
    if ( matrixD_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    // UInt size = matrixA_->GetNumRows(); // TODO: Check if this is still needed

    // Now generate a new complex matrix
    BaseMatrix::EntryType eType = BaseMatrix::COMPLEX;
    BaseMatrix::StorageType sType = matrixA_->GetStorageType();
    Integer nrows = matrixA_->GetNumRows();
    Integer ncols = matrixA_->GetNumCols();
    Integer fill = matrixA_->GetNnz(); // not really needed for
    StdMatrix * newComplexMat = GenerateStdMatrixObject( eType, sType,
                                                         nrows, ncols, fill );
    // Copy sparsity pattern
    matrixA_->CopySparsityPattern( *newComplexMat );

    // Now "add" the double matrix (e.g. matrixA_) to the entries of
    // newComplexMat
    newComplexMat->Add( 1.0, *matrixA_ );

//    // Print original double matrix into .las file
//    *cla << matrixA_->ToString();
//
//    // Print complex matrix into .las file
//    *cla << newComplexMat->ToString();


    exit(-1);

    // ... TO BE IMPLEMENTED ...

 }

  UInt ArpackEigenSolver::CalcEigenFrequencies( BaseVector &sol,
                                                BaseVector &err) {

    // Find the eigenvalues and calculate the eigenvectors
    UInt numEVs = arpackSolver_->FindEigenvalues();



    // Save eigenvalues
    // if non-negative eigenvalue, convert to eigenfrequency

    if( !isQuadratic_ ) {
      // case1: generalized problem
      Vector<Double> & solConverted = dynamic_cast<Vector<Double>&>(sol);
      solConverted.Resize( numEVs );
      for (UInt i = 0; i < numEVs; i++ ) {
        solConverted[i] = arpackSolver_->Eigenvalue(i);
        if (solConverted[i] >= 0.0)
          solConverted[i] = sqrt(solConverted[i])/(8.0*atan(1.0));
      }
    } else {
      // case2: quadratic problem

      // ... TO BE IMPLEMENTED ...
      EXCEPTION( "Not implemented yet" );
      
    }

    // Save error norms
    Vector<Double> & errVec = dynamic_cast<Vector<Double>&>(err);
    errVec.Resize( numEVs );
    for (UInt i = 0; i < numEVs; i++ ) {
        errVec[i] = arpackSolver_->Tolerance(i);
    }

    // Return number of converged eigenvalues
    return numEVs;
  }

  void ArpackEigenSolver::
  CalcConditionNumber( const BaseMatrix& mat, 
                       Double& condNumber, Vector<Double>& evs,
                       Vector<Double>& err ) {
    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;
    isGeneralized_ = false;
    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(mat);
    if ( matrixA_ == NULL ) {
      EXCEPTION( WRONG_CAST_MSG );
    }

    UInt size = matrixA_->GetNumRows();

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( matrixA_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver();

    // Set additional parameters for tolerance, number of Arnoldi vectors and
    // number of iterations
    Double tol = myParams_->GetDoubleValue("ARPACK_tolerance");
    Integer maxIt = myParams_->GetIntValue("ARPACK_maxIt");
    Integer numVec = myParams_->GetIntValue("ARPACK_numVec");

    // set mode: we look at both ends of the spectrum for eigenvalues
    std::string whichString = "BE";
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Save frequency parameters
    numFreq_ = 10; 
    freqShift_ = 0.0;

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
    Setup( interface_, size, numFreq_, freqShift_, which_ , "I", shiftAndInvert_);


    if (tol > 0.0)
      arpackSolver_->SetTolerance(tol);
    if (maxIt > 0)
      arpackSolver_->SetIterations(maxIt);
    if (numVec > 0)
      arpackSolver_->SetNumVectors(numVec);

    // Check trace settings
    if ( myParams_->GetBoolValue( "ARPACK_logging" ) == true ) {
      arpackSolver_->DebugOn();
    }

    // Print log-info about EigenSolver
    PrintInfo();

    // Create standard solver
    SolverType solver;
    myParams_->GetEnumValue( "Solver", solver );
    solver_ = GenerateSolverObject( *matrixA_, solver, xml_, myParams_,
                                    myReport_ );

    // Create preconditioner
    PrecondType precond;
    myParams_->GetEnumValue("Precond", precond);

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixA_ );
      precond_ = GenerateStdPrecondObject( mat, precond,
                                           myParams_, myReport_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->Setup( solver_, precond_ );

    //Find the eigenvalues and calculate the eigenvectors
    UInt numEVs = arpackSolver_->FindEigenvalues();
    if( numEVs == 0) {
      EXCEPTION( "No convergence in search for eigenvalues");
    }

    evs.Resize( numEVs );
    err.Resize( numEVs );
    for (UInt i = 0; i < numEVs; i++ ) {
      evs[i] = arpackSolver_->Eigenvalue(i);
      err[i] = arpackSolver_->Tolerance(i);
    }

    condNumber = evs[evs.GetSize()-1] / evs[0];
  }

  void ArpackEigenSolver::CalcEigenMode( UInt modeNr, Vector<Double> & mode ) {

    UInt size = matrixA_->GetNumRows();
    mode.Resize( size );
    mode.Init();

    for ( UInt i = 0; i < size; i++ ) {
        mode[i] = (arpackSolver_->GetEigenvector( modeNr ))[i];
    }
  }

  void ArpackEigenSolver::PrintInfo() {

    (*cla) << " -------------------------------------------------------"
           << "-----------------------\n"
           << " ARPACK Eigenvalue Solver - Information\n\n";

    // Type of problem to be solved
    (*cla) << "Problem type: ";
    if ( isQuadratic_ ) {
      (*cla) << "QUADRATIC";
    } else {
      (*cla) << "GENERALIZED";
    }
    (*cla) << "\n";

    // Number of eigenvalues
    (*cla) << " Calculating " << arpackSolver_->GetNev()
           << " eigenvalues\n";

    // 'which'-settings
    (*cla) << " Searching for eigenvalues with '"
           << arpackSolver_->GetWhich() << "'\n";

    // Shift-And-Invert info
    if ( shiftAndInvert_ == true ) {
      (*cla) << " Applying Shift-And-Invert mode with shift of "
             << arpackSolver_->GetShift() << "\n";
    } else {
      (*cla) << " Applying regular mode\n";
    }

    // Tolerance, number of arnoldi vectors, number of Iterations
    (*cla) << " Using tolerance of "
            << arpackSolver_->GetTol() << "\n";
    (*cla) << " Maximum number of iterations: "
           << arpackSolver_->GetMaxit() << "\n";
    (*cla) << " Maximum number of Arnoldi vectors: "
           << arpackSolver_->GetNcv() << "\n";

    // Trace settings
    if ( myParams_->GetBoolValue( "ARPACK_logging" ) == true ) {
      (*cla) << " Logging activated\n";
    }
    (*cla) << " -------------------------------------------------------"
           << "-----------------------\n\n";
  }

}
