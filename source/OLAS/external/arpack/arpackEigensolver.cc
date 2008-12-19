// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "arpackEigensolver.hh"
#include "arpackSolver.cc"

#include <limits>
#include "solver/solver.hh"
#include "precond/precond.hh"
#include "matvec/matvec.hh"

namespace OLAS {

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
    
    shiftAndInvert_ = FALSE;
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

  void ArpackEigenSolver::Setup(const  BaseMatrix & stiffMat,
                                const  BaseMatrix & massMat,
                                UInt numFreq, Double freqShift ) {

    // Set flag for indicating a non-quadratic problem
    isQuadratic_ = false;

    // Copy matrix references and determine size of system
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);
    if ( matrixA_ == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
    
    matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);
    if ( matrixB_ == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

    UInt size = matrixA_->GetNrows();
    
    // Save frequency parameters
    numFreq_ = numFreq;
    freqShift_ = freqShift;
    
    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Create matrix interface for arpack
    interface_ = New ArpackMatInterface( matrixA_, matrixB_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = New ArpackSolver();
    
    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = (myParams_->GetStringValue( "ARPACK_which" ));
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
        Setup( interface_, size, numFreq_, freqShift_, which_ , shiftAndInvert_);

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
    if ( myParams_->GetBoolValue( "ARPACK_logging" ) == TRUE ) {
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
    if ( matrixA_->GetStructureType() == STDMATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixB_ );
      precond_ = GenerateStdPrecondObject( mat, precond,
                                           myParams_, myReport_ );
    } else {
      (*error) << "No preconditioner available for SBM-matrices!";
      Error( __FILE__, __LINE__ );      
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

    // Temporary hard coded debug information
    std::cerr << "--- Setup of QUADRATIC EV problem! ---\n";

    // Copy matrix references and convert them to StdMatrices
    matrixA_ = & dynamic_cast<const StdMatrix&>(stiffMat);
    if ( matrixA_ == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
    
    matrixB_ = & dynamic_cast<const StdMatrix&>(massMat);
    if ( matrixB_ == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
    
    matrixD_ = & dynamic_cast<const StdMatrix&>(dampMat);    
    if ( matrixD_ == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

    // COMPWARNING unused variable UInt size = matrixA_->GetNrows();

    // Now generate a new complex matrix 
    MatrixEntryType eType = COMPLEX;
    MatrixStorageType sType = matrixA_->GetStorageType();
    Integer dof = matrixA_->GetBlockSize();
    Integer nrows = matrixA_->GetNrows();
    Integer ncols = matrixA_->GetNcols();
    Integer fill = matrixA_->GetNnz(); // not really needed for 
    StdMatrix * newComplexMat = GenerateStdMatrixObject( eType, sType, dof,
                                                         nrows, ncols, fill );
    // Copy sparsity pattern
    matrixA_->CopySparsityPattern( *newComplexMat );

    // Now "add" the double matrix (e.g. matrixA_) to the entries of 
    // newComplexMat
    newComplexMat->Add( 1.0, *matrixA_ );

    // Print original double matrix into .las file
    matrixA_->Print( *cla );

    // Print complex matrix into .las file
    newComplexMat->Print( *cla );
    

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
        solConverted[i+1] = arpackSolver_->Eigenvalue(i);
        if (solConverted[i+1] >= 0.0)
          solConverted[i+1] = sqrt(solConverted[i+1])/(8.0*atan(1.0));
      }
    } else {
      // case2: quadratic problem

      // ... TO BE IMPLEMENTED ...
    }

    // Save error norms
    Vector<Double> & errVec = dynamic_cast<Vector<Double>&>(err);
    errVec.Resize( numEVs );
    for (UInt i = 0; i < numEVs; i++ ) {
        errVec[i+1] = arpackSolver_->Tolerance(i);
    }

    // Return number of converged eigenvalues
    return numEVs;
  }
 
  void ArpackEigenSolver::CalcEigenMode( UInt modeNr, Vector<Double> & mode ) {
    
    UInt size = matrixA_->GetNrows();
    mode.Resize( size );
    mode.Init();

    for ( UInt i = 0; i < size; i++ ) {
        mode[i+1] = (arpackSolver_->GetEigenvector( modeNr ))[i];
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
    if ( shiftAndInvert_ == TRUE ) {
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
    if ( myParams_->GetBoolValue( "ARPACK_logging" ) == TRUE ) {
      (*cla) << " Logging activated\n";
    }
    (*cla) << " -------------------------------------------------------"
           << "-----------------------\n\n";
  }
  
}
