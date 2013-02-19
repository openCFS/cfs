#include <limits>
#include <string.h>

#include "MatVec/StdMatrix.hh"
#include "MatVec/generatematvec.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "OLAS/algsys/SolStrategy.hh"

#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"

#include "ArpackEigenSolver.hh"

namespace CoupledField {

  ArpackEigenSolver::ArpackEigenSolver( shared_ptr<SolStrategy> strat,
                                        PtrParamNode xml,
                                        PtrParamNode solverList,
                                        PtrParamNode precondList,
                                        PtrParamNode eigenInfo)
    : BaseEigenSolver( strat, xml, solverList, precondList, eigenInfo ){

    interface_    = NULL;
    matrixA_      = NULL;
    matrixB_      = NULL;
    matrixD_      = NULL;
    zStiff_       = NULL;
    zMass_        = NULL;
    zDamp_        = NULL;
    arpackSolver_ = NULL;
    solver_       = NULL;
    precond_      = NULL;
    which_        = NULL;
    xml_          = xml;
    
    isGeneralized_ = false;
    shiftAndInvert_ = false;
    logging_ = false;
  }

  ArpackEigenSolver::~ArpackEigenSolver() {

    delete solver_;
    delete precond_;
    delete interface_;
    delete arpackSolver_;

    delete zStiff_;
    delete zDamp_;
    delete zMass_;

    delete[] which_;

    matrixA_ = NULL;
    matrixB_ = NULL;
    matrixD_ = NULL;
    zStiff_ = NULL;
    zMass_  = NULL;
    zDamp_  = NULL;
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
    std::string whichString = "LM";
    xml_->GetValue("which", whichString, ParamNode::INSERT);
    
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
    Setup( interface_, size, numFreq_, freqShift_, which_ , (char*) "I", shiftAndInvert_);

    // Set additional parameters for tolerance, number of Arnoldi vectors and
    // number of iterations
    Double tol = -1.0;
    xml_->GetValue("tolerance", tol, ParamNode::INSERT );
    Integer maxIt = -1;
    
    xml_->GetValue("maxIt", maxIt, ParamNode::INSERT );
    Integer numVec = -1;
    
    xml_->GetValue("numVec", numVec, ParamNode::INSERT );

    if (tol > 0.0)
      arpackSolver_->SetTolerance(tol);
    if (maxIt > 0)
      arpackSolver_->SetIterations(maxIt);
    if (numVec > 0)
      arpackSolver_->SetNumVectors(numVec);

    // Check trace settings
    arpackSolver_->DebugOff();
    if( xml_->Has("logging") ) {
      xml_->GetValue("logging", logging_);
      if (logging_) {
        arpackSolver_->DebugOn();
        logging_ = true;
      }
    }
    
    // Print log-info about EigenSolver
    PrintInfo();

    // Create solver
    solver_ = GenerateSolverObject( *matrixA_, solStrat_,
                                    solverList_, eigenInfo_ );

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixA_ );

      // Create preconditioner
      std::string precondId = solStrat_->GetPrecondId();
      precond_ = GenerateStdPrecondObject( mat, precondId, 
                                           precondList_, eigenInfo_ );
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
    std::string whichString = "LM";
    xml_->GetValue("which", whichString, ParamNode::INSERT );
    
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
        Setup( interface_, size, numFreq_, freqShift_, which_ , (char*) "G", shiftAndInvert_);

    // Set additional parameters for tolerance, number of Arnoldi vectors and
    // number of iterations
    Double tol = -1.0;
    xml_->GetValue("tolerance", tol, ParamNode::INSERT);
    
    Integer maxIt = -1;
    xml_->GetValue("maxIt", maxIt, ParamNode::INSERT);
    
    Integer numVec = -1;
    xml_->GetValue("numVec", numVec, ParamNode::INSERT);

    if (tol > 0.0)
        arpackSolver_->SetTolerance(tol);
    if (maxIt > 0)
        arpackSolver_->SetIterations(maxIt);
    if (numVec > 0)
        arpackSolver_->SetNumVectors(numVec);

    // Check trace settings
    arpackSolver_->DebugOff();
    if( xml_->Has("logging") ) {
      xml_->GetValue("logging", logging_);
      if (logging_) {
        arpackSolver_->DebugOn();
        logging_ = true;
      }
    }
    

    // Print log-info about EigenSolver
    PrintInfo();

    // Create solver
    solver_ = GenerateSolverObject( *matrixB_, solStrat_, 
                                    solverList_, eigenInfo_ );

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixB_ );

      // Create preconditioner
      std::string precondId = solStrat_->GetPrecondId();
      precond_ = GenerateStdPrecondObject( mat, precondId, 
                                           precondList_, eigenInfo_ );
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

    // need twice the original system size and twice the number of frequencies
    UInt size = 2*matrixA_->GetNumRows();
    //numFreq_  = 2*numFreq;
    numFreq_  = numFreq;

    // currently no shifts are allowed
    freqShift_ = freqShift;
    
    // NOTE: Hard coded as true!!!
    shiftAndInvert_ = true;

    // Now generate a new complex matrix
    BaseMatrix::EntryType eType = BaseMatrix::COMPLEX;
    BaseMatrix::StorageType sType = matrixA_->GetStorageType();
    Integer nrows = matrixA_->GetNumRows();
    Integer ncols = matrixA_->GetNumCols();
    Integer fill = matrixA_->GetNnz(); // not really needed for

    zStiff_ = GenerateStdMatrixObject( eType, sType,
                                       nrows, ncols, fill );
    zStiff_->SetSparsityPattern( *matrixA_ );
    zStiff_->Add( 1.0, *matrixA_ );

    sType = matrixB_->GetStorageType();
    nrows = matrixB_->GetNumRows();
    ncols = matrixB_->GetNumCols();
    fill = matrixB_->GetNnz();
    zMass_ = GenerateStdMatrixObject( eType, sType,
                                      nrows, ncols, fill );
    zMass_->SetSparsityPattern( *matrixB_ );
    zMass_->Add( 1.0, *matrixB_ );

    sType = matrixD_->GetStorageType();
    nrows = matrixD_->GetNumRows();
    ncols = matrixD_->GetNumCols();
    fill = matrixD_->GetNnz();
    zDamp_ = GenerateStdMatrixObject( eType, sType,
                                      nrows, ncols, fill );
    zDamp_->SetSparsityPattern( *matrixD_ );
    zDamp_->Add( 1.0, *matrixD_ );

    // Create matrix interface for arpack
    interface_ = new ArpackMatInterface( zStiff_, zMass_, zDamp_, shiftAndInvert_, freqShift_ );

    // Create solver class
    arpackSolver_ = new ArpackSolver();
    
    // Check 'which'-settings regarding the type of eigenvalues searched for
    std::string whichString = "LM";
    xml_->GetValue("which", whichString, ParamNode::INSERT );
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
        QuadSetup( interface_, size, numFreq_, freqShift_, which_ , shiftAndInvert_);

    // Set additional parameters for tolerance, number of Arnoldi vectors and
    // number of iterations
    Double tol = -1.0;
    xml_->GetValue("tolerance", tol, ParamNode::INSERT);
    
    Integer maxIt = -1;
    xml_->GetValue("maxIt", maxIt, ParamNode::INSERT);
    
    Integer numVec = -1;
    xml_->GetValue("numVec", numVec, ParamNode::INSERT);

    if (tol > 0.0)
        arpackSolver_->SetTolerance(tol);
    if (maxIt > 0)
        arpackSolver_->SetIterations(maxIt);
    if (numVec > 0)
        arpackSolver_->SetNumVectors(numVec);

    // Set additional scaling parameter for matrix interface
    Double scaleFac = 1.0;
    xml_->GetValue("scaleFac", scaleFac, ParamNode::INSERT);
    if (scaleFac > 0.0)
      interface_->SetDiagScaling(scaleFac);

    // Check trace settings
    arpackSolver_->DebugOff();
    if( xml_->Has("logging") ) {
      xml_->GetValue("logging", logging_);
      if (logging_) {
        arpackSolver_->DebugOn();
        logging_ = true;
      }
    }

    // Print log-info about EigenSolver
    PrintInfo();

    // also report to terminal
    std::cout << "++ Running quadratic eigenvalue problem with the following settings\n" 
	      << "       number of frequencies       :  " << arpackSolver_->GetNev() << "\n"
	      << "       type of eigenvalues (which) :  " << arpackSolver_->GetWhich() << "\n"
	      << "       applied shift               :  " << arpackSolver_->GetShift() << "\n"
	      << "       convergence tolerance       :  " << arpackSolver_->GetTol() << "\n"
	      << "       maximum number of iterations:  " << arpackSolver_->GetMaxit() << "\n"
	      << "       number of Arnoldi vectors   :  " << arpackSolver_->GetNcv() << "\n";

    if (scaleFac != 1.0)
      std::cout << "       applied scaling factor      :  " << interface_->GetDiagScaling() << "\n";

    std::cout << std::endl;

    // Create solver
    solver_ = GenerateSolverObject( *zStiff_, solStrat_, 
                                    solverList_, eigenInfo_ );


    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixB_ );

      // Create preconditioner
      std::string precondId = solStrat_->GetPrecondId();

      precond_ = GenerateStdPrecondObject( mat, precondId, 
                                           precondList_, eigenInfo_ );
    } else {
      EXCEPTION( "No preconditioner available for SBM-matrices!" );
    }

    // Setup matrixinterface
    interface_->QuadSetup( solver_, precond_ );
  }

  UInt ArpackEigenSolver::CalcEigenFrequencies( BaseVector &sol,
                                                BaseVector &err) {


    UInt numEVs = 0;

    // Save eigenvalues
    // if non-negative eigenvalue, convert to eigenfrequency

    if( !isQuadratic_ ) {
      // Find the eigenvalues and calculate the eigenvectors
      numEVs = arpackSolver_->FindEigenvalues();

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
      numEVs = arpackSolver_->FindQuadEigenvalues();

      Vector<Complex> & solConverted = dynamic_cast<Vector<Complex>&>(sol);
      solConverted.Resize( numEVs );
      for (UInt i = 0; i < numEVs; i++ ) {
        solConverted[i] = arpackSolver_->CmplxEigenvalue(i);
      }
      
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
    Double tol = -1.0;
    xml_->GetValue("tolerance", tol, ParamNode::INSERT);

    Integer maxIt = -1;
    xml_->GetValue("maxIt", maxIt, ParamNode::INSERT);

    Integer numVec = -1;
    xml_->GetValue("numVec", numVec, ParamNode::INSERT);
    
    // set mode: we look at both ends of the spectrum for eigenvalues
    std::string whichString = "BE";
    which_ = new char[whichString.size()+1];
    strncpy(which_, whichString.c_str(), whichString.size()+1 );

    // Save frequency parameters
    numFreq_ = 10; 
    freqShift_ = 0.0;

    // Create the solver dependent on the problem type ( regular / shiftAndInvert )
    arpackSolver_->
    Setup( interface_, size, numFreq_, freqShift_, which_ , (char*) "I", shiftAndInvert_);


    if (tol > 0.0)
      arpackSolver_->SetTolerance(tol);
    if (maxIt > 0)
      arpackSolver_->SetIterations(maxIt);
    if (numVec > 0)
      arpackSolver_->SetNumVectors(numVec);


    // Check trace settings
    arpackSolver_->DebugOff();
    if( xml_->Has("logging") ) {
      xml_->GetValue("logging", logging_);
      if (logging_) {
        arpackSolver_->DebugOn();
        logging_ = true;
      }
    }
    
    // Print log-info about EigenSolver
    PrintInfo();
    
    // Create standard solver
    solver_ = GenerateSolverObject( *matrixA_, solStrat_, 
                                    solverList_, eigenInfo_ );

    // Perform check, if matrix is std or sbm
    if ( matrixA_->GetStructureType() == BaseMatrix::SPARSE_MATRIX ) {
      const StdMatrix & mat = dynamic_cast< const StdMatrix &>( *matrixA_ );

      // Create preconditioner
      std::string precondId = solStrat_->GetPrecondId();
      precond_ = GenerateStdPrecondObject( mat, precondId, 
                                           precondList_, eigenInfo_ );
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

  void ArpackEigenSolver::CalcEigenMode( UInt modeNr, Vector<Complex> & mode ) {

    UInt size = matrixA_->GetNumRows();
    mode.Resize( size );
    mode.Init();

    for ( UInt i = 0; i < size; i++ ) {
        mode[i] = Complex((arpackSolver_->GetEigenvector( modeNr ))[i],0);
    }
  }

  void ArpackEigenSolver::CalcQuadEigenMode( UInt modeNr, Vector<Complex> & mode ) {
    
    UInt size = matrixA_->GetNumRows();
    mode.Resize( size );
    mode.Init();

//    Double  rad2deg = 45.0/atan(1.0);
    Complex *zMode = new Complex [2*size];

//    (*cla) << "Mode no " << modeNr << std::endl;
    for ( UInt i = 0; i < 2*size; i++ ) {
      zMode[i] = (arpackSolver_->GetQuadEigenvector( modeNr ))[i];
//      (*cla) << modeNr << ", " << i << ":  " << std::abs(zMode[i]);
//      (*cla) << ", " << std::atan2(zMode[i].imag(),zMode[i].real())*rad2deg << std::endl;
    }
    
    for ( UInt i = 0; i < size; i++ ) {
      mode[i] = zMode[i+size];
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
    if ( logging_ ) {
      (*cla) << " Logging activated\n";
    }
    (*cla) << " -------------------------------------------------------"
           << "-----------------------\n\n";
  }

}
