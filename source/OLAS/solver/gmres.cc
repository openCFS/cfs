// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "algsys/algsys.hh"

#include "solver/gmres.hh"
#include "General/exception.hh"
using CoupledField::Exception;

namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  GMRESSolver<T>::GMRESSolver( OLAS_Params *myParams, OLAS_Report *myReport ){

    ENTER_FCN( "GMRESSolver::GMRESSolver" );

    // Set pointers to communication objects
    myParams_ = myParams;
    myReport_ = myReport;

    // Prepare remaining internal attributes
    c_              = NULL;
    s_              = NULL;
    hMat_           = NULL;
    pVec_           = NULL;
    maxKrylovDim_   = 0;
    problemDim_     = 0;
    hMatNumRows_    = 0;

    // Generate object for performing Givens rotations
    givens_ = new GivensRotation( GivensRotation::OLAS );
  }


  // **********************
  //   Default destructor
  // **********************
  template<typename T>
  GMRESSolver<T>::~GMRESSolver() {

    ENTER_FCN( "GMRESSolver::~GMRESSolver" );

    // ---------------
    //   Free memory
    // ---------------

    // Givens coefficients
    DeleteArray( c_ );
    DeleteArray( s_ );

    // Upper Hessenberg matrix
    DeleteHessenbergMatrix();

    // Krylov base vectors
    DeleteKrylovBasis();

    // Givens rotation object
    delete givens_;

    // Auxilliary BaseVector
    delete pVec_;

  }


  // **************************
  //   Setup (public version)
  // **************************
  template<typename T>
  void GMRESSolver<T>::Setup( BaseMatrix &sysMat ) {
    ENTER_FCN( "GMRESSolver::Setup" );
    PrivateSetup( sysMat );
  }


  // ***************************
  //   Setup (private version)
  // ***************************
  template<typename T>
  void GMRESSolver<T>::PrivateSetup( const BaseMatrix &sysMat ) {

    ENTER_FCN( "GMRESSolver::PrivateSetup" );

    // ----------------------------------------
    //   Preform some consistency checks
    // ----------------------------------------

    // Obtain relevant dimension information
    UInt nCols = sysMat.GetNcolsScalar();
    UInt nRows = sysMat.GetNrowsScalar();
    UInt newKrylovDim = (UInt)myParams_->GetIntValue( "GMRES_maxKrylovDim" );

    // Check that matrix is square
    if ( nCols != nRows ) 
      EXCEPTION("system matrix is not square but " << nRows << " x " << nCols);

    // Check that the maximal dimension of the Krylov subspace is not
    // larger than the dimension of the problem
    if ( newKrylovDim > nCols ) {
      (*warning) << "maxKrylovDim = " << newKrylovDim << " exceeds matrix "
                 << " dimension = " << nCols << ". Using " << nCols
                 << " instead.";
      Warning( __FILE__, __LINE__ );
      newKrylovDim = nCols;
      myParams_->SetValue( "GMRES_maxKrylovDim", (Integer)nCols );
    }

    // -------------------------------------
    //   (Re-)allocate Givens coefficients
    // -------------------------------------
    if ( newKrylovDim > maxKrylovDim_ ) {
      DeleteArray( c_ );
      DeleteArray( s_ );
      NewArray( c_, Double, newKrylovDim );
      NewArray( s_, T, newKrylovDim );
    }

    // -------------------------------------
    //   (Re-)allocate right-hand side for
    //   the least squares problem and
    //   also auxilliary BaseVector
    // -------------------------------------
    if ( nCols != problemDim_ ) {
      bVec_.Resize( nCols );
      delete pVec_;
      pVec_ = GenerateVectorObject( sysMat );
    }

    // ----------------------------------------
    //   Must set the sizes to the 'new'
    //   values before, calling the next
    //   two methods
    // ----------------------------------------
    maxKrylovDim_ = newKrylovDim;
    problemDim_ = nCols;

    // ----------------------------------------
    //   (Re-)allocate upper HessenbergMatrix
    // ----------------------------------------
    AllocateHessenbergMatrix();

    // ----------------------------------------
    //   (Re-)allocate Krylov base vectors
    // ----------------------------------------
    AllocateKrylovBasis( sysMat );

  }


  // *********
  //   Solve
  // *********
  template<typename T>
  void GMRESSolver<T>::Solve( const BaseMatrix &sysMat,
                              const BasePrecond &precond,
                              const BaseVector &rhs, BaseVector &sol ) {

    ENTER_FCN( "GMRESSolver::Solve" );

    bool logging = myParams_->GetBoolValue( "GMRES_logging" );


    // ------------------------------------------------------
    //   We call PrivateSetup to (re)-allocate the internal
    //   data structures, if this is necessary
    // ------------------------------------------------------
    PrivateSetup( sysMat );


    // ----------------------
    //   A simple safeguard
    // ----------------------
    if ( hMat_ == NULL ) {
      Error( "GMRESSolver::Solve() Detected NULL pointer!", __FILE__,
             __LINE__ );
    }


    // ---------------------------------------
    //   Determine residual of initial guess
    // ---------------------------------------
    sysMat.CompRes( *(vMat_[1]), sol, rhs );


    // -------------------------------------------------------------
    //   Prepare threshold for stopping criterion. We use a
    //   relative test. An approximation x^k is accepted as good
    //   enough, if 
    //
    //                  || r^k || <= tau
    //
    //   where tau is computed from a user specified threshold tol.
    //   For absNorm tau = tol, for relNormRHS tau = tol * || b ||_2
    //   and for relNormRes0, if b = 0, or the penalty formulation
    //   is used to handle inhomogeneous Dirichlet boundary
    //   conditions, we replace ||b||_2 by || r^0 ||_2
    //
    //   We let ComputeThreshold do this, and also obtain norm of
    //   initial residual in this way
    // -------------------------------------------------------------
    Double resNorm = 0;
    Double eps = myParams_->GetDoubleValue( "GMRES_epsilon" );
    Double tolerance = ComputeThreshold( eps, rhs, *(vMat_[1]), resNorm,
                                         logging );
    if ( logging == true ) {
      LogConvergence( resNorm, 0, true );
    }


    // ----------------------------------
    //   Perform the GMRES(m) iteration
    // ----------------------------------
    Integer maxIter = myParams_->GetIntValue( "GMRES_maxIter" );
    Integer loopsDone = maxIter;
    UInt stepCount = 0;
    UInt stepCountGlobal = 0;
    bool approxIsGood = false;

    for ( Integer k = 1; k <= maxIter; k++ ) {

      // Perform one inner loop
      InnerLoop( sysMat, precond, rhs, sol, resNorm, tolerance,
                 approxIsGood, stepCount, stepCountGlobal );

      // Increase global iteration counter
      stepCountGlobal += stepCount;

      // Compute residual of new approximation and its norm
      // (Will serve in constructing base vector for next inner loop)
      sysMat.CompRes( *(vMat_[1]), sol, rhs );
      resNorm = vMat_[1]->NormEuclid();

      // If InnerLoop says approximation is fine, test for false
      // convergence. If everything is fine, stop iteration.
      if ( approxIsGood == true ) {
        if ( resNorm <= tolerance ) {
          loopsDone = k;
          if ( logging == true ) {
            (*cla) << "GMRES: Approximate satisfies stopping rule"
                   << std::endl;
          }
          break;
        }
        else {
          (*cla) << "# GMRES: False convergence!\n"
		 << "# Real resnorm = " << resNorm << std::endl;
          approxIsGood = false;
        }
      }
    }

    // Make clear, if GMRES failed to converge
    if ( !approxIsGood ) {
      (*cla) << "\n GMRES: Solution fails to satisfy stopping test!"
	     << std::endl;
    }


    // ----------------------------
    //   Generate solution report
    // ----------------------------

    // Number of iterations: Depends on GMRES(m) -> Full GMRES
    if ( maxIter == 1 ) {
      myReport_->SetValue( "numIter", (Integer)stepCount );
      myReport_->SetValue( "numGlobalIter", (Integer)stepCount );
    }
    else {
      myReport_->SetValue( "numIter", loopsDone );
      myReport_->SetValue( "numGlobalIter", (Integer)stepCountGlobal );
    }

    // Final relative residual norm
    myReport_->SetValue( "finalNorm", resNorm );

    if ( logging == true ) {
      (*cla) << "\n --> Final norm = " << resNorm << std::endl;
    }

    // Status of solution
    myReport_->SetValue( "solutionIsOkay", approxIsGood );

    // Report to standard logfile
    (*cla) << "\n GMRESSolver done\n" << std::endl;

  }


  // *************
  //   InnerLoop
  // *************
  template<typename T>
  void GMRESSolver<T>::InnerLoop( const BaseMatrix &sysMat,
                                  const BasePrecond &precond,
                                  const BaseVector &rhs,
                                  BaseVector &sol,
                                  Double beta,
                                  Double threshold,
                                  bool &approxIsGood,
                                  UInt &stepCount,
                                  UInt globNum ) {

    ENTER_IFCN( "GMRESSolver::InnerLoop" );

    UInt i, k;
    T aux = 0;
    Double resNormOld = 0.0;
    Double resNormNew = beta;
    approxIsGood = false;
    bool logging = myParams_->GetBoolValue( "GMRES_logging" );


    // ------------------
    //   Initialisation
    // ------------------
    
    // Setup first basis vector
    vMat_[1]->ScalarDiv( beta );

    // Setup right-hand side of least squares problem
    bVec_.Init();
    bVec_[1] = beta;


    // ------------------
    //   The inner loop
    // ------------------
    stepCount = maxKrylovDim_ - 1;
    for ( i = 1; i < maxKrylovDim_; i++ ) {


      // ------------------------------------
      //   Computation of next basis vector
      // ------------------------------------

      // Apply preconditioner to previous basis vector
      precond.Apply( sysMat, *(vMat_[i]), *pVec_ );

      // Candidate for next basis vector
      sysMat.Mult( *pVec_, *(vMat_[i+1]) );

      // Now orthonormalise the candidate using modified Gram-Schmidt
      for ( k = 1; k <= i; k++ ) {
        vMat_[i+1]->Inner( *(vMat_[k]), hMat_[k][i] );
        vMat_[i+1]->Add( - hMat_[k][i], *(vMat_[k]) );
      }
      hMat_[i+1][i] = vMat_[i+1]->NormEuclid();
      vMat_[i+1]->ScalarDiv( hMat_[i+1][i] );


      // -------------------------------------------------------------------
      //   The new vector increases the number of columns in the H matrix.
      //   We have to make up in this new last column with the Givens
      //   rotations already computed for the upper rows 2 ... i
      // -------------------------------------------------------------------
      for ( k = 2; k <= i; k++ ) {
        aux = hMat_[k-1][i];
        hMat_[k-1][i] =       c_[k-1]  * aux + s_[k-1] * hMat_[k][i];
        hMat_[ k ][i] = -Conj(s_[k-1]) * aux + c_[k-1] * hMat_[k][i];
      }


      // ----------------------------------------------------------------
      //   Eliminate new subdiagonal entry H_i+1,i by a Givens rotation
      // ----------------------------------------------------------------
      givens_->gRot( hMat_[i][i], hMat_[i+1][i], c_[i], s_[i],
                     hMat_[i][i] );


      // ---------------------------------------------------------------
      //   Compute the effect of the Givens rotation on the right-hand
      //   side vector in the least-squares problem
      // ---------------------------------------------------------------
      bVec_[i+1] = -Conj(s_[i]) * bVec_[i];
      bVec_[ i ] =     c_[i]    * bVec_[i];


      // -----------------------------------------------------------------
      //   The (i+1)-th entry of bVec_ now contains the current residual
      //   norm, that is the minimal value of the solution of the least-
      //   squares problem. If this is small enough, we can terminate
      //   this loop.
      // -----------------------------------------------------------------

      resNormOld = resNormNew;
      resNormNew = std::abs(bVec_[i+1]); // needed to obtain a Double

      if ( logging == true ) {
        LogConvergence( resNormNew, i + globNum );
      }

      if ( std::abs(bVec_[i+1]) <= threshold ) {
        (*cla) << "# Inner loop: Approximate is fine" << std::endl;
        approxIsGood = true;
        stepCount = i;
        break;
      }
    }
    (*cla) << "# Inner loop: finished" << std::endl;

    // --------------------------------------------
    //   Compute actual approximation to solution
    // --------------------------------------------

    // Determine minimiser of least squares problem
    // by backward substitution
    bVec_[stepCount] /= hMat_[stepCount][stepCount];
    for ( k = stepCount - 1; k >= 1; k-- ) {
      for ( i = k + 1; i <= stepCount; i++ ) {
        bVec_[k] -= hMat_[k][i] * bVec_[i];
      }
      bVec_[k] /= hMat_[k][k];
    }

    // Compute update vector for new preconditioned
    // solution approximation vector y
    vMat_[1]->ScalarMult( bVec_[1] );
    for ( i = 2; i <= stepCount; i++ ) {
      vMat_[1]->Add( bVec_[i], *(vMat_[i]) );
    }

    // Apply preconditioner to update vector
    precond.Apply( sysMat, *(vMat_[1]), *pVec_ );

    // Update approximate solution (initial guess)
    sol.Add( *pVec_ );

  }


  // ****************************
  //   AllocateHessenbergMatrix
  // ****************************
  template<typename T>
  void GMRESSolver<T>::AllocateHessenbergMatrix() {

    ENTER_FCN( "GMRESSolver::AllocateHessenbergMatrix" );

    bool logging = myParams_->GetBoolValue( "GMRES_logging" );

    // First test, if it is really necessary to allocate
    // a new upper Hessenberg matrix. This is the case,
    // if no matrix has been allocated up to now, or if
    // the maximal dimension of the Krylov subspace has
    // increased.
    if ( hMat_ == NULL || ( maxKrylovDim_ + 1 > hMatNumRows_ ) ) {

      // Report to standard logfile
      if ( logging == true ) {
        (*cla) << "GMRESSolver:\n"
               << " --> Re-allocating upper Hessenberg matrix, since maximal"
               << "\n --> dimension of Krylov subspace has increased to "
               << maxKrylovDim_
               << std::endl;
      }

      // Delete old matrix (safe, even if hMat_ is NULL )
      DeleteHessenbergMatrix();

      // Generate space for new upper Hessenberg matrix
      hMatNumRows_ = maxKrylovDim_ + 1;
      NewArray( hMat_, T*, hMatNumRows_ );
      NewArray( hMat_[1], T, maxKrylovDim_ );
      UInt rowSize = 0;
      UInt offset = 0;
      for ( UInt i = 2; i <= hMatNumRows_; i++ ) {
        rowSize = maxKrylovDim_ - i + 2;
        offset = i - 1;
        hMat_[i] = New T[rowSize];
        AssertMem( hMat_[i], sizeof(T) * rowSize );
        hMat_[i] -= offset;
      }
    }
  }


  // **************************
  //   DeleteHessenbergMatrix
  // **************************
  template<typename T>
  void GMRESSolver<T>::DeleteHessenbergMatrix() {

    ENTER_FCN( "GMRESSolver::DeleteHessenbergMatrix" );

    if ( hMat_ != NULL ) {
      delete[] ( hMat_[1] + 1 );
      for ( UInt i = 2; i <= hMatNumRows_; i++ ) {
        if ( hMat_[i] != NULL ) {
          delete[] ( hMat_[i] + i - 1 );
        }
      }
      DeleteArray( hMat_ );
    }
  }


  // ***********************
  //   AllocateKrylovBasis
  // ***********************
  template<typename T>
  void GMRESSolver<T>::AllocateKrylovBasis( const BaseMatrix &sysMat ) {

    ENTER_FCN( "GMRESSolver::AllocateKrylovBasis" );

    bool logging = myParams_->GetBoolValue( "GMRES_logging" );

    // Determine length of base vectors
    UInt oldLength = 0;
    if ( vMat_.size() > 0 ) {
      oldLength = vMat_[1]->GetSize();
    }

    // First test, if the number of columns of the system matrix
    // or the maximal dimension of the Krylov subspace has increased.
    // Only if this is the case, we must re-allocate the basis.
    if ( vMat_.size() == 0 || problemDim_ != oldLength ||
         maxKrylovDim_ > vMat_.size() ) {

      // Report to standard logfile
      if ( logging == true ) {
        if ( vMat_.size() == 0 ) {
          (*cla) << "GMRESSolver:\n"
                 << " --> Allocating " << maxKrylovDim_ << " Krylov base "
                 << "vectors"
                 << std::endl;
        }
        else if ( problemDim_ != oldLength ) {
          (*cla) << "GMRESSolver:\n"
                 << " --> Re-allocating Krylov base vectors, since column "
                 << "number\n of system matrix has changed. Now is "
                 << problemDim_
                 << std::endl;
        }
        else if ( maxKrylovDim_ > vMat_.size() ) {
          (*cla) << "GMRESSolver:\n"
                 << " --> Re-allocating Krylov base vectors, since maximal\n"
                 << " --> dimension of Krylov subspace has increased to "
                 << maxKrylovDim_
                 << std::endl;
        }
      }

      // Delete old base vectors ( safe for vMat_.size() = 0 )
      DeleteKrylovBasis();

      // Generate space for new basis vectors
      vMat_.reserve( maxKrylovDim_ + 1 );
      vMat_.push_back( NULL );
      for ( UInt i = 1; i <= maxKrylovDim_; i++ ) {
        vMat_.push_back( GenerateVectorObject( sysMat ) );
      }
    }
  }


  // *********************
  //   DeleteKrylovBasis
  // *********************
  template<typename T>
  void GMRESSolver<T>::DeleteKrylovBasis() {

    ENTER_FCN( "GMRESSolver::DeleteKrylovBasis" );

    for ( UInt k = 0; k < vMat_.size(); k++ ) {
      delete vMat_[k];
    }
    vMat_.clear();
  }

}
