#include <cassert>

#include "MatVec/opdefs.hh"
#include "MatVec/generatematvec.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "Utils/ToolsFull.hh"
#include "General/Exception.hh"


#include "OLAS/solver/GMRESSolver.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  GMRESSolver<T>::GMRESSolver( PtrParamNode solverNode, PtrParamNode olasInfo ){
    xml_ = solverNode;

    // Set pointers to communication objects
    infoNode_ = olasInfo->Get("gmres");
    
    // Prepare remaining internal attributes
    c_              = NULL;
    s_              = NULL;
    hMat_           = NULL;
    pVec_           = NULL;
    maxKrylovDim_   = 0;
    problemDim_     = 0;
    hMatNumRows_    = 0;
    consoleConvergence_ = false;

    // Generate object for performing Givens rotations
    givens_ = new GivensRotation( GivensRotation::OLAS );
  }


  // **********************
  //   Default destructor
  // **********************
  template<typename T>
  GMRESSolver<T>::~GMRESSolver() {


    // ---------------
    //   Free memory
    // ---------------

    // Givens coefficients
    delete [] ( c_ );  c_  = NULL;
    delete [] ( s_ );  s_  = NULL;

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
    PrivateSetup( sysMat );
  }


  // ***************************
  //   Setup (private version)
  // ***************************
  template<typename T>
  void GMRESSolver<T>::PrivateSetup( const BaseMatrix &sysMat ) {


    // ----------------------------------------
    //   Preform some consistency checks
    // ----------------------------------------

    // Obtain relevant dimension information
    UInt nCols = sysMat.GetNumCols();
#ifndef NDEBUG
    UInt nRows = sysMat.GetNumRows();
#endif
    UInt newKrylovDim = 50;
    xml_->GetValue("maxKrylovDim", newKrylovDim, ParamNode::INSERT);
    xml_->GetValue("consoleConvergence", consoleConvergence_, ParamNode::INSERT);

    // Check that matrix is square
    assert(nCols == nRows); 

    // Check that the maximal dimension of the Krylov subspace is not
    // larger than the dimension of the problem
    if ( newKrylovDim > nCols ) {
      WARN("maxKrylovDim = " << newKrylovDim << " exceeds matrix "
           << " dimension = " << nCols << ". Using " << nCols
           << " instead.");
      newKrylovDim = nCols;
    }

    // -------------------------------------
    //   (Re-)allocate Givens coefficients
    // -------------------------------------
    if ( newKrylovDim > maxKrylovDim_ ) {
      delete [] ( c_ );  c_  = NULL;
      delete [] ( s_ );  s_  = NULL;
      NEWARRAY( c_, Double, newKrylovDim );
      NEWARRAY( s_, T, newKrylovDim );
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
  void GMRESSolver<T>::Solve( const BaseMatrix &sysMat, const BaseVector &rhs, BaseVector &sol) {

    // ------------------------------------------------------
    //   We call PrivateSetup to (re)-allocate the internal
    //   data structures, if this is necessary
    // ------------------------------------------------------
    PrivateSetup( sysMat );


    // ----------------------
    //   A simple safeguard
    // ----------------------
    if ( hMat_ == NULL ) {
      EXCEPTION( "GMRESSolver::Solve() Detected NULL pointer!" );
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
    Double eps = 1e-6;
    xml_->GetValue("tol", eps, ParamNode::INSERT);
    
    Double tolerance = ComputeThreshold( eps, rhs, *(vMat_[1]), resNorm,
                                         false );

    // ----------------------------------
    //   Perform the GMRES(m) iteration
    // ----------------------------------
    Integer maxIter = 1;
    xml_->GetValue("maxIter", maxIter, ParamNode::INSERT);

    Integer loopsDone = maxIter;
    UInt stepCount = 0;
    UInt stepCountGlobal = 0;
    bool approxIsGood = false;

    for ( Integer k = 1; k <= maxIter; k++ ) {

      // Perform one inner loop
      InnerLoop( sysMat, *ptPrecond_, rhs, sol, resNorm, tolerance,
                 approxIsGood, stepCount, stepCountGlobal );

      // Increase global iteration counter
      stepCountGlobal += stepCount;

      // Compute residual of new approximation and its norm
      // (Will serve in constructing base vector for next inner loop)
      sysMat.CompRes( *(vMat_[1]), sol, rhs );
      resNorm = vMat_[1]->NormL2();

      if(consoleConvergence_ == true){
        std::cout<<"Residual L2 norm of iteration "<<stepCountGlobal<<" = "<<resNorm<<std::endl;
      }

      // If InnerLoop says approximation is fine, test for false
      // convergence. If everything is fine, stop iteration.
      if ( approxIsGood == true ) {
        if ( resNorm <= tolerance ) {
          loopsDone = k;
          break;
        }
        else {
          approxIsGood = false;
        }
      }
    }

    // Make clear, if GMRES failed to converge
    if ( !approxIsGood ) { } // removed logging


    // ----------------------------
    //   Generate solution report
    // ----------------------------

    // Number of iterations: Depends on GMRES(m) -> Full GMRES
    PtrParamNode out = infoNode_->Get(ParamNode::PROCESS)->Get("solver", ParamNode::APPEND);
    
    if ( maxIter == 1 ) {
      out->Get("numIter")->SetValue( (Integer)stepCount );
      out->Get("numGlobalIter")->SetValue( (Integer)stepCount );
    }
    else {
      out->Get("numIter")->SetValue( loopsDone );
      out->Get("numGlobalIter")->SetValue( (Integer)stepCountGlobal );
    }

    // Final relative residual norm
    out->Get("finalNorm")->SetValue(resNorm);

    // Status of solution
    out->Get("solutionIsOkay")->SetValue(resNorm);

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


    UInt i, k;
    T aux = 0;
//    Double resNormNew = beta;
    approxIsGood = false;

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
      ptPrecond_->Apply( sysMat, *(vMat_[i]), *pVec_ );

      // Candidate for next basis vector
      sysMat.Mult( *pVec_, *(vMat_[i+1]) );

      // Now orthonormalise the candidate using modified Gram-Schmidt
      for ( k = 1; k <= i; k++ ) {
        vMat_[i+1]->Inner( *(vMat_[k]), hMat_[k][i] );
        vMat_[i+1]->Add( - hMat_[k][i], *(vMat_[k]) );
      }
      hMat_[i+1][i] = vMat_[i+1]->NormL2();
      vMat_[i+1]->ScalarDiv( hMat_[i+1][i] );


      // -------------------------------------------------------------------
      //   The new vector increases the number of columns in the H matrix.
      //   We have to make up in this new last column with the Givens
      //   rotations already computed for the upper rows 2 ... i
      // -------------------------------------------------------------------
      for ( k = 2; k <= i; k++ ) {
        aux = hMat_[k-1][i];
        hMat_[k-1][i] =       c_[k-1]  * aux + s_[k-1] * hMat_[k][i];
        hMat_[ k ][i] = -conj(s_[k-1]) * aux + c_[k-1] * hMat_[k][i];
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
      bVec_[i+1] = -conj(s_[i]) * bVec_[i];
      bVec_[ i ] =       c_[i]  * bVec_[i];


      // -----------------------------------------------------------------
      //   The (i+1)-th entry of bVec_ now contains the current residual
      //   norm, that is the minimal value of the solution of the least-
      //   squares problem. If this is small enough, we can terminate
      //   this loop.
      // -----------------------------------------------------------------

//      resNormNew = std::abs(bVec_[i+1]); // needed to obtain a Double

      if ( std::abs(bVec_[i+1]) <= threshold ) {
        approxIsGood = true;
        stepCount = i;
        break;
      }
    }

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
    ptPrecond_->Apply( sysMat, *(vMat_[1]), *pVec_ );

    // Update approximate solution (initial guess)
    sol.Add( *pVec_ );

  }


  // ****************************
  //   AllocateHessenbergMatrix
  // ****************************
  template<typename T>
  void GMRESSolver<T>::AllocateHessenbergMatrix() {

    // First test, if it is really necessary to allocate
    // a new upper Hessenberg matrix. This is the case,
    // if no matrix has been allocated up to now, or if
    // the maximal dimension of the Krylov subspace has
    // increased.
    if ( hMat_ == NULL || ( maxKrylovDim_ + 1 > hMatNumRows_ ) ) {

      // Delete old matrix (safe, even if hMat_ is NULL )
      DeleteHessenbergMatrix();

      // Generate space for new upper Hessenberg matrix
      hMatNumRows_ = maxKrylovDim_ + 1;
      hMat_ = new T*[hMatNumRows_];
      assert(hMatNumRows_ > 1);
      hMat_[1] = new T[hMatNumRows_];
      UInt rowSize = 0;
      UInt offset = 0;
      for ( UInt i = 2; i < hMatNumRows_; i++ ) {
        rowSize = maxKrylovDim_ - i + 2;
        offset = i - 1;
        hMat_[i] = new T[rowSize];
        ASSERTMEM( hMat_[i], sizeof(T) * rowSize );
        hMat_[i] -= offset;
      }
    }
  }


#pragma GCC diagnostic push
// GMRESSolver.cc:452:11: warning: ‘void operator delete [](void*)’ called on pointer ‘*prephitmp_110 + ivtmp.619_119’ with nonzero offset [8, 34359738344] [-Wfree-nonheap-object]
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"

  // **************************
  //   DeleteHessenbergMatrix
  // **************************
  template<typename T>
  void GMRESSolver<T>::DeleteHessenbergMatrix()
  {
    if ( hMat_ != NULL )
    {
      delete [] ( hMat_[1]  );
      for ( UInt i = 2; i < hMatNumRows_; i++ )
      {
        if ( hMat_[i] != NULL )
          delete [] ( hMat_[i] + i - 1 ); // that appears to be wrong (compared to hMat_[i] = new T[rowSize] but chaning the line causes segfaults?!
      }
      delete [] ( hMat_ );  hMat_  = NULL;
    }
  }

#pragma GCC diagnostic pop

  // ***********************
  //   AllocateKrylovBasis
  // ***********************
  template<typename T>
  void GMRESSolver<T>::AllocateKrylovBasis( const BaseMatrix &sysMat ) {

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


    for ( UInt k = 0; k < vMat_.size(); k++ ) {
      delete vMat_[k];
    }
    vMat_.clear();
  }

  // Explicit template instantiation
  template class GMRESSolver<Double>;
  template class GMRESSolver<Complex>;
  
}
