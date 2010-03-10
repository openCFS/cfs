// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/opdefs.hh"
#include "MatVec/generatematvec.hh"

#include "OLAS/solver/minres.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  MINRESSolver<T>::MINRESSolver( ParamNode* solverNode, InfoNode *olasInfo){
    xml_ = solverNode;

    // Set pointers to communication objects
    solverInfo_ = olasInfo->Get("minres");

    // Prepare remaining internal attributes
    w0_ = NULL;
    w1_ = NULL;
    w2_ = NULL;
    qN_ = NULL;
    q0_ = NULL;
    q1_ = NULL;
    pV_ = NULL;
    bV_ = NULL;
    vectorLength_ = 0;

    // Generate object for performing Givens rotations
    givens_ = new GivensRotation( GivensRotation::OLAS );

    // Generate object for computing Lanczos vectors
    lanczos_ = new LanczosMethod;
  }


  // **********************
  //   Default destructor
  // **********************
  template<typename T>
  MINRESSolver<T>::~MINRESSolver() {


    // Delete auxillary math objects
    delete givens_;
    delete lanczos_;

    // Delete vectors for 3-term recurrences
    delete w0_;
    delete w1_;
    delete w2_;
    delete qN_;
    delete q0_;
    delete q1_;
    delete pV_;
    delete bV_;
  }


  // **************************
  //   Setup (public version)
  // **************************
  template<typename T>
  void MINRESSolver<T>::Setup( BaseMatrix &sysMat, InfoNode* analysis_step) {
    PrivateSetup( sysMat );
  }


  // ***************************
  //   Setup (private version)
  // ***************************
  template<typename T>
  void MINRESSolver<T>::PrivateSetup( const BaseMatrix &sysMat) {


    // Test, whether the problem dimension has changed
    UInt newDim_ = sysMat.GetNumCols();

    if ( newDim_ > vectorLength_ ) {

      bool allocate = false;

      // Check, if there are already vector objects
      if ( w0_ != NULL ) {

        // Maybe, we can simply re-size them
        Vector<T> *ww0 = dynamic_cast<Vector<T>*>( w0_ );
        Vector<T> *ww1 = dynamic_cast<Vector<T>*>( w1_ );
        Vector<T> *ww2 = dynamic_cast<Vector<T>*>( w2_ );
        Vector<T> *qqN = dynamic_cast<Vector<T>*>( qN_ );
        Vector<T> *qq0 = dynamic_cast<Vector<T>*>( q0_ );
        Vector<T> *qq1 = dynamic_cast<Vector<T>*>( q1_ );
        Vector<T> *ppV = dynamic_cast<Vector<T>*>( pV_ );
        Vector<T> *bbV = dynamic_cast<Vector<T>*>( bV_ );
        if ( ww0 != NULL ) {
          ww0->Resize( newDim_ );
          ww1->Resize( newDim_ );
          ww2->Resize( newDim_ );
          qqN->Resize( newDim_ );
          qq0->Resize( newDim_ );
          qq1->Resize( newDim_ );
          ppV->Resize( newDim_ );
          bbV->Resize( newDim_ );
        }

        // Else delete them and allocate new ones
        else {
          delete w0_;
          delete w1_;
          delete w2_;
          delete qN_;
          delete q0_;
          delete q1_;
          delete pV_;
          delete bV_;
          allocate = true;
        }
      }

      // If not, we must generate vector objects
      else {
        allocate = true;
      }

      // Generate vector objects, if requried
      if ( allocate == true ) {
        w0_ = GenerateVectorObject( sysMat );
        w1_ = GenerateVectorObject( sysMat );
        w2_ = GenerateVectorObject( sysMat );
        qN_ = GenerateVectorObject( sysMat );
        q0_ = GenerateVectorObject( sysMat );
        q1_ = GenerateVectorObject( sysMat );
        pV_ = GenerateVectorObject( sysMat );
        bV_ = GenerateVectorObject( sysMat );
      }
    }


    // Now initialise auxilliary vectors to zero
    //
    // NOTE: This is probably not necessary, but
    //       before removing it carefully check algorithm
    w0_->Init();
    w1_->Init();
    w2_->Init();
    qN_->Init();
    q0_->Init();
    q1_->Init();
    pV_->Init();
    bV_->Init();

  }


  // *********
  //   Solve
  // *********
  template<typename T>
  void MINRESSolver<T>::Solve( const BaseMatrix &sysMat,
                               const BasePrecond &precond,
                               const BaseVector &rhs, BaseVector &sol, InfoNode* analysis_step ) {


    // ----------------------------------------
    //   Let private setup test the type of
    //   the system matrix and determine,
    //   whether the auxilliary vectors must
    //   be re-sized.
    // ---------------------------------------
    PrivateSetup( sysMat );


    // ------------------------
    //   Auxilliary variables
    // ------------------------
    T beta2, beta1, beta0;
    T l2, l1, l0;
    T alpha;
    double cOld, cNew;
    T sOld, sNew;
    T aux, tmp;
    double rho;
    BaseVector *tmpVec = NULL;


    // ---------------------------------
    //   Prepare start-up of algorithm
    // ---------------------------------

    // Determine norm of preconditioned right hand side
    precond.Apply( sysMat, rhs, *pV_ );
    Double rhsNorm = pV_->NormL2();

    // Compute residual of initial guess
    sysMat.CompRes( *pV_, sol, rhs );

    // Apply preconditioner
    precond.Apply( sysMat, *pV_, *q0_ );

    // Compute norm of residual of preconditioned system
    rho = q0_->NormL2();

    // Be verbose
    (*cla) << "\n MINRESSolver:\n"
           << "\n Norm of rhs              = " << rhs.NormL2()
	   << "\n Norm of precond rhs      = " << rhsNorm
	   << "\n Norm of residual         = " << pV_->NormL2()
	   << "\n Norm of precond residual = " << rho
	   << std::endl;

    // Determine first base and update vector
    q0_->ScalarDiv( rho );

    // Initialise right-hand side of least-squares problem
    bV_->SetEntry( 1, rho );

    // For start of recurrence, set last but one vectors
    // and some coefficients
    q1_->Init();
    w2_->Init();
    w1_->Init();
    beta2 = 0;
    beta1 = 0;
    cOld  = 0;
    sOld  = 0;
    cNew  = 1;
    sNew  = 0;

    // Compute stopping threshold for loop
    double threshold = 1e-6;

    ParamNode *sNode = NULL;
    sNode = xml_->Get("minres", false);
    if(sNode) {
      sNode->Get("tol", threshold, false);
    }
    
    threshold *= rhsNorm;


    // ----------------------
    //   Initialise logging
    // ----------------------
    bool logging = false;
    if ( logging == true ) {
      LogConvergence( rho, 0, true );
    }


    // ------------------------------
    //   Main loop of the algorithm
    // ------------------------------
    UInt maxIter = 50;
    if(sNode) {
      sNode->Get("maxIter", maxIter, false);
    }
    UInt k = 1;
    bool loop = true;
    while ( k <= maxIter && loop == true ) {


      // --------------------------------
      //   Compute basis vector q_(k+1)
      //   of the Krylov subspace
      // --------------------------------
      lanczos_->CompNextVector( sysMat, precond, *q1_, *q0_, *qN_, *pV_,
				alpha, beta1, beta0 );


#ifdef DEBUG_MINRES
      std::cerr << "\n Norm of basis vector q[" << k-1 << "] = "
		<< q1_->NormL2();
      std::cerr << "\n Norm of basis vector q[" <<  k  << "] = "
		<< q0_->NormL2();
      std::cerr << "\n Norm of basis vector q[" << k+1 << "] = "
		<< qN_->NormL2();
      std::cerr << "\n beta2 = " << beta2
		<< "\n beta1 = " << beta1
		<< "\n beta0 = " << beta0
		<< "\n alpha = " << alpha
		<< std::endl;
#endif

      // ---------------------------
      //   Perform Givens rotation
      // ---------------------------

      // Determine l2 by replaying previous Givens
      // rotation from step (k-2)
      l2 = cOld * beta2;

      // Replay Givens rotation from step (k-1)
      // This yields the final l1 and alters l0
      l1 = cNew * beta0 +      sNew  * alpha;
      l0 = cNew * alpha - Conj(sNew) * beta0;

      // Shift Givens coefficients
      cOld = cNew;
      sOld = sNew;

      // Compute Givens rotation for step k
      // This yields the new coefficients and
      // the final l0
      givens_->gRot( l0, beta0, cNew, sNew, l0 );

      // Compute the effect of the Givens rotation on the
      // right-hand side vector of the least-squares problem
      bV_->GetEntry( k, aux );

      tmp = -Conj(sNew) * aux;
      bV_->SetEntry( k+1, tmp );

      tmp =       cNew  * aux;
      bV_->SetEntry(  k , tmp );

#ifdef DEBUG_MINRES
      std::cerr << "\n cNew = " << cNew
		<< "\n sNew = " << sNew
		<< "\n b[" << k << "] = " << aux << std::endl;
      bV_->GetEntry( k+1, aux );
      std::cerr << "\n b[" << k+1 << "] = " << aux
		<< std::endl;
#endif


      // --------------------------------------
      //   Determine update vector by 3-term
      //   recurrence and compute iterate x_k
      // --------------------------------------
      w0_->Add( -l1, *w1_, -l2, *w2_ );
      w0_->Add( *q0_ );
      w0_->ScalarDiv( l0 );
      bV_->GetEntry( k, aux );
      sol.Add( cNew * aux, *w0_ );


      // -------------------------
      //   Cycle Lanczos vectors
      // -------------------------
      tmpVec = q1_;
      q1_ = q0_;
      q0_ = qN_;
      qN_ = tmpVec;


      // ---------------------------
      //   Test stopping threshold
      // ---------------------------
      bV_->GetEntry( k+1, aux );
      if ( Abs(aux) <= threshold ) {

        // Test for false convergence
        sysMat.CompRes( *pV_, sol, rhs );
        precond.Apply( sysMat, *pV_, *q0_ );
        rho = q0_->NormL2();
        if ( rho > threshold ) {
          WARN(" MINRESSolver::Solve\n"
		           << " False convergence detected.\n"
      		     << " Predicted res.norm = " << Abs(aux) << '\n'
		           << " Actual res.norm = " << rho);
        }
	else {
	  loop = false;
	}
      }

      // Log convergence
      if ( logging == true ) {
	LogConvergence( Abs(aux), k );
      }


      // Increase loop counter
      k++;

    }

    // Compute real residual of preconditioned system
    sysMat.CompRes( *pV_, sol, rhs );
    precond.Apply( sysMat, *pV_, *q0_ );
    rho = q0_->NormL2();

    // Compose report
    InfoNode* out = solverInfo_->Get(InfoNode::PROCESS)->Get("solver", InfoNode::APPEND);
    out->Get("finalNorm")->SetValue(rho);
    out->Get("numIter")->SetValue((Integer)(k-1));

  }


  // ******************************
  //   Test Type of Input Matrix
  // ******************************
  template<typename T>
  void MINRESSolver<T>::TestMatrixType( const BaseMatrix &sysMat ) const {


    if ( sysMat.GetStructureType() == BaseMatrix::SBM_MATRIX ) {
      WARN( "MINRESSolver expects matrix entries to be scalars!"
            " We do not test this for SBM_Matrix class" );
    }
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class MINRESSolver<Double>;
  template class MINRESSolver<Complex>;
#endif

}
