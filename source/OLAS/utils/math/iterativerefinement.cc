#include "iterativerefinement.hh"
#include "solver/solver.hh"

namespace OLAS {


  // *****************************
  //   GenerateAuxilliaryVectors
  // *****************************
  void IterativeRefinement::
  GenerateAuxilliaryVectors( const BaseMatrix &sysMat ) {

    ENTER_FCN( "IterativeRefinement::GenerateAuxilliaryVectors" );

    TRY_CAST {

      ConstRefCast( sysMat, StdMatrix, stdMat );

      bool genNewVectors = false;

      // Check whether we must generate new auxilliary vectors
      if ( residual_ == NULL ) {
        genNewVectors = true;
      }
      else if ( residual_->GetSize() != stdMat.GetNcols() ) {
        genNewVectors = true;
      }

      // Generate new vectors (discarding old ones) if necessary
      if ( genNewVectors == true ) {

        delete residual_;
        delete update_;

        residual_ = GenerateStdVectorObject( stdMat.GetStorageType(),
                                             stdMat.GetEntryType(),
                                             stdMat.GetBlockSize(),
                                             stdMat.GetNcols() );

        update_   = GenerateStdVectorObject( stdMat.GetStorageType(),
                                             stdMat.GetEntryType(),
                                             stdMat.GetBlockSize(),
                                             stdMat.GetNcols() );
      }
    } CATCH_CAST;
  }


  // **********
  //   Refine
  // **********
  void IterativeRefinement::Refine( BaseDirectSolver &mySolver,
                                    const BaseMatrix &sysMat,
                                    BaseVector &sol,
                                    const BaseVector &rhs,
                                    UInt numSteps, UInt logLevel ) {

    ENTER_FCN( "IterativeRefinement::Refine" );

    // Be verbose
    if ( logLevel > 0 ) {
      (*cla) << " Performing " << numSteps << " steps of iterative"
             << " refinement:" << std::endl;
    }

    // Check, if we must generate auxilliary vectors
    // and do it if necessary
    GenerateAuxilliaryVectors( sysMat );


    for ( UInt i = 1; i <= numSteps; i++ ) {

      // STEP 1: Compute Residual (r = b - A * x)
      sysMat.CompRes( *residual_, sol, rhs );

      // Report residual norm
      if ( logLevel > 1 ) {
        (*cla) << " step " << i-1 << ": residual norm = "
               << std::scientific
               << residual_->NormEuclid()
               << std::endl;
        cla->unsetf( std::ios::scientific | std::ios::fixed );
      }

      // STEP 2: Solve update equation ( A * dx = r )
      mySolver.Solve( sysMat, *dummyPrecond_, *residual_, *update_ );

      // STEP 3: Perform upate ( x <- x + dx )
      sol.Add( *update_ );
    }

    // For report compute final residual norm
    if ( logLevel > 1 ) {
      sysMat.CompRes( *residual_, sol, rhs );
      (*cla) << " step " << numSteps << ": residual norm = "
             << std::scientific
             << residual_->NormEuclid()
             << std::endl;
    }

  }

}
