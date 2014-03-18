#include "MatVec/BaseMatrix.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/BaseVector.hh"
#include "MatVec/generatematvec.hh"
#include "OLAS/precond/IdPrecondStd.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "IterativeRefinement.hh"

namespace CoupledField {


  IterativeRefinement::IterativeRefinement() :
    residual_(NULL),
    update_(NULL)
  {
  }
  
  IterativeRefinement::~IterativeRefinement() {
    delete residual_;
    delete update_;
  }

  // *****************************
  //   GenerateAuxilliaryVectors
  // *****************************
  void IterativeRefinement::
  GenerateAuxilliaryVectors( const BaseMatrix &sysMat ) {

    const StdMatrix& stdMat = dynamic_cast<const StdMatrix&>(sysMat);

    bool genNewVectors = false;

    // Check whether we must generate new auxilliary vectors
    if ( residual_ == NULL || residual_->GetSize() != stdMat.GetNumCols() ) {
      genNewVectors = true;
    }

    // Generate new vectors (discarding old ones) if necessary
    if ( genNewVectors == true ) {

      delete residual_;
      delete update_;

      residual_ = GenerateSingleVectorObject( stdMat.GetEntryType(),
                                              stdMat.GetNumCols() );

      update_   = GenerateSingleVectorObject( stdMat.GetEntryType(),
                                              stdMat.GetNumCols() );
    }
  }


  // **********
  //   Refine
  // **********
  void IterativeRefinement::Refine( BaseDirectSolver &mySolver,
                                    const BaseMatrix &sysMat,
                                    BaseVector &sol,
                                    const BaseVector &rhs,
                                    UInt& numSteps, UInt logLevel ) {


//    if(!numSteps)
//      return;
      
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
               << residual_->NormL2()
               << std::endl;
        cla->unsetf( std::ios::scientific | std::ios::fixed );
      }

      // STEP 2: Solve update equation ( A * dx = r )
      mySolver.Solve( sysMat, *residual_, *update_ );

      // STEP 3: Perform upate ( x <- x + dx )
      sol.Add( *update_ );
    }

    // For report compute final residual norm
    if ( logLevel > 1 ) {
      sysMat.CompRes( *residual_, sol, rhs );
      (*cla) << " step " << numSteps << ": residual norm = "
             << std::scientific
             << residual_->NormL2()
             << std::endl;
    }
  }

}
