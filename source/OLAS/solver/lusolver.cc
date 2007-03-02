#include "solver/lusolver.hh"

// Include source code of CroutLU class for template instantiation
// Note: Might lead to double instantiation, since CroutLU is also
// used in ILUTP_Precond. Going to implement better concept as soon
// as time permits.
#include "utils/math/croutlu.cc"

namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  LUSolver<T>::LUSolver( OLAS_Params *myParams, OLAS_Report *myReport ) {

    ENTER_FCN( "LUSolver::LUSolver" );

    // Set pointers to communication objects
    myParams_ = myParams;
    myReport_ = myReport;

    // No factorisation was performed yet
    amFactorised_ = false;

    // With CFS we only treat structurally symmetric problems (in principle)
    myParams_->SetValue( "LUSOLVER_symPattern", true );

  }

  // **************
  //   Destructor
  // **************
  template<typename T>
  LUSolver<T>::~LUSolver() {

    ENTER_FCN( "LUSolver::~LUSolver" );

  }

  // ************************
  //   Forced Instantiation
  // ************************
  template<typename T>
  void LUSolver<T>::InstantiateAdditionalPublicMethods( BaseMatrix &sysMat ) {
    ENTER_FCN( "LUSolver::InstantiateAdditionalPublicMethods" );
    this->ExportILUFactorisation( "dummy.mtx" );
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void LUSolver<T>::Setup( BaseMatrix &sysMat ) {

    ENTER_FCN( "LUSolver::Setup" );

    // Check that we have a StdMatrix
    if ( sysMat.GetStructureType() != STDMATRIX ) {
      Error( "LUSolver cannot deal with matrices other than StdMatrix",
	     __FILE__, __LINE__ );
    }

    TRY_CAST {

      // Down-cast to StdMatrix
      RefCast( sysMat, StdMatrix, stdMat );

      // Now test the storage layout
      MatrixStorageType sType = stdMat.GetStorageType();
      if ( sType != SPARSE_NONSYM ) {
	(*error) << "LUSolver::Setup: The LUSolver requires the system matrix"
		 << " to be a CRS_Matrix, i.e. sparseNonSym. The system matrix"
		 << " you supplied is a matrix in " << Enum2String( sType )
		 << " format.";
	Error( __FILE__, __LINE__ );
      }

      // Down-cast to CRS_Matrix
      RefCast( sysMat, CRS_Matrix<T>, crsMat );

      // Logging
      (*cla) << " -----------------------------------------------\n"
	     << " LUSolver: Factorisation of a "
	     << crsMat.GetNrows() << " x " << crsMat.GetNcols()
	     << " matrix (nnz = " << crsMat.GetNnz() << ")"
	     << std::endl;

      // Set estimate for growth of sparsity pattern
      // Being generous here improves performance
      this->memGrowthEstimate_ = 25;

      // Perform the factorisation
      Factorise( crsMat );
      amFactorised_ = true;

      // Logging
      (*cla) << " -----------------------------------------------"
	     << std::endl;

    } CATCH_CAST;

    // If the user wishes, we can export the LU factorisation to a file
    if ( myParams_->GetBoolValue( "CROUT_saveFacToFile" ) ) {
      std::string filename;
      filename = myParams_->GetStringValue( "CROUT_facFileName" );
      this->ExportILUFactorisation( filename.c_str() );
    }

  }

  // *********
  //   Solve
  // *********
  template<typename T>
  void LUSolver<T>::Solve( const BaseMatrix &sysMat,
			   const BasePrecond &precond,
			   const BaseVector &rhs, BaseVector &sol ) {

    ENTER_FCN( "LUSolver::Solve" );

    // Test that a factorisation is available, if not issue a warning.
    // Note: We cannot initiate the factorisation from here, since we
    // only have a const reference, but the factorisation requires a
    // non-const reference (currently only for sorting the CRS matrix)
    if ( amFactorised_ == false ) {
      (*error) << "LUSolver::Solve: No factorisation available. Call Setup() "
	       << "first!";
      Error( __FILE__, __LINE__ );
    }

    // Solve the problem
    TRY_CAST {
      ConstRefCast( rhs, Vector<T>, myRHS );
      RefCast( sol, Vector<T>, mySol );

      // Logging
      bool logging = myParams_->GetBoolValue( "LUSOLVER_logging" );
      if ( logging ) {
        (*cla) << " -----------------------------------------------\n"
               << " LUSolver: Solving a problem with "
               << sysMat.GetNcols() << " unknowns" << std::endl;
      }

      // Actual solve is done by CroutLU class
      CroutLU<T>::Solve( myRHS, mySol );

      // If desired perform iterative refinement
      UInt numSteps = (UInt)myParams_->GetIntValue( "LUSOLVER_itRefSteps" );
      UInt logLevel =
        (UInt)myParams_->GetIntValue( "LUSOLVER_itRefVerbosity" );

      if ( numSteps > 0 ) {

        // Avoid recursion, we do not want to do refinement on the
        // solution in the refinement step
        myParams_->SetValue( "LUSOLVER_itRefSteps", (Integer)0 );
        myParams_->SetValue( "LUSOLVER_logging", false );

        // Refine
        iterativeRefiner_.Refine( (*this), sysMat, sol, rhs, numSteps,
                                  logLevel );

        // Re-set parameter object
        myParams_->SetValue( "LUSOLVER_itRefSteps", (Integer)numSteps );
        myParams_->SetValue( "LUSOLVER_logging", logging );

      }

      // Logging
      if ( logging ) {
        (*cla) << " -----------------------------------------------"
               << std::endl;
      }

    } CATCH_CAST;


    // Generate Report

    // Now this currently is of dubious value, since the two things queried
    // from olasReport are actually meaningless in the context of a direct
    // solver. Nevertheless we supply some values for consistency
    if ( myReport_ != NULL ) {
      myReport_->SetValue( "numIter", -1 );
      myReport_->SetValue( "finalNorm", -1.0 );
    }

  }

  // *************************
  //   ComputeRequiredMemory
  // *************************
  template<typename T>
  void LUSolver<T>::ComputeRequiredMemory( UInt &maxColLengthL,
					   UInt &maxRowLengthU,
					   UInt &maxEntriesL,
					   UInt &maxEntriesU,
					   CRS_Matrix<T> &crsMat ) {

    ENTER_FCN( "LUSolver:ComputeRequiredMemory:" );


    // If matrix is unsorted, change to LEX format
    if ( crsMat.GetCurrentLayout() == CRS_Matrix<T>::UNSORTED ) {
      crsMat.ChangeLayout( CRS_Matrix<T>::LEX );
    }


    // ================================
    //  Problem structurally symmetric
    // ================================
    if ( myParams_->GetBoolValue( "LUSOLVER_symPattern" ) == true ) {

      UInt i;
      UInt bwGlobal = 0;
      UInt bwLocal  = 0;
      UInt profile  = 0;
      Integer aux;

      // Get hold of index arrays
      const Integer *cidx = crsMat.GetColPointer();
      const Integer *rptr = crsMat.GetRowPointer();
      const Integer *diag = crsMat.GetDiagPointer();

      // Implementation for LEX_DIAG_FIRST sub-format
      if ( crsMat.GetCurrentLayout() == CRS_Matrix<T>::LEX_DIAG_FIRST ) {

        for ( i = 1; i <= crsMat.GetNrows(); i++ ) {

          // Determine bandwidth of this row
          aux = cidx[ rptr[i+1] - 1 ] - i;
          bwLocal = aux > 0 ? (UInt)aux : 0;

          // Add to profile
          profile += bwLocal;

          // Compare to global bandwidth
          bwGlobal = bwLocal > bwGlobal ? bwLocal : bwGlobal;
        }
      }
      
      // Implementation for LEX sub-format
      else if ( crsMat.GetCurrentLayout() == CRS_Matrix<T>::LEX ) {

        for ( i = 1; i <= crsMat.GetNrows(); i++ ) {

          // Determine bandwidth of this row
          aux = cidx[rptr[i]] - i;
          bwLocal = aux > 0 ? (UInt)aux : 0;

          // Add to profile
          profile += bwLocal;

          // Compare to global bandwidth
          bwGlobal = bwLocal > bwGlobal ? bwLocal : bwGlobal;
        }
      }

      // Should not happen
      else {
        (*error) << "Reached a branch that should not exist!";
        Error( __FILE__, __LINE__ );
      }

      // For the L factor we do not need to store the ones
      // on the diagonal
      maxColLengthL = bwGlobal;
      maxEntriesL = profile;
  
      // For the U factor we need the diagonal entries
      maxRowLengthU = bwGlobal + 1;
      maxEntriesU = profile + crsMat.GetNrows();
    }

    // ====================================
    //  Problem not structurally symmetric
    // ====================================
    else {
      Error( "Missing implementation! Go find a code monkey ;-)",
	     __FILE__, __LINE__ );
    }

  }

}
