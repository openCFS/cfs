// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/crs_matrix.hh"

// Include source code of CroutLU class for template instantiation
// Note: Might lead to double instantiation, since CroutLU is also
// used in ILUTP_Precond. Going to implement better concept as soon
// as time permits.
#include "OLAS/utils/math/croutlu.hh"

#include "OLAS/solver/lusolver.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  LUSolver<T>::LUSolver( ParamNode* solverNode, InfoNode *olasInfo ) {


    // Set pointers to communication objects
    xml_ = solverNode;
    solverInfo_ = olasInfo->Get("directLU");

    // No factorisation was performed yet
    amFactorised_ = false;
    
    itRefSteps_ = 2;
    ParamNode *sNode = NULL;
    sNode = xml_->Get("directLU", false);
    sNode->Get("itRefSteps", itRefSteps_, false);
  }

  // **************
  //   Destructor
  // **************
  template<typename T>
  LUSolver<T>::~LUSolver() {


  }

  // *********
  //   Setup
  // *********
  template<typename T>
  void LUSolver<T>::Setup( BaseMatrix &sysMat, InfoNode* analysis_step ) {


    // Check that we have a StdMatrix
    if ( sysMat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "LUSolver cannot deal with matrices other than StdMatrix" );
    }

    StdMatrix& stdMat = dynamic_cast<StdMatrix&>(sysMat);

      // Now test the storage layout
      BaseMatrix::StorageType sType = stdMat.GetStorageType();
      if ( sType != BaseMatrix::SPARSE_NONSYM ) {
        EXCEPTION( "LUSolver::Setup: The LUSolver requires the system matrix"
            << " to be a CRS_Matrix, i.e. sparseNonSym. The system matrix"
            << " you supplied is a matrix in " << BaseMatrix::storageType.ToString( sType )
            << " format." );
      }

      // Down-cast to CRS_Matrix
      CRS_Matrix<T>& crsMat = dynamic_cast<CRS_Matrix<T>&>(sysMat);

      // Logging
      (*cla) << " -----------------------------------------------\n"
	     << " LUSolver: Factorisation of a "
	     << crsMat.GetNumRows() << " x " << crsMat.GetNumCols()
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

    // If the user wishes, we can export the LU factorisation to a file
    bool saveFacToFile = false;
    std::string facFileName = "fac.out";
    
    ParamNode *sNode = NULL;
    sNode = xml_->Get("directLU", false);
    if(sNode) {
      if(sNode->Has("saveFacFile")) {
        saveFacToFile = true;
        sNode->Get("saveFacFile", facFileName, false);
      }
    }
    
    if ( saveFacToFile ) {
      std::string filename;
      this->ExportILUFactorisation( facFileName.c_str() );
    }

  }

  // *********
  //   Solve
  // *********
  template<typename T>
  void LUSolver<T>::Solve( const BaseMatrix &sysMat,
			   const BasePrecond &precond,
			   const BaseVector &rhs, BaseVector &sol, InfoNode* analysis_step ) {


    // Test that a factorisation is available, if not issue a warning.
    // Note: We cannot initiate the factorisation from here, since we
    // only have a const reference, but the factorisation requires a
    // non-const reference (currently only for sorting the CRS matrix)
    if ( amFactorised_ == false ) {
      EXCEPTION( "LUSolver::Solve: No factorisation available. Call Setup() "
	       << "first!" );
    }

    // Solve the problem
    const Vector<T>& myRHS = dynamic_cast<const Vector<T>&>(rhs);
    Vector<T>& mySol = dynamic_cast<Vector<T>&>(sol);

    // Logging
    bool logging = false;
    if ( logging ) {
      (*cla) << " -----------------------------------------------\n"
             << " LUSolver: Solving a problem with "
             << sysMat.GetNumCols() << " unknowns" << std::endl;
    }

    // Actual solve is done by CroutLU class
    CroutLU<T>::Solve( myRHS, mySol );

    // If desired perform iterative refinement
    UInt logLevel = 2;
    UInt numSteps = itRefSteps_;
      
    ParamNode *sNode = NULL;
    sNode = xml_->Get("directLU", false);
    sNode->Get("itRefVerbosity", logLevel, false);
      

    if ( numSteps > 0 ) {

        // Avoid recursion, we do not want to do refinement on the
        // solution in the refinement step
        itRefSteps_ = 0;

        // Refine
        iterativeRefiner_.Refine( (*this), sysMat, sol, rhs, numSteps,
                                  logLevel );

        // Re-set parameter object
        itRefSteps_ = numSteps;

    }

      // Logging
      if ( logging ) {
        (*cla) << " -----------------------------------------------"
               << std::endl;
      }

    // Generate Report

    // Now this currently is of dubious value, since the two things queried
    // from olasReport are actually meaningless in the context of a direct
    // solver. Nevertheless we supply some values for consistency
    InfoNode* out = solverInfo_->Get(InfoNode::PROCESS)->Get("solver", InfoNode::APPEND);
    out->Get("numIter")->SetValue(-1);
    out->Get("finalNorm")->SetValue(-1.0);

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



    // If matrix is unsorted, change to LEX format
    if ( crsMat.GetCurrentLayout() == CRS_Matrix<T>::UNSORTED ) {
      crsMat.ChangeLayout( CRS_Matrix<T>::LEX );
    }


    // ================================
    //  Problem structurally symmetric
    // ================================

    // With CFS we only treat structurally symmetric problems (in principle)
    bool symPattern = true;
    
    if ( symPattern ) {

      UInt i;
      UInt bwGlobal = 0;
      UInt bwLocal  = 0;
      UInt profile  = 0;
      Integer aux;

      // Get hold of index arrays
      const UInt *cidx = crsMat.GetColPointer();
      const UInt *rptr = crsMat.GetRowPointer();
      // TODO: Check if this is still needed      
      // const UInt *diag = crsMat.GetDiagPointer();

      // Implementation for LEX_DIAG_FIRST sub-format
      if ( crsMat.GetCurrentLayout() == CRS_Matrix<T>::LEX_DIAG_FIRST ) {

        for ( i = 0; i < crsMat.GetNumRows(); i++ ) {

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

        for ( i = 0; i < crsMat.GetNumRows(); i++ ) {

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
        EXCEPTION( "Reached a branch that should not exist!" );
      }

      // For the L factor we do not need to store the ones
      // on the diagonal
      maxColLengthL = bwGlobal;
      maxEntriesL = profile;
  
      // For the U factor we need the diagonal entries
      maxRowLengthU = bwGlobal + 1;
      maxEntriesU = profile + crsMat.GetNumRows();
    }

    // ====================================
    //  Problem not structurally symmetric
    // ====================================
    else {
      EXCEPTION( "Missing implementation! Go find a code monkey ;-)" );
    }

  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class LUSolver<Double>;
  template class LUSolver<Complex>;
#endif
  
}
