// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <algorithm>

#include <complex>

#include "MatVec/crs_matrix.hh"
#include "OLAS/precond/ilutpprecond.hh"

// Include source code of CroutLU class for template instantiation
// Note: Might lead to double instantiation, since CroutLU is also
// used in LUSolver. Going to implement better concept as soon as
// time permits.
#include "OLAS/utils/math/croutlu.hh"

namespace CoupledField {

  // =====================================================
  //   Constructor (for use in GenerateStdPrecondObject)
  // =====================================================
  template <typename T>
  ILUTP_Precond<T>::ILUTP_Precond( const StdMatrix& stdMat, 
                                   PtrParamNode solverNode,
				   PtrParamNode olasInfo ) {


    // Set pointers to communication objects
    this->xml_ = solverNode;
    this->infoNode_ = olasInfo;

    // No factorisation was performed yet
    this->readyToUse_ = false;

    // We will set this in the Setup() method
    // right before the factorisation
    this->memGrowthEstimate_ = 1;

    // Initialise attributes
    maxFill_ = 0;
    tau_ = 0.0;

  }


  // ==============
  //   Destructor
  // ==============
  template <typename T>
  ILUTP_Precond<T>::~ILUTP_Precond() {


  }


  // =================================
  //   Application of Preconditioner
  // =================================
  template <typename T>
  void ILUTP_Precond<T>::Apply( const CRS_Matrix<T> &sysMat,
				const Vector<T> &res, Vector<T> &sol ) {


    // Test that a factorisation is available, if not issue an error.
    if ( amFactorised_ == false ) {
      EXCEPTION( "ILUTP_Precond::Apply: No factorisation available. "
          << "Call Setup() first!" );
    }

    // Solve the problem
    CroutLU<T>::Solve( res, sol );
  }


  // ===========================
  //   Setup of Preconditioner
  // ===========================
  template <typename T>
  void ILUTP_Precond<T>::Setup( CRS_Matrix<T> &sysMat,
                                PtrParamNode analysis_id ) {


    // Query parameter object for factorisation parameters
    tau_ = 0.01;
    Integer aux = -2;
    PtrParamNode pNode = this->xml_->Get("ILUTP", ParamNode::INSERT);
    pNode->GetValue("tau", tau_, ParamNode::INSERT);
    pNode->GetValue("aux", aux, ParamNode::INSERT);
    
    if ( aux >= 0 ) {
      maxFill_ = (UInt)aux;
    }
    else {
      maxFill_ = (-aux) * ( sysMat.GetNnz() / sysMat.GetNumCols() - 1 );
    }

    // Report parameters to standard log stream
    (*cla) << "ILUTP_Precond: Performing an ILU( " << tau_ << " , " << maxFill_
      	   << " ) factorisation" << std::endl;

    // Perform the factorisation
    this->Factorise( sysMat );
    amFactorised_ = true;

    // If the user wishes, we can export the LU factorisation to a file
    std::string saveFacFile = "crout_fac.out";
    if(pNode->Has("saveFacFile")) {
      pNode->GetValue("saveFacFile", saveFacFile, ParamNode::INSERT);
      this->ExportILUFactorisation( saveFacFile.c_str() );
    }

  }


  // ===============
  //   DropEntries
  // ===============
  template <typename T>
  void ILUTP_Precond<T>::DropEntries( UInt k, std::vector<T> &vecZ,
				      std::vector<UInt> &vecZFill,
				      std::vector<T> &vecW,
				      std::vector<UInt> &vecWFill ) {


    EXCEPTION( "ILUTP_Precond currently not operational" );
    // COMPWARNING: unused variable UInt i, j;
    // COMPWARNING: unused variable Double curTau;
    // COMPWARNING: unused variable Double norm = 0.0;

    // ****************************
    //   Drop entries in row of U
    // ****************************

    // Remove index of diagonal entry from the set, since this one is always
    // kept and does not enter norm computation
    // vecZFill.clear( k );

    // PrintContainer( vecZFill, "vecZFill" );

    // Compute 1-norm of row
    // norm = 0.0;
    //     if ( !vecZFill.empty() ) {
    //       for ( std::vector<UInt>::iterator it = vecZFill.begin();
    // 	    it != vecZFill.end(); it++ ) {
    // 	norm += abs( vecZ[ *it ] );
    //       }
    // }

    // Determine current threshold
    //     curTau = tau_ * norm;

    //     UInt numUDropTau = 0;

    // Drop all entries smaller than threshold
    //     for ( std::vector<UInt>::iterator it = vecZFill.begin();
    // 	  it != vecZFill.end(); it++ ) {
    //       if ( abs(vecZ[ *it ]) < curTau ) {

    // 	// Order is important, since erasing *it
    // 	// invalidates the iterator it!
    // 	vecZ[ *it ] = 0.0;
    // 	vecZFill.erase( it );
    // 	numUDropTau++;
    //       }
    //     }

    //     std::cerr << " Row " << k << ": Dropped " << numUDropTau << " entries < "
    // 	      << curTau << std::endl;

    // If there are more than maxFill_ entries left,
    // we must find the maxFill_ largest ones (plus
    // the diagonal entry)
    //     if ( vecZFill.size() > maxFill_ + 1 ) {
    
    //       std::cerr << " --> Searching " << maxFill_ << "largest entries left"
    // 		<< std::endl;

    // Step 1: Remove index of diagonal entry from the set,
    //         since this one is always keep
    // vecZFill.erase( k );

    // Step 2: Copy all indices in auxilliary vector (omit diagonal)
    //       indexVec_.clear();
    //       for ( std::set<UInt>::iterator it = vecZFill.begin();
    // 	    it != vecZFill.end(); it++ ) {
    // 	indexVec_.push_back( *it );
    //       }
    
    // Step 3: Sort the index vector, so that the maxFill_ largest
    //         entries are up front
    //       nth_element( indexVec_.begin(),
    // 		   indexVec_.begin() + maxFill_,
    // 		   indexVec_.end(),
    // 		   ILUTP_Precond< T >::FindMaxEntries( vecZ ) );
    
    // Step 4: Eliminate all drop entries from the index set and nullify
    //         them in the dense vector
    //       for ( std::vector<UInt>::iterator it = indexVec_.begin() + maxFill_;
    // 	    it != indexVec_.end(); it++ ) {
    // 	vecZ[ *it ] = 0.0;
    // 	vecZFill.erase( *it );
    //       }

    // Step 5: Re-add index of diagonal entry to set
    // vecZFill.insert( k );
    // }

    // Re-add index of diagonal entry to set
    //     vecZFill.insert( k );


    // ****************************
    //   Drop entries in col of L
    // ****************************

    //     PrintContainer( vecWFill, "vecWFill" );

    // Compute 1-norm of column
    //     norm = 0.0;
    //     for ( i = 1; i <= sysMatDim_; i++ ) {
    //       norm += abs( vecW[i] );
    //     }
    //     norm /= (Double)sysMatDim_;

    // Determine current threshold
    //     curTau = tau_ * norm;

    //     UInt numLDropTau = 0;

    // Drop all entries smaller than threshold
    //     for ( std::set<UInt>::iterator it = vecWFill.begin();
    // 	  it != vecWFill.end(); it++ ) {
    //       if ( abs(vecW[ *it ]) < curTau ) {

    // 	// Order is important, since erasing *it
    // 	// invalidates the iterator it!
    // 	vecW[ *it ] = 0.0;
    // 	vecWFill.erase( it );
    // 	numLDropTau++;
    //       }
    //     }

    //     std::cerr << " Col " << k << ": Dropped " << numLDropTau << " entries < "
    // 	      << curTau << std::endl;

    // If there are more than maxFill_ entries left,
    // we must find the maxFill_ largest ones (plus
    // the diagonal entry)
    //     if ( vecWFill.size() > maxFill_ + 1 ) {

    // Step 1: Remove index of diagonal entry from the set,
    //         since this one is always keep
    // vecWFill.erase( k );

    // Step 2: Copy all indices in auxilliary vector (omit diagonal)
    //       indexVec_.clear();
    //       for ( std::set<UInt>::iterator it = vecWFill.begin();
    // 	    it != vecWFill.end(); it++ ) {
    // 	indexVec_.push_back( *it );
    //       }

    // Step 3: Sort the index vector, so that the maxFill_ largest
    //         entries are up front
    //       nth_element( indexVec_.begin(),
    // 		   indexVec_.begin() + maxFill_,
    // 		   indexVec_.end(),
    // 		   ILUTP_Precond< T >::FindMaxEntries( vecW ) );

    // Step 4: Eliminate all drop entries from the index set and nullify
    //         them in the dense vector
    //       for ( std::vector<UInt>::iterator it = indexVec_.begin() + maxFill_;
    // 	    it != indexVec_.end(); it++ ) {
    // 	vecW[ *it ] = 0.0;
    // 	vecWFill.erase( *it );
    //       }

    // Step 5: Re-add index of diagonal entry to set
    // vecWFill.insert( k );

  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class ILUTP_Precond<Double>;
  template class ILUTP_Precond<Complex>;
#endif
  
}
