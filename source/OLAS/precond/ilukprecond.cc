// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <algorithm>
#include <iterator>

#include "MatVec/crs_matrix.hh"

#include "OLAS/precond/ilukprecond.hh"

// Include source code of CroutLU class for template instantiation
// Note: Might lead to double instantiation, since CroutLU is also
// used in LUSolver. Going to implement better concept as soon as
// time permits.
#include "OLAS/utils/math/croutlu.hh"

namespace CoupledField {


  // *****************************************************
  //   Constructor (for use in GenerateStdPrecondObject)
  // *****************************************************
  template <typename T>
  ILUK_Precond<T>::ILUK_Precond( const StdMatrix& stdMat, 
                                 PtrParamNode solverNode,
                                 PtrParamNode olasInfo ) {

    // Set pointers to communication objects
    this->xml_ = solverNode;
    this->infoNode_ = olasInfo;

    // No factorisation was performed yet
    this->readyToUse_ = false;

    // No fill level known yet
    maxLevel_ = 0;
    
    // Deactivate logging by default
    logging_ = false;

    // Problem dimension not known yet
    this->sysMatDim_ = 0;

    // We will set this in the Setup() method
    // right before the factorisation
    this->memGrowthEstimate_ = 1;

    // Must set this to true in order to force CroutLU
    // to compute the fill-level information
    this->compFillLevels_ = true;

  }


  // **************
  //   Destructor
  // **************
  template <typename T>
  ILUK_Precond<T>::~ILUK_Precond() {


  }


  // *********************************
  //   Application of Preconditioner
  // *********************************
  template <typename T>
  void ILUK_Precond<T>::Apply( const CRS_Matrix<T> &sysMat,
				const Vector<T> &res, Vector<T> &sol ) {


    // Test that a factorisation is available, if not issue an error.

    if ( this->readyToUse_ == false ) {
      EXCEPTION( "ILUK_Precond::Apply: No factorisation available. "
	       << "Call Setup() first!" );
    }

    // Solve the problem
    CroutLU<T>::Solve( res, sol );

  }


  // ***************************
  //   Setup of Preconditioner
  // ***************************
  template <typename T>
  void ILUK_Precond<T>::Setup( CRS_Matrix<T> &sysMat,
                               PtrParamNode infoNode ) {



    // Query parameter object for factorisation parameter
    maxLevel_ = 1;
    PtrParamNode pNode = this->xml_->Get("ILUK", ParamNode::INSERT );
    pNode->GetValue("level", maxLevel_, ParamNode::INSERT );
    pNode->GetValue("logging", logging_, ParamNode::INSERT ) ;
    
    // Obtain and check dimensions of matrix
    this->sysMatDim_ = sysMat.GetNumCols();
    if ( this->sysMatDim_ != sysMat.GetNumRows() ) {
      EXCEPTION( "ILUK_Precond: Input matrix is "
          << sysMat.GetNumRows() << " x " << sysMat.GetNumCols()
          << ", but needs to be square" );
    }

    // Report parameters to standard log stream
    if ( logging_ == true ) {
      (*cla) << " -----------------------------------------------\n"
	     << " ILUK_Precond: Performing an ILU( " << maxLevel_
	     << " ) factorisation\n of a "
	     << this->sysMatDim_ << " x " << this->sysMatDim_
	     << " matrix (nnz = " << sysMat.GetNnz() << ")"
	     << std::endl;
    }

    // Try to figure out an approximate factor for memory
    // pre-allocation
    if ( maxLevel_ <= 3 ) {
      this->memGrowthEstimate_ = 2;
    }
    else {
      this->memGrowthEstimate_ = maxLevel_ / 2;
    }

    // Perform the factorisation
    this->Factorise( sysMat );
    this->readyToUse_ = true;

    // If the user wishes, we can export the LU factorisation to a file
    std::string saveFacFile = "crout_fac.out";
    if(pNode->Has("saveFacFile")) {
      pNode->GetValue("saveFacFile", saveFacFile, ParamNode::INSERT);

      this->ExportILUFactorisation( saveFacFile.c_str() );
      if ( logging_ == true ) {
	(*cla) << " Exported factor matrix to file '" << saveFacFile << "'"
	       << std::endl;
      }
    }

    // Close log section
    if ( logging_ == true ) {
      (*cla) << " -----------------------------------------------"
	     << std::endl;
    }

  }


  // ***************
  //   DropEntries
  // ***************
  template <typename T>
  void ILUK_Precond<T>::DropEntries( UInt k, std::vector<T> &vecZ,
				      std::vector<UInt> &vecZFill,
				      std::vector<T> &vecW,
				      std::vector<UInt> &vecWFill ) {


    std::vector<UInt>::iterator it;

#ifdef DEBUG_ILUKPRECOND
    (*debug) << " -> length vecZFill = " << vecZFill.size() << '\n'
	     << " -> length vecWFill = " << vecWFill.size() << std::endl;
    PrintVector( vecZLevel_, "vecZLevel_" );
    PrintVector( vecWLevel_, "vecWLevel_" );
#endif

    // Drop entries in row of U
    it = std::remove_if( vecZFill.begin(), vecZFill.end(),
			 TestFillLevel( maxLevel_, this->vecZLevel_,
                                        vecZ, this->sysMatDim_ + 1 ) );
    vecZFill.erase( it, vecZFill.end() );
    
    // Drop entries in column of L
    it = std::remove_if( vecWFill.begin(), vecWFill.end(),
			 TestFillLevel( maxLevel_, this->vecWLevel_,
                                        vecW, this->sysMatDim_ + 1 ) );
    vecWFill.erase( it, vecWFill.end() );

#ifdef DEBUG_ILUKPRECOND
    (*debug) << " -> length vecZFill = " << vecZFill.size() << '\n'
	     << " -> length vecWFill = " << vecWFill.size() << std::endl;
#endif

  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class ILUK_Precond<Double>;
  template class ILUK_Precond<Complex>;
#endif
  
}
