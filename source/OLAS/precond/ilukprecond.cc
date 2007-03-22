// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <algorithm>
#include <iterator>

#include "precond/ilukprecond.hh"

// Include source code of CroutLU class for template instantiation
// Note: Might lead to double instantiation, since CroutLU is also
// used in LUSolver. Going to implement better concept as soon as
// time permits.
#include "utils/math/croutlu.cc"

namespace OLAS {


  // *****************************************************
  //   Constructor (for use in GenerateStdPrecondObject)
  // *****************************************************
  template <typename T>
  ILUK_Precond<T>::ILUK_Precond( const StdMatrix& stdMat, 
                                 OLAS_Params *myParams,
                                 OLAS_Report *myReport ) {

    ENTER_FCN( "ILUK_Precond::ILUK_Precond" );

    // Set pointers to communication objects
    this->myParams_ = myParams;
    this->myReport_ = myReport;

    // No factorisation was performed yet
    this->readyToUse_ = false;

    // No fill level known yet
    maxLevel_ = 0;

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

    ENTER_FCN( "ILUK_Precond::~ILUK_Precond" );

  }


  // *********************************
  //   Application of Preconditioner
  // *********************************
  template <typename T>
  void ILUK_Precond<T>::Apply( const CRS_Matrix<T> &sysMat,
				const Vector<T> &res, Vector<T> &sol ) const {

    ENTER_FCN( "ILUK_Precond::Apply" );

    // Test that a factorisation is available, if not issue an error.

    if ( this->readyToUse_ == false ) {
      (*error) << "ILUK_Precond::Apply: No factorisation available. "
	       << "Call Setup() first!";
      Error( __FILE__, __LINE__ );
    }

    // Solve the problem
    CroutLU<T>::Solve( res, sol );

  }


  // ***************************
  //   Setup of Preconditioner
  // ***************************
  template <typename T>
  void ILUK_Precond<T>::Setup( CRS_Matrix<T> &sysMat ) {

    ENTER_FCN( "ILUK_Precond::Setup" );

    bool logging = this->myParams_->GetBoolValue( "ILUK_logging" );

    // Query parameter object for factorisation parameter
    maxLevel_ = this->myParams_->GetIntValue( "ILUK_level" );

    // Obtain and check dimensions of matrix
    this->sysMatDim_ = sysMat.GetNcols();
    if ( this->sysMatDim_ != sysMat.GetNrows() ) {
      (*error) << "ILUK_Precond: Input matrix is "
	       << sysMat.GetNrows() << " x " << sysMat.GetNcols()
	       << ", but needs to be square";
      Error( __FILE__, __LINE__ );
    }

    // Report parameters to standard log stream
    if ( logging == true ) {
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
    Factorise( sysMat );
    this->readyToUse_ = true;

    // If the user wishes, we can export the LU factorisation to a file
    if ( this->myParams_->GetBoolValue( "CROUT_saveFacToFile" ) ) {
      std::string filename;
      filename = this->myParams_->GetStringValue( "CROUT_facFileName" );
      this->ExportILUFactorisation( filename.c_str() );
      if ( logging == true ) {
	(*cla) << " Exported factor matrix to file '" << filename << "'"
	       << std::endl;
      }
    }

    // Close log section
    if ( logging == true ) {
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

    ENTER_IFCN( "ILUK_Precond::DropEntries" );

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


  // ************************
  //   Forced Instantiation
  // ************************
  template<typename T>
  void ILUK_Precond<T>::
  InstantiateAdditionalPublicMethods( BaseMatrix &sysMat ) {
    ENTER_FCN( "ILUK_Precond::InstantiateAdditionalPublicMethods" );
    this->ExportILUFactorisation( "dummy.mtx" );
  }

}
