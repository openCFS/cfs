// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// #include <string>
// #include <iostream>
// #include <fstream>

#include "OLAS/precond/baseprecond.hh"
#include "OLAS/precond/sbmdiagprecond.hh"

namespace CoupledField {

  void SBMDiagPrecond::Setup( SBM_Matrix &mat ) {


    // Obtain size information from matrix
    Ncols_ = mat.GetNumCols();
    Nrows_ = mat.GetNumRows();

    // Check that matrix is square
    if ( Ncols_ != Nrows_ ) {
      EXCEPTION("Expected a square SBMMatrix");
    }

    // Construct the individual preconditioners
    if ( precond_ != NULL ) {
      delete[] precond_;
    }
    precond_ = new BaseStdPrecond[ Ncols_ ];
    for ( int i = 1; i <= Ncols_; i++ ) {
      // precond_[i] = Call to Preconditioner creation routine???
      // Get types of preconditioners from parameter object???
    }
  
    // Trigger setup of the individual preconditioners for the diagonal sub-
    // matrices
    for ( int i = 1; i <= Ncols_; i++ ) {
      precond_[i].Setup( mat(i,i) );
    }

    // Now the preconditioner is ready to use
    readyToUse_ = true;

  }


  void SBMDiagPrecond::Apply( const SBM_Matrix mat, const SBM_Vector r,
			      SBM_Vector z ) const {


    for ( int i = 1; i <= Ncols_; i++ ) {
      precond_[i].Setup( mat(i,i) );
    }



  };

}
