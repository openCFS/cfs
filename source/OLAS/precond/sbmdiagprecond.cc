// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// #include <string>
// #include <iostream>
// #include <fstream>

#include "matvec/matvec.hh"
#include "precond/baseprecond.hh"
#include "precond/sbmdiagprecond.hh"

namespace OLAS {

  void SBMDiagPrecond::Setup( SBM_Matrix &mat ) {

    ENTER_FCN( "SBMDiagPrecond::Setup" );

    // Obtain size information from matrix
    Ncols_ = mat.GetNcols();
    Nrows_ = mat.GetNrows();

    // Check that matrix is square
    if ( Ncols_ != Nrows_ ) {
      Error( "Expected a square SBMMatrix", __FILE__, __LINE__ );
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
    readyToUse_ = TRUE;

  }


  void SBMDiagPrecond::Apply( const SBM_Matrix mat, const SBM_Vector r,
			      SBM_Vector z ) const {

    ENTER_FCN( "SBMDiagPrecond::Apply" );

    for ( int i = 1; i <= Ncols_; i++ ) {
      precond_[i].Setup( mat(i,i) );
    }



  };

}
