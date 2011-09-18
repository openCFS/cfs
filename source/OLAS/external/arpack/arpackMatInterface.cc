// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>

#include "MatVec/stdmatrix.hh"
#include "MatVec/vector.hh"
#include "MatVec/generatematvec.hh"
#include "OLAS/precond/baseprecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/basesolver.hh"

#include "arpackMatInterface.hh"

namespace CoupledField {


  ArpackMatInterface::ArpackMatInterface( const BaseMatrix * matA, 
                                          bool shiftMode, Double shift) {

    matrixA_ = matA;
    matrixB_ = NULL;
    shift_ = pow(shift*8.0*atan(1.0),2);
    size_ = matA->GetNumRows();
    isGeneralized_ = false;
    shiftAndInvert_ = shiftMode;
  }

  ArpackMatInterface::ArpackMatInterface( const BaseMatrix * matA, 
                                          const BaseMatrix * matB, 
                                          bool shiftMode, Double shift) {
    
    matrixA_ = matA;
    matrixB_ = matB;
    shift_ = pow(shift*8.0*atan(1.0),2);
    size_ = matA->GetNumRows();
    isGeneralized_ = true;
    shiftAndInvert_ = shiftMode;
  }
    
  
  ArpackMatInterface::~ArpackMatInterface() {

  }
  
  
  void ArpackMatInterface::Setup( BaseSolver* solver, BasePrecond* precond ) {

    // Copy references
    solver_ = solver;
    precond_ = precond;
    solver_->SetPrecond(precond_);

    // Note: At this point I am not really sure, if we have to copy the
    // matrix into a new one
    // Copy matrix b to matrix c

    if( isGeneralized_ ) {

      // Depending on calculation mode (regular / shift and invert)
      // copy the correct matrix into matrixC_
      if ( shiftAndInvert_ == false ) {
        const StdMatrix & matB = dynamic_cast<const StdMatrix &>(*matrixB_);
        matrixC_ = CopyStdMatrixObject( matB );
      } else {
        const StdMatrix & matA = dynamic_cast<const StdMatrix &>(*matrixA_);
        const StdMatrix & matB = dynamic_cast<const StdMatrix &>(*matrixB_);
        matrixC_ = CopyStdMatrixObject( matB );

        matrixC_->Scale( -shift_ );
        matrixC_->Add( 1.0, matA );
      } 
    } else {
      // Standard EV problem
      if ( shiftAndInvert_ == false ) {
      EXCEPTION("Non-Shift-and-Invert mode not implemented for standard eigenvalue problem");
      } else {
        const StdMatrix & matA = dynamic_cast<const StdMatrix &>(*matrixA_);
        matrixC_ = CopyStdMatrixObject( matA );
        StdMatrix & matC = dynamic_cast<StdMatrix &>(*matrixC_);
        
        // calculate C = ( A - shift * I) with I being the unit matrix
        Double diag = 0.0;
        for( UInt i = 0; i < matC.GetNumRows(); i++) {
          matC.GetDiagEntry( i, diag );
          diag = (diag - shift_ * 1.0 );
          matC.SetDiagEntry( i, diag );
        }
      }
    }
    //matrixC_->Export("cmat.mat");
      
    // Setup solver and precond-object
    precond_->Setup( *matrixC_, shared_ptr<ParamNode>() );
    solver_->Setup( *matrixC_, shared_ptr<ParamNode>() );
  }
    
  void ArpackMatInterface::MultShiftOpV( Double* x, Double* y ) {

    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
//     x1--;
//     y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Solve system 
    solver_->Solve( *matrixC_, vecX, vecY );

  }

  void ArpackMatInterface::MultRegularOpV( Double* x, Double* y ) {

    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
//     x1--;
//     y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Create temporary vector
    Vector<Double> ax(size_);
    matrixA_->Mult( vecX, ax );

    // Solve system 
    solver_->Solve( *matrixC_,   ax, vecY );

  }
  
  void ArpackMatInterface::MultBV( Double* x, Double* y ) {
    
    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
//     x1--;
//     y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Perform matrix-vector multiplication
    if( isGeneralized_ ) {
      matrixB_->Mult( vecX, vecY );
    } else {
      vecY = vecX;
    }
    
  }
  
  void ArpackMatInterface::MultAV( Double* x, Double* y ) {
    
    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
//     x1--;
//     y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Perform matrix-vector multiplication
    matrixA_->Mult( vecX, vecY );
    
  }
  
}
