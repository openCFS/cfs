// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>

#include "matvec/matvec.hh"
#include "precond/precond.hh"
#include "solver/solver.hh"

#include "arpackMatInterface.hh"

namespace OLAS {


  ArpackMatInterface::ArpackMatInterface( const BaseMatrix * matA, 
                                          const BaseMatrix * matB, 
                                          bool shiftMode, Double shift) {
    ENTER_FCN( "ArpackMatInterface::ArpackMatInterface" );
    
    matrixA_ = matA;
    matrixB_ = matB;
    shift_ = pow(shift*8.0*atan(1.0),2);
    size_ = matA->GetNrows();
    shiftAndInvert_ = shiftMode;
  }
    
  
  ArpackMatInterface::~ArpackMatInterface() {
    ENTER_FCN( "ArpackMatInterface::~ArpackMatInterface" );

  }
  
  
  void ArpackMatInterface::Setup( BaseSolver* solver, BasePrecond* precond ) {
    ENTER_FCN( "ArpackMatInterface::Setup" );

    // Copy references
    solver_ = solver;
    precond_ = precond;

    // Note: At this point I am not really sure, if we have to copy the
    // matrix into a new one
    // Copy matrix b to matrix c

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
    //matrixC_->Export("cmat.mat");
      
    // Setup solver and precond-object
    precond_->Setup( *matrixC_ );
    solver_->Setup( *matrixC_ );
  }
    
  void ArpackMatInterface::MultShiftOpV( Double* x, Double* y ) {
    ENTER_FCN( "ArpackMatInterface::MultShiftOpV" );

    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
    x1--;
    y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Solve system 
    solver_->Solve( *matrixC_, *precond_, vecX, vecY );

  }

  void ArpackMatInterface::MultRegularOpV( Double* x, Double* y ) {
    ENTER_FCN( "ArpackMatInterface::MultRegularOpV" );

    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
    x1--;
    y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Create temporary vector
    Vector<Double> ax(size_);
    matrixA_->Mult( vecX, ax );

    // Solve system 
    solver_->Solve( *matrixC_, *precond_, ax, vecY );

  }
  
  void ArpackMatInterface::MultBV( Double* x, Double* y ) {
    ENTER_FCN( "ArpackMatInterface::MultBV" );
    
    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
    x1--;
    y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Perform matrix-vector multiplication
    matrixB_->Mult( vecX, vecY );
    
  }
  
  void ArpackMatInterface::MultAV( Double* x, Double* y ) {
    ENTER_FCN( "ArpackMatInterface::MultAV" );
    
    // Create two temporary vectors as wrappers for x and y
    Double *x1 = x;
    Double *y1 = y;
    x1--;
    y1--;
    Vector<Double> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Perform matrix-vector multiplication
    matrixA_->Mult( vecX, vecY );
    
  }
  
}
