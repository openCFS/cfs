// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_ARPACKMATINTERFACE_HH
#define OLAS_ARPACKMATINTERFACE_HH

#include "General/defs.hh"

namespace CoupledField {

  class BaseSolver;
  class BasePrecond;
  class BaseMatrix;
  
  //! Helper class for defining an interface to our matrices for the case of
  //! the "shift-and-invert" mode of arpack
  class ArpackMatInterface {
    
  public:
    
    //! Constructor
    ArpackMatInterface( const BaseMatrix * matA, const BaseMatrix * matB,
                        bool shiftMode, Double shift );
    
    //! Destructor
    ~ArpackMatInterface();
    
    //! Setup the internal methods vor solving
    void Setup( BaseSolver* solver, BasePrecond* precond );
    
    //! Performs the operation ( A - sigma * B ) y = x
    void MultShiftOpV( Double* x, Double* y );

    //! Performs the operation ( B * y = A * x )
    void MultRegularOpV( Double* x, Double* y );
    
    //! Performs the operation y=A*x
    void MultAV( Double* x, Double* y );

    //! Performs the operation y=B*x
    void MultBV( Double* x, Double* y );
    
    
  private:
    
    //! size of rows / columns
    UInt size_;
    
    //! Value of frequency shift
    Double shift_;
    
    //! Pointer to matrix which gets interfaced
    const BaseMatrix * matrixA_, * matrixB_;
    
    //! Pointer to matrix (A - sigma * B);
    BaseMatrix * matrixC_;
    
    //! Pointer to solver
    BaseSolver * solver_;
    
    //! Pointer to preconditionier
    BasePrecond * precond_;

    //! Flag for shift-and-invert mode
    bool shiftAndInvert_;
    
  };
}

#endif
