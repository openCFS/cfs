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
    
    //! Constructor (standard EV problem)
    ArpackMatInterface( const BaseMatrix * matA,
                        bool shiftMode, Double shift );
    
    //! Constructor (generalized EV problem)
    ArpackMatInterface( const BaseMatrix * matA, const BaseMatrix * matB,
                        bool shiftMode, Double shift );
    ArpackMatInterface( BaseMatrix * matA, BaseMatrix * matB,
                        BaseMatrix * matD, bool shiftMode, Double shift );    
    //! Destructor
    ~ArpackMatInterface();
    
    //! Setup the internal methods for solving
    void Setup( BaseSolver* solver, BasePrecond* precond );
    
    //! Setup the internal methods for solving, quadratic case
    void QuadSetup( BaseSolver* solver, BasePrecond* precond );

    //! Method to set diagonal scaling factor
    void SetDiagScaling( Double scaleFac );

    //! Method to get diagonal scaling factor
    Double GetDiagScaling( );


    //! use stiffnes matrix in nonsingular N matrix in quadratic solver
    bool UseStiffInNMat ( ) { return useStiffInNMat_; } ;

    //! Performs the operation ( A - sigma * B ) y = x
    void MultShiftOpV( Double* x, Double* y );

    //! Performs the operation ( B * y = A * x )
    void MultRegularOpV( Double* x, Double* y );
    
    //! Performs the operation y=A*x
    void MultAV( Double* x, Double* y );

    //! Performs the operation y=B*x
    void MultBV( Double* x, Double* y );

    //! Performs the operation ( A - sigma * B ) y = x
    void MultShiftOpV( Complex* x, Complex* y );

    //! Performs the operation ( B * y = A * x )
    void MultRegularOpV( Complex* x, Complex* y );
    
    //! Performs the operation y=A*x
    void MultAV( Complex* x, Complex* y );

    //! Performs the operation y=B*x
    void MultBV( Complex* x, Complex* y );
    
    //! Perform the operation y=complex(Stiff|Mass|Damp)*x
    void MultZStiffV( Complex* x, Complex* y );
    void MultZMassV( Complex* x, Complex* y );
    void MultZDampV( Complex* x, Complex* y );
    void ScaleDiag( Complex* x );
    void InvScaleDiag( Complex* x );
    
  private:
    
    //! size of rows / columns
    UInt size_;
    
    //! Value of frequency shift
    Double shift_;
    
    //! Pointer to matrix which gets interfaced
    const BaseMatrix * matrixA_, * matrixB_, *matrixD_;
    
    //! Pointer to matrix (A - sigma * B);
    BaseMatrix * matrixC_;

    //! scaling factor for N matrix in quad eigenvalue problem
    Double diagScale_;
    
    //! Pointer to solver
    BaseSolver * solver_;
    
    //! Pointer to preconditionier
    BasePrecond * precond_;

    //! Flag for shift-and-invert mode
    bool shiftAndInvert_;
    
    //! Flag for generalized eigenvalue problem
    bool isGeneralized_;

    //! Flag for shift-and-invert mode
    bool useStiffInNMat_;
  };
}

#endif
