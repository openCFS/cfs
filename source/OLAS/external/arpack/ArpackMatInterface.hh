#ifndef OLAS_ARPACKMATINTERFACE_HH
#define OLAS_ARPACKMATINTERFACE_HH

#include <boost/shared_ptr.hpp>
#include "General/defs.hh"

namespace CoupledField {

  class BaseSolver;
  class BasePrecond;
  class BaseMatrix;
  class Timer;
  
  //! Helper class for defining an interface to our matrices for the case of
  //! the "shift-and-invert" mode of arpack
  class ArpackMatInterface {
    
  public:
    
    //! Computational mode / problem transformation used for solving
    //! see e.g. https://www.caam.rice.edu/software/ARPACK/UG/node39.html
    typedef enum {REGULAR, SHIFT_INVERT, BUCKLING} ComputeMode;

    //! Constructor (standard EV problem)
    ArpackMatInterface( const BaseMatrix * matA, ComputeMode computeMode);
    
    //! Constructor (generalized EV problem)
    ArpackMatInterface( const BaseMatrix * matA, const BaseMatrix * matB, ComputeMode computeMode);
    ArpackMatInterface( BaseMatrix * matA, BaseMatrix * matB, BaseMatrix * matD, ComputeMode computeMode);

    //! Destructor
    ~ArpackMatInterface();
    
    //! Setup the internal methods for solving. setupTimer and solveTimer are the main solver times and used here as parents
    template <class TYPE>
    void Setup( BaseSolver* solver, BasePrecond* precond, TYPE shift, boost::shared_ptr<Timer> parentSetupTimer, boost::shared_ptr<Timer> parentSolveTimer );

    //! Setup the internal methods for solving, quadratic case
    template <class TYPE>
    void QuadSetup( BaseSolver* solver, BasePrecond* precond, TYPE shift);


    //! Method to set diagonal scaling factor
    void SetDiagScaling( Double scaleFac ) { diagScale_ = scaleFac; };

    //! Method to get diagonal scaling factor
    Double GetDiagScaling( ) { return diagScale_; };


    //! use stiffnes matrix in nonsingular N matrix in quadratic solver
    bool UseStiffInNMat ( ) { return useStiffInNMat_; } ;

    //! Performs the operation ( A - sigma * B ) y = x
    template <class TYPE>
    void MultShiftOpV(TYPE* x, TYPE* y);

    //! Performs the operation ( B * y = A * x )
    template <class TYPE>
    void MultRegularOpV(TYPE* x, TYPE* y);
    
    //! Performs the operation y=A*x
    template <class TYPE>
    void MultAV(TYPE* x, TYPE* y);

    //! Performs the operation y=B*x
    template <class TYPE>
    void MultBV(TYPE* x, TYPE* y);


    // TODO Bloch rename complex methods to Quad and implement complex methods (with old name) for bloch

    //! Performs the operation ( A - sigma * B ) y = x
    void MultShiftOpVQuad( Complex* x, Complex* y );

    //! Performs the operation ( B * y = A * x )
    void MultRegularOpVQuad( Complex* x, Complex* y );
    
    //! Performs the operation y=A*x
    void MultAVQuad( Complex* x, Complex* y );

    //! Performs the operation y=B*x
    void MultBVQuad( Complex* x, Complex* y );
    
    //! Perform the operation y=complex(Stiff|Mass|Damp)*x
    void MultZStiffV( Complex* x, Complex* y );
    void MultZMassV( Complex* x, Complex* y );
    void MultZDampV( Complex* x, Complex* y );
    void ScaleDiag( Complex* x );
    void InvScaleDiag( Complex* x );
    
  private:
    
    //! size of rows / columns
    UInt size_;
    bool complex_shift_ = false;
    
    //! Value of eigenvalue shift
    Complex shift_ = Complex(0.0,0.0);
    
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

    //! Computational mode / problem transformation used for solving
    ComputeMode computeMode_;

    //! Flag for generalized eigenvalue problem
    bool isGeneralized_;

    //! Flag for shift-and-invert mode
    bool useStiffInNMat_;
  };
}



#endif
