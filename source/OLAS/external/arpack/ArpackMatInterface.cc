#include <iostream>

#include "MatVec/StdMatrix.hh"
#include "MatVec/Vector.hh"
#include "MatVec/generatematvec.hh"
#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "Utils/Timer.hh"
#include "ArpackMatInterface.hh"

namespace CoupledField {


  ArpackMatInterface::ArpackMatInterface( const BaseMatrix * matA, ComputeMode computeMode) {

    matrixA_ = matA;
    matrixB_ = NULL;
    matrixC_ = NULL;
    matrixD_ = NULL;
    
    shift_ = 0;
    size_ = matA->GetNumRows();
    isGeneralized_ = false;
    computeMode_ = computeMode;
    precond_ = NULL;
    solver_ = NULL;
    diagScale_ = 1.0;
    useStiffInNMat_ = false;
  }

  ArpackMatInterface::ArpackMatInterface( const BaseMatrix * matA, const BaseMatrix * matB, ComputeMode computeMode) {
    
    matrixA_ = matA;
    matrixB_ = matB;
    matrixC_ = NULL;
    matrixD_ = NULL;
    
    shift_ = 0;
    size_ = matA->GetNumRows();
    isGeneralized_ = true;
    computeMode_ = computeMode;
    precond_ = NULL;
    solver_ = NULL;
    diagScale_ = 1.0;
    useStiffInNMat_ = false;

  }

  ArpackMatInterface::ArpackMatInterface( BaseMatrix * matA, BaseMatrix * matB, BaseMatrix * matD, ComputeMode computeMode) {
    
    matrixA_ = matA;
    matrixB_ = matB;
    matrixC_ = NULL;
    matrixD_ = matD;
    
    shift_ = 0;
    size_    = 2*matA->GetNumRows();
    computeMode_ = computeMode;
    useStiffInNMat_ = false;
    isGeneralized_ = true;
    precond_ = NULL;
    solver_ = NULL;
    diagScale_ = 1.0;
  }

  ArpackMatInterface::~ArpackMatInterface() {
    if( matrixC_ ) {
      delete matrixC_;
      matrixC_ = NULL;
    }
    solver_  = NULL;
    precond_ = NULL;
  }
  
  template<class TYPE>
  void ArpackMatInterface::Setup( BaseSolver* solver, BasePrecond* precond, TYPE shift, shared_ptr<Timer> parentSetupTimer, shared_ptr<Timer> parentSolveTimer ) {

    complex_shift_ = boost::is_complex<TYPE>::value;
    // Copy references
    solver_ = solver;
    precond_ = precond;
    solver_->SetPrecond(precond_);
    shift_= (Complex) shift;

    // Setup() might be called multiple times for Bloch mode analysis
    if(matrixC_) {
      delete matrixC_; // ReRegister from PatternPool
      matrixC_ = NULL;
    }

    // Note: At this point I am not really sure, if we have to copy the
    // matrix into a new one
    // Copy matrix b to matrix c
    if( isGeneralized_ ) {
      // Depending on calculation mode (regular | shift and invert | buckling)
      // copy the correct matrix into matrixC_
      if ( computeMode_ == ComputeMode::REGULAR ) {
        const StdMatrix & matB = dynamic_cast<const StdMatrix &>(*matrixB_);
        matrixC_ = CopyStdMatrixObject( matB );
      } else {
        const StdMatrix & matA = dynamic_cast<const StdMatrix &>(*matrixA_);
        const StdMatrix & matB = dynamic_cast<const StdMatrix &>(*matrixB_);
        matrixC_ = CopyStdMatrixObject( matB );

        if(complex_shift_)
        	matrixC_->Scale( -shift_ );
        else
        	matrixC_->Scale( -shift_.real());
        matrixC_->Add( 1.0, matA );
      }
    } else {
      // Standard EV problem
      if ( computeMode_ == ComputeMode::REGULAR ) {
        EXCEPTION("Non-Shift-and-Invert mode not implemented for standard eigenvalue problem");
      } else {
        const StdMatrix & matA = dynamic_cast<const StdMatrix &>(*matrixA_);
        matrixC_ = CopyStdMatrixObject( matA );
        StdMatrix & matC = dynamic_cast<StdMatrix &>(*matrixC_);

        // calculate C = ( A - shift * I) with I being the unit matrix
        TYPE diag = 0.0;
        for( UInt i = 0; i < matC.GetNumRows(); i++) {
          matC.GetDiagEntry( i, diag );
          diag = (diag - shift * 1.0 );
          matC.SetDiagEntry( i, diag );
        }
      }
    }

    // Setup solver and precond-object
    precond_->GetSetupTimer()->SetSub(parentSetupTimer); // is in the service of arpack
    precond_->GetPrecondTimer()->SetSub(parentSolveTimer);
    precond_->GetSetupTimer()->Start();
    precond_->Setup( *matrixC_ );
    precond_->GetSetupTimer()->Stop();

    solver_->GetSetupTimer()->SetSub(parentSetupTimer);
    solver_->GetSolveTimer()->SetSub(parentSolveTimer);
    solver_->GetSetupTimer()->Start();
    solver_->Setup( *matrixC_ );
    solver_->GetSetupTimer()->Stop();
  }

  template<class TYPE>
  void ArpackMatInterface::QuadSetup( BaseSolver* solver, BasePrecond* precond, TYPE shift) {

    // Copy references
    solver_ = solver;
    precond_ = precond;
    complex_shift_ = boost::is_complex<TYPE>::value;
    if (matrixC_ == NULL || shift_ != shift) {
      shift_= (Complex) shift;

      const StdMatrix & matA = dynamic_cast<const StdMatrix &>(*matrixA_);
      const StdMatrix & matB = dynamic_cast<const StdMatrix &>(*matrixB_);
      const StdMatrix & matD = dynamic_cast<const StdMatrix &>(*matrixD_);

      // form (B*shift + D)*shift + A) = A + sigma*D + sigma**2*B
      if(complex_shift_)
      {
        matrixC_ = CopyStdMatrixObject( matB );
        matrixC_->Scale( shift_ );
        matrixC_->Add( 1.0, matD );
        matrixC_->Scale( shift_ );
        matrixC_->Add( 1.0, matA );
      }
      else
      {
        matrixC_ = CopyStdMatrixObject( matB );
        matrixC_->Scale( shift_.real() );
        matrixC_->Add( 1.0, matD );
        matrixC_->Scale( shift_.real() );
        matrixC_->Add( 1.0, matA );
      }

      // set diagonal scaling entry (hard coded = 1)
      diagScale_ = 1.0;

      // Setup solver and precond-object
      // we must do the whole timer management manually as we do not go via a non-overloaded BaseSolver function
      // For the standard FEA this is done in AlgebraicSys::SetupPrecond()
      precond_->GetSetupTimer()->Start();
      precond_->Setup( *matrixC_ );
      precond_->GetSetupTimer()->Stop();

      solver_->GetSetupTimer()->Start();
      solver_->Setup( *matrixC_ );
      solver_->GetSetupTimer()->Stop();
    }
  }

  template <class TYPE>
  void ArpackMatInterface::MultShiftOpV(TYPE* x, TYPE* y) {

    // Create two temporary vectors as wrappers for x and y
    TYPE* x1 = x;
    TYPE* y1 = y;
    Vector<TYPE> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Solve system 
    solver_->GetSolveTimer()->Start();
    solver_->Solve(*matrixC_, vecX, vecY);
    solver_->GetSolveTimer()->Stop();
  }

  template <class TYPE>
  void ArpackMatInterface::MultRegularOpV(TYPE* x, TYPE* y) {

    // Create two temporary vectors as wrappers for x and y
    TYPE* x1 = x;
    TYPE* y1 = y;
    Vector<TYPE> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Create temporary vector
    Vector<TYPE> ax(size_);
    matrixA_->Mult(vecX, ax);

    // Solve system
    solver_->GetSolveTimer()->Start();
    solver_->Solve(*matrixC_, ax, vecY);
    solver_->GetSolveTimer()->Stop();
  }
  
  template <class TYPE>
  void ArpackMatInterface::MultBV(TYPE* x, TYPE* y) {
    
    // Create two temporary vectors as wrappers for x and y
    TYPE* x1 = x;
    TYPE* y1 = y;
    Vector<TYPE> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Perform matrix-vector multiplication
    if( isGeneralized_ ) {
      matrixB_->Mult( vecX, vecY );
    } else {
      vecY = vecX;
    }
    
  }
  
  template <class TYPE>
  void ArpackMatInterface::MultAV(TYPE* x, TYPE* y ) {
    
    // Create two temporary vectors as wrappers for x and y
    TYPE* x1 = x;
    TYPE* y1 = y;
    Vector<TYPE> vecX, vecY;
    vecX.Replace( size_, x1, false );
    vecY.Replace( size_, y1, false );
    
    // Perform matrix-vector multiplication
    matrixA_->Mult( vecX, vecY );
  }
  
  
    
  void ArpackMatInterface::MultShiftOpVQuad( Complex* x, Complex* y ) {

    // do a solve operation A*y = x for the system matrix of the
    // quadratic eigenvalue problem
    // important notes: 
    // x and y are vectors of length size_, which is twice the 
    //       number of unknowns in this case x = (b1,b2), y=(a1,a2)
    // the system matrix is (as in MultAV) given by
    //   ( -C -K | N 0 )
    // where K and C are the system stiffness and damping matrices, 
    // and N is an arbitrary, nonsingular matrix. 
    // we use a multiple of the identiy here (or, alternatively, of K).
    // Solution is done in a 2 step procedure
    // solve
    //      a1 = sigma*a2 + N**(-1)*b2
    // and
    //     (K +sigma*C + sigma**2*M) a2 = -b1 - (C+sigma*M)*a1

    // Create two temporary vectors as wrappers for x and y
    Complex *b1 = x, *b2;
    Complex *a1 = y, *a2;
    Integer neq = size_/2;

    b2 = b1 + neq;
    a2 = a1 + neq;
    Vector<Complex> vecB1, vecA1, vecB2, vecA2;

    vecA1.Replace( neq, a1, false );
    vecB1.Replace( neq, b1, false );
    vecA2.Replace( neq, a2, false );
    vecB2.Replace( neq, b2, false );
    
    Vector<Complex> cx(neq),mx(neq),nInvB2(neq);

    // solve N*nInvB2 = b2
    nInvB2.Init();
    if ( useStiffInNMat_ ) {
      solver_->GetSolveTimer()->Start();
      solver_->Solve( *matrixC_, vecB2, nInvB2);
      solver_->GetSolveTimer()->Stop();
    } else {
      nInvB2.Add(vecB2);
    }
    nInvB2.ScalarMult(1./diagScale_);

    // upper row
    // mx = (C+sigma*M)*N**(-1)*b2
    matrixD_->Mult(nInvB2,cx);
    matrixB_->Mult(nInvB2,mx);
    if(complex_shift_)
    	mx.Axpy(shift_,cx);
    else
        mx.Axpy(shift_.real(),cx);
    vecB1.Add(mx);
    // full rhs
    vecB1.ScalarMult(-1.0);
    solver_->GetSolveTimer()->Start();
    solver_->Solve( *matrixC_, vecB1, vecA2 );
    solver_->GetSolveTimer()->Stop();
    // solve for a1 (lower part of equation system)
    if(complex_shift_)
    	vecA1.Add(shift_, vecA2, 1., nInvB2);
    else
    	vecA1.Add(shift_.real(), vecA2, 1., nInvB2);
  }

  void ArpackMatInterface::MultBVQuad( Complex* x, Complex* y ) {

    // do a Matrix times vector operation for the solution of the
    // quadratic eigenvalue problem
    // important notes: 
    // x and y are vectors of length size_, which is twice the 
    //       number of unknowns in this case
    // B is a complex matrix which is of the form
    //   ( M 0 | 0 N )
    // where M is the system "mass" matrix, and N is an arbitrary,
    // nonsingular matrix. we just take N to be a multiple of the identity
    // (or, alternatively K).

    // Create temporary vectors as wrappers for x and y
    // note: x2,y2, vecX2, vecY2 are only required for N = K
    Complex *x1 = x, *x2;
    Complex *y1 = y, *y2;

    x2 = x1 + size_/2;
    y2 = y1 + size_/2;

    Vector<Complex> vecX1, vecY1, vecX2, vecY2;
    vecX1.Replace( size_/2, x1, false );
    vecY1.Replace( size_/2, y1, false );
    vecX2.Replace( size_/2, x2, false );
    vecY2.Replace( size_/2, y2, false );

    // Perform matrix-vector multiplication
    matrixB_->Mult( vecX1, vecY1 );

    //for (Integer i=size_/2; i < size_; i++) {
    //  y[i] = diagScale_*x[i];
    //}
    // Perform matrix-vector multiplication
    vecY2.Init();
    if ( useStiffInNMat_ ) {
      matrixC_->Mult( vecX2, vecY2 );
    } else {
      vecY2.Add(vecX2);
    }
    vecY2.ScalarMult(diagScale_);
    
  }
  
  void ArpackMatInterface::MultAVQuad( Complex* x, Complex* y ) {
    
    // do a Matrix times vector operation for the solution of the
    // quadratic eigenvalue problem
    // important notes: 
    // x and y are vectors of length size_, which is twice the 
    //       number of unknowns in this case
    // A is the complex matrix which has of the form
    //   ( -C -K | N 0 )
    // where K and C are the system stiffness and damping matrices, 
    // and N is an arbitrary, nonsingular matrix. (as in MultBV)
    // We just take N to be scaled identity matrix.

    // Create two temporary vectors as wrappers for x and y
    Complex *x1 = x, *x2;
    Complex *y1 = y, *y2;
    Integer neq = size_/2;

    x2 = x1 + neq;
    y2 = y1 + neq;
    Vector<Complex> vecX1, vecY1, vecX2, vecY2, vecTemp(neq);

    vecX1.Replace( neq, x1, false );
    vecY1.Replace( neq, y1, false );
    vecX2.Replace( neq, x2, false );
    vecY2.Replace( neq, y2, false );
    
    // Perform matrix-vector multiplication - upper part
    matrixD_->Mult( vecX1, vecY1 );
    matrixA_->Mult( vecX2, vecTemp );
    vecY1.Axpy(1.0,vecTemp);
    vecY1.ScalarMult(-1.0);

    // Perform matrix-vector multiplication - lower part
    vecY2.Init();
    vecY2.Axpy(diagScale_,vecX1);
    
  }

  void ArpackMatInterface::MultZStiffV( Complex* x, Complex* y ) {
    
    // Create two temporary vectors as wrappers for x and y
    Complex *x1 = x;
    Complex *y1 = y;
    Integer neq = size_/2;

    Vector<Complex> vecX, vecY;
    vecX.Replace( neq, x1, false );
    vecY.Replace( neq, y1, false );
    
    // Perform matrix-vector multiplication
    matrixA_->Mult( vecX, vecY );
    
  }

  void ArpackMatInterface::MultZMassV( Complex* x, Complex* y ) {
    
    // Create two temporary vectors as wrappers for x and y
    Complex *x1 = x;
    Complex *y1 = y;
    Integer neq = size_/2;

    Vector<Complex> vecX, vecY;
    vecX.Replace( neq, x1, false );
    vecY.Replace( neq, y1, false );
    
    // Perform matrix-vector multiplication
    matrixB_->Mult( vecX, vecY );
    
  }

  void ArpackMatInterface::MultZDampV( Complex* x, Complex* y ) {
    
    // Create two temporary vectors as wrappers for x and y
    Complex *x1 = x;
    Complex *y1 = y;
    Integer neq = size_/2;

    Vector<Complex> vecX, vecY;
    vecX.Replace( neq, x1, false );
    vecY.Replace( neq, y1, false );
    
    // Perform matrix-vector multiplication
    matrixD_->Mult( vecX, vecY );
    
  }

  void ArpackMatInterface::ScaleDiag( Complex* x ) {
    
    for (UInt i=0; i<size_/2; i++) {
      x[i] = x[i]*diagScale_;
    }
    
  }

  void ArpackMatInterface::InvScaleDiag( Complex* x ) {
    
    for (UInt i=0; i<size_/2; i++) {
      x[i] = x[i]/diagScale_;
    }
    
  }

  template void ArpackMatInterface::MultShiftOpV<double>(double*, double*);
  template void ArpackMatInterface::MultRegularOpV<double>(double*, double*);
  template void ArpackMatInterface::MultBV<double>(double*, double*);
  template void ArpackMatInterface::MultAV<double>(double*, double*);

  template void ArpackMatInterface::MultShiftOpV<Complex>(Complex*, Complex*);
  template void ArpackMatInterface::MultRegularOpV<Complex>(Complex*, Complex*);
  template void ArpackMatInterface::MultBV<Complex>(Complex*, Complex*);
  template void ArpackMatInterface::MultAV<Complex>(Complex*, Complex*);

  template void ArpackMatInterface::Setup<Double>( BaseSolver* , BasePrecond*, Double, shared_ptr<Timer>, shared_ptr<Timer>);
  template void ArpackMatInterface::Setup<Complex>( BaseSolver* , BasePrecond*, Complex, shared_ptr<Timer>, shared_ptr<Timer>);

  template void ArpackMatInterface::QuadSetup<Double>( BaseSolver*, BasePrecond*, Double );
  template void ArpackMatInterface::QuadSetup<Complex>( BaseSolver*, BasePrecond*, Complex );

}
