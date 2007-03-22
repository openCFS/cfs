// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "precond/ssorprecond.hh"

namespace OLAS {

  // ***************
  //   Constructor
  // ***************
  template <typename T>
  SSORPrecond<T>::SSORPrecond( const StdMatrix& mat, OLAS_Params *myParams,
			       OLAS_Report *myReport ) {
    ENTER_FCN( "SSORPrecond::SSORPrecond" );
    this->myParams_ = myParams;
    this->myReport_ = myReport;
    size_     = mat.GetNrows();
    NewArray( diagInv_, T, size_ );
  }


  // **************
  //   Destructor
  // **************
  template <typename T>
  SSORPrecond<T>::~SSORPrecond() {
    ENTER_FCN( "SSORPrecond::~SSORPrecond" );
    DeleteArray( diagInv_ );
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void SSORPrecond<T>::Setup( CRS_Matrix<T> &sysmat ) {
    ENTER_FCN( "SSORPrecond::Setup" );
    for ( UInt i=1; i<=size_; i++ ) {
      diagInv_[i] = opType<T>::invert(sysmat.GetDiag(i));
    }
  }


  // *******************************
  //   Apply: apply preconditioner
  // *******************************
  template <typename T> 
  void SSORPrecond<T>::Apply( const CRS_Matrix<T> &sysmat,
			      const Vector<T> &r, Vector<T> &z ) const {

    ENTER_FCN( "SSORPrecond::Apply" );

    UInt i, j;
    T_Vtype sum;

    // Query relaxation parameter
    Double omega = this->myParams_->GetDoubleValue( "SSOR_Omega" );

    // Extract matrix info/data
    const T *valA = sysmat.GetDataPointer();
    const Integer *rowP = sysmat.GetRowPointer();
    const Integer *colP = sysmat.GetColPointer();

    // Set initial guess to zero
    z.Init();

    // =============================
    //   Perform forward SOR sweep
    // =============================
    for ( i = 1; i <= size_; i++ ) {

      // Compute current local residual: sum = b - (L*xnew + U*xold)
      sum = r[i];
      for ( j = rowP[i]; j < rowP[i+1]; j++ ) {
	sum -= valA[j] * z[colP[j]];
      }

      // Compute new approximation for this unknown (starting from zero)
      z[i] = omega * diagInv_[i] * sum;
    }

    // =============================
    //   Perform backward SOR sweep
    // =============================
    for( i = size_; i >= 1; i-- ) {

      // Compute current local residual: sum = b - (L*xnew + U*xold)
      sum = r[i];
      for ( j = rowP[i]; j < rowP[i+1]; j++ ) {
	sum -= valA[j] * z[colP[j]];
      }

      // Compute new approximation for this unknown
      z[i] += omega * diagInv_[i] * sum;
    }
  }

} // namespace
