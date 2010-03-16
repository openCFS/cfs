// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/stdmatrix.hh"
#include "MatVec/crs_matrix.hh"
#include "MatVec/opdefs.hh"

#include "OLAS/precond/ssorprecond.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  template <typename T>
  SSORPrecond<T>::SSORPrecond( const StdMatrix& mat, PtrParamNode solverNode,
			       PtrParamNode olasInfo ) {
    this->xml_ = solverNode;
    this->olasInfo_ = olasInfo;
    size_     = mat.GetNumRows();
    NEWARRAY( diagInv_, T, size_ );
  }


  // **************
  //   Destructor
  // **************
  template <typename T>
  SSORPrecond<T>::~SSORPrecond() {
    delete [] ( diagInv_ );
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void SSORPrecond<T>::Setup( CRS_Matrix<T> &sysmat ) {
    for ( UInt i=1; i<=size_; i++ ) {
      diagInv_[i] = OpType<T>::invert(sysmat.GetDiag(i));
    }
  }


  // *******************************
  //   Apply: apply preconditioner
  // *******************************
  template <typename T> 
  void SSORPrecond<T>::Apply( const CRS_Matrix<T> &sysmat,
			      const Vector<T> &r, Vector<T> &z ) const {


    UInt i;
    UInt j;
    T_Vtype sum;

    // Query relaxation parameter
    Double omega = 1.7;
    PtrParamNode pNode = this->xml_->Get("SSOR", ParamNode::INSERT);
    pNode->GetValue("omega", omega, ParamNode::INSERT);

    // Extract matrix info/data
    const T *valA = sysmat.GetDataPointer();
    const UInt *rowP = sysmat.GetRowPointer();
    const UInt *colP = sysmat.GetColPointer();

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

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class SSORPrecond<Double>;
  template class SSORPrecond<Complex>;
#endif
  
} // namespace
