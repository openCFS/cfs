#include "MatVec/StdMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/opdefs.hh"
#include "DataInOut/ProgramOptions.hh"
#include "OLAS/precond/SSORPrecond.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  template <typename T>
  SSORPrecond<T>::SSORPrecond( const StdMatrix& mat, PtrParamNode precondNode,
			       PtrParamNode olasInfo ) {
    
    // The SSOR preconditioner is currently not working
    WARN("The SSOR-Preconditioner is not yet ported from 1 to 0 based"
        "numbering, i.e. it will NOT work correctly!" );
    
    this->xml_ = precondNode;
    this->infoNode_ = olasInfo->Get("ssor", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);
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
  void SSORPrecond<T>::Setup( CRS_Matrix<T> &sysmat) {
    for ( UInt i=1; i<=size_; i++ ) {
      diagInv_[i] = OpType<T>::invert(sysmat.GetDiag(i));
    }
  }


  // *******************************
  //   Apply: apply preconditioner
  // *******************************
  template <typename T> 
  void SSORPrecond<T>::Apply( const CRS_Matrix<T> &sysmat,
			      const Vector<T> &r, Vector<T> &z ) {


    UInt i;
    UInt j;
    T sum;

    // Query relaxation parameter
    Double omega = 1.7;
    this->xml_->GetValue("omega", omega, ParamNode::INSERT);

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
