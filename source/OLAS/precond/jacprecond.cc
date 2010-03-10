// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/opdefs.hh"
#include "MatVec/crs_matrix.hh"
#include "MatVec/diag_matrix.hh"
#include "MatVec/scrs_matrix.hh"

#include "jacprecond.hh"

namespace CoupledField {

  // ****************
  //   Constructor
  // ****************
  // template<class T_storage,typename T>
  // JacPrecond<T_storage,T>::JacPrecond(Integer asize ) {
  //   size_   = asize;
  //   NEWARRAY(diagInv_,T,size_);
  // }


  // ***********************
  //   Another Constructor
  // ***********************
  template<class T_storage,typename T>
  JacPrecond<T_storage,T>::JacPrecond( const StdMatrix &mat, 
                                       ParamNode *solverNode,
				                               InfoNode *olasInfo )  {
    this->xml_ = solverNode;
    this->olasInfo_ = olasInfo;
    size_     = mat.GetNumRows();
    NEWARRAY( diagInv_, T, size_ );
  }

  
  // **************
  //   Destructor
  // **************
  template<class T_storage,typename T>
  JacPrecond<T_storage,T>::~JacPrecond() {
    delete [] ( diagInv_ );
  }


  // ************************
  //   Apply preconditioner
  // ************************
  template<class T_storage,typename T>
  void JacPrecond<T_storage,T>::Apply( const T_storage &sysmat,	
				       const Vector<T> &r,
				       Vector<T> &z ) const {

#pragma omp parallel for 
    for ( UInt i = 0; i < size_; i++ ) {
      z[i] = diagInv_[i] * r[i];
    }
  }


  // ***************************
  //   Setup of preconditioner
  // ***************************
  template<class T_storage,typename T>
  void JacPrecond<T_storage,T>::Setup( T_storage &sysmat ) {

#pragma omp parallel for 
    for ( UInt i = 0; i < size_; i++ ) {
      diagInv_[i] = OpType<T>::invert( sysmat.GetDiag(i) );
    }

#ifdef DEBUG_JACPRECOND
    (*debug) << "JacPrecond calculated:" << std::endl;
    for ( UInt i = 0; i < size_; i++ ) {
      (*debug) << " " << diagInv_[i];
    }
    (*debug) << std::endl;
#endif

  }
  
// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class JacPrecond< CRS_Matrix<Double>, Double>;
  template class JacPrecond< CRS_Matrix<Complex>, Complex>;
  template class JacPrecond< Diag_Matrix<Double>, Double>;
  template class JacPrecond< Diag_Matrix<Complex>, Complex>;
  template class JacPrecond< SCRS_Matrix<Double>, Double>;
  template class JacPrecond< SCRS_Matrix<Complex>, Complex>;
#endif
  
}
