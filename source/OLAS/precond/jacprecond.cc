// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/opdefs.hh"
#include "MatVec/crs_matrix.hh"
#include "MatVec/diag_matrix.hh"
#include "MatVec/scrs_matrix.hh"
#include "MatVec/vbr_matrix.hh"
#include "MatVec/matrixBLASSupport.hh"

#include "jacprecond.hh"

namespace CoupledField {



  // ========================================================================
  //  S C A L A R   J A C O B I   P R E C O N D I T I O N E R
  // ========================================================================
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
                                       PtrParamNode solverNode,
				                               PtrParamNode olasInfo )  {
    this->xml_ = solverNode;
    this->olasInfo_ = olasInfo->Get("jacobi");
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
  
  


  // ========================================================================
  //   B L O C K   J A C O B I   P R E C O N D I T I O N E R
  // ========================================================================
   // ***********************
   template<class T_storage,typename T>
   BlockJacPrecond<T_storage,T>::BlockJacPrecond( const StdMatrix &mat, 
                                                  PtrParamNode solverNode,
                                                  PtrParamNode olasInfo )  {
     this->xml_ = solverNode;
     this->olasInfo_ = olasInfo->Get("blockJacobi");
     numRows_ = 0;
     
   }

   
   // **************
   //   Destructor
   // **************
   template<class T_storage,typename T>
   BlockJacPrecond<T_storage,T>::~BlockJacPrecond() {
   }


   // ************************
   //   Apply preconditioner
   // ************************
   template<class T_storage,typename T>
   void BlockJacPrecond<T_storage,T>::Apply( const T_storage &sysmat, 
                                             const Vector<T> &rhs,
                                             Vector<T> &sol ) const {
     EXCEPTION("No general implementation available");
   }
       
   template<>
   void BlockJacPrecond<VBR_Matrix<Double>,Double>::Apply( const VBR_Matrix<Double> &sysmat, 
                                                           const Vector<Double> &rhs,
                                                           Vector<Double> &sol ) const {

     UInt rStart = 0;
     char trans = 'T';
     
     // BLAS related stuff
     Double alpha = 1.0;
     Double beta = 0.0;
     Integer inc = 1;
     sol.Init();
     
     // loop over all diagonals
     for( UInt i = 0; i < numRows_; ++i ) {

       // perform matrix-vector multiplication of given block
       const Matrix<Double> & inv = factors_[i];
       Integer size = inv.GetNumRows();
       DGEMV( &trans, &size, &size, &alpha, inv[0], &size, 
              &rhs[rStart], &inc, &beta, &sol[rStart], &inc);

       // sum up current offset
       rStart += size;
     }
   }


   // ***************************
   //   Setup of preconditioner
   // ***************************
   template<class T_storage,typename T>
     void BlockJacPrecond<T_storage,T>::Setup( T_storage &sysmat ) {
     EXCEPTION("No general implementation available");
   }
   
   template<>
   void BlockJacPrecond<VBR_Matrix<Double>,Double>::Setup( VBR_Matrix<Double> &sysmat ) {

     
     // get block information
     UInt nbRow, nbCol, nBlocks;
     sysmat.GetNumBlocks( nbRow, nbCol, nBlocks );
     
     numRows_ = nbRow;
     factors_.Resize( numRows_ );
     
     // make sure we have a square matrix
     if( nbRow != nbCol ) {
       EXCEPTION( "Block Jacobi preconditioner only works with"
                  " square matrices");
     }

     Double * ptDiag = NULL;
     UInt numEntries = 0;
     
     // loop over all diagonal blocks
     UInt bs = 0, rStart = 0;
     for( UInt i = 0; i < numRows_; ++i ) {

       // obtain block information
       ptDiag = sysmat.GetDiagBlock( i,bs, rStart);
       
       // get factorization matrix for this block and copy entries
       Matrix<Double> & inv = factors_[i];
       inv.Resize(bs);
       std::copy(ptDiag,ptDiag+bs*bs, inv[0]);
       numEntries += (bs * bs);
       
       // ===================================================================
       // VERSION 1: LU BASED INVERSION (also for non-symmetric or  non SPD)
       // ===================================================================
       int *ipiv = new int[bs+1];
       int n = bs;
       int lwork = bs*bs;
       double *work = new double[lwork];
       int info;

       // calculate LU-factorization of block
       LP_DGETRF(&n,&n,inv[0],&n,ipiv,&info);
       if( info != 0 ) {
           EXCEPTION("Error during LU-factorization of block #" << i
                     << ". Error value is " << info );
       }
       // invert matrix using previous LU factorization
       LP_DGETRI(&n,inv[0],&n,ipiv,work,&lwork,&info);
       if( info != 0 ) {
           EXCEPTION("Error during inversion of block #" << i
                     << ". Error value is " << info );
       }
       
       delete[] ipiv;
       delete[] work;
     }
     

     // log some information to info
     olasInfo_->Get("numDiagBlocks")->SetValue(numRows_);
     olasInfo_->Get("numEntries")->SetValue(numEntries);
   }
   
 // Explicit template instantiation
 #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
   template class BlockJacPrecond< VBR_Matrix<Double>, Double>;
   template class BlockJacPrecond< VBR_Matrix<Complex>, Complex>;
 #endif
}
