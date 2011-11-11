#include "MatVec/opdefs.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/Diag_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/VBR_Matrix.hh"
#include "MatVec/matrixBLASSupport.hh"
#include "MatVec/generatematvec.hh"

#include "JacPrecond.hh"

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
    this->infoNode_ = olasInfo->Get("jacobi",ParamNode::APPEND);
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
				       Vector<T> &z ) {

#pragma omp parallel for 
    for ( UInt i = 0; i < size_; i++ ) {
      z[i] = diagInv_[i] * r[i];
    }
  }


  // ***************************
  //   Setup of preconditioner
  // ***************************
  template<class T_storage,typename T>
  void JacPrecond<T_storage,T>::Setup( T_storage &sysmat, PtrParamNode analysis_id ) {

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
  
  // ***********************************
  //   Export of preconditioned matrix
  // ***********************************
  template<class T_storage,typename T>
  void JacPrecond<T_storage,T>::
  GetPrecondSysMat( BaseMatrix& sysMat ) {
    
    // generate copy of system matrix
    StdMatrix & std =  dynamic_cast<StdMatrix&>(sysMat);
    std.Init();
    
    // loop over diagonals and set it to 1
    for( UInt i = 0; i < size_; ++i ) {
      std.SetDiagEntry(i, diagInv_[i]);
   }
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
     this->infoNode_ = olasInfo->Get("blockJacobi",ParamNode::APPEND);
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
                                             Vector<T> &sol ) {
     EXCEPTION("No general implementation available");
   }
       
   template<>
   void BlockJacPrecond<VBR_Matrix<Double>,Double>::Apply( const VBR_Matrix<Double> &sysmat, 
                                                           const Vector<Double> &rhs,
                                                           Vector<Double> &sol ){

     UInt rStart = 0;
     char trans = 'T';
     
     // BLAS related stuff
     Double alpha = 1.0;
     Double beta = 0.0;
     Integer inc = 1;
     sol.Init();
     //std::cerr << "rhs is \n" << rhs << std::endl;
     
     // loop over all diagonals
     for( UInt i = 0; i < numRows_; ++i ) {

       // perform matrix-vector multiplication of given block
       const Matrix<Double> & inv = factors_[i];
       Integer size = inv.GetNumRows();
       DGEMV( &trans, &size, &size, &alpha, inv[0], &size, 
              &rhs[rStart], &inc, &beta, &sol[rStart], &inc);
       
//       std::cerr << "\nBlockPrecond: Row "<< i << ", Diag is \n" << inv.ToString() << std::endl;
//       std::cerr << "rhs is\n";
//       for( int j = 0; j < size; ++j  ) {
//         std::cerr << rhs[j+rStart] << ";";
//       }
//       std::cerr << "\nsol is\n";
//              for( int j = 0; j < size; ++j  ) {
//                std::cerr << sol[j+rStart] << ";";
//              }

       // sum up current offset
       rStart += size;
     }
//     std::cerr << "sol is \n" << sol.ToString() << std::endl;
   }


   // ***************************
   //   Setup of preconditioner
   // ***************************
   template<class T_storage,typename T>
     void BlockJacPrecond<T_storage,T>::Setup( T_storage &sysmat,
                                               PtrParamNode analysis_id ) {
     EXCEPTION("No general implementation available");
   }
   
   template<>
   void BlockJacPrecond<VBR_Matrix<Double>,Double>::Setup( VBR_Matrix<Double> &sysmat,
                                                           PtrParamNode infoNode ) {

     
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
//#define DEBUG_JAC_PRECOND
#ifdef DEBUG_JAC_PRECOND
       std::cerr << "\n original block #" << i <<"\n";
       for( UInt k = 0; k < bs; ++k ) {
         for( UInt l = 0; l < bs; ++l ) {
           std::cerr << ptDiag[k*bs+l] << ", ";
         }
         std::cerr << "\n";
       }
#endif
       
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

#ifdef DEBUG_JAC_PRECOND       
       std::cerr << "\n inverted block #" << i <<"\n";
       std::cerr << inv << std::endl;
//       for( UInt k = 0; k < bs; ++k ) {
//         for( UInt l = 0; l < bs; ++l ) {
//           std::cerr << inv[k*bs+l] << ", ";
//         }
//         std::cerr << "\n";
//       }
#endif
     }
     

     // log some information to info
     infoNode_->Get("numDiagBlocks")->SetValue(numRows_);
     infoNode_->Get("numEntries")->SetValue(numEntries);
   }

   template<class T_storage,typename T>
   void BlockJacPrecond<T_storage,T>::
   GetPrecondSysMat( BaseMatrix& sysMat ) {

     VBR_Matrix<T>& ret = dynamic_cast<VBR_Matrix<T>& >(sysMat);

     T * ptDiag = NULL;
     sysMat.Init();

     // loop over all diagonal blocks
     UInt bs = 0, rStart = 0;
     for( UInt i = 0; i < numRows_; ++i ) {
       
       Matrix<T> & inv = factors_[i];

       // obtain block information
       ptDiag = ret.GetDiagBlock( i,bs, rStart);

       // set b
       for( UInt i = 0; i < bs ; ++i ) {
         for( UInt j = 0; j < bs ; ++j ) {
         ptDiag[i*bs +j] = inv[i][j];
         }
       }

     }

   }

 // Explicit template instantiation
 #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
   template class BlockJacPrecond< VBR_Matrix<Double>, Double>;
   template class BlockJacPrecond< VBR_Matrix<Complex>, Complex>;
 #endif
}
