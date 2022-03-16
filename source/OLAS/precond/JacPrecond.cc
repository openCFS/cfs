#include <def_expl_templ_inst.hh>

#include "MatVec/opdefs.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/Diag_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/VBR_Matrix.hh"
#include "MatVec/BLASLAPACKInterface.hh"
#include "MatVec/generatematvec.hh"
#include "DataInOut/ProgramOptions.hh"
#include "JacPrecond.hh"

// in case of debugging
//#define DEBUG_JAC_PRECOND

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
                                       PtrParamNode precondNode,
				                               PtrParamNode olasInfo )  {
    this->xml_ = precondNode;
    this->infoNode_ = olasInfo->Get("jacobi",progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);
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
    for ( Integer i = 0; i < (Integer) size_; i++ ) {
      z[i] = diagInv_[i] * r[i];
    }
  }


  // ***************************
  //   Setup of preconditioner
  // ***************************
  template<class T_storage,typename T>
  void JacPrecond<T_storage,T>::Setup( T_storage &sysmat ) {

#pragma omp parallel for 
    for ( Integer i = 0; i < (Integer) size_; i++ ) {
      diagInv_[i] = OpType<T>::invert( sysmat.GetDiagEntry(i) );
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
  //   I M P L E M E N T A T I O N    P A R T
  // ========================================================================
  template<class T>
  BlockJacPrecondImpl<T>::BlockJacPrecondImpl( const StdMatrix &mat ) {
    // nothing to be done here
  }
  
  template<class T>
   BlockJacPrecondImpl<T>::~BlockJacPrecondImpl( ) {
    numRows_ = 0;
    factors_.Clear();
  }
  
  
  template<class T>
  void BlockJacPrecondImpl<T>::Setup( StdMatrix &mat ) {
    EXCEPTION("Not imeplemented for the general case");
  }
  template<>
  void BlockJacPrecondImpl<Double>::Setup( StdMatrix &sysmat ) {
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

    // loop over all diagonal blocks
    UInt bs = 0;
    for( UInt i = 0; i < numRows_; ++i ) {

      // get factorization matrix for this block and copy entries
      Matrix<Double> & inv = factors_[i];
      sysmat.GetDiagBlock(i, inv);
      bs = inv.GetNumRows();
#ifdef DEBUG_JAC_PRECOND
      std::cerr << "\n original block #" << i <<"\n";
      for( UInt k = 0; k < bs; ++k ) {
        for( UInt l = 0; l < bs; ++l ) {
          std::cerr << inv[k][l] << ", ";
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
      dgetrf(&n,&n,inv[0],&n,ipiv,&info);
      if( info != 0 ) {
        EXCEPTION("Error during LU-factorization of block #" << i
                  << ". Error value is " << info );
      }
      // invert matrix using previous LU factorization
      dgetri(&n,inv[0],&n,ipiv,work,&lwork,&info);
      if( info != 0 ) {
        EXCEPTION("Error during inversion of block #" << i
                  << ". Error value is " << info );
      }

      delete[] ipiv;
      delete[] work;

#ifdef DEBUG_JAC_PRECOND       
      std::cerr << "\n inverted block #" << i <<"\n";
      std::cerr << inv << std::endl;
#endif
    }
  }
       
      
  
  
  template<class T>
  void BlockJacPrecondImpl<T>::Apply( const Vector<T> &r, Vector<T> &z ) {
    EXCEPTION("Not implemented for the general case");
  }
  
  template<>
  void BlockJacPrecondImpl<Double>::Apply( const Vector<Double> &rhs,
                                           Vector<Double> &sol ) {
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
      dgemv( &trans, &size, &size, &alpha, inv[0], &size, 
             &rhs[rStart], &inc, &beta, &sol[rStart], &inc);

      // sum up current offset
      rStart += size;
    }
  }
  
  template<class T>
    void BlockJacPrecondImpl<T>::GetPrecondSysMat( BaseMatrix& sysMat ) {
    EXCEPTION( "Not implemeneted for the general case" )
  }
  
  template<>
  void BlockJacPrecondImpl<Double>::
     GetPrecondSysMat( BaseMatrix& sysMat ) {

       VBR_Matrix<Double>& ret = dynamic_cast<VBR_Matrix<Double>& >(sysMat);

       sysMat.Init();
       // loop over all diagonal blocks
       for( UInt i = 0; i < numRows_; ++i ) {
         
         Matrix<Double> & inv = factors_[i];
         ret.SetDiagBlock(i, inv);
       }

     }

  
  // ***********************
  // Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class BlockJacPrecondImpl<Double>;
  template class BlockJacPrecondImpl<Complex>;
#endif

  // ========================================================================
  //   B L O C K   J A C O B I   P R E C O N D I T I O N E R
  // ========================================================================
  
    
   template<class T_storage,typename T>
   BlockJacPrecond<T_storage,T>::BlockJacPrecond( const StdMatrix &mat, 
                                                  PtrParamNode precondNode,
                                                  PtrParamNode olasInfo )  {
     this->xml_ = precondNode;
     this->infoNode_ = olasInfo->Get("blockJacobi",progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);
     //numRows_ = 0;
     this->pimpl_ = new BlockJacPrecondImpl<T>(mat );
     
   }

   
   // **************
   //   Destructor
   // **************
   template<class T_storage,typename T>
   BlockJacPrecond<T_storage,T>::~BlockJacPrecond() {
     delete this->pimpl_;
   }


   // ************************
   //   Apply preconditioner
   // ************************
   template<class T_storage,typename T>
      void BlockJacPrecond<T_storage,T>::Apply( const T_storage &sysmat,
                                                const Vector<T> &rhs,
                                                Vector<T> &sol ) {
    this->pimpl_->Apply( rhs, sol );
   }
   


   // ***************************
   //   Setup of preconditioner
   // ***************************
   template<class T_storage,typename T>
   void BlockJacPrecond<T_storage,T>::Setup( T_storage &sysmat) {
     this->pimpl_->Setup( sysmat );
   }

   template<class T_storage,typename T>
   void BlockJacPrecond<T_storage,T>::
   GetPrecondSysMat( BaseMatrix& sysMat ) {

     this->pimpl_->GetPrecondSysMat( sysMat );
   }
   // ***********************
   // Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
   template class BlockJacPrecond< VBR_Matrix<Double>, Double>;
   template class BlockJacPrecond< VBR_Matrix<Complex>, Complex>;
   template class BlockJacPrecond< SCRS_Matrix<Double> , Double>;
   template class BlockJacPrecond< SCRS_Matrix<Complex>, Complex>;
#endif

}
