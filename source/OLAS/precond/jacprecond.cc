// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "precond/jacprecond.hh"

namespace OLAS {

  // ****************
  //   Constructor
  // ****************
  // template<class T_storage,typename T>
  // JacPrecond<T_storage,T>::JacPrecond(Integer asize ) {
  //   ENTER_FCN("JacPrecond::JacPrecond");
  //   size_   = asize;
  //   NewArray(diagInv_,T,size_);
  // }


  // ***********************
  //   Another Constructor
  // ***********************
  template<class T_storage,typename T>
  JacPrecond<T_storage,T>::JacPrecond( const StdMatrix &mat, 
                                       OLAS_Params *myParams,
				       OLAS_Report *myReport )  {
    ENTER_FCN( "JacPrecond::JacPrecond" );
    this->myParams_ = myParams;
    this->myReport_ = myReport;
    size_     = mat.GetNrows();
    NewArray( diagInv_, T, size_ );
  }

  
  // **************
  //   Destructor
  // **************
  template<class T_storage,typename T>
  JacPrecond<T_storage,T>::~JacPrecond() {
    ENTER_FCN( "JacPrecond::~JacPrecond" );
    DeleteArray( diagInv_ );
  }


  // ************************
  //   Apply preconditioner
  // ************************
  template<class T_storage,typename T>
  void JacPrecond<T_storage,T>::Apply( const T_storage &sysmat,	
				       const Vector<T> &r,
				       Vector<T> &z ) const {
    ENTER_FCN( "JacPrecond::Apply" );

    PROFILE( "JacPrecond::Apply",
             size_ * BlockSize<T>::size * BlockSize<T>::size );

#pragma omp parallel for 
    for ( Integer i = 1; i <= size_; i++ ) {
      z[i] = diagInv_[i] * r[i];
    }
  }


  // ***************************
  //   Setup of preconditioner
  // ***************************
  template<class T_storage,typename T>
  void JacPrecond<T_storage,T>::Setup( T_storage &sysmat ) {
    ENTER_FCN( "JacPrecond::Setup" );  

#pragma omp parallel for 
    for ( Integer i = 1; i <= size_; i++ ) {
      diagInv_[i] = opType<T>::invert( sysmat.GetDiag(i) );
    }

#ifdef DEBUG_JACPRECOND
    (*debug) << "JacPrecond calculated:" << std::endl;
    for ( UInt i = 1; i <= size_; i++ ) {
      (*debug) << " " << diagInv_[i];
    }
    (*debug) << std::endl;
#endif

  }

}
