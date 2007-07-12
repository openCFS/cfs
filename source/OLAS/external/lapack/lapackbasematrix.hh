// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_LAPACKBASEMATRIX_HH
#define OLAS_LAPACKBASEMATRIX_HH

#include "matvec/matvec.hh"
#include "utils/utils.hh"

namespace OLAS {
  
  //! This is the base class for LAPACK matrices in OLAS.

  //! This class represents the base class for all LAPACK matrices in OLAS.
  //! This base class will perform a dummy implementation of all the pure
  //! virtual methods of the BaseMatrix and StdMatrix classes, which are not
  //! required for a direct solution of the linear system using LAPACK
  //! routines. Besides this it also implements some "small" methods that
  //! are independent of the storage layout and entry types.
  //! \htmlonly
  //! <center><img src="../AddDoc/lapack.jpg"></center>
  //! \endhtmlonly
  //! For further information see the
  //! <a href="http://www.netlib.org/lapack">LAPACK</a> pages at Netlib.
  class LapackBaseMatrix : public StdMatrix {

  public:

    //! Return the block size of the matrix

    //! The method returns the blocksize (= degrees of freedom) of the matrix.
    //! Currently the LAPACK matrices only support scalar entries. Thus, the
    //! answer will always be one.
    Integer GetBlockSize() const {
      ENTER_IFCN( "LapackBaseMatrix::GetBlockSize" );
      return 1;
    };


    // ************************************************************************
    //   Dummy interface/implementation of methods which are not useful in the
    //   context of a LAPACK Matrix.
    // ************************************************************************

    //@{
    //! Dummy implementation

    //! This method is currently not implemented for any type of LAPACK Matrix.
    //! Since these types of matrices are intended to be used in conjunction
    //! with the direct methods LAPACK offers it is unclear, whether this
    //! functionality will be needed in OLAS.
    void Mult(const SparseVector& mvec, SparseVector& rvec) const {
      Error( "Mult not implemented by any LAPACK Matrix", __FILE__, __LINE__ );
    }
    void MultT(const SparseVector& mvec, SparseVector& rvec) const {
      Error( "MultT not implemented by any LAPACK Matrix", __FILE__,
	     __LINE__ );
    }
    void MultAdd( const SparseVector& mvec, SparseVector& rvec ) const {
      Error( "MultAdd not implemented by any LAPACK Matrix", __FILE__,
	     __LINE__ );
    }
    void MultTAdd( const SparseVector& mvec, SparseVector& rvec ) const {
      Error( "MultTAdd not implemented by any LAPACK Matrix", __FILE__,
	     __LINE__ );
    }
    void MultSub( const SparseVector& mvec, SparseVector& rvec ) const {
      Error( "MultSub not implemented by any LAPACK Matrix", __FILE__,
	     __LINE__ );
    }
    void CompRes( SparseVector &r, const SparseVector &x, const SparseVector& b ) const{
      Error( "CompRes not implemented by any LAPACK Matrix", __FILE__,
	     __LINE__ );
    }
    Double GetMaxDiag() const {
      Error( "GetMaxDiag not implemented by any LAPACK Matrix", __FILE__,
	     __LINE__ );
      return 0;
    }
    //@}

    //@{
    //! Not yet implemented

    //! This method is not yet implemented. Due to its nature is must be
    //! implemented in the corresponding derived class!
    // void Export( char *fname, char *comment = NULL ) const {
    //   Error( "Export not yet implemented for any LAPACK Matrix", __FILE__,
    //          __LINE__);
    // }
    void Print( std::ostream &os ) const {
      Warning( "Print not yet implemented for any LAPACK Matrix", __FILE__,
	       __LINE__);
    }
    void Add( const Double a, const StdMatrix& mat ) {
      Error( "Add not yet implemented for any LAPACK Matrix", __FILE__,
	     __LINE__);
    }
    //@}

  };

}

#endif
