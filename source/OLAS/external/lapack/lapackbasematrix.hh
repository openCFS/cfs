// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef OLAS_LAPACKBASEMATRIX_HH
#define OLAS_LAPACKBASEMATRIX_HH
#include <string>
#include "General/defs.hh"
#include "General/exception.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"
namespace CoupledField {
 
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
class SingleVector;
  class LapackBaseMatrix : public StdMatrix {
  public:
   
    virtual ~LapackBaseMatrix() {}
    //! Return the block size of the matrix
    //! The method returns the blocksize (= degrees of freedom) of the matrix.
    //! Currently the LAPACK matrices only support scalar entries. Thus, the
    //! answer will always be one.
    Integer GetBlockSize() const {
      return 1;
    };
    // ************************************************************************
    //   Dummy interface/implementation of methods which are not useful in the
    //   context of a LAPACK Matrix.
    // ************************************************************************
    using BaseMatrix::Mult;
    using BaseMatrix::MultT;
    using BaseMatrix::MultAdd;
    using BaseMatrix::MultTAdd;
    using BaseMatrix::MultSub;
    using BaseMatrix::CompRes;
    using BaseMatrix::Add;
    //@{
    //! Dummy implementation
    //! This method is currently not implemented for any type of LAPACK Matrix.
    //! Since these types of matrices are intended to be used in conjunction
    //! with the direct methods LAPACK offers it is unclear, whether this
    //! functionality will be needed in OLAS.
    void Mult(const SingleVector& mvec, SingleVector& rvec) const {
      EXCEPTION( "Mult not implemented by any LAPACK Matrix" );
    }
    void MultT(const SingleVector& mvec, SingleVector& rvec) const {
      EXCEPTION( "MultT not implemented by any LAPACK Matrix" );
    }
    void MultAdd( const SingleVector& mvec, SingleVector& rvec ) const {
      EXCEPTION( "MultAdd not implemented by any LAPACK Matrix" );
    }
    void MultTAdd( const SingleVector& mvec, SingleVector& rvec ) const {
      EXCEPTION( "MultTAdd not implemented by any LAPACK Matrix" );
    }
    void MultSub( const SingleVector& mvec, SingleVector& rvec ) const {
      EXCEPTION( "MultSub not implemented by any LAPACK Matrix" );
    }
    void CompRes( SingleVector &r, const SingleVector &x, const SingleVector& b ) const{
      EXCEPTION( "CompRes not implemented by any LAPACK Matrix" );
    }
    Double GetMaxDiag() const {
      EXCEPTION( "GetMaxDiag not implemented by any LAPACK Matrix" );
      return 0;
    }
    //@}
    //@{
    //! Not yet implemented
    //! This method is not yet implemented. Due to its nature is must be
    //! implemented in the corresponding derived class!
    // void Export( char *fname, char *comment = NULL ) const {
    //   EXCEPTION( "Export not yet implemented for any LAPACK Matrix");
    // }
    std::string ToString( char colSeparator = ' ',
                          char rowSeparator = '\n' ) const {
      return "LapackBaseMatrix: Print not yet implemented for any LAPACK Matrix";
    }
    void Add( const Double a, const StdMatrix& mat ) {
      EXCEPTION( "Add not yet implemented for any LAPACK Matrix" );
    }
    //@}
  };
}
#endif