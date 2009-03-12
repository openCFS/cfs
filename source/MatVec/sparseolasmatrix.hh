// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef SPARSE_MATRIX_HH
#define SPARSE_MATRIX_HH

#include "basematrix.hh"
#include "stdmatrix.hh"
#include "vector.hh"

namespace CoupledField {


  //! This class constitutes an intermediate level of inheritance between
  //! the StdMatrix class and the proprietory sparse matrix classes of OLAS.
  //! It serves only two purposes. The first one is to provide all the
  //! interfaces for arthimetic methods required by the StdMatrix. All that the
  //! respective implementations of these methods do is to down-cast the input
  //! vectors to the Vector<T> class and, if this succeeds, to call the
  //! respective method of the derived matrix class. The second purpose is to
  //! implement some methods that do not rely on the storage layout of the
  //! derived matrix.
  template<typename T>
  class SparseOLASMatrix : public StdMatrix {

  public:

    using BaseMatrix::Mult;
    using BaseMatrix::MultT;
    using BaseMatrix::MultAdd;
    using BaseMatrix::MultTAdd;
    using BaseMatrix::MultSub;
    using BaseMatrix::CompRes;
    using BaseMatrix::Add;

    // ========================================================================
    // ARITHMETIC METHODS
    // ========================================================================

    //! \name Methods for arithmetic operations
    //! These methods try to down-cast the input vectors to the respective
    //! Vector<T> class and, if this succeeds, call the respective method
    //! of the derived class, the latter is forced to implement it, since
    //! we specify a purely virtual interface in this class.

    //@{
    //! Compute the residual \f$r=b-Ax\f$
    void CompRes( SingleVector &r, const SingleVector &x, const SingleVector& b ) const{
      TRY_CAST {
        REFCAST     ( r, Vector<T>, tr );
        CONSTREFCAST( x, Vector<T>, tx );
        CONSTREFCAST( b, Vector<T>, tb );
        CompRes( tr, tx ,tb );
      } CATCH_CAST;
    }

    virtual void CompRes( Vector<T> &r, const Vector<T> &x,
                          const Vector<T>& b ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec = this*mvec
    void Mult( const SingleVector& mvec, SingleVector& rvec ) const 
    {
      TRY_CAST {
        CONSTREFCAST( mvec, Vector<T>, tmvec );
        REFCAST     ( rvec, Vector<T>, trvec );
        Mult( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void Mult( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
    void MultT( const SingleVector& mvec, SingleVector& rvec ) const {
      TRY_CAST {
        CONSTREFCAST( mvec, Vector<T>, tmvec );
        REFCAST     ( rvec, Vector<T>, trvec );
        MultT( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void MultT( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec += this*mvec
    void MultAdd( const SingleVector& mvec, SingleVector& rvec ) const {
      TRY_CAST {
        CONSTREFCAST( mvec, Vector<T>, tmvec );
        REFCAST     ( rvec, Vector<T>, trvec );
        MultAdd( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void MultAdd( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
    void MultTAdd( const SingleVector& mvec, SingleVector& rvec ) const {
      TRY_CAST {
        CONSTREFCAST( mvec, Vector<T>, tmvec );
        REFCAST     ( rvec, Vector<T>, trvec );
        MultTAdd( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void MultTAdd( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec -= this*mvec
    void MultSub( const SingleVector& mvec, SingleVector& rvec ) const {
      TRY_CAST {
        CONSTREFCAST( mvec, Vector<T>, tmvec );
        REFCAST     ( rvec, Vector<T>, trvec );
        MultSub( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void MultSub( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    // ========================================================================
    // QUERY GENERAL MATRIX INFORMATION
    // ========================================================================

    //@{

    //! \name Query general matrix format

    //! Return the Entry type of the matrix
  
    //! The method returns the entry type of the matrix (i.e. Double, Complex,
    //! or whatever it is). This is encoded as a value of the enumeration data
    //! type MatrixEntryType.
    BaseMatrix::EntryType GetEntryType() const {
      return  CoupledField::EntryType<T>::M_EntryType;
    }

    //! Get the number of non-zero entries

    //! The method returns the number of non-zero entries. Note that this
    //! value is based upon the stored entry type. In the case of block entries
    //! the latter are counted and the number of entries on the scalar level
    //! will actually be higher depending on the block size.
    UInt GetNnz() const {
      return nnz_;
    }

    //! Trasforms SCRS matrix into vector containing all upper
    // triangle elements further usage in CFS++
    virtual void CopySCRSMatrix2Vec(Complex* &A){;};

    //@}

  protected:

    //! Number of (potentially) non-zero matrix entries

    //! The attribute stores the number of potentially non-zero matrix
    //! entries. Potentially means that the entries belong to the non-zero
    //! pattern of the matrix, but are not guaranteed to be non-zero on
    //! the numerical level. Note that this count is based upon the stored
    //! entry type. In the case of block entries the latter are counted and the
    //! number of entries on the scalar level is actually higher depending on
    //! the block size.
    UInt nnz_;

  };

}

#endif
