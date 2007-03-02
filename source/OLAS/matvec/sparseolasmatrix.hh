#ifndef SPARSE_OLAS_MATRIX_HH
#define SPARSE_OLAS_MATRIX_HH

namespace OLAS {


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
    void CompRes( StdVector &r, const StdVector &x, const StdVector& b ) const{
      ENTER_FCN( "SparseOLASMatrix::CompRes" );
      TRY_CAST {
	RefCast     ( r, Vector<T>, tr );
	ConstRefCast( x, Vector<T>, tx );
	ConstRefCast( b, Vector<T>, tb );
	CompRes( tr, tx ,tb );
      } CATCH_CAST;
    }

    virtual void CompRes( Vector<T> &r, const Vector<T> &x,
			  const Vector<T>& b ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec = this*mvec
    void Mult( const StdVector& mvec, StdVector& rvec ) const {
      ENTER_FCN( "SparseOLASMatrix::Mult" );
      TRY_CAST {
	ConstRefCast( mvec, Vector<T>, tmvec );
	RefCast     ( rvec, Vector<T>, trvec );
	Mult( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void Mult( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
    void MultT( const StdVector& mvec, StdVector& rvec ) const {
      ENTER_FCN( "SparseOLASMatrix::MultT" );
      TRY_CAST {
	ConstRefCast( mvec, Vector<T>, tmvec );
	RefCast     ( rvec, Vector<T>, trvec );
	MultT( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void MultT( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec += this*mvec
    void MultAdd( const StdVector& mvec, StdVector& rvec ) const {
      ENTER_FCN( "SparseOLASMatrix::MultAdd" );
      TRY_CAST {
	ConstRefCast( mvec, Vector<T>, tmvec );
	RefCast     ( rvec, Vector<T>, trvec );
	MultAdd( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void MultAdd( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
    void MultTAdd( const StdVector& mvec, StdVector& rvec ) const {
      ENTER_FCN( "SparseOLASMatrix::MultTAdd" );
      TRY_CAST {
	ConstRefCast( mvec, Vector<T>, tmvec );
	RefCast     ( rvec, Vector<T>, trvec );
	MultTAdd( tmvec, trvec );
      } CATCH_CAST;
    }

    virtual void MultTAdd( const Vector<T>& mvec, Vector<T>& rvec ) const = 0;
    //@}

    //@{
    //! Perform a matrix-vector multiplication rvec -= this*mvec
    void MultSub( const StdVector& mvec, StdVector& rvec ) const {
      ENTER_FCN( "SparseOLASMatrix::MultSub" );
      TRY_CAST {
	ConstRefCast( mvec, Vector<T>, tmvec );
	RefCast     ( rvec, Vector<T>, trvec );
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
    MatrixEntryType GetEntryType() const {
      return  EntryType<T>::M_EntryType;
    }

    //! Return the block size of the matrix
    
    //! The method returns the blocksize (= degrees of freedom)
    //! of the matrix.
    Integer GetBlockSize() const {
      return BlockSize<T>::size;
    };

    //! Get the number of non-zero entries

    //! The method returns the number of non-zero entries. Note that this
    //! value is based upon the stored entry type. In the case of block entries
    //! the latter are counted and the number of entries on the scalar level
    //! will actually be higher depending on the block size.
    Integer GetNnz() const {
      return nnz_;
    }

    //! Get the number of matrix rows on scalar level

    //! The method returns the number of rows of the matrix on the scalar
    //! level, i.e. the count is not based upon the stored entry type but
    //! accumulated over the scalar real/complex entries.
    //! This is different from the GetNrows() method.
    //! \return The same value as GetNrows(), since HyperMatrices only store
    //!         scalar entries anyhow.
    Integer GetNrowsScalar() const {
      ENTER_IFCN( "SparseOLASMatrix::GetNrowsScalar" );
      return nrows_ * BlockSize<T>::size;
    }

    //! Get the number of matrix columns on scalar level

    //! The method returns the number of columns of the matrix on the scalar
    //! level, i.e. the count is not based upon the stored entry type but
    //! accumulated over the scalar real/complex entries.
    //! This is different from the GetNcols() method.
    //! \return The same value as GetNcols(), since HyperMatrices only store
    //!         scalar entries anyhow.
    Integer GetNcolsScalar() const {
      ENTER_IFCN( "SparseOLASMatrix::GetNcolsScalar" );
      return ncols_ * BlockSize<T>::size;
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
    Integer nnz_;

  };

}

#endif
