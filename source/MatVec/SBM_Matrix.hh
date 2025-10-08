#ifndef OLAS_SBM_MATRIX_HH
#define OLAS_SBM_MATRIX_HH

#include <iostream>

#include "OLAS/graph/GraphManager.hh"

#include "BaseMatrix.hh"
#include "StdMatrix.hh"
#include "Vector.hh"
#include "SBM_Vector.hh"
#include "generatematvec.hh"

namespace CoupledField {

  //! This class implements the super-block-matrix concept

  //! A super-block-matrix is a dense matrix whose entries are sparse matrices
  //! that can be composed of scalar or block entries.\n
  //! The following describes the storage scheme of the %SBM_Matrix class.
  //! Assume that we have a \f$2\times 2\f$ block-matrix with the sub-matrices
  //! \f$A\f$, \f$B\f$, \f$C\f$ and \f$D\f$ where \f$B\f$ and \f$C\f$ are not
  //! necessarily square.
  //! - If all four sub-matrices are not identical to zero, then
  //! \f$\displaystyle
  //! \left(\begin{array}{cc}
  //! \displaystyle A & \displaystyle B \\
  //! \displaystyle C & \displaystyle D
  //! \end{array}\right)
  //! \f$ is stored as
  //! \f$\displaystyle
  //! \left[\begin{array}{cc}
  //! \displaystyle  A  & \displaystyle B \\
  //! \displaystyle  C  & \displaystyle D
  //! \end{array}\right]
  //! \f$
  //!
  //! - If the matrix is symmetric, then the constructor can be informed on
  //!   this during instantiation and only the upper triangle of the matrix
  //!   needs to be stored, i.e.
  //! \f$\displaystyle
  //! \left(\begin{array}{cc}
  //! \displaystyle A   & \displaystyle B \\
  //! \displaystyle B^T & \displaystyle D
  //! \end{array}\right)
  //! \f$ is stored as
  //! \f$\displaystyle
  //! \left[\begin{array}{cc}
  //! \displaystyle  A  & \displaystyle B \\
  //! \displaystyle  -  & \displaystyle D
  //! \end{array}\right]
  //! \f$
  //! Actually a symmetric %SBM_Matrix object will discard all attempts to
  //! set a matrix in the lower triangular part and issue a warning.
  //! The arithmetic method will handle the symmetry of a super-block matrix
  //! automatically if the symmetry flag was set during construction.
  //!
  //! - Also, if a sub-matrix is identical to zero, then it is not required to
  //!   store it explicitely, i.e.
  //! \f$\displaystyle
  //! \left(\begin{array}{cc}
  //! \displaystyle A & \displaystyle 0 \\
  //! \displaystyle C & \displaystyle D
  //! \end{array}\right)
  //! \f$ can be stored as
  //! \f$\displaystyle
  //! \left[\begin{array}{cc}
  //! \displaystyle  A  & \displaystyle - \\
  //! \displaystyle  C  & \displaystyle D
  //! \end{array}\right]
  //! \f$
  //! The arithmetic operations treat non-existent sub-matrices as zero
  //! automatically.
  //!
  //! - It is safe to add %SBM_Matrices with different block-sparsity
  //!   patterns, i.e.
  //!   \f[
  //!   \left[\begin{array}{cc}
  //!   \displaystyle  A_1 & \displaystyle B_1 \\
  //!   \displaystyle   -  & \displaystyle D_2
  //!   \end{array}\right]
  //!   + \beta
  //!   \left[\begin{array}{cc}
  //!   \displaystyle   -   & \displaystyle  - \\
  //!   \displaystyle   C_2 & \displaystyle D_2
  //!   \end{array}\right]
  //!   =
  //!   \left[\begin{array}{cc}
  //!   \displaystyle     A_1     & \displaystyle  B_1 \\
  //!   \displaystyle   \beta C_2 & \displaystyle D_1 + \beta D_2
  //!   \end{array}\right]
  //!   \f]
  //!   will be treated correctly and new sub-matrices (in this example
  //!   \f$\beta D_2\f$ will automatically be generated, if necessary.
  //!   However, it is not allowed to add symmetric and unsymmetric matrices
  //!   and non-zero sub-matrices must have identical storage layout,
  //!   dimension and sparsity pattern.
  //!
  //! - The %SBM_Matrix class does not support sub-matrices with different
  //!   entry types and/or different block sizes on the scalar level.
  class SBM_Matrix : public BaseMatrix {

  public:

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION AND SETUP
    // =======================================================================

    //! \name Construction, Destruction and Setup
    //! The following methods are related to construction and destruction of
    //! SBM_Matrix objects and to the initialisation of sub-matrices.
    //@{

    //! This is the allowed constructor. It requires the matrix dimensions

    //! The default constructor is not allowed. Instead this constructor
    //! should be used. It expects to receive the number nrows of rows and
    //! ncols of columns to setup the internal data array.
    //! \param nrows    number of sub-matrix blocks per matrix row
    //! \param ncols    number of sub-matrix blocks per matrix column
    //! \param amSymm   flag setting symmetry of the SBM_Matrix
    SBM_Matrix( UInt nrows, UInt ncols, bool amSymm = false );

    //! Destructor

    //! The destructor is deep, i.e. it deletes all dynamically allocated
    //! objects as well as the sub-matrices
    ~SBM_Matrix() {
      if( ownSubMatrices_ ) {
        for ( UInt k = 0; k < subMat_.GetSize(); k++ ) {
          delete subMat_[k];
        }
      }
    }
    
    //! Copy Constructor
    
    //! This method creates a deep copy of the original matrix, i.e. a copy of
    //! all the sub-matrices.
    SBM_Matrix( const SBM_Matrix& origMat );

    
    //! Create a weak copy of the SBM-matrix
    
    //! This method creates a weak copy of #numRows and #numCols of the
    //! original matrix.
    //! This weak copy will not get ownership of the sub-matrices, but just
    //! get hold of the pointers. Thus, the weak copy is also not allowed
    //! to set new sub-matrices.
    //! \param origMat matrix to be copied
    //! \param numRows number of rows of the sub-matix
    //! \param numCols number of cols of the sub-matrx
    SBM_Matrix( const SBM_Matrix& origMat, UInt numRows, UInt numCols );
    
    //! Insert a StdMatrix into the SBM_Matrix

    //! The method inserts a standard StdMatrix into the super-block-matrix
    //! as the entry in row i and column j. The sub-matrix is newly generated
    //! by calling the GenerateStdMatrixObject() method of the matvec factory.
    //! If there exists already an entry at the position (i,j), then the
    //! method will report an error. Besides this, the method performs several
    //! other consistency checks to avoid generating an inconsistent matrix,
    //! it ensures that
    //! - all sub-matrices have identical block sizes
    //! - all sub-matrices have identical entry types
    //! - all sub-matrices in a block row/column have the same row/column
    //!   dimension
    //! - discards attempts to set sub-matrices in the lower triangle of a
    //!   symmetric %SBM_Matrix issuing a warning
    //!
    //! \param i                row index for inserting sub-matrix
    //! \param j                column index for inserting sub-matrix
    //! \param eType            entry type of sub-matrix
    //! \param sType            storage format of sub-matrix
    //! \param nrows            number of rows of sub-matrix
    //! \param ncols            number of columns of sub-matrix
    //! \param nnz              number of non-zeros in sub-matrix
    //! \param forceOverwrite   forces a matrix overwrite - use carefully!!!
    //! \param silentOVerwrite  no warning will be issued when the matrix is overwritten - use even more carefully!!!
    void SetSubMatrix( UInt i, UInt j, BaseMatrix::EntryType eType,
                       BaseMatrix::StorageType sType, 
                       UInt nrows, UInt ncols, UInt nnz, 
                       bool forceOverwrite = false, bool silentOverWrite = false );

    //! Initialize matrix to zero
    
    //! This matrix sets all matrix entries to zero. Note that the sparsity
    //! structure of the sparse sub-matrices is not changed by this method.
    void Init();

    //@}


    // =======================================================================
    // QUERY GENERAL MATRIX INFORMATION
    // =======================================================================

    //@{
    //! \name Query general matrix format information

    //! Get the number of matrix columns

    //! The method returns the number of columns of the matrix. Note that this
    //! count is based upon the stored entry type. So e.g. for a 2x3 block
    //! matrix the returned value will be 3 independent of the size of the
    //! blocks.
    UInt GetNumCols() const {
      return ncols_;
    }

    //! Get the number of matrix rows

    //! The method returns the number of rows of the matrix. Note that this
    //! count is based upon the stored entry type. So e.g. for a 2x3 block
    //! matrix the returned value will be 2 independent of the size of the
    //! blocks.
    UInt GetNumRows() const {
      return nrows_;
    }

    //! Returns the structural type of the matrix

    //! The method returns the structural type of the matrix. This is encoded
    //! as a value of the enumeration data type StructureType. For an
    //! SBM_Matrix the return value is of course SBM_MATRIX.
    BaseMatrix::StructureType GetStructureType() const {
      return SBM_MATRIX;
    }

    BaseMatrix::StorageType GetStorageType() const {
        WARN("what should be the storage type of a SBM_Matrix?") // FIXME: define StorageType of SBM_Matrix
        return NOSTORAGETYPE;
    }

    //! Return the Entry type of the matrix

    //! The method returns the entry type of the matrix on the scalar level.
    //! This is encoded as a value of the enumeration data type
    //! MatrixEntryType. If no sub-matrix has yet been set, NOENTRYTYPE is
    //! returned.
    BaseMatrix::EntryType GetEntryType() const {
      return myEntryType_;
    }
    
    //! \copydoc BaseMatrix::GetMemoryUsage()
    Double GetMemoryUsage() const;

    //! Return symmetry flag of matrix object

    //! This method can be used to query the symmetry flag of the matrix.
    //! The latter is set during instantiation of the %SBM_Matrix object
    //! by the corresponding constructor argument.
    bool IsSymmetric() {
      return amSymm_;
    }

    //! Determine maximum absolute value of diagonal entries

    //! This method determines the the maximal absolute value (on the scalar)
    //! level of the entries on the main diagonal of the matrix.
    Double GetMaxDiag() const;
    //@}


    // =======================================================================
    // I/O operations
    // =======================================================================

    //! \name I/O operations
    //@{
    //! Return matrix as separated string
    std::string ToString( char colSeparator = ' ',
    char rowSeparator = '\n' ) const {

      std::stringstream os;
      for( UInt i = 0; i < nrows_; i++ ) {
        for( UInt j = 0; j < ncols_; j++ ) {

          if( subMat_[ComputeIndex(i,j)] != NULL ) {
            os <<   "sub-matrix #" << i << " - " << j << "\n";
            os << "\n--------------\n";
            os <<  subMat_[ComputeIndex(i,j)]->ToString(colSeparator, rowSeparator );
          } else {
            os <<   "sub-matrix #" << i << " - " << j << "\n";
            os << "\n-EMPTY-\n";
          }
        }
      }
      return os.str();

      //      WARN( "ToString not implemented for SBM_Matrix");
      //
      //      return "SBM_Matrix";
    }

    /** stupid debug info for developers to understand the stuff better */
    std::string Dump() const;

    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. For details of the specification see
    //! http://math.nist.gov/MatrixMarket
    //! \param fname used for building name of output files (see below)
    //! \param comment string to be inserted into file header (optional)
    //! \note The method does currently export each sub-matrix to a separate
    //!       file. The file names are constructed by appending the suffix
    //!       <em>_i_j.mtx</em> to fname with (i,j) replaced by the index
    //!       tuple for the sub-matrix.
    virtual void Export( const std::string& fname,
                         BaseMatrix::OutputFormat format,
                         const std::string& comment = "" ) const;

    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. For details of the specification see
    //! http://math.nist.gov/MatrixMarket
    //! \param fname used for building name of output files (see below)
    //! \param N     number of harmonics
    //! \param sbmInd  indices of sbm-Blocks to be exported
    //! \param comment string to be inserted into file header (optional)
    //! \note The method does currently export each sub-matrix to a separate
    //!       file. The file names are constructed by appending the suffix
    //!       <em>_i_j.mtx</em> to fname with (i,j) replaced by the index
    //!       tuple for the sub-matrix.
    virtual void Export_MultHarm( const std::string& fname,
                                  BaseMatrix::OutputFormat format,
                                  const UInt& N,
                                  const StdVector<UInt>& sbmInd,
                                  const std::string& comment = "") const;


    //@}


    // =======================================================================
    // ACCESS TO SUB-MATRICES
    // =======================================================================

    //! \name Methods to access sub-matrices

    //@{
    //! Overload () operator to retrieve submatrices

    //! The () operator is overloaded to allow to retrieve individual
    //! submatrices. The access follows the one-based indexing convention
    //! of OLAS.
    StdMatrix& operator()( UInt i, UInt j ) const {
      return *(subMat_[ComputeIndex(i,j)]);
    }

    //! Obtain pointer to a submatrix

    //! The method returns a pointer to a sub-matrix. The access follows the
    //! one-based indexing convention of OLAS.
    const StdMatrix* GetPointer( UInt i, UInt j ) const {
      return (subMat_[ComputeIndex(i,j)]);
    }

    //! Obtain pointer to a submatrix

    //! The method returns a pointer to a sub-matrix. The access follows the
    //! one-based indexing convention of OLAS.
    StdMatrix* GetPointer( UInt i, UInt j ) {
      return (subMat_[ComputeIndex(i,j)]);
    }

    //@}


    // =======================================================================
    // MATRIX-VECTOR PRODUCTS
    // =======================================================================

    //! \name Matrix-Vector products 

    //! The following methods perform arithmetic operations of the form
    //! \f$ \vec{y} = \alpha \vec{y} + \beta A \vec{x}\f$ or similar

    //@{

    //! Perform a matrix-vector multiplication

    //! This method will multiply the vector mvec with this matrix object
    //! and return the result in rvec. For this to work both vectors must
    //! correspond to the layout of the super-block-matrix.
    void Mult( const BaseVector& mvec, BaseVector& rvec ) const;

    //! Perform a matrix-vector multiplication

    //! This method will multiply the vector mvec with this matrix object
    //! and add the result to rvec. For this to work both vectors must
    //! correspond to the layout of the super-block-matrix.
    void MultAdd( const BaseVector& mvec, BaseVector& rvec ) const;

    //! Perform a matrix-vector multiplication

    //! This method will multiply the vector mvec with this matrix object
    //! and subtract the result from rvec. For this to work both vectors must
    //! correspond to the layout of the super-block-matrix.
    void MultSub( const BaseVector& mvec, BaseVector& rvec ) const;

    //! Compute the residual for a linear system with this matrix

    //! The method computes the residual \f$r=b-Ax\f$ where \f$A\f$ is the
    //! matrix represented by this matrix object.
    void CompRes( BaseVector &r, const BaseVector &x,
			  const BaseVector& b ) const;

    //@}


    // =======================================================================
    // MATRIX-VECTOR-PRODUCT WITH TRANSPOSE
    // =======================================================================

    //! \name Matrix-vector multiplication with transpose of matrix

    //! The following methods perform arithmetic operations of the form
    //! \f$ \vec{y} = \alpha \vec{y} + \beta A^T \vec{x}\f$ or similar

    //@{

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object rvec = transpose(this)*mvec.
    void MultT(const BaseVector& mvec, BaseVector& rvec) const;

    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object followed by an addition:
    //! rvec += transpose(this)*mvec.
    void MultTAdd(const BaseVector& mvec, BaseVector& rvec) const;

    //@}


    // =======================================================================
    // MISCELLANEOUS ARITHMETIC OPERATIONS
    // =======================================================================

    //! \name Miscellaneous arithmetic matrix operations
    //@{

    //! Add a scalar multiple of another matrix to this matrix object

    //! This method will add the matrix mat multiplied by the factor fac to
    //! this matrix object. It is safe to add %SBM_Matrices with different
    //! block-sparsity patterns, i.e.
    //!   \f[
    //!   \left[\begin{array}{cc}
    //!   \displaystyle  A_1 & \displaystyle B_1 \\
    //!   \displaystyle   -  & \displaystyle D_2
    //!   \end{array}\right]
    //!   + \beta
    //!   \left[\begin{array}{cc}
    //!   \displaystyle   -   & \displaystyle  - \\
    //!   \displaystyle   C_2 & \displaystyle D_2
    //!   \end{array}\right]
    //!   =
    //!   \left[\begin{array}{cc}
    //!   \displaystyle     A_1     & \displaystyle  B_1 \\
    //!   \displaystyle   \beta C_2 & \displaystyle D_1 + \beta D_2
    //!   \end{array}\right]
    //!   \f]
    //!   will be treated correctly and new sub-matrices (in this example
    //!   \f$\beta D_2\f$ will automatically be generated, if necessary.
    //!   However, it is not allowed to add symmetric and unsymmetric matrices
    //!   and non-zero sub-matrices must have identical storage layout,
    //!   dimension and sparsity pattern.
    void Add( const Double fac, const BaseMatrix& mat );

    
    //! Add a scalar multiple of entries of another matrix
    
    //! This method is similar to the Add(Double, mat) method, but only some
    //! entries are considered, specified by SBM-matrix block and an index set.
    //! If an entry can not be found or the matrix is
    //! not quadratic, an exception is thrown.
    //! If the row / col index set is empty, all rows / cols are considered.
    template <typename U>
    void Add( const U fac, const BaseMatrix& mat,
              std::map<UInt, std::set<UInt> >& rowIndPerBlock,
              std::map<UInt, std::set<UInt> >& colIndPerBlock );
    //@}

    //! This array contains pointers to the entries, i.e. the standard
    //! matrices, of the dense super-block-matrix. The storage follows the
    //! standard C layout, but the submat_ vector is 1-based.
    StdVector<StdMatrix *> subMat_;
  private:

    //! Default constructor

    //! The default constructor is not allowed. Pass the dimension of the
    //! super-block-matrix to the constructor instead.
    SBM_Matrix();

    //! Number of rows in the super-block matrix
    UInt nrows_;

    //! Number of columns in the super-block matrix
    UInt ncols_;

    //! This array contains pointers to the sub-matrices



    //! Computes the index of sub-matrix in 1D array submat_

    //! The class stores the pointers to the sub-matrices of the super block
    //! matrix in a 'virtual' 2D matrix. The latter is stored as 1D array.
    //! This method offers a way to compute for a given sub-matrix index
    //! pair (i,j) the corresponding index into the one-based 1D array.
    inline UInt ComputeIndex( UInt i, UInt j ) const {
      return ncols_ * i + j;
    }


    //! Computes the row and column of a given index
    inline StdVector<UInt> DeflattenIndex(UInt ind) const {
      UInt nr = (nrows_ - 1)/2;
      std::vector<UInt> a = {ind / (2*nr+1), ind % (2*nr+1)};
      return a;
    }

    //! Flag signalling symmetry of the matrix

    //! This boolean attribute is used as a flag to signal the symmetry of
    //! the matrix stored in super-block matrix format. In this case only
    //! the upper-triangular sub-matrix blocks are stored and the diagonal
    //! sub-matrices are assumed to be symmetric (though they need not be
    //! stored in a symmetric format).
    bool amSymm_;
    
    //! Flag indicating if object is responsible for deletion of submatrices
    bool ownSubMatrices_;

    //! Attribute storing the entry type of the matrix

    //! Currently the %SBM_Matrix class only allows for storing sub-matrices
    //! that share a common entry type. This attribute stores the entry type
    //! of the first sub-matrix that is set and uses this as the over-all
    //! entry type of the %SBM_Matrix object.
    BaseMatrix::EntryType myEntryType_;

  };

}

#endif
