#ifndef OLAS_SPARSE_MATRIX_HH
#define OLAS_SPARSE_MATRIX_HH

#include <iostream>
#include <set>

#include "TypeDefs.hh"
#include "BaseMatrix.hh"
#include "SingleVector.hh"
#include "OLAS/graph/BaseGraph.hh"


namespace CoupledField {


  //! Forward declarations
  class PatternPool;
  class DenseMatrix;

  //! base class for unstructured matrices
  class StdMatrix : public BaseMatrix {

  public:
    
    virtual ~StdMatrix() {}

    // ========================================================================
    // QUERY METHODS
    // ========================================================================

    //@{ \name Query methods to obtain information on the matrix

    //! Returns the structural type of the matrix

    //! The method returns the structural type of the matrix. This is encoded
    //! as a value of the enumeration data type StructureType. For an
    //! StdMatrix the return value is of course SPARSE_MATRIX.
    BaseMatrix::StructureType GetStructureType() const {
      return BaseMatrix::SPARSE_MATRIX;
    };


    //! Returns the parallelisation type of the matrix

    //! The method returns the parallelisation type of the matrix. This is
    //! encoded as a value of the enumeration data type ParMatrixType. The
    //! default value is to return SEQMATRIX, which describes a sequential
    //! matrix object. Thus, this method must be over-written by all parallel
    //! matrix classes. So that they can give back a meaningful value.
    virtual ParMatrixType GetParMatrixType() const {
      return SEQMATRIX;
    };


    //! Return the Entry type of the matrix

    //! The method returns the entry type of the matrix (i.e. Double,
    //! Complex ..). This is encoded as a value of the enumeration data type
    //! MatrixEntryType.
    virtual BaseMatrix::EntryType GetEntryType() const {
      EXCEPTION( "Not implemented here" );
      return BaseMatrix::NOENTRYTYPE;
    };


    //! Return the storage type of the matrix
    
    //! The method returns the storagfe type of the matrix. This is encoded
    //! as a value of the enumeration data type MatrixStorageType.
    virtual BaseMatrix::StorageType GetStorageType() const {
      EXCEPTION("Not implemented here");
      return NOSTORAGETYPE;
    };

    //! Get the number of matrix rows

    //! The method returns the number of rows of the matrix. Note that this
    //! count is based upon the stored entry type.
    UInt GetNumRows() const {
      return nrows_;
    };

    //! Get the number of matrix columns

    //! The method returns the number of columns of the matrix. Note that this
    //! count is based upon the stored entry type.
    UInt GetNumCols() const {
      return ncols_;
    };

    //! Get the number of non-zero entries

    //! The method returns the number of non-zero entries. Note that this
    //! value is based upon the stored entry type. In the case of block entries
    //! the latter are counted and the number of entries on the scalar level
    //! will actually be higher depending on the block size.
    //! \note The implementation of this method in the derivated class
    //!       ParSSCRS_Matrix might cause dead locks in parallel programs
    //!       and must be handled with care. See ParSSCRS_Matrix::GetNnz
    virtual UInt GetNnz() const {
      EXCEPTION( "Method not implemented by derived class" );
      return 0;
    }

    //! Get pointer to sparsity pattern pool 

    //! The method returns the pointer to the pattern pool, from which the
    //! matrix got its sparsity pattern. If the matrix was created by using
    //! the graph object, the return pointer will be NULL.
    PatternPool* GetPatternPool() const {
      EXCEPTION( "Method not implemented by derived class" );
      return NULL;
    }

    //! Get id of the sparsity pattern of the pattern pool

    //! This method return the id for the sparsity pattern which is used
    //! for communication with the pattern pool. If no pattern pool was used
    //! for the creation of the matrix, NO_PATTERN_ID is returned.
    PatternIdType GetPatternId() const {
      EXCEPTION( "Method not implemented by derived class" );
      return NO_PATTERN_ID;
    }

    //! Export the matrix to a file in MatrixMarket or Harwell-Boeing format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket or Harwell-Boeing specifications. This method will
    //! redirect the calls to the corresponding implementations in the sub
    //! classes.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! and http://people.sc.fsu.edu/~jburkardt/data/hb/hb.html
    //! \param fname base name of output file
    //! \param format matrix output format 1 = MatrixMarket, 2 = Harwell-Boeing (optional)
    //! \param comment string to be inserted into file header (optional)
    virtual void Export( const std::string& fname,
                         OutputFormat format,
                         const std::string& comment = "") const;

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket and is only implemented in sub classes of StdMatrix.
    virtual void ExportMatrixMarket( const std::string& fname, const std::string& comment = "") const = 0;

    /** Exports in Harwell-Boeing format. This format is used for the
   * ILUPACK frontend. A benefit is the included RHS, inital guesses, ...
     * are not implemented here. */  
    void ExportHarwellBoeing(const std::string& file, BaseVector* rhs) const; 

    //@}

    /** Helper to export the columns for a Compressed row Storage (CRS -> wikipedia).
     * Alternatively use e.g. SCRS_Matrix::GetColPointer()
     * @param out will be resized to nnz. Int because this is standard for Fortran?! Size
     * @param base is either 0 or 1. 0 keeps this 0-based structure, 1 makes it 1-based */
    void ExportCRSColumns(StdVector<int>& out, int base) const;

    /** @see ExportCRSColumns()
     * @param out resized to num-rows. If tailing_nnz the nnz are added
     * @param base see ExportCRSColumns()
     * @param tailing_nnz shall nnz be added? */
    void ExportCRSRows(StdVector<int>& out, int base, bool tailing_nnz) const;
      

    // ========================================================================
    // ARITHMETIC METHODS
    // ========================================================================

    //@{

    //! \name Methods for arithmetic operations
    //! These methods try to down-cast the input vector / matrix to the
    //! SingleVector / StdMatrix class and, if this succeeds, call the respective
    //! method of the derived class, the latter is forced to implement it,
    //! since we specify a purely virtual interface in this class.

    //! Compute the residual for a linear system with this matrix

    //! The method computes the residual \f$r=b-Ax\f$ where \f$A\f$ is the
    //! matrix represented by this matrix object. This method in fact only
    //! downcasts the BaseVectors to SingleVectors and delegates the work to
    //! the method with the appropriate interface. It implements the method
    //! defined in the BaseMatrix class.
    void CompRes( BaseVector &r, const BaseVector &x,
                          const BaseVector& b ) const {
      const SingleVector& std_x = dynamic_cast<const SingleVector&>(x);
      const SingleVector& std_b = dynamic_cast<const SingleVector&>(b);
      SingleVector& std_r = dynamic_cast<SingleVector&>(r);

      CompRes( std_r, std_x, std_b );
    }

    //! Compute the residual for a linear system with this matrix

    //! The method computes the residual \f$r=b-Ax\f$ where \f$A\f$ is the
    //! matrix represented by this matrix object. This is the variant of the
    //! method with the appropriate interface for SingleVectors.
    virtual void CompRes( SingleVector &r, const SingleVector &x,
                          const SingleVector& b ) const = 0;

    //! Add the multiple of a matrix to this matrix.

    //! The method adds the multiple of a matrix to this matrix. In doing so
    //! the sparsity structure of the matrix mat is assumed to be identical to
    //! this matrix' structure. However, all it does is downcast the
    //! BaseMatrices to StdMatrices and delegate the work to the method with
    //! the appropriate interface. It implements the method defined in the
    //! BaseMatrix class.
    void Add( const Double a, const BaseMatrix& mat )
    {
      Add( a, dynamic_cast<const StdMatrix&>(mat) );
    };

    //! Add the multiple of a matrix to this matrix.

    //! The method adds the multiple of a matrix to this matrix. In doing so
    //! the sparsity structure of the matrix mat is assumed to be identical to
    //! this matrix' structure. It offers an appropriate interface for
    //! StdMatrices.
    virtual void Add( const Double a, const StdMatrix& mat) = 0;
    
    //! Add the multiple of a matrix to this matrix (entries given by index)

    //! The method adds the multiple of a subset of the of a matrix to 
    //! this matrix. The entries are specified in a row and column index set.
    //! If a set is empty, all rows / cols are considered. 
    //! If an entry can not be found in the matrix, an exception is thrown.
    //! The sparsity structure of the matrix mat is assumed to be identical to
    //! this matrix' structure. It offers an appropriate interface for
    //! StdMatrices.
    //! \param a scalar factor the matrix #mat gets weighted with upon addition 
    //! \param mat matrix which gets added
    //! \param rowIndices set of row indices to be considered. If the set is
    //!                   empty, all rows will be taken.
    //! \param colIndices set of row indices to be considered. If the set is
    //!                   empty, all rows will be taken.
    virtual void Add( const Double a, const StdMatrix& mat,
                      const std::set<UInt>& rowIndices,
                      const std::set<UInt>& colIndices ) = 0;

    //! Perform a matrix-vector multiplication rvec = this*mvec

    //! This method performs a matrix-vector multiplication rvec = this*mvec.
    //! However, all it does is downcast the BaseVectors to SingleVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void Mult(const BaseVector& mvec, BaseVector& rvec) const
    {
      Mult(dynamic_cast<const SingleVector&>(mvec), 
           dynamic_cast<SingleVector&>(rvec));
    }

    //! Perform a matrix-vector multiplication rvec = this*mvec

    //! This method performs a matrix-vector multiplication rvec = this*mvec.
    //! It offers an appropriate interface for SingleVectors.
    virtual void Mult(const SingleVector& mvec, SingleVector& rvec) const = 0;

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object rvec = transpose(this)*mvec.
    //! However, all it does is downcast the BaseVectors to SingleVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void MultT(const BaseVector& mvec, BaseVector& rvec) const
    {
      MultT(dynamic_cast<const SingleVector&>(mvec), dynamic_cast<SingleVector&>(rvec));
    }

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object rvec = transpose(this)*mvec.
    //! It offers an appropriate interface for SingleVectors.
    virtual void MultT(const SingleVector& mvec, SingleVector& rvec) const {
      EXCEPTION( "Function not re-implemented in derived class!" );
    }

    //! Perform a matrix-vector multiplication rvec += this*mvec

    //! This method performs a matrix-vector multiplication rvec += this*mvec.
    //! However, all it does is downcast the BaseVectors to SingleVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    virtual void MultAdd(const BaseVector& mvec, BaseVector& rvec) const
    {
      MultAdd(dynamic_cast<const SingleVector&>(mvec), dynamic_cast<SingleVector&>(rvec));
    }
        
    //! Perform a matrix-vector multiplication rvec += this*mvec

    //! This method performs a matrix-vector multiplication rvec += this*mvec.
    //! It offers an appropriate interface for SingleVectors.
    virtual void MultAdd( const SingleVector& mvec, SingleVector& rvec ) const = 0;

    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object followed by an addition:
    //! rvec += transpose(this)*mvec.
    //! However, all it does is downcast the BaseVectors to SingleVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void MultTAdd(const BaseVector& mvec, BaseVector& rvec) const
    {
      MultTAdd(dynamic_cast<const SingleVector&>(mvec), dynamic_cast<SingleVector&>(rvec));
    }

    virtual void MultTAdd(const SingleVector& mvec, SingleVector& rvec) const {
      EXCEPTION( "Function MultTAdd not re-implemented in derived class!" );
    }

    //! Perform a matrix-vector multiplication rvec -= this*mvec

    //! This method performs a matrix-vector multiplication rvec -= this*mvec.
    //! However, all it does is downcast the BaseVectors to SingleVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void MultSub( const BaseVector& mvec, BaseVector& rvec ) const
    {
      MultSub(dynamic_cast<const SingleVector&>(mvec), dynamic_cast<SingleVector&>(rvec));
    }
        
    //! Perform a matrix-vector multiplication rvec -= this*mvec

    //! This method performs a matrix-vector multiplication rvec -= this*mvec.
    //! It offers an appropriate interface for SingleVectors.
    virtual void MultSub( const SingleVector& mvec, SingleVector& rvec ) const = 0;

    
    //! Perform a matrix-vector multiplication rvec -= transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object followed by a subtraction:
    //! rvec -= transpose(this)*mvec.
    //! However, all it does is downcast the BaseVectors to SingleVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void MultTSub(const BaseVector& mvec, BaseVector& rvec) const
    {
      MultTSub(dynamic_cast<const SingleVector&>(mvec), dynamic_cast<SingleVector&>(rvec));
    }

    virtual void MultTSub(const SingleVector& mvec, SingleVector& rvec) const {
      EXCEPTION( "Function MultTSub not re-implemented in derived class!" );
       }
    
    //! \copydoc BaseMatrix::Scale(Double)
    virtual void Scale( Double factor ) {
      EXCEPTION("BaseMatrix::Scale: Method must be implemented by derived " \
                "class, but is not!");
    };
    
    //! \copydoc BaseMatrix::Scale(Double)
    virtual void Scale( Complex factor ) {
      EXCEPTION("BaseMatrix::Scale: Method must be implemented by derived " \
                "class, but is not!");
    };

    //! Scale matrix entries - given by an index set - by a constant factor

    //! This method allows to scale all non-zero matrix entries, which are
    //! contained in an index set, by a constant real valued factor. 
    //! The operation assumes to work
    //! on a quadratic matrix, i.e. the indices denote the  and column
    //! entries of the matrix. If an entry can not be found or the matrix is
    //! not quadratic, an exception is 
    //! thrown. It is currently pseudo-implemented in this class,
    //! since not all derived classes implement it, yet. The pseudo
    //! implementation will issue an error, if not over-written.
    //! \param factor scalar the matrix gets multiplied with
    //! \param rowIndices set of row indices to be considered. If the set is
    //!                   empty, all rows will be taken.
    //! \param colIndices set of row indices to be considered. If the set is
    //!                   empty, all rows will be taken.
    virtual void Scale( Double factor, 
                        const std::set<UInt>& rowIndices,
                        const std::set<UInt>& colIndices) {
      EXCEPTION("StdMatrix::Scale: Method must be implemented by derived " \
                "class, but is not!");
    };
    //@}
    
    // ========================================================================
    // BLOCK STRUCTURE RELATED METHODS
    // ========================================================================
    //@{
    //! \name Block Structure Related Methods
    
    //! Get dimension of matrix blocks

    //! This method returns the dimension of the matrix in terms of blocks.
    //! \param nRowBlocks number of row blocks
    //! \param nColBlocks number of column blocks
    //! \param numBlocks total number of (nonzero) blocks
    virtual void GetNumBlocks( UInt& nRowBlocks, UInt& nColBlocks, 
                               UInt& numBlocks ) const {
      EXCEPTION( "Not implemented in the base class" );
    }
        
    
    //! Return the diagonal block
    
    //! This method returns the dense diagonal matrix block in case block 
    //! information is supported by the matrix structure and information is 
    //! provided by the graph during setup.
    virtual void GetDiagBlock( UInt blockRow, DenseMatrix& diagBlock ) const {
      EXCEPTION( "Not implemented in the base class" );
    }
    
    //! Set the diagonal block

    //! This method sets the dense diagonal matrix block in case block 
    //! information is supported by the matrix structure and information is 
    //! provided by the graph during setup.
    virtual void SetDiagBlock( UInt blockRow, const DenseMatrix& diagBlock ) {
      EXCEPTION( "Not implemented in the base class" );
    }

    //@}
    
    // *****************************************************
    //   The following two methods can probably be removed
    // *****************************************************

    //! Setup the sparsity pattern of the matrix

    //! This method uses the graph to generate the internal information needed
    //! for representing the sparsity pattern of the matrix correctly. Before
    //! this function has not been called the matrix is not ready to be used.
    virtual void SetSparsityPattern( BaseGraph &graph ) = 0;

    //! Setup the sparsity pattern of the matrix
  
    //! This method provides an alternative form to set the sparsity pattern
    //! of the matrix. If this method is used a fitting sparsity pattern must
    //! already have been generated, e.g. by another instance of this class
    //! that has the same pattern, and been added to the specified pool. In
    //! this case we can simply obtain copy the address information from the
    //! pool object.
    //! \param patternPool  pointer to a PatternPool object that stores the
    //!                     desired sparsity pattern
    //! \param patternID    identifier required to obtain the sparsity
    //!                     pattern from the pool of patterns
    virtual void SetSparsityPattern( PatternPool *patternPool,
                                     PatternIdType patternID ) {
      EXCEPTION("Sorry, but SetSparsityPattern is only implemented "
               << "for the SCRS_Matrix class currently!");
    }

    //! Transfer ownership of sparsity pattern from matrix to pool

    //! This method can be used to instruct the matrix to transfer the
    //! ownership of its sparsity pattern to a PatternPool object such
    //! that the sparsity pattern can be re-used by other instances of
    //! the class. After a call to the method the matrix object will no
    //! longer attempt to de-allocate memory for the pattern in its
    //! destructor, but instead only de-register itself from the pool.
    //! \param patternPool pointer to the PatternPool object to which the
    //!        ownership of the sparsity pattern is to be transfered.
    //! \return An identifier that can be used to identify the matrix
    //!         pattern when communicating with the PatternPool object.
    virtual PatternIdType TransferPatternToPool( PatternPool *patternPool ) {
      EXCEPTION("Sorry, but TransferPatternToPool is only implemented "
               << "for the SCRS_Matrix class currently!");
      return NO_PATTERN_ID;
    }

    //! Copy the sparsity pattern of another matrix to this matrix 
    
    //! This method copies the sparsity pattern of the source matrix to the
    //! current matrix
    virtual void SetSparsityPattern ( const StdMatrix &srcMat ) {
      EXCEPTION( "StdMatrix::SetSparsityPattern(): Not implemented here" );
    }

    //! Set a specific matrix entry

    /*! this functions sets the Matrix entry in row i and column j
      to the value v. Depending on the entry type of the matrix,
      v might be a Double, Complex or a tiny Matrix of either type
    */
    virtual void SetMatrixEntry( UInt i, UInt j, const Double &v, 
                                 bool setCounterPart ) {
      EXCEPTION( "StdMatrix::SetMatrixEntry(): Not implemented here" );
    }

    //! Get a specific matrix entry

    /*! this function returns a reference to the Matrix entry in 
      row i and column j. Depending on the entry type of the matrix,
      v might be a Double, Complex or a tiny Matrix of either type
    */
    virtual void GetMatrixEntry( UInt i, UInt j, Double &v ) const{
      EXCEPTION( "StdMatrix::GetMatrixEntry(): Not implemented here" );
    }

    //! gets entries of system matrix as a vector

    /*! copies the system matrix from SCRS format into a vector (complex)
      containing all (even all zeros) upper triangle elements
      is called by method GetFullSystemMatrixAsVec in standardsys which 
      provides an interface to openCFS
    */
    virtual void CopySCRSMatrix2Vec(Complex* &vec){;};
 

    //! Add value to a matrix entry

    //! This method adds the value v to the matrix entry at position (i,j).
    //! The method is documented here for the case of Double scalar entries,
    //! but is defined for all entry types.
    virtual void AddToMatrixEntry(UInt i, UInt j, const Double &v ) {
      EXCEPTION( "StdMatrix::AddToMatrixEntry(): Not implemented here");
    }

    //!Set the diagonal entry of row i to v
    virtual void SetDiagEntry(UInt i, const Double &v) {
      EXCEPTION( "StdMatrix::SetDiagEntry(): Not implemented here" );
    }

    //!Set the diagonal entry of row i to v     
    virtual void GetDiagEntry(UInt i, Double &v) const {
      EXCEPTION("StdMatrix::GetDiagEntry(): Not implemented here");
    }

    //! slow direct matrix access
#define DECL_SPARSE_MATRIX_FCN(TYPE) \
virtual void SetMatrixEntry(UInt i, UInt j, const TYPE &v,\
bool setCounterPart){\
EXCEPTION("StdMatrix::SetMatrixEntry(): Not implemented here");}\
\
virtual void GetMatrixEntry(UInt i, UInt j, TYPE &v) const {\
EXCEPTION("StdMatrix::GetMatrixEntry(): Not implemented here");}\
\
virtual void AddToMatrixEntry(UInt i, UInt j, const TYPE &v ) { \
EXCEPTION("StdMatrix::AddToMatrixEntry(): Not implemented here"); }\
\
virtual void SetDiagEntry(UInt i, const TYPE &v){ \
EXCEPTION("StdMatrix::SetDiagEntry(): Not implemented here");}\
\
virtual void GetDiagEntry(UInt i, TYPE &v) const { \
EXCEPTION("StdMatrix::GetDiagEntry(): Not implemented here");}

DECL_SPARSE_MATRIX_FCN(Complex);

  protected:

    //! Number of matrix rows
    UInt nrows_;

    //! Number of matrix columns
    UInt ncols_;

   private:
     /** A small harwell boint export class, templated to handel double and complex
      * easier */
    template<typename T>
    class HarwellBoeing
    {
        public:
        HarwellBoeing(const StdMatrix* matrix);
        
        void Export(const std::string& file, BaseVector* rhs);
        
        private:
        void Fill(const T* values, const UInt* ints, UInt length, std::vector<std::string>& out);
        
        bool is_complex_;
        const StdMatrix* matrix_;
    };   

  };


} // namespace



#endif // OLAS_BASEMATRIX_HH
