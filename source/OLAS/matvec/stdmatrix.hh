// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_STDMATRIX_HH
#define OLAS_STDMATRIX_HH


#include <iostream>

#include "matvec/typedefs.hh"
#include "utils/utils.hh"
#include "matvec/basematrix.hh"
#include "matvec/SparseVector.hh"
#include "graph/basegraph.hh"


namespace OLAS {


  //! Forward declarations
  class PatternPool;


  //! base class for unstructured matrices
  class StdMatrix : public BaseMatrix {

  public:

    // ========================================================================
    // QUERY METHODS
    // ========================================================================

    //@{ \name Query methods to obtain information on the matrix

    //! Returns the structural type of the matrix

    //! The method returns the structural type of the matrix. This is encoded
    //! as a value of the enumeration data type StructureType. For an
    //! StdMatrix the return value is of course STDMATRIX.
    MatrixStructureType GetStructureType() const {
      return STDMATRIX;
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
    virtual MatrixEntryType GetEntryType() const {
      Error( "Not implemented here", __FILE__, __LINE__ );
      return NOENTRYTYPE;
    };


    //! Return the storage type of the matrix
    
    //! The method returns the storagfe type of the matrix. This is encoded
    //! as a value of the enumeration data type MatrixStorageType.
    virtual MatrixStorageType GetStorageType() const {
      Error("Not implemented here", __FILE__, __LINE__);
      return NOSTORAGETYPE;
    };

    //! Get the number of matrix rows

    //! The method returns the number of rows of the matrix. Note that this
    //! count is based upon the stored entry type.
    Integer GetNrows() const {
      return nrows_;
    };

    //! Get the number of matrix columns

    //! The method returns the number of columns of the matrix. Note that this
    //! count is based upon the stored entry type.
    Integer GetNcols() const {
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
    virtual Integer GetNnz() const {
      (*error) << "Method not implemented by derived class";
      Error( __FILE__, __LINE__ );
      return -1;
    }

    //! Get pointer to sparsity pattern pool 

    //! The method returns the pointer to the pattern pool, from which the
    //! matrix got its sparsity pattern. If the matrix was created by using
    //! the graph object, the return pointer will be NULL.
    PatternPool* GetPatternPool() const {
      (*error) << "Method not implemented by derived class";
      Error( __FILE__, __LINE__ );
      return NULL;
    }

    //! Get id of the sparsity pattern of the pattern pool

    //! This method return the id for the sparsity pattern which is used
    //! for communication with the pattern pool. If no pattern pool was used
    //! for the creation of the matrix, NO_PATTERN_ID is returned.
    PatternIdType GetPatternId() const {
     (*error) << "Method not implemented by derived class";
      Error( __FILE__, __LINE__ );
      return NO_PATTERN_ID;
    }

    /** Exports in Harwell-Boeing format. This format is used for the
     * ILUPACK frontend. A benefit is the included RHS, inital guesses, ...
     * are not implemented here. */  
    void ExportHarwellBoeing(const std::string& file, const BaseVector& rhs); 

    //@}

    //! Method to force instantiation of members in templated matrices

    //! This method is used by the GenerateStdMatrix function to force the
    //! instantiation of the public member functions of the templated matrix
    //! classes.
    virtual void InstantiatePublicMethods() {
      Error( "Not implemented here", __FILE__, __LINE__ );
    };


      

    // ========================================================================
    // ARITHMETIC METHODS
    // ========================================================================

    //@{

    //! \name Methods for arithmetic operations
    //! These methods try to down-cast the input vector / matrix to the
    //! SparseVector / StdMatrix class and, if this succeeds, call the respective
    //! method of the derived class, the latter is forced to implement it,
    //! since we specify a purely virtual interface in this class.

    //! Compute the residual for a linear system with this matrix

    //! The method computes the residual \f$r=b-Ax\f$ where \f$A\f$ is the
    //! matrix represented by this matrix object. This method in fact only
    //! downcasts the BaseVectors to SparseVectors and delegates the work to
    //! the method with the appropriate interface. It implements the method
    //! defined in the BaseMatrix class.
    void CompRes( BaseVector &r, const BaseVector &x,
                          const BaseVector& b ) const {
      TRY_CAST
      ConstRefCast( x, SparseVector, std_x );
      ConstRefCast( b, SparseVector, std_b );
      RefCast( r, SparseVector, std_r );
      CompRes( std_r, std_x, std_b );
      CATCH_CAST
    }

    //! Compute the residual for a linear system with this matrix

    //! The method computes the residual \f$r=b-Ax\f$ where \f$A\f$ is the
    //! matrix represented by this matrix object. This is the variant of the
    //! method with the appropriate interface for SparseVectors.
    virtual void CompRes( SparseVector &r, const SparseVector &x,
                          const SparseVector& b ) const = 0;

    //! Add the multiple of a matrix to this matrix.

    //! The method adds the multiple of a matrix to this matrix. In doing so
    //! the sparsity structure of the matrix mat is assumed to be identical to
    //! this matrix' structure. However, all it does is downcast the
    //! BaseMatrices to StdMatrices and delegate the work to the method with
    //! the appropriate interface. It implements the method defined in the
    //! BaseMatrix class.
    void Add( const Double a, const BaseMatrix& mat ){
      TRY_CAST
      ConstRefCast( mat, StdMatrix, stdmat );
      Add( a, stdmat );
      CATCH_CAST
    };

    //! Add the multiple of a matrix to this matrix.

    //! The method adds the multiple of a matrix to this matrix. In doing so
    //! the sparsity structure of the matrix mat is assumed to be identical to
    //! this matrix' structure. It offers an appropriate interface for
    //! StdMatrices.
    virtual void Add( const Double a, const StdMatrix& mat) = 0;

    //! Perform a matrix-vector multiplication rvec = this*mvec

    //! This method performs a matrix-vector multiplication rvec = this*mvec.
    //! However, all it does is downcast the BaseVectors to SparseVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void Mult(const BaseVector& mvec, BaseVector& rvec) const {
      TRY_CAST
      ConstRefCast(mvec,SparseVector,stdmvec);
      RefCast(rvec,SparseVector,stdrvec);
      Mult(stdmvec,stdrvec);
      CATCH_CAST
    };

    //! Perform a matrix-vector multiplication rvec = this*mvec

    //! This method performs a matrix-vector multiplication rvec = this*mvec.
    //! It offers an appropriate interface for SparseVectors.
    virtual void Mult(const SparseVector& mvec, SparseVector& rvec) const = 0;

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object rvec = transpose(this)*mvec.
    //! However, all it does is downcast the BaseVectors to SparseVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void MultT(const BaseVector& mvec, BaseVector& rvec) const {
      TRY_CAST
      ConstRefCast(mvec,SparseVector,stdmvec);
      RefCast(rvec,SparseVector,stdrvec);
      MultT(stdmvec,stdrvec);
      CATCH_CAST
    };

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object rvec = transpose(this)*mvec.
    //! It offers an appropriate interface for SparseVectors.
    virtual void MultT(const SparseVector& mvec, SparseVector& rvec) const {
      Error( "Function not re-implemented in derived class!", __FILE__,
             __LINE__);
    }

    //! Perform a matrix-vector multiplication rvec += this*mvec

    //! This method performs a matrix-vector multiplication rvec += this*mvec.
    //! However, all it does is downcast the BaseVectors to SparseVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    virtual void MultAdd(const BaseVector& mvec, BaseVector& rvec) const {
      TRY_CAST
        ConstRefCast(mvec,SparseVector,stdmvec);
      RefCast(rvec,SparseVector,stdrvec);
      MultAdd(stdmvec,stdrvec);
      CATCH_CAST
    };
        
    //! Perform a matrix-vector multiplication rvec += this*mvec

    //! This method performs a matrix-vector multiplication rvec += this*mvec.
    //! It offers an appropriate interface for SparseVectors.
    virtual void MultAdd( const SparseVector& mvec, SparseVector& rvec ) const = 0;

    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object followed by an addition:
    //! rvec += transpose(this)*mvec.
    //! However, all it does is downcast the BaseVectors to SparseVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void MultTAdd(const BaseVector& mvec, BaseVector& rvec) const {
      TRY_CAST
      ConstRefCast(mvec,SparseVector,stdmvec);
      RefCast(rvec,SparseVector,stdrvec);
      MultTAdd(stdmvec,stdrvec);
      CATCH_CAST
    };

    virtual void MultTAdd(const SparseVector& mvec, SparseVector& rvec) const {
      Error( "Function MultTAdd not re-implemented in derived class!",
	     __FILE__, __LINE__);
    }

    //! Perform a matrix-vector multiplication rvec -= this*mvec

    //! This method performs a matrix-vector multiplication rvec -= this*mvec.
    //! However, all it does is downcast the BaseVectors to SparseVectors and
    //! delegate the work to the method with the appropriate interface. It
    //! implements the method defined in the BaseMatrix class.
    void MultSub( const BaseVector& mvec, BaseVector& rvec ) const {
      TRY_CAST
        ConstRefCast(mvec,SparseVector,stdmvec);
      RefCast(rvec,SparseVector,stdrvec);
      MultSub(stdmvec,stdrvec);
      CATCH_CAST
    };
        
    //! Perform a matrix-vector multiplication rvec -= this*mvec

    //! This method performs a matrix-vector multiplication rvec -= this*mvec.
    //! It offers an appropriate interface for SparseVectors.
    virtual void MultSub( const SparseVector& mvec, SparseVector& rvec ) const = 0;

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
      (*error) << "Sorry, but SetSparsityPattern is only implemented "
               << "for the SCRS_Matrix class currently!";
      Error( __FILE__, __LINE__ );
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
      (*error) << "Sorry, but TransferPatternToPool is only implemented "
               << "for the SCRS_Matrix class currently!";
      Error( __FILE__, __LINE__ );
      return NO_PATTERN_ID;
    }

    //! Copy the sparsity pattern of this matrix to another 
    
    //! This method copies the sparsity pattern of this matrix and
    //! sets it in another matrix
    virtual void CopySparsityPattern ( StdMatrix &mat ) const {
      Error( "StdMatrix::CopySparsityPattern(): Not implemented here",
             __FILE__, __LINE__ );
    }

    //! Set a specific matrix entry

    /*! this functions sets the Matrix entry in row i and column j
      to the value v. Depending on the entry type of the matrix,
      v might be a Double, Complex or a tiny Matrix of either type
    */
    virtual void SetMatrixEntry( Integer i, Integer j, Double &v, 
                                 bool setCounterPart ) {
      Error( "StdMatrix::SetMatrixEntry(): Not implemented here",
            __FILE__, __LINE__);
    }

    //! Get a specific matrix entry

    /*! this function returns a reference to the Matrix entry in 
      row i and column j. Depending on the entry type of the matrix,
      v might be a Double, Complex or a tiny Matrix of either type
    */
    virtual void GetMatrixEntry( Integer i, Integer j, Double &v ) const{
      Error( "StdMatrix::GetMatrixEntry(): Not implemented here",
            __FILE__, __LINE__);
    }

    //! gets entries of system matrix as a vector

    /*! copies the system matrix from SCRS format into a vector (complex)
      containing all (even all zeros) upper triangle elements
      is called by method GetFullSystemMatrixAsVec in standardsys which 
      provides an interface to CFS++
    */
    virtual void CopySCRSMatrix2Vec(Complex* &vec){;};
 

    //! Add value to a matrix entry

    //! This method adds the value v to the matrix entry at position (i,j).
    //! The method is documented here for the case of Double scalar entries,
    //! but is defined for all entry types.
    virtual void AddToMatrixEntry(Integer i, Integer j, Double &v ) {
      Error( "StdMatrix::AddToMatrixEntry(): Not implemented here",
            __FILE__, __LINE__);
    }

    //!Set the diagonal entry of row i to v
    virtual void SetDiagEntry(Integer i, Double &v) {
      Error( "StdMatrix::SetDiagEntry(): Not implemented here",
             __FILE__, __LINE__ );
    }

    //!Set the diagonal entry of row i to v     
    virtual void GetDiagEntry(Integer i, Double &v) const {
      Error("StdMatrix::GetDiagEntry(): Not implemented here",
            __FILE__, __LINE__ );
    }

    //! slow direct matrix access
#define DECL_STDMATRIX_FCN(TYPE) \
virtual void SetMatrixEntry(Integer i, Integer j, TYPE &v,\
bool setCounterPart){\
Error("StdMatrix::SetMatrixEntry(): Not implemented here",__FILE__,__LINE__);}\
\
virtual void GetMatrixEntry(Integer i, Integer j, TYPE &v) const {\
Error("StdMatrix::GetMatrixEntry(): Not implemented here",__FILE__,__LINE__);}\
\
virtual void AddToMatrixEntry(Integer i, Integer j, TYPE &v ) { \
Error("StdMatrix::AddToMatrixEntry(): Not implemented here", __FILE__, \
__LINE__); }\
\
virtual void SetDiagEntry(Integer i, TYPE &v){ \
Error("StdMatrix::SetDiagEntry(): Not implemented here",__FILE__,__LINE__);}\
\
virtual void GetDiagEntry(Integer i, TYPE &v) const { \
Error("StdMatrix::GetDiagEntry(): Not implemented here",__FILE__,__LINE__);}

DECL_STDMATRIX_FCN(Complex);


     /** for debugging purpose it is cool to have a name to identify the matrix */
     std::string name;
     
     /** Do some debug output */
     std::string ToString();
  protected:

    //! Number of matrix rows
    Integer nrows_;

    //! Number of matrix columns
    Integer ncols_;

   private:
     /** A small harwell boint export class, templated to handel double and complex
      * easier */
    template<typename T>
    class HarwellBoeing
    {
        public:
        HarwellBoeing(const StdMatrix* matrix);
        
        void Export(const std::string& file, const BaseVector& rhs);
        
        private:
        void Fill(const T* values, const Integer* ints, UInt length, std::vector<std::string>& out);
        
        bool is_complex_;
        const StdMatrix* matrix_;
    };   

  };


} // namespace



#endif // OLAS_BASEMATRIX_HH
