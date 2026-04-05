#ifndef OLAS_SBM_VECTOR_HH
#define OLAS_SBM_VECTOR_HH

#include <iostream>

#include "BaseVector.hh"
#include "SingleVector.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {


  //! This is the vector class associated to the SBM_Matrix class
  class SBM_Vector : public BaseVector {

  public:

    //! Default Constructor

    //! The default constructor generates an empty vector, so size is set to
    //! 0 and subVec_ to NULL. In addition the entryType of this vector is set,
    //! so any subsequent call to Resize will create sub-vectors of the 
    //! given type. 
    explicit SBM_Vector( BaseMatrix::EntryType entryType = 
                         BaseMatrix::NOENTRYTYPE );
    
    //! Constructor for empty vectors of specified length

    //! The constructor will generate a SBM_Vector of length \a size and
    //! dynamically allocate memory to store the sub-vectors.
    //! If the \a entryType is not given explicitly, the entries will be
    //! NULL-pointers. Otherwise 0-size sub-vectors of the given type 
    //! will be created. 
    explicit SBM_Vector( UInt size, BaseMatrix::EntryType entryType = 
                         BaseMatrix::NOENTRYTYPE );
    
    //! Create a weak copy of the SBM-vector
    
    //! This method creates a weak copy of size #numRows of the original 
    //! vector. This weak copy will not get ownership of the sub-vectors,
    //! but just get hold of the pointers. Thus, the weak copy is not
    //! allowed to set new sub-vectors.
    //! \param origVector vector to be copied
    //! \param numRows number of rows to be copied
    SBM_Vector( const SBM_Vector& origVector, UInt numRows );
    
    //! Destructor

    //! The destructor is deep, i.e. it will free all dynamically
    //! allocated memory.
    virtual ~SBM_Vector();

    //! Method to re-set the entry type of a SBM_Vector

    //! This method checks, if the SBM_Vector was created with the first constructor,
    //! where only the EntryType is set AND if no further actions were performed.
    //! If these requirements are met, we change the entrytype to the given one.
    void ResetEntryType(BaseMatrix::EntryType entryType);

    //! Return the Entry type of the vector

    //! The method returns the entry type of the vector on the scalar level.
    //! This is encoded as a value of the enumeration data type
    //! MatrixEntryType. If no sub-vector has yet been set, NOENTRYTYPE is
    //! returned.
    BaseMatrix::EntryType GetEntryType() const {
      return myEntryType_;
    }

    //! Set number of vector entries without creation of sub-vectors

    //! This sets the length of the vector and allocates memory accordingly.
    //! In the case the size was not zero previously the internal array subVec_
    //! will be re-sized. The new vector will have NULL pointers as entries.
    void SetSize( UInt size );
    
    //! Obtain size of vector, i.e. number of sub-vector entries
    UInt GetSize() const {
      return size_;
    }
    
    //! Set ownership status for sub-vectors
    
    //! If set to true, the SBM vector will not delete the sub-vectors upon
    //! destruction. If set to yes, the sub-vectors will be deleted,
    //! independent of the original setting.
    void SetOwnership( bool ownsVectors ) {
      ownSubVectors_ = ownsVectors;
    }

    //! Resize the vector to new size with optional creation of sub-vectors
    
    //! This sets the length of the vector to the specified size.
    //! This method tries to keep as much data as possible, i.e. if the new
    //! size is smaller than the current one, smaller number of entries is
    //! maintained. If the new size is larger than the old one, all old
    //! entries are copied into the new Vector. Newly created entries
    //! sill be initialized with 0-size vectors of the given #entryType_,
    //! i.e. after the call of this methods, no new NULL-pointers are 
    //! created.
    virtual void Resize( UInt newSize );

    //! Insert a sub-vector

    //! A call to this method will set the i-th entry of the SBM_Vector to
    //! be the given SingleVector subvec.
    void SetSubVector( SingleVector *subvec, UInt i );
   
    //! Overload assignment operator
    virtual BaseVector& operator= ( const BaseVector &bVec ) {
      const SBM_Vector& sVec = dynamic_cast<const SBM_Vector&>(bVec);
      *this = sVec;
      return *this;
    }

    //! Overload assignment operator
    virtual SBM_Vector& operator= ( const SBM_Vector &bVec );
    

    //! Add to sub-vector

    //! A call to this method will add the given standard vector to the i-th
    //! entry of the SBM_Vector.
    void AddToSubVector( const SingleVector &vec, const UInt i );

    //! Initialise entries to zero

    //! The method sets all entries of the vector to zero. Here entry does not
    //! denote the sub-vectors but the individual entries of the SBM_Vector
    //! if not seen in a block fashion.
    void Init() {
      for ( UInt i = 0; i < size_; i++ ) {
        if ( subVec_[i] != NULL ) {
          subVec_[i]->Init();
        }
      }
    }

    //! Add the given vector to the this vector object

    //! The method will try to add the given vector to this vector object.
    //! This, however, can only work, if both vectors are of the same type
    //! and so are all their respective sub-vectors.
    void Add( const BaseVector &vec );


    /** Return vector as separated string
     * @see Vector::ToString() */
    std::string ToString(ToStringFormat format = TS_PYTHON, const std::string& sep="", int digits=-1) const;

    //! Export vector to file

    //! This method can be used to export the non-zero sub-vectors to ascii
    //! files. The format is extremely simple. If the sub-vector is of
    //! dimension \f$n\times 1\f$, then the output file will contain
    //! \f$n+1\f$ rows.
    //! The first row contains the dimension \f$n\f$, while the remaining
    //! rows contain the vector's entries, so row (k+1) contains entry
    //! \f$v_k\f$.
    //! \note The filenames are determined as fname_<em>i</em>.vec where
    //!       <em>i</em> is replaced by the sub-vector's index.
    void Export(const std::string& fname, BaseMatrix::OutputFormat format ) const;

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x.
    void Axpy( const Double alpha, const BaseVector &y );
    void Axpy( const Complex alpha, const BaseVector &y );
    //@}

    //@{
    //! Compute inner product

    //! The method computes the value of the inner product between this vector
    //! and the input vector vec. The value is returned in retval.
    void Inner( const BaseVector& vec, Double &retval ) const;
    void Inner( const BaseVector& vec, Complex &retval ) const;
    //@}

    //@{
    //! Add a scaled version of a vector to this vector object

    //! The method takes this vector object \f$x\f$ and replaces it with
    //! the result of \f$x + \alpha v\f$.
    void Add( Double alpha, const BaseVector& v );
    void Add( Complex alpha, const BaseVector& v );
    //@}

    //@{
    //! Replace this vector object by the sum of two scaled vectors

    //! This method replaces this vector object by the sum
    //! \f$\alpha x +\beta y\f$.
    void Add( Double alpha, const BaseVector& y,
	      Double beta, const BaseVector& z );
    void Add( Complex alpha, const BaseVector& y,
	      Complex beta, const BaseVector& z );
    //! This should add only something for the N-th entry
    void Add( Complex alpha, const BaseVector& y,
          Complex beta, const BaseVector& z , UInt h);
    //@}

    //! Compute the Euclidean norm of this vector object
    Double NormL2() const;

    //! Compute the Euclidean norm with respect to the other object
    Double NormL2(const SBM_Vector& other) const;

    //! Divide each vector entry by specified real-valued scalar
    void ScalarDiv( const Double factor );

    //! Multiply each vector entry by specified real-valued scalar
    void ScalarMult( const Double factor );

    //! Divide each vector entry by specified complex-valued scalar
    void ScalarDiv( const Complex factor );

    //! Multiply each vector entry by specified complex-valued scalar
    void ScalarMult( const Complex factor );


    // =======================================================================
    // ACCESS TO SUB-VECTORS
    // =======================================================================

    //! \name Methods to access sub-vectors

    //@{
    //! Overload () operator to retrieve sub-vectors

    //! The () operator is overloaded to allow to retrieve individual
    //! sub-vectors. 
    SingleVector& operator()( UInt i ) const {
      if ( subVec_[i] == NULL ) {
        EXCEPTION( "Cannot return reference to non-existing sub-vector "
                 << "with index " << i);
      }
      return *(subVec_[i]);
    }

    //! Obtain pointer to a sub-vector

    //! The method returns a pointer to a sub-vector. The access follows the
    //! one-based indexing convention of OLAS.
    SingleVector* GetPointer( UInt i ) {
      return subVec_[i];
    }

    //! Obtain pointer to a sub-vector

    //! The method returns a pointer to a sub-vector. The access follows the
    //! one-based indexing convention of OLAS.
    const SingleVector* GetPointer( UInt i ) const {
      return subVec_[i];
    }

    //@}

  private:

    //! Array containing the sub-vectors

    //! Array containing pointers to the sub-vectors fitting to the
    //! sub-matrices in the associated SBM_Matrix
    StdVector<SingleVector *> subVec_;

    //! Number of sub-vectors
    UInt size_;
    
    //! Flag indicating if object is is responsible for deletion of subvectors
    bool ownSubVectors_;

    //! Attribute storing the entry type of the vector

    //! Currently the %SBM_Vector class only allows for storing sub-matrices
    //! that share a common entry type. This attribute stores the entry type
    //! of the first sub-vector that is set and uses this as the over-all
    //! entry type of the %SBM_Vector object.
    BaseMatrix::EntryType myEntryType_;

  };

}

#endif
