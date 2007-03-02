#ifndef OLAS_SBM_VECTOR_HH
#define OLAS_SBM_VECTOR_HH

#include <iostream>
#include "utils/utils.hh"
#include "matvec/basevector.hh"
#include "matvec/stdvector.hh"

namespace OLAS {

  //! This is the vector class associated to the SBM_Matrix class
  class SBM_Vector : public BaseVector {

  public:

    //! Default Constructor

    //! The default constructor generates an empty vector, so size is set to
    //! 0 and subVec_ to NULL
    SBM_Vector() : subVec_(NULL), size_(0) {};

    //! Constructor for empty vectors of specified length

    //! The constructor will generate a SBM_Vector of length nlength and
    //! dynamically allocate memory to store the sub-vectors. In this initial
    //! state the entries will all be NULL pointers.
    SBM_Vector( UInt size );

    //! Destructor

    //! The destructor is deep, i.e. it will free all dynamically
    //! allocated memory.
    virtual ~SBM_Vector();

    //! Set number of vector entries

    //! This sets the length of the vector and allocates memory accordingly.
    //! In the case the size was not zero previously the internal array subVec_
    //! will be re-sized. The new vector will have NULL pointers as entries.
    void SetSize( Integer size );

    //! Obtain size of vector, i.e. number of sub-vector entries
    UInt GetSize() const {
      return size_;
    };

    //! Insert a sub-vector

    //! A call to this method will set the i-th entry of the SBM_Vector to
    //! be the given StdVector subvec.
    void SetSubVector( StdVector *subvec, Integer i ) {
      ENTER_FCN( "SBM_Vector::SetSubVector" );
      delete subVec_[i];
      subVec_[i] = subvec;
    }

    //! Add to sub-vector

    //! A call to this method will add the given standard vector to the i-th
    //! entry of the SBM_Vector.
    void AddToSubVector( const StdVector &vec, const Integer i );

    //! Initialise entries to zero

    //! The method sets all entries of the vector to zero. Here entry does not
    //! denote the sub-vectors but the individual entries of the SBM_Vector
    //! if not seen in a block fashion.
    void Init() {
      ENTER_FCN( "SBM_Vector::Init" );
      for ( UInt i = 1; i <= size_; i++ ) {
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

    //! Print vector to output stream os
    void Print ( std::ostream &os ) const {
      Warning( "Print not implemented for SBM_Vector", __FILE__, __LINE__ );
    }

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
    void Export( const Char *fname ) const;

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
    //@}

    //! Compute the Euclidean norm of this vector object
    Double NormEuclid() const;

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
    //! sub-vectors. The access follows the one-based indexing convention
    //! of OLAS.
    StdVector& operator()( UInt i ) const {
      if ( subVec_[i] == NULL ) {
        (*error) << "Cannot return reference to non-existing sub-vector "
                 << "with index " << i;
        Error( __FILE__, __LINE__ );
      }
      return *(subVec_[i]);
    }

    //! Obtain pointer to a sub-vector

    //! The method returns a pointer to a sub-vector. The access follows the
    //! one-based indexing convention of OLAS.
    StdVector* GetPointer( UInt i ) {
      return subVec_[i];
    };

    //! Obtain pointer to a sub-vector

    //! The method returns a pointer to a sub-vector. The access follows the
    //! one-based indexing convention of OLAS.
    const StdVector* GetPointer( UInt i ) const {
      return subVec_[i];
    };

    //@}

  private:

    //! Array containing the sub-vectors

    //! Array containing pointers to the sub-vectors fitting to the
    //! sub-matrices in the associated SBM_Matrix
    StdVector **subVec_;

    //! Number of sub-vectors
    UInt size_;

  };

}

#endif
