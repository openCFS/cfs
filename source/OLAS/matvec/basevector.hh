#ifndef OLAS_BASEVECTOR_HH
#define OLAS_BASEVECTOR_HH

#include "matvec/typedefs.hh"
#include "utils/utils.hh"

namespace OLAS {

  //! Generic algebraic vector

  //! This is a pure umbrella class for passing references/pointers around.
  class BaseVector {
  
  public:	

    //! Defaul Constructor
    BaseVector() {
      ENTER_FCN( "BaseVector::BaseVector" );
    }

    //! Destructor
    virtual ~BaseVector() {
      ENTER_FCN( "BaseVector::~BaseVector" );
    }

    //! Get number of vector entries
    virtual UInt GetSize() const = 0;

    //! Initialize entries with 0
    virtual void Init() = 0;
	
    //! Print vector to output stream os
    virtual void Print(std::ostream &os) const = 0;
	
    //! Print vector to output stream *cla
    void Print(void) const {
      Print(*cla);
    }

    //! Export vector to file

    //! This method can be used to export the vector to an ascii file. The
    //! format is extremely simple. If the vector is of dimension
    //! \f$n\times 1\f$, then the output file will contain \f$n+1\f$ rows.
    //! The first row contains the dimension \f$n\f$, while the remaining
    //! rows contain the vector's entries, so row (k+1) contains entry
    //! \f$a_k\f$.
    virtual void Export( const Char *fname ) const = 0;

    //! Add vec to this vector
    virtual void Add(const BaseVector& vec) = 0;

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x.
    virtual void Axpy( const Double alpha, const BaseVector &y ) = 0;
    virtual void Axpy( const Complex alpha, const BaseVector &y ) = 0;
    //@}

    //@{
    //! Compute inner product

    //! The method computes the value of the inner product between this vector
    //! and the input vector vec. The value is returned in retval.
    virtual void Inner( const BaseVector& vec, Double &retval ) const = 0;
    virtual void Inner( const BaseVector& vec, Complex &retval ) const = 0;
    //@}

    //@{
    //! Add a scaled version of a vector to this vector object

    //! The method takes this vector object \f$x\f$ and replaces it with
    //! the result of \f$x + \alpha v\f$.
    virtual void Add( Double alpha, const BaseVector& v ) = 0;
    virtual void Add( Complex alpha, const BaseVector& v ) = 0;
    //@}

    //@{
    //! Replace this vector object by the sum of two scaled vectors

    //! This method replaces this vector object by the sum
    //! \f$\alpha x +\beta y\f$.
    virtual void Add( Double alpha, const BaseVector& y,
		      Double beta, const BaseVector& z ) = 0;
    virtual void Add( Complex alpha, const BaseVector& y,
		      Complex beta, const BaseVector& z ) = 0;
    //@}

    //! Compute the Euclidean norm of this vector object
    virtual Double NormEuclid() const = 0;
    
    //! Overload assignment operator
    virtual BaseVector &operator= ( const BaseVector &bvec ) {
      Error( "BaseVector::operator= not over-written by derived class!",
	     __FILE__, __LINE__ );
      return *this;
    }

    //! Query wether the vector is parallel or sequential
    virtual bool IsParallel() const {
      return false;
    }

    //! Divide each vector entry by specified real-valued scalar
    virtual void ScalarDiv( const Double factor ) = 0;

    //! Multiply each vector entry by specified real-valued scalar
    virtual void ScalarMult( const Double factor ) = 0;

    //! Divide each vector entry by specified complex-valued scalar
    virtual void ScalarDiv( const Complex factor ) = 0;

    //! Multiply each vector entry by specified complex-valued scalar
    virtual void ScalarMult( const Complex factor ) = 0;

    //! Set the value of a vector entry

    //! This method sets the entry of the vector at position i to the
    //! specified value. The method is documented here for the case of Double
    //! scalar entries, but is defined for all entry types.
    virtual void SetVectorEntry( Integer i, const Double &val ) {
      Error( "BaseVector::SetVectorEntry: Not over-written by derived class",
	     __FILE__, __LINE__);
    }

    //! Add val to the value of a vector entry

    //! This method adds val to the entry of the vector at position i.
    //! The method is documented here for the case of Double
    //! scalar entries, but is defined for all entry types.
    virtual void AddToVectorEntry( Integer i, const Double &val ) {
      Error( "BaseVector::AddToVectorEntry: Not over-written by derived class",
	     __FILE__, __LINE__);
    }


    // SetEntry versions for non-Double entry types

#define DECL_SETENTRY_FCN(TYPE) \
virtual void SetVectorEntry( Integer i, const TYPE &val ){\
Error( "BaseVector::SetVectorEntry: Not over-written by derived class",\
__FILE__, __LINE__);}\
virtual void AddToVectorEntry( Integer i, const TYPE &val ){\
Error( "BaseVector::AddToVectorEntry: Not over-written by derived class",\
__FILE__, __LINE__);}

    DECL_SETENTRY_FCN(Complex);

    //! Query the value of a vector entry

    //! This method queries the entry of the vector at position i. The method
    //! is documented here for the case of Double scalar entries, but is
    //! defined for all entry types.
    virtual void GetVectorEntry( Integer i, Double &val ) const {
      Error( "BaseVector::GetVectorEntry: Not over-written by derived class",
	     __FILE__, __LINE__);
    }

    // GetEntry versions for non-Double entry types

#define DECL_GETENTRY_FCN(TYPE) \
virtual void GetVectorEntry( Integer i, TYPE &val ) const {\
Error( "BaseVector::GetVectorEntry: Not over-written by derived class",\
__FILE__,__LINE__);}

    DECL_GETENTRY_FCN(Complex);
  };



} // namespace

#endif // OLAS_BASEVECTOR_HH
