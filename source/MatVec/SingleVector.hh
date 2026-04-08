#ifndef SingleVector_HH
#define SingleVector_HH

#include "General/defs.hh"
#include "BaseVector.hh"

namespace CoupledField {

  //! Base class for all vectors with one continuous block of data
  class SingleVector : public BaseVector {

  public:

    //! Default constructor
    SingleVector() {};

    //! Default destructor
    virtual ~SingleVector() {};

    //! Return the length of the vector
    UInt GetSize() const {
      return size_;
    }
    
    //! Resize the vector to new size 
    virtual void Resize( UInt newSize ) = 0;
    
    //! Export vector to file

    //! This method can be used to export the vector to an ascii file. The
    //! format is extremely simple. If the vector is of dimension
    //! \f$n\times 1\f$, then the output file will contain \f$n+1\f$ rows.
    //! The first row contains the dimension \f$n\f$, while the remaining
    //! rows contain the vector's entries, so row (k+1) contains entry
    //! \f$a_k\f$.
    virtual void Export( const std::string& fname,
                         BaseMatrix::OutputFormat format ) const {
      EXCEPTION( "Class not implemented by derived class" );
    }

    //! Add vec to this vector object

    //! This method adds vec to this vector object. However, it does not
    //! perform this task itself, but only performs a downcast and calls
    //! the method with an SingleVector interface.
    void Add(const BaseVector& vec){
      Add(dynamic_cast<const SingleVector&>(vec));
    }

    //! Add vec to this vector object

    //! This method adds vec to this vector object. It offers an interface
    //! for SingleVectors.
    virtual void Add(const SingleVector &vec) = 0;

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x. \n
    //! This implementation only performs a downcast and calls the method
    //! with an interface for SingleVectors.
    virtual void Axpy( const Double alpha, const BaseVector &y ) {
      Axpy( alpha, dynamic_cast<const SingleVector&>(y) );
    }
    
    virtual void Axpy( const Complex alpha, const BaseVector &y ) {
      Axpy( alpha, dynamic_cast<const SingleVector&>(y) );
    }
    //@}

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x. This method offers an interface for SingleVectors.
    virtual void Axpy( const Double alpha, const SingleVector &y ) {
      EXCEPTION( "Axpy not over-written by derived class" );
    }
    virtual void Axpy( const Complex alpha, const SingleVector &y ) {
      EXCEPTION( "Axpy not over-written by derived class" );
    }
    //@}

    //! Compute Euclidean norm of this vector object
    Double NormL2() const {
      EXCEPTION("Method NormL2 not implemented");
    }

    virtual Double SignedMax() const {
      EXCEPTION("Method SignedMax() not implemented");
    }

    //! Overload assignment operator

    //! Following our standard policy this implementation of the overloaded
    //! assignment operator does simply down-cast the base vector and call
    //! the implementation for standard vectors.
    virtual BaseVector &operator= ( const BaseVector &bvec ) {
      const SingleVector& svec = dynamic_cast<const SingleVector&>(bvec);
      *this = svec;
      return *this;
    }
    
    //! Overload assignment operator

    //! This is the implementation of the overloaded assignment operator with
    //! an interface for standard vectors.
    virtual SingleVector &operator= ( const SingleVector &stdvec ) {
      EXCEPTION( "SingleVector::operator= not over-written by derived class!" );
      return *this;
    }

    void Add(Double a, const BaseVector& vec)
    {
      Add(a, dynamic_cast<const SingleVector&>(vec));
    }

    virtual void Add(Double a, const SingleVector &vec)
    {
      EXCEPTION("SingleVector::Add(): Not implemented here");
    }

    void Add(Double a, const BaseVector& vec1,
             Double b, const BaseVector& vec2)
    {
      const SingleVector& stdvec1 = dynamic_cast<const SingleVector&>(vec1);
      const SingleVector& stdvec2 = dynamic_cast<const SingleVector&>(vec2);
      Add(a,stdvec1,b,stdvec2);
    }

    virtual void Add(Double a, const SingleVector& vec,
                     Double b, const SingleVector& vec2)
    {
      EXCEPTION("SingleVector::Add(): Not implemented here");
    }

    void Inner(const BaseVector& vec, Double& s) const
    {
      Inner(dynamic_cast<const SingleVector&>(vec), s);
    }

    virtual void Inner(const SingleVector& vec, Double& s) const
    {
      EXCEPTION("SingleVector::Inner(): Not implemented here");
    }
    
    void Add(Complex a, const BaseVector& vec)
    {
      Add(a, dynamic_cast<const SingleVector&>(vec));
    }

    virtual void Add(Complex a,const SingleVector &vec)
    {
      EXCEPTION("SingleVector::Add(): Not implemented here");
    }

    void Add(Complex a, const BaseVector& vec1,
             Complex b, const BaseVector& vec2)
    {
      const SingleVector& stdvec1 = dynamic_cast<const SingleVector&>(vec1);
      const SingleVector& stdvec2 = dynamic_cast<const SingleVector&>(vec2);
      Add(a,stdvec1,b,stdvec2);
    }

    virtual void Add(Complex a, const SingleVector& vec,
        Complex b, const SingleVector& vec2)
    {
      EXCEPTION("SingleVector::Add(): Not implemented here");
    }

    void Inner(const BaseVector& vec, Complex& s) const
    {
      Inner(dynamic_cast<const SingleVector&>(vec), s);
    }

    virtual void Inner(const SingleVector& vec, Complex& s) const
    {
      EXCEPTION("SingleVector::Inner(): Not implemented here");
    }   
    

  protected:

    //! Length of the vector object
    UInt size_ = 0;
  };

} // namespace


#endif // OLAS_BASEVECTOR_HH
