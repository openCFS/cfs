// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef SingleVector_HH
#define SingleVector_HH

#include "General/defs.hh"
#include "MatVec/basevector.hh"

namespace CoupledField {

  //! Base class for all vectors with one continuous block of data
  class SingleVector : public BaseVector {

  public:

    //! Default constructor
    SingleVector();

    //! Default destructor
    virtual ~SingleVector();

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
    virtual void Export( const char *fname ) const {
      EXCEPTION( "Class not implemented by derived class" );
    }

    //! Add vec to this vector object

    //! This method adds vec to this vector object. However, it does not
    //! perform this task itself, but only performs a downcast and calls
    //! the method with an SingleVector interface.
    void Add(const BaseVector& vec){
      TRY_CAST
      CONSTREFCAST(vec,SingleVector,stdvec);
      Add(stdvec);
      CATCH_CAST
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
      TRY_CAST
      CONSTREFCAST(y,SingleVector,std_y);
      Axpy( alpha, std_y );
      CATCH_CAST
    }
    virtual void Axpy( const Complex alpha, const BaseVector &y ) {
      TRY_CAST
      CONSTREFCAST(y,SingleVector,std_y);
      Axpy( alpha, std_y );
      CATCH_CAST
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
      EXCEPTION( " Method NormL2 not implemented " );
      return 0;
    }

    //! Overload assignment operator

    //! Following our standard policy this implementation of the overloaded
    //! assignment operator does simply down-cast the base vector and call
    //! the implementation for standard vectors.
    virtual BaseVector &operator= ( const BaseVector &bvec ) {
      TRY_CAST {
        CONSTREFCAST( bvec, SingleVector, svec );
        *this = svec;
      } CATCH_CAST;
      return *this;
    }
    
    //! Overload assignment operator

    //! This is the implementation of the overloaded assignment operator with
    //! an interface for standard vectors.
    virtual SingleVector &operator= ( const SingleVector &stdvec ) {
      EXCEPTION( "SingleVector::operator= not over-written by derived class!" );
      return *this;
    }

#define DECL_SingleVector_FCN(TYPE)\
	void Add(TYPE a, const BaseVector& vec){\
	TRY_CAST\
	CONSTREFCAST(vec,SingleVector,stdvec);\
	Add(a,stdvec);\
	CATCH_CAST\
	}\
\
        virtual void Add(TYPE a,const SingleVector &vec)\
        {EXCEPTION("SingleVector::Add(): Not implemented here");};\
\
	void Add(TYPE a, const BaseVector& vec1,\
		TYPE b,const BaseVector& vec2){\
	TRY_CAST\
	CONSTREFCAST(vec1,SingleVector,stdvec1);\
	CONSTREFCAST(vec2,SingleVector,stdvec2);\
	Add(a,stdvec1,b,stdvec2);\
	CATCH_CAST\
	}\
\
	virtual void Add(TYPE a, const SingleVector& vec,\
                         TYPE b, const SingleVector& vec2)\
	 {EXCEPTION("SingleVector::Add(): Not implemented here");};\
\
	void Inner(const BaseVector& vec,TYPE& s) const {\
	TRY_CAST\
	CONSTREFCAST(vec,SingleVector,stdvec);\
	Inner(stdvec,s);\
	CATCH_CAST\
	}\
\
	virtual void Inner(const SingleVector& vec,TYPE& s) const\
	 {EXCEPTION("SingleVector::Inner(): Not implemented here");}

    DECL_SingleVector_FCN(Double);
    DECL_SingleVector_FCN(Complex);

  protected:

    //! Length of the vector object
    UInt size_;

  };

} // namespace


#endif // OLAS_BASEVECTOR_HH
