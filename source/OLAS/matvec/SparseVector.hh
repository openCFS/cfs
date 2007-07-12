// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_SparseVector_HH
#define OLAS_SparseVector_HH

#include "utils/utils.hh"
#include "matvec/basevector.hh"

namespace OLAS {

  //! Base class for algebraic vector
  class SparseVector : public BaseVector {

  public:

    //! Default constructor
    SparseVector();

    //! Default destructor
    virtual ~SparseVector();

    //! Return the length of the vector
    UInt GetSize() const {
      return size_;
    }

    //! Get number of degrees of freedom

    //! The method returns the number of degrees of freedom. For a vector this
    //! means the size of the actual entries.
    UInt GetNumDoF() const {
      return dof_;
    }

    //! Method to force instantiation of members in templated vectors

    //! This method is used by the GenerateVectorObject function to force the
    //! instantiation of the public member functions of the templated m
    //! classes.
    virtual void InstantiatePublicMethods() {
      Error( "Not implemented here", __FILE__, __LINE__ );
    };

    //! Export vector to file

    //! This method can be used to export the vector to an ascii file. The
    //! format is extremely simple. If the vector is of dimension
    //! \f$n\times 1\f$, then the output file will contain \f$n+1\f$ rows.
    //! The first row contains the dimension \f$n\f$, while the remaining
    //! rows contain the vector's entries, so row (k+1) contains entry
    //! \f$a_k\f$.
    virtual void Export( const Char *fname ) const {
      Error( "Class not implemented by derived class", __FILE__, __LINE__ );
    }

    //! Add vec to this vector object

    //! This method adds vec to this vector object. However, it does not
    //! perform this task itself, but only performs a downcast and calls
    //! the method with an SparseVector interface.
    void Add(const BaseVector& vec){
      ENTER_FCN("SparseVector::Add");
      TRY_CAST
      ConstRefCast(vec,SparseVector,stdvec);
      Add(stdvec);
      CATCH_CAST
    }

    //! Add vec to this vector object

    //! This method adds vec to this vector object. It offers an interface
    //! for SparseVectors.
    virtual void Add(const SparseVector &vec) = 0;

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x. \n
    //! This implementation only performs a downcast and calls the method
    //! with an interface for SparseVectors.
    virtual void Axpy( const Double alpha, const BaseVector &y ) {
      ENTER_FCN("SparseVector::Axpy");
      TRY_CAST
      ConstRefCast(y,SparseVector,std_y);
      Axpy( alpha, std_y );
      CATCH_CAST
    }
    virtual void Axpy( const Complex alpha, const BaseVector &y ) {
      ENTER_FCN("SparseVector::Axpy");
      TRY_CAST
      ConstRefCast(y,SparseVector,std_y);
      Axpy( alpha, std_y );
      CATCH_CAST
    }
    //@}

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x. This method offers an interface for SparseVectors.
    virtual void Axpy( const Double alpha, const SparseVector &y ) {
      Error( "Axpy not over-written by derived class", __FILE__, __LINE__ );
    }
    virtual void Axpy( const Complex alpha, const SparseVector &y ) {
      Error( "Axpy not over-written by derived class", __FILE__, __LINE__ );
    }
    //@}

    //! Compute Euclidean norm of this vector object
    Double NormEuclid() const {
      Error( " Method NormEuclid not implemented ", __FILE__, __LINE__ );
      return 0;
    }

    //! Overload assignment operator

    //! Following our standard policy this implementation of the overloaded
    //! assignment operator does simply down-cast the base vector and call
    //! the implementation for standard vectors.
    virtual BaseVector &operator= ( const BaseVector &bvec ) {
      ENTER_IFCN( "SparseVector::operator=" );
      TRY_CAST {
	ConstRefCast( bvec, SparseVector, svec );
	*this = svec;
      } CATCH_CAST;
      return *this;
    }
    
    //! Overload assignment operator

    //! This is the implementation of the overloaded assignment operator with
    //! an interface for standard vectors.
    virtual SparseVector &operator= ( const SparseVector &stdvec ) {
      Error( "SparseVector::operator= not over-written by derived class!",
	     __FILE__, __LINE__ );
      return *this;
    }

#define DECL_SparseVector_FCN(TYPE)\
	void Add(TYPE a, const BaseVector& vec){\
	ENTER_FCN("SparseVector::Add");\
	TRY_CAST\
	ConstRefCast(vec,SparseVector,stdvec);\
	Add(a,stdvec);\
	CATCH_CAST\
	}\
\
        virtual void Add(TYPE a,const SparseVector &vec)\
        {Error("SparseVector::Add(): Not implemented here",__FILE__,__LINE__);};\
\
	void Add(TYPE a, const BaseVector& vec1,\
		TYPE b,const BaseVector& vec2){\
	ENTER_FCN("SparseVector::Add");\
	TRY_CAST\
	ConstRefCast(vec1,SparseVector,stdvec1);\
	ConstRefCast(vec2,SparseVector,stdvec2);\
	Add(a,stdvec1,b,stdvec2);\
	CATCH_CAST\
	}\
\
	virtual void Add(TYPE a, const SparseVector& vec,\
                         TYPE b, const SparseVector& vec2)\
	 {Error("SparseVector::Add(): Not implemented here",__FILE__,__LINE__);};\
\
	void Inner(const BaseVector& vec,TYPE& s) const {\
	ENTER_FCN("SparseVector::Inner");\
	TRY_CAST\
	ConstRefCast(vec,SparseVector,stdvec);\
	Inner(stdvec,s);\
	CATCH_CAST\
	}\
\
	virtual void Inner(const SparseVector& vec,TYPE& s) const\
	 {Error("SparseVector::Inner(): Not implemented here",__FILE__,__LINE__);}

    DECL_SparseVector_FCN(Double);
    DECL_SparseVector_FCN(Complex);

  protected:

    //! Length of the vector object
    UInt size_;

    //! Number of degrees of freedom, i.e. block size of vector entries.
    UInt dof_;

  };

} // namespace


#endif // OLAS_BASEVECTOR_HH
