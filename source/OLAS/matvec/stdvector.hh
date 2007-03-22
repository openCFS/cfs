// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_STDVECTOR_HH
#define OLAS_STDVECTOR_HH

#include "utils/utils.hh"
#include "matvec/basevector.hh"

namespace OLAS {

  //! Base class for algebraic vector
  class StdVector : public BaseVector {

  public:

    //! Default constructor
    StdVector();

    //! Default destructor
    virtual ~StdVector();

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
    //! the method with an StdVector interface.
    void Add(const BaseVector& vec){
      ENTER_FCN("StdVector::Add");
      TRY_CAST
      ConstRefCast(vec,StdVector,stdvec);
      Add(stdvec);
      CATCH_CAST
    }

    //! Add vec to this vector object

    //! This method adds vec to this vector object. It offers an interface
    //! for StdVectors.
    virtual void Add(const StdVector &vec) = 0;

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x. \n
    //! This implementation only performs a downcast and calls the method
    //! with an interface for StdVectors.
    virtual void Axpy( const Double alpha, const BaseVector &y ) {
      ENTER_FCN("StdVector::Axpy");
      TRY_CAST
      ConstRefCast(y,StdVector,std_y);
      Axpy( alpha, std_y );
      CATCH_CAST
    }
    virtual void Axpy( const Complex alpha, const BaseVector &y ) {
      ENTER_FCN("StdVector::Axpy");
      TRY_CAST
      ConstRefCast(y,StdVector,std_y);
      Axpy( alpha, std_y );
      CATCH_CAST
    }
    //@}

    //@{
    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the factor
    //! alpha and adds the vector y to it. The result will over-write the
    //! vector x. This method offers an interface for StdVectors.
    virtual void Axpy( const Double alpha, const StdVector &y ) {
      Error( "Axpy not over-written by derived class", __FILE__, __LINE__ );
    }
    virtual void Axpy( const Complex alpha, const StdVector &y ) {
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
      ENTER_IFCN( "StdVector::operator=" );
      TRY_CAST {
	ConstRefCast( bvec, StdVector, svec );
	*this = svec;
      } CATCH_CAST;
      return *this;
    }
    
    //! Overload assignment operator

    //! This is the implementation of the overloaded assignment operator with
    //! an interface for standard vectors.
    virtual StdVector &operator= ( const StdVector &stdvec ) {
      Error( "StdVector::operator= not over-written by derived class!",
	     __FILE__, __LINE__ );
      return *this;
    }

#define DECL_STDVECTOR_FCN(TYPE)\
	void Add(TYPE a, const BaseVector& vec){\
	ENTER_FCN("StdVector::Add");\
	TRY_CAST\
	ConstRefCast(vec,StdVector,stdvec);\
	Add(a,stdvec);\
	CATCH_CAST\
	}\
\
        virtual void Add(TYPE a,const StdVector &vec)\
        {Error("StdVector::Add(): Not implemented here",__FILE__,__LINE__);};\
\
	void Add(TYPE a, const BaseVector& vec1,\
		TYPE b,const BaseVector& vec2){\
	ENTER_FCN("StdVector::Add");\
	TRY_CAST\
	ConstRefCast(vec1,StdVector,stdvec1);\
	ConstRefCast(vec2,StdVector,stdvec2);\
	Add(a,stdvec1,b,stdvec2);\
	CATCH_CAST\
	}\
\
	virtual void Add(TYPE a, const StdVector& vec,\
                         TYPE b, const StdVector& vec2)\
	 {Error("StdVector::Add(): Not implemented here",__FILE__,__LINE__);};\
\
	void Inner(const BaseVector& vec,TYPE& s) const {\
	ENTER_FCN("StdVector::Inner");\
	TRY_CAST\
	ConstRefCast(vec,StdVector,stdvec);\
	Inner(stdvec,s);\
	CATCH_CAST\
	}\
\
	virtual void Inner(const StdVector& vec,TYPE& s) const\
	 {Error("StdVector::Inner(): Not implemented here",__FILE__,__LINE__);}

    DECL_STDVECTOR_FCN(Double);
    DECL_STDVECTOR_FCN(Complex);

  protected:

    //! Length of the vector object
    UInt size_;

    //! Number of degrees of freedom, i.e. block size of vector entries.
    UInt dof_;

  };

} // namespace


#endif // OLAS_BASEVECTOR_HH
