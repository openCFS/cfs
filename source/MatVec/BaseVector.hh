// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEVECTOR_HH
#define OLAS_BASEVECTOR_HH

#include <iostream>
#include "MatVec/BaseMatrix.hh"
#include "General/EnvironmentTypes.hh"
#include "General/Exception.hh"

//#include "MatVec/typedefs.hh"


namespace CoupledField {

  //! Generic algebraic vector

  //! This is a pure umbrella class for passing references/pointers around.
  class BaseVector {

  public:

    //! Default Constructor
    BaseVector() {
    }

    //! Destructor
    virtual ~BaseVector() {
    }

    //! Get number of vector entries
    virtual UInt GetSize() const = 0;

    //! Resize the vector to new size 
    virtual void Resize( UInt newSize ) = 0;

    //! Initialize entries with 0
    virtual void Init() = 0;

    /** @see Vector::ToString() */
    virtual std::string ToString(ToStringFormat format = TS_PYTHON, const std::string& sep="", int digits=-1) const = 0;

    //! Return the Entry type of the vector

    //! The method returns the entry type of the vector on the scalar level.
    //! This is encoded as a value of the enumeration data type
    //! MatrixEntryType.
    virtual BaseMatrix::EntryType GetEntryType() const = 0;


    /** is the entry type a constant one? */
    bool IsComplex() const { return GetEntryType() == BaseMatrix::COMPLEX || GetEntryType() == BaseMatrix::F77COMPLEX16; }

    //! Export vector to file

    //! This method can be used to export the vector to an ascii file. The
    //! format is extremely simple. If the vector is of dimension
    //! \f$n\times 1\f$, then the output file will contain \f$n+1\f$ rows.
    //! The first row contains the dimension \f$n\f$, while the remaining
    //! rows contain the vector's entries, so row (k+1) contains entry
    //! \f$a_k\f$.
    virtual void Export(const std::string& fname,
                         BaseMatrix::OutputFormat format = 
                           BaseMatrix::MATRIX_MARKET ) const = 0;

    /** Import ascii file to vector
     * This method can be used to import the vector from an ascii file.
     * Only MatrixMarket formated ascii files can be imported.
     * Information about the MatrixMarket format can be found here:
     * https://math.nist.gov/MatrixMarket/formats.html
     */
    virtual void Import(const std::string& fname, bool checkSize)
    {
      EXCEPTION( "Class not implemented by derived class" );
    }

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
    virtual Double NormL2() const = 0;

    //! Overload assignment operator
    virtual BaseVector &operator= ( const BaseVector &bvec ) {
      EXCEPTION( "BaseVector::operator= not over-written by derived class!");
      return *this;
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
    virtual void SetEntry( UInt i, const Double &val ) {
      EXCEPTION( "BaseVector::SetEntry: Not over-written by derived class");
    }

    //! Add val to the value of a vector entry

    //! This method adds val to the entry of the vector at position i.
    //! The method is documented here for the case of Double
    //! scalar entries, but is defined for all entry types.
    virtual void AddToEntry( UInt i, const Double &val ) {
      EXCEPTION( "BaseVector::AddToEntry: Not over-written by derived class");
    }

    //! Query the value of a vector entry

    //! This method queries the entry of the vector at position i. The method
    //! is documented here for the case of Double scalar entries, but is
    //! defined for all entry types.
    virtual void GetEntry( UInt i, Double &val ) const {
      EXCEPTION( "BaseVector::GetEntry: Not over-written by derived class" );
    }

    /** shortcuts */
    double GetDoubleEntry(unsigned int i) const {
      double v;
      GetEntry(i, v);
      return v;
    }

    Complex GetComplexEntry(unsigned int i) const {
      Complex v;
      GetEntry(i, v);
      return v;
    }


    // Get/Set/AddToEntry versions for non-Double entry 
#define DECL_ENTRY_FCN(TYPE) \
virtual void SetEntry( UInt i, const TYPE &val ){                \
EXCEPTION( "BaseVector::SetEntry: Not over-written by derived class");  } \
virtual void AddToEntry( UInt i, const TYPE &val ){              \
EXCEPTION( "BaseVector::AddToEntry: Not over-written by derived class");} \
virtual void GetEntry( UInt i, TYPE &val ) const {               \
EXCEPTION( "BaseVector::GetEntry: Not over-written by derived class" );}

DECL_ENTRY_FCN(Complex);
#undef DECL_ENTRY_FCN



  private:    
  };

} // namespace

#endif // OLAS_BASEVECTOR_HH
