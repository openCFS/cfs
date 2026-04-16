#ifndef FILE_CFS_POLYNOMIAL_HH
#define FILE_CFS_POLYNOMIAL_HH

#include <map>
#include <algorithm>
#include <array>

#include "General/Environment.hh"
#include "MatVec/Vector.hh"

namespace CoupledField {

//! Class representing a multidimensional polynomial with constant coefficients

/** This class represents a polynomial with several variables (see template
 parameter DIM) and real or complex valued coefficients (see template 
 parameter T). A polynomial depending on the variables \f$ x,y,z\f$ can be 
 written as
 \f[ P(x,y,z) = \sum_{i=0}^{p_x}  \sum_{j=0}^{p_y} \sum_{k=0}^{p_z}
                C_{ijk} \cdot x^i y^j z^k \f]
 The type of the coefficients \f$ C_{ijk} \f$ can be determined with the
 template parameter T.
 In this class the non-zero coefficients  \f$ C_{ijk} \neq 0 \f$ are stored
 in a stl-map, with the exponents \f$(i,j,k\f$ stored in a std::array 
 used as key to the map.
 The internal coefficients can be queried using the method 
 Polynomial::Coefs. A string representation can be accessed using 
 the method Polynomial::ToString(), where the variable names are taken
 as 'x', 'y' and 'z'. In addition, the polynomial can be evaluated
 using the Polynomial::Eval() method. This implementation is currently 
 very basic and not optimized. 

 \tparam T Type of the coefficients of the polynomial (Double or Complex)
 \tparam DIM Number of variables / dimension of the polynomial.

 \note The implementation of the polynomial is currently limited to 3
       variables due to the explicit generation of variable names. **/
template<typename T, UInt DIM>
class Polynomial {

public:
  
  //! Typedef for specifying polynomial instance
  typedef Polynomial<T, DIM> POL;   
  
  //! Typedef for order type
  
  //! This typedef encapsulates the data type used to represent the 
  //! tuple of variable coefficients, e.g. for the expression
  //! \f$ x^i y^j z^k \f$ the entries \f$ i,j,k \f$ are stored in 
  //! a fixed-size std::array instance.
  typedef typename std::array<UInt, DIM> ORD;
  
  //! Typedef for coefficient map
  typedef typename std::map<ORD, T> COEFMAP;
  
  //! Typedef for coefficient map iterator
  typedef typename COEFMAP::iterator COEF_IT;
  
  //! Typedef for coefficient map iterator (const)
  typedef typename COEFMAP::const_iterator COEF_CONST_IT;

  // ========================================================================
  //  Constructors & Assignment
  // ========================================================================
  //@{ \name Constructors / Assignment
  
  //! Default constructor 
  Polynomial();
  
  //! Constructor for constant value 
  //! \param val Constant value use to initialize this polynomial
  Polynomial(const T& val);
  
  //! Constructor for monomial
  
  //! This constructor generates a polynomial instance of a monomial
  //! \f$ coef \cdot v^{order} \f$, where the variable \f$ v \f$
  //! is determined by its index.
  //! \param coef Constant coefficient
  //! \param index Index to the variable vector
  //! \param order Exponent of the monomial
  explicit Polynomial( const T& coef, UInt index, UInt order = 1 );
  
  //! Copy constructor
  Polynomial( const Polynomial& a );
  
  //! Destructror
  ~Polynomial();

  //! Assignment from constant value
  Polynomial& operator=( T val);
  //@}
    
  // ========================================================================
  //  Access methods
  // ========================================================================
  //@{ \name Access methods

  //! Return dimension of polynomial
  UInt GetDim() const;
  
  //! Return maximum order for given dimension
  
  //! Return the maximum order of the polynomial in a given variable.
  //! \param  index Index to the variable vector
  UInt GetMaxOrder( UInt index );
  
  //! Evaluate polynomial at with given variable values
  
  //! This method evaluates the polynomial at a given location, denoted by
  //! the vars-vector.
  //! \param vars Vector containing the values for all variables
  T Eval( const Vector<T>& vars ) const;
  
  //! Get all coefficients
  COEFMAP& Coefs();

  //! Get all coefficients (const.)
  const COEFMAP& Coefs() const;
  
  //! Remove all entries, which evaluate to 0 (= epsilon value)
  
  //! This method removes all coefficients, which are numerically zero but
  //! still are stored explicitly. This happens e.g. by addition / 
  //! subtraction of several polynomials. 
  //! It is safe to call this method at any time. 
  //! \param tol Tolerance to determine, if a coefficient is assumed to be zero
  void Cleanup( Double tol = 1e-16);

  //! Generate a string representation of the polynomial
  
  //! This method returns a string representation of the polynomial, where
  //! the variable names are taken as 'x', 'y' and 'z'. 
  //! Additionally the output of numerically zero elements can be prevented.
  //! \param omitZeros Flag, to suppress the output of entries, which are
  //!                  numerically zero
  std::string ToString(bool omitZeros = true) const;
  
  //@}
  
  // ========================================================================
  //  Operators
  // ========================================================================
  //@{ \name Operators

  //! operator+= Polynomial
  POL& operator+=( const POL& a );

  //! operator+= Scalar
  POL& operator+=( const T& a );

  //! operator-= Polynomial
  POL& operator-=( const POL& a );

  //! operator-= Scalar
  POL& operator-=( const T& a );

  //! operator*= Polynomial
  POL& operator*=( const POL& a );
  
  //! operator*= Scalar
  POL& operator*=( const T& a );

  //! operator/= Polynomial -> not implemented

  //! operator/= Scalar
  POL& operator/=( const T& a );
  //@}
  
private:
  
  //! Coefficients of the polynomial
  
  //! This map stores all non-zero coefficients of the polynomial (value).
  //! The coefficients are associated to the exponents, described by a
  //! boost array.
  COEFMAP coefs_;
};

// ========================================================================
//  Free Operators
// ========================================================================
//@{ \name Polynomial helper functions
//! -Polynomial
template<typename T, UInt DIM>
Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a );

//! Polynomial + Polynomial
template<typename T, UInt DIM>
Polynomial<T,DIM> operator+( const Polynomial<T,DIM>& a,
                             const Polynomial<T,DIM>& b );

//! Polynomial + Scalar
template<typename T, UInt DIM>
Polynomial<T,DIM> operator+( const Polynomial<T,DIM>& a,
                             const T b );

//! Scalar + Polynomial
template<typename T, UInt DIM>
Polynomial<T,DIM> operator+( const T a,
                             const Polynomial<T,DIM>& b );

//! Polynomial - Polynomial
template<typename T, UInt DIM>
Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a,
                             const Polynomial<T,DIM>& b );

//! Polynomial - Scalar
template<typename T, UInt DIM>
Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a,
                             const T b );

//! Scalar - Polynomial
template<typename T, UInt DIM>
Polynomial<T,DIM> operator-( const T a,
                             const Polynomial<T,DIM>& b );

//! Polynomial * Polynomial
template<typename T, UInt DIM>
Polynomial<T,DIM> operator*( const Polynomial<T,DIM>& a,
                             const Polynomial<T,DIM>& b );

//! Polynomial * Scalar
template<typename T, UInt DIM>
Polynomial<T,DIM> operator*( const Polynomial<T,DIM>& a,
                             const T b );

//! Scalar * Polynomial
template<typename T, UInt DIM>
Polynomial<T,DIM> operator*( const T a,
                             const Polynomial<T,DIM>& b );

//! Polynomial / Polynomial -> not implemented
//! Polynomial / Scalar -> not implemented
//! Scalar / Polynomial -> not implemented

//! Inverse
//! sqr
//! sqrt
//! exp
//! fabs

//! Miscellaneous stuff
//! ostream operator
template<typename T, UInt DIM> 
 std::ostream& operator<<( std::ostream& out, 
                              const Polynomial<T,DIM> a );
//@}
} // namespace CoupledField
#endif
