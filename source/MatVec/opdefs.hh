#ifndef OLAS_OPDEFS_HH
#define OLAS_OPDEFS_HH

#include <stdio.h>

#include "TypeDefs.hh"
#include "General/Exception.hh"

//! \file opdefs.hh

//! This file contains specialized implementations for the 
//! primitive data types in OLAS. Primitive data types are
//! scalar real values (Double), scalar complex values
//! (Complex) and tiny vectors/matrices of these data types.
//! Using specialized functions allows
//! to implement algorithms for all data types, i.e. by using the
//! optype::Invert function instead of writing 1/x.
//! 
//! \note This file has to be read bottom up in order to be understood
//! since it contains some meta loops which are terminated at
//! the bottom.
//!       
//! \note Be aware that the tvmet matrices are 0-based. This file is
//! therefore an exception to the general rule that OLAS data
//! structures are 1-based

//! Data type for integer template arguments in tvmet

//! This type definition must match the data type used in tvmet
typedef std::size_t _T_int_;


namespace CoupledField {


  // ========================================================================
  //
  // AUXILLIARY FUNCTIONS
  //
  // ========================================================================


  //! Computation of absolute value for non-scalar types

  //! Computation of absolute value for non-scalar types. This functionality
  //! is currently not implemented.
  template <typename T>
  inline Double Abs(T a){ 
    EXCEPTION("Abs not implemented for this data type ");
  }

  template <typename T>
  inline T Sqrt(T a){ 
    EXCEPTION("Sqrt not implemented for this data type ");
  }
  template<>
  inline Double Sqrt<Double> (Double a){
    return sqrt(a);
  }
#ifndef USE_ADOLC
   template<>
  inline Complex Sqrt<Complex> (Complex a){
    return sqrt(a);
  }
#else
  template<>
  inline Complex Sqrt<Complex> (Complex a){
    EXCEPTION("Sqrt not implemented for this data type ");
  }
#endif

  //! Specialize Abs to compute absolute value of a Double variable
  template<>
  inline Double Abs<Double> (Double a){
    return fabs(a);
  }

  template<>
  inline Double Abs<Integer> (Integer a){
    return abs(a);
  }

  /** this is the reason why we cannot use std::abs() as there is no definition and the compiler complaines :( */
  template<>
  inline Double Abs<UInt> (UInt a){
    return a;
  }

  //! Specialize Abs to compute modulus of a Complex variable
  template<>
  inline Double Abs<Complex>(Complex a) {
    return sqrt((a.real())*(a.real())+(a.imag())*(a.imag()));
  }

  //! Compute \f$z\cdot\bar{z}\f$
  template <typename T>
  inline Double Abs2( const T &arg ) {
    return arg * arg;
  }

  //! Compute \f$z\cdot\bar{z}\f$
  template<>
  inline Double Abs2<Complex>( const Complex &arg ) {
    return (arg.real()) * (arg.real()) + (arg.imag()) * (arg.imag());
  }

  /** this is defined only in C++11 - how can this idiots take so long?
   * defining the same for complex leads to overloading problems with std::conj with icc 15 even w/o using std::conj */
  inline static Double conj(const Double &a1) {
    return a1;
  }

  template <typename T>
  inline Double Real(T a) {  return ((std::complex<Double>) a).real(); }

  template <>
  inline Double Real<Double>(Double a) { return a; }

  template <>
  inline Double Real<std::complex<Double> >(std::complex<Double> a) { return a.real(); }

  template <typename T>
  inline Double Imag(T a) { return ((std::complex<Double>) a).imag(); }

  template <>
  inline Double Imag<Double>(Double a) { return 0; }

  template <>
  inline Double Imag<std::complex<Double> >(std::complex<Double> a) { return a.imag(); }
  
  template <typename T>
  inline bool IsZero(T a) { return a == 0; }

  template <>
  inline bool IsZero<Double>(Double a) { return a == 0.0; }

  template <>
  inline bool IsZero<std::complex<Double> >(std::complex<Double> a) { return a.real() == 0.0 && a.imag() == 0.0; }

  template <typename T>
  inline void PrintSingleEntry( T val, FILE *fp ) {
    EXCEPTION("PrintSingleEntry not implemented for this data type ");
  }

  //! Print a single Double value to a file
  template<>
  inline void PrintSingleEntry<Double>( Double val, FILE *fp ) {
#ifndef USE_ADOLC
    fprintf( fp, " % 22.16E", val );
#else
    fprintf( fp, " % 22.16E", val.getValue() );
#endif
}
  
  //! Print a single complex value to a file
  template<>
  inline void PrintSingleEntry<Complex>( Complex val, FILE *fp ) {
#ifndef USE_ADOLC
    fprintf( fp, " % 22.16E % 22.16E", (val.real()), (val.imag()) );
#else
     fprintf( fp, " % 22.16E % 22.16E", (val.real()).getValue(), (val.imag()).getValue() );
#endif
  }

  //! Dummy implementation of multiplication of a real value and complex
  //! scalar. We cannot scale a real value by a complex one and get a real
  //! result.
  template <typename T>
  inline void MultScalarWithComplex( T &multiplicant,
                                     const Complex factor ) {
    EXCEPTION("We cannot multiply a real-valued scalar with a "
             << "complex-valued one! The result would not be Double!");
  }

  //! Multiplication of two complex scalars
  template<>
  inline void MultScalarWithComplex<Complex>( Complex &multiplicant,
                                              const Complex factor ) {
    multiplicant *= factor;
  }

  //! Division of real-value scalar by complex-valued one

  //! Dummy implementation of division of a real-valued scalar by a
  //! complex-valued one. We cannot perform this operation in the realm
  //! of real-numbers, since the result would be complex.
  template <typename T>
  inline void DivScalarByComplex( T &nominator,
                                  const Complex denominator ) {
    EXCEPTION("We cannot divide a real-valued scalar by a "
             << "complex-valued one! The result would not be Double!");
  }

  //! Division of one complex scalar by another
  template<>
  inline void DivScalarByComplex<Complex>( Complex &nominator,
                                           const Complex denominator ) {
    nominator /= denominator;
  }


  // ========================================================================
  //
  // ARITHMETIC OPERATIONS
  //
  // ========================================================================




  // *************************************************************************
  //   OPTYPE for DOUBLE/COMPLEX SCALARS
  // *************************************************************************

  //! Struct containing standard arithmetic operations

  //! This struct is a container for standard arithmetic operations like e.g.
  //! inversion or the dot product. This version is intended for scalar
  //! variables of type Double or Complex.
  template<typename T>
  struct OpType{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //bjurgel AD conversion comment:                                                                             //
    // this is a bad idea. These operators return  Double or Complex etc. However, I've seen a few cases         //
    // where you are actually using the type T after this call.                                                  //
    // This means you are relying on (automatic) implicit conversion of Double to T. Bad style and bad for AD.   //
    // You probably want to return the type T instead but do you really? How can we be sure?                     //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

      
    //! Compute the inverse of the input argument.

    //! This method computes the inverse of the input argument. It is intended
    //! for variables of Double and Complex type, where the inverse of arg is
    //! given as 1/arg.
    inline static T invert(const T& arg) {
      return static_cast<T>(T(1.0)/arg);
    }

    inline static Double dotProduct(const Double &a1, const Double &a2 ) {
      return a1 * a2;
    }

    inline static Complex dotProduct(const Complex &a1, const Complex &a2 ) {
      return a1 * std::conj(a2);
    }

    //! Computation of "Euclidean Norm" of Double or Complex scalar
    inline static Double NormL2( const T& a ) {
      return Abs( a );
    }

    //! Compute \f$z\cdot\bar{z}\f$
    inline static Double zConjz( const T &z ) {
      return Abs2( z );
    }

    /** this implementation does noting als than mutiplying two scalars */
    inline static T multT( const T& A, const T& x) {
      return A*x;
    }

    //! Method needed for exporting matrices to a file
    inline static void ExportEntry( const T &val, Integer i, Integer j,
				    FILE *fp ) {
      PrintSingleEntry( val, fp );
    }

    //! Method needed for exporting vectors to a file
    inline static void ExportEntry( const T &val, Integer i, FILE *fp ) {
      PrintSingleEntry( val, fp );
    }

    //! this method returns the absolute maximum diagonal entry 
    //! for tiny matrices or the absolute value of scalars
    inline static Double MaxDiag(const T &v) {
      return Abs(v);
    }

    //! This method multiplies the scalar with a complex factor.

    //! This method multiplies the scalar with a complex factor. This will
    //! only work for a complex-valued scalar and not for a real-valued one.
    inline static void MultWithComplex( T &multiplicant,
                                        const Complex factor ) {
      MultScalarWithComplex( multiplicant, factor );
    }

    //! This method multiplies the scalar with a complex factor.

    //! This method multiplies the scalar with a complex factor. This will
    //! only work for a complex-valued scalar and not for a real-valued one.
    inline static void DivByComplex( T &nominator,
                                     const Complex denominator ) {
      DivScalarByComplex( nominator, denominator );
    }

  };

}


#endif // OLAS_OPDEFS_HH
