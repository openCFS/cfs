#ifndef OLAS_OPDEFS_HH
#define OLAS_OPDEFS_HH

#include <stdio.h>

#include "utils/utils.hh"
#include "matvec/typedefs.hh"

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


namespace OLAS {




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
  struct opType{

    //! Compute the inverse of the input argument.

    //! This method computes the inverse of the input argument. It is intended
    //! for variables of Double and Complex type, where the inverse of arg is
    //! given as 1/arg.
    inline static T invert(const T& arg) {
      return ((T)1.0)/arg;
    }

    //! Inner product for real- and complex-valued scalars
    inline static T dotProduct( const T &a1, const T &a2 ) {
      return a1 * Conj(a2);
    }

    //! Computation of "Euclidean Norm" of Double or Complex scalar
    inline static Double NormEuclid( const T& a ) {
      return Abs( a );
    }

    //! Compute \f$z\cdot\bar{z}\f$
    inline static Double zConjz( const T &z ) {
      return Abs2( z );
    }

    //! this function takes the transpose of the first argument A
    //! and multiplies it by the second argument x. This is the same
    //! as the standard product defined on the two types T_Mtype and
    //! T_Vtype unless T_Mtype is an actual tiny matrix in which case
    //! the transpose of the matrix is used. The result is stored in y.
    inline static T multT( const T &A,
                           const typename assocType<T>::T_Vtype &x ) {
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




  // ========================================================================
  //
  // AUXILLIARY FUNCTIONS
  //
  // ========================================================================


  //! Computation of absolute value for non-scalar types

  //! Computation of absolute value for non-scalar types. This functionality
  //! is currently not implemented.
  template <typename T>
  static T Abs(T a){ 
    Error("Abs not implemented for this data type ",__FILE__,__LINE__);
  }

  //! Compute absolute value of a double variable
  static Double Abs(Double a){
    return fabs(a);
  }

  //! Compute modulus of a Complex variable
  static Double Abs(Complex a) {
    return sqrt((a.real())*(a.real())+(a.imag())*(a.imag()));
  }

  //! Compute \f$z\cdot\bar{z}\f$
  inline static Double Abs2( const Double &arg ) {
    return arg * arg;
  }

  //! Compute \f$z\cdot\bar{z}\f$
  inline static Double Abs2( const Complex &arg ) {
    return (arg.real()) * (arg.real()) + (arg.imag()) * (arg.imag());
  }

  //! Compute the conjugate of a real number
  static Double Conj( Double a ) {
    return a;
  }

  //! Compute the conjugate of a complex number
  static Complex Conj( Complex a ) {
    return std::conj(a);
  }
  
  //! Print a single double value to a file
  static void PrintSingleEntry( Double val, FILE *fp ) {
    fprintf( fp, "% 22.16e", val );
  }
  
  //! Print a single complex value to a file
  static void PrintSingleEntry( Complex val, FILE *fp ) {
    fprintf( fp, "% 22.16e\t% 22.16e", (val.real()), (val.imag()) );
  }

  //! Multiplication of two complex scalars
  static void MultScalarWithComplex( Complex &multiplicant,
				     const Complex factor ) {
    multiplicant *= factor;
  }

  //! Dummy implementation of multiplication of a real value and complex
  //! scalar. We cannot scale a real value by a complex one and get a real
  //! result.
  static void MultScalarWithComplex( Double &multiplicant,
				     const Complex factor ) {
    (*error) << "We cannot multiply a real-valued scalar with a "
             << "complex-valued one! The result would not be Double!";
    Error( __FILE__, __LINE__ );
  }

  //! Division of one complex scalar by another
  static void DivScalarByComplex( Complex &nominator,
                                  const Complex denominator ) {
    nominator /= denominator;
  }

  //! Division of real-value scalar by complex-valued one

  //! Dummy implementation of division of a real-valued scalar by a
  //! complex-valued one. We cannot perform this operation in the realm
  //! of real-numbers, since the result would be complex.
  static void DivScalarByComplex( Double &nominator,
                                  const Complex denominator ) {
    (*error) << "We cannot divide a real-valued scalar by a "
             << "complex-valued one! The result would not be Double!";
    Error( __FILE__, __LINE__ );
  }

}


#endif // OLAS_OPDEFS_HH
