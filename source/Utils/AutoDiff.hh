#ifndef FILE_CFS_AUTODIFF_HH
#define FILE_CFS_AUTODIFF_HH

#include "General/Environment.hh"

namespace CoupledField {

// helper function
#define AD_LOOP(WORK) \
  for( UInt i = 0; i < DIM; ++i) { WORK; }
 
//! Class for automatic differentiation

//! This class can calculate automatic n-dimensional derivatives.
template<typename T, UInt DIM>
class AutoDiff {
  
public:
  
  //! Typedef for specific AutoDiff instance
  typedef AutoDiff<T,DIM> TA;
  
  // ========================================================================
  //  Constructors & Assignment
  // ========================================================================
  //@{ \name Constructors / Assignment
  
  //! Default constructor (initialize to zero)
  AutoDiff() {
    // initialize object to constant 0
    val_ = 0.0;
    AD_LOOP( deriv_[i] = 0.0; )
  };
  
  //! Constructor for constant value
  AutoDiff( T val ) {
    val_ = val;
    AD_LOOP( deriv_[i] = 0.0; )
  }

  //! Constructor for defined derivative with index
  AutoDiff( T val, UInt derivIndex ) {
    val_ = val;
    AD_LOOP( deriv_[i] = 0.0 );
    deriv_[derivIndex] = 1.0;
  }

  //! Copy constructor
  AutoDiff( const AutoDiff& a ) {
    val_ = a.val_;
    AD_LOOP( deriv_[i] = a.deriv_[i] );
  }

  //! Assignment from constant value
  AutoDiff& operator=( T val) {
    val_ = val;
    AD_LOOP( deriv_[i] = 0.0 );
    return *this;
  }
  //@}
    
  // ========================================================================
  //  Access methods
  // ========================================================================
  //@{ \name Access methods
  
  //! Obtain value (const)
  T Val() const { return val_;}
  
  //! Obtain value
  T& Val() {return val_;}
  
  //! Obtain derivative (const)
  T DVal( UInt i) const {return deriv_[i];} 
  
  //! Obtain derivative
  T& DVal( UInt i) {return deriv_[i];}
  //@}
  
  // ========================================================================
  //  Operators
  // ========================================================================
  //@{ \name Operators

  //! operator+= AutoDiff
  TA& operator+=( const TA& a ) {
    val_ += a.val_;
    AD_LOOP( deriv_[i] += a.deriv_[i] );
    return *this;
  }
  //! operator+= Scalar
  TA& operator+=( const T& a ) {
    val_ += a;
    // derivative of constant is 0
    return *this;
  }

  //! operator-= AutoDiff
  TA& operator-=( const TA& a ) {
    val_ -= a.val_;
    AD_LOOP( deriv_[i] -= a.deriv_[i] );
    return *this;
  }
  //! operator-= Scalar
  TA& operator-=( const T& a ) {
    val_ -= a;
    // derivative of constant is 0
    return *this;
  }

  //! operator*= AutoDiff
  TA& operator*=( const TA& a ) {
    val_ *= a.val_;
    AD_LOOP( deriv_[i] = val_ * a.deriv_[i] + deriv_[i] * a.val_ );
    return *this;
  }
  //! operator*= Scalar
  TA& operator*=( const T& a ) {
    val_ *= a;
    AD_LOOP( deriv_[i] *= a);
    return *this;
  }

  //! operator/= AutoDiff -> not implemented

  //! operator/= Scalar
  TA& operator/=( const T& a ) {
    T inv = 1.0 / a;
    val_ *= inv;
    AD_LOOP( deriv_[i] *= inv);
    return *this;
  }

  //@}
  
private:
  
  //! Value of the variable
  T val_;
  
  //! Derivatives of the variable
  T deriv_[DIM];
};

// ========================================================================
//  Free Operators
// ========================================================================
//@{ \name AutoDiff helper functions
//! -AutoDiff
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator-( const AutoDiff<T,DIM>& a ) {
  AutoDiff<T,DIM> ret;
  ret.Val() = -a.Val();
  AD_LOOP( ret.DVal(i) = -a.DVal(i) ); 
  return ret;
}

 

//! AutoDiff + AutoDiff
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator+( const AutoDiff<T,DIM>& a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  ret.Val() = a.Val() + b.Val();
  AD_LOOP( ret.DVal(i) = a.DVal(i) + b.DVal(i) ); 
  return ret;
}

//! AutoDiff + Scalar
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator+( const AutoDiff<T,DIM>& a,
                                  const T b ) {
  AutoDiff<T,DIM> ret;
  ret.val_ = a.Val() + b;
  AD_LOOP( ret.DVal(i) = a.DVal(i) ); 
  return ret;
}

//! Scalar + AutoDiff
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator+( const T a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  ret.Val() = a + b.Val();
  AD_LOOP( ret.DVal(i) = b.DVal(i) ); 
  return ret;
}

//! AutoDiff - AutoDiff
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator-( const AutoDiff<T,DIM>& a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  ret.Val() = a.Val() - b.Val();
  AD_LOOP( ret.DVal(i) = a.DVal(i) - b.DVal(i) ); 
  return ret;
}

//! AutoDiff - Scalar
template<typename T, UInt DIM>
inline  AutoDiff<T,DIM> operator-( const AutoDiff<T,DIM>& a,
                                   const T b ) {
  AutoDiff<T,DIM> ret;
  ret.Val() = a.Val() - b;
  AD_LOOP( ret.DVal(i) = a.DVal(i) ); 
  return ret;
}

//! Scalar - AutoDiff
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator-( const T a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  ret.Val() = a - b.Val();
  AD_LOOP( ret.DVal(i) = -b.DVal(i) ); 
  return ret;
}

//! AutoDiff * AutoDiff
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator*( const AutoDiff<T,DIM>& a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  const T av = a.Val();
  const T bv = b.Val();
  ret.Val() = av * bv;
  AD_LOOP( ret.DVal(i) = a.DVal(i)* bv + b.DVal(i) * av ); 
  return ret;
}

//! AutoDiff * Scalar
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator*( const AutoDiff<T,DIM>& a,
                                  const T b ) {
  AutoDiff<T,DIM> ret;
  const T av = a.Val();
  ret.Val() = av * b;
  AD_LOOP( ret.DVal(i) = a.DVal(i)* b ); 
  return ret;
}

//! Scalar * AutoDiff
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator*( const T a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  const T bv = b.Val();
  ret.Val() = a * bv;
  AD_LOOP( ret.DVal(i) = a * b.DVal(i) ); 
  return ret;
}

//! AutoDiff / AutoDiff 
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator/( const AutoDiff<T,DIM>& a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  const T av = a.Val();
  const T bv = b.Val();
  ret.Val() = av / bv;
  AD_LOOP( ret.DVal(i) = (a.DVal(i) *  bv - b.DVal(i) * av ) / (bv * bv) ); 
  return ret;
}

//! AutoDiff / Scalar 
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator/( const AutoDiff<T,DIM>& a,
                                  const T b ) {
  AutoDiff<T,DIM> ret;
  const T av = a.Val();
  ret.Val() = av / b;
  AD_LOOP( ret.DVal(i) = a.DVal(i) / b ); 
  return ret;
}

//! Scalar / AutoDiff 
template<typename T, UInt DIM>
inline AutoDiff<T,DIM> operator/( const T a,
                                  const AutoDiff<T,DIM>& b ) {
  AutoDiff<T,DIM> ret;
  const T bv = b.Val();
  ret.Val() = a / bv;
  AD_LOOP( ret.DVal(i) = -( b.DVal(i) * a ) / (bv*bv) ); 
  return ret;
}


//! Inverse
//! sqr
//! sqrt
//! exp
//! fabs


  //! Miscellaneous stuff
  //! ostream operator
template<typename T, UInt DIM> 
inline std::ostream& operator<<( std::ostream& out, 
                                const AutoDiff<T,DIM> a ) {
  out << a.Val();
  out << ", D = ";
  AD_LOOP(out << a.DVal(i) << " ")
  return out;
}
//@}





} // namespace NACS

#undef AD_LOOP
#endif
