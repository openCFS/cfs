#include "Polynomial.hh"

namespace CoupledField {


// ============================================================================
//  H E L P E R    M E T H O D S 
// ============================================================================

//! Helper method for adding up two boost arrays
template<UInt DIM>
std::array<UInt,DIM> Add(const std::array<UInt,DIM>&a, 
                           const std::array<UInt,DIM>&b ) {
  std::array<UInt,DIM> ret;
  for( UInt i=0; i < DIM; ++i ) {
    ret[i] = a[i] + b[i];
  }
  return ret;
}

//! Helper method for subtracting two boost arrays
template<UInt DIM>
std::array<UInt,DIM> Sub(const std::array<UInt,DIM>&a, 
                           const std::array<UInt,DIM>&b ) {
  std::array<UInt,DIM> ret;
  for( UInt i=0; i < DIM; ++i ) {
    ret[i] = a[i] - b[i];
  }
  return ret;
}

//! Helper method to print the exponents
template<UInt DIM>
std::string OrdToString( const std::array<UInt, DIM>& a ) {
  std::stringstream out;
   for( UInt i = 0; i < DIM; ++i ) {
     out << a[i];
   }
   return out.str();
}

//! Helper method to zero an array
template<UInt DIM>
void Zero( std::array<UInt, DIM>& a) {
  for( UInt i = 0; i < DIM; ++i ) {
    a[i] = 0;
  }
}

// ============================================================================
//  P O L Y N O M I A L    C L A S S 
// ============================================================================

template<typename T, UInt DIM> 
Polynomial<T,DIM>::Polynomial() {
  }
  
template<typename T, UInt DIM> 
Polynomial<T,DIM>::Polynomial(const T& val) {
  ORD ord;
  Zero<DIM>(ord);
  coefs_[ord] = val;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>:: Polynomial( const T& val, UInt index, UInt order ) {
  ORD ord;
  Zero<DIM>(ord);
  ord[index] = order;
  coefs_[ord] = val;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>::Polynomial( const Polynomial& a ) {
    coefs_ = a.coefs_;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>::~Polynomial() {
  coefs_.clear();
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator=( T val) {
    coefs_.clear();
    ORD ord;
    Zero<DIM>(ord);
    coefs_[ord] = val;
    return *this;
  }
    
// ========================================================================
//  Access methods
// ========================================================================

template<typename T, UInt DIM> 
UInt Polynomial<T,DIM>::GetMaxOrder( UInt i) {
  assert(i < DIM);
  // make sure to clean up the polynomial before
  Cleanup();

  UInt maxOrder = 0;
  COEF_IT it = coefs_.begin();
  for( ; it != coefs_.end(); ++it ) {
    const ORD & ord = it->first;
    maxOrder = std::max( maxOrder, ord[i]);
  }
  return maxOrder;
}
  

template<typename T, UInt DIM> 
T Polynomial<T,DIM>::Eval( const Vector<T>& vars ) const {
  EXCEPTION( "Not implemented yet.");
  return T();
}

template<typename T, UInt DIM> 
typename Polynomial<T,DIM>::COEFMAP& Polynomial<T,DIM>::Coefs() {
  return coefs_;
}

template<typename T, UInt DIM> 
const typename Polynomial<T,DIM>::COEFMAP& Polynomial<T,DIM>::Coefs() const {
  return coefs_;
}

template<typename T, UInt DIM>
void Polynomial<T,DIM>::Cleanup( Double tol ) {
  std::set<std::array<UInt,DIM> > del;
  COEF_IT it = coefs_.begin();
  for( ; it != coefs_.end(); ++it ) {
    if(  std::abs(it->second) < tol ) {
      del.insert(it->first);
    }
  }
  typename std::set<std::array<UInt,DIM> >::const_iterator it2 = del.begin();
  for( ; it2 != del.end(); ++it2 ) {
    coefs_.erase(*it2);
  }
}

template<typename T, UInt DIM> 
std::string Polynomial<T,DIM>::ToString(bool omitZeros ) const {
  std::stringstream out;

  // use monomial names 'x', 'y', 'z' to describe the polynomial
  StdVector<std::string> names;
  names = "x", "y", "z";
  COEF_CONST_IT it = coefs_.begin();
  for( ; it != coefs_.end(); ++it ) {

    if( std::abs(it->second) < 1e-16 && omitZeros) 
      continue;

    // print sign: + or -
    if( it->second > 0.0 ) {
      out << " + ";
    } else {
      out << " - ";
    }

    // print coefficient value
    out << std::abs(it->second);
    for( UInt i = 0; i < DIM; ++i ) {
      UInt order = it->first[i];

      // print monomial
      if (order > 0 ) {
        out << names[i];
      }
      if( order > 1 ) {
        out << "^" << order;
      }
    }
  }

  return out.str();
}

// ========================================================================
//  Operators
// ========================================================================
//@{ \name Operators

template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator+=( const Polynomial<T,DIM>& a ) {
  COEF_CONST_IT it = a.coefs_.begin();
  for( ; it != a.coefs_.end(); ++it ) {
    coefs_[it->first] += it->second;
  }
  return *this;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator+=( const T& a ) {
  std::array<UInt,DIM> zero;
  Zero<DIM>(zero);
  coefs_[zero] += a;
  return *this;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator-=( const Polynomial<T,DIM>& a ) {
  COEF_CONST_IT it = a.coefs_.begin();
  for( ; it != a.coefs_.end(); ++it ) {
    coefs_[it->first] -= it->second;
  }
  return *this;
}
template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator-=( const T& a ) {
  std::array<UInt,DIM> zero;
  Zero<DIM>(zero);
  coefs_[zero] -= a;
  return *this;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator*=( const Polynomial<T,DIM>& a ) {
  COEF_IT itSelf = coefs_.begin();
  COEF_CONST_IT itA = a.coefs_.begin(); 
  COEFMAP cp;
  ORD ord;
  for( ; itSelf != coefs_.end(); ++itSelf ) {
    for( itA=a.coefs_.begin(); itA != a.coefs_.end(); ++itA ) {
      ord = Add<DIM>(itSelf->first, itA->first);
      cp[ord] += itSelf->second * itA->second; 
    }
  }
  coefs_ = cp;
  return *this;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator*=( const T& a ) {
  typename Polynomial<T,DIM>::COEF_IT itSelf = coefs_.begin();
  for( ; itSelf != coefs_.end(); ++itSelf ) {
    itSelf->second *= a;
  }
  return *this;
}

template<typename T, UInt DIM> 
Polynomial<T,DIM>& Polynomial<T,DIM>::operator/=( const T& a ) {
  typename Polynomial<T,DIM>::COEF_IT itSelf = coefs_.begin();
  for( ; itSelf != coefs_.end(); ++itSelf ) {
    itSelf->second /= a;
  }
  return *this;
  }


// ========================================================================
//  Free Operators
// ========================================================================
template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a ) {
  Polynomial<T,DIM> ret(a);
  typename Polynomial<T,DIM>::COEF_CONST_IT it = a.Coefs().begin();
  for( ; it != a.Coefs().end(); ++it ) {
    ret.Coefs()[it->first] = -(it->second);
  }
  return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator+( const Polynomial<T,DIM>& a,
                                    const Polynomial<T,DIM>& b ) {
  Polynomial<T,DIM> ret(a);
  typename Polynomial<T,DIM>::COEF_CONST_IT itB = b.Coefs().begin();
  for( ; itB != b.Coefs().end(); ++itB ) {
    ret.Coefs()[itB->first] += (itB->second);
  }
  return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator+( const Polynomial<T,DIM>& a,
                                    const T b ) {
  Polynomial<T,DIM> ret(a);
  
  // add constant part
  typename Polynomial<T,DIM>::ORD ord;
  Zero<DIM>( ord );
  ret.Coefs()[ord] += b;
  return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator+( const T a,
                                    const Polynomial<T,DIM>& b ) {
  Polynomial<T,DIM> ret(b);
  
  // add constant part
  typename Polynomial<T,DIM>::ORD ord;
  Zero<DIM>(ord);
  ret.Coefs()[ord] += a;
  return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a,
                                    const Polynomial<T,DIM>& b ) {
  Polynomial<T,DIM> ret(a);
  typename Polynomial<T,DIM>::COEF_CONST_IT itB = b.Coefs().begin();
  for( ; itB != b.Coefs().end(); ++itB ) {
    ret.Coefs()[itB->first] -= (itB->second);
  }
  return ret;
}

template<typename T, UInt DIM>
inline  Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a,
                                     const T b ) {
  Polynomial<T,DIM> ret(a);

  // subtract constant part
  typename Polynomial<T,DIM>::ORD ord;
  Zero<DIM>(ord);
  ret.Coefs()[ord] -= b;
  return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator-( const T a,
                                    const Polynomial<T,DIM>& b ) {
  Polynomial<T,DIM> ret;
   typename Polynomial<T,DIM>::COEF_CONST_IT itB = b.Coefs().begin();

   // first, reverse sign
   for( ; itB != b.Coefs().end(); ++itB ) {
     ret.Coefs()[itB->first] = -(itB->second);
   }
   // than add constant part
   typename Polynomial<T,DIM>::ORD ord;
   Zero<DIM>(ord);
   ret.Coefs()[ord] += a;
   
   return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator*( const Polynomial<T,DIM>& a,
                                    const Polynomial<T,DIM>& b ) {
  Polynomial<T,DIM> ret;
  typename Polynomial<T,DIM>::COEF_CONST_IT itA = a.Coefs().begin();
  typename Polynomial<T,DIM>::COEF_CONST_IT itB = b.Coefs().begin();
  typename Polynomial<T,DIM>::ORD ord;
  for( ; itA != a.Coefs().end(); ++itA) {
    for( itB=b.Coefs().begin(); itB != b.Coefs().end(); ++itB ) {
      ord = Add<DIM>(itA->first, itB->first);
      ret.Coefs()[ord] += itA->second * itB->second; 
    }
  }
  return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator*( const Polynomial<T,DIM>& a,
                                    const T b ) {
  Polynomial<T,DIM> ret(a);
  typename Polynomial<T,DIM>::COEF_CONST_IT itA = a.Coefs().begin();
  for( ; itA != a.Coefs().end(); ++itA) {
    ret.Coefs()[itA->first] = itA->second * b;
  }
  return ret;
}

template<typename T, UInt DIM>
inline Polynomial<T,DIM> operator*( const T a,
                                    const Polynomial<T,DIM>& b ) {
  Polynomial<T,DIM> ret(b);
  typename Polynomial<T,DIM>::COEF_CONST_IT itB = b.Coefs().begin();
  for( ; itB != b.Coefs().end(); ++itB) {
    ret.Coefs()[itB->first] = a* itB->second;
  }
  return ret;
}

template<typename T, UInt DIM> 
std::ostream& operator<<( std::ostream& out, 
                              const Polynomial<T,DIM> a ) {
   out << a.ToString(true);
   return out;
}

// ========================================================================
//  EXPLICIT TEMPLATE INSTANTIATION
// ========================================================================
#define INSTANTIATE_POLY( T, DIM)                                   \
template class Polynomial<T, DIM>;                                  \
template Polynomial<T,DIM> operator-( const Polynomial<T,DIM>&);    \
template Polynomial<T,DIM> operator+( const Polynomial<T,DIM>& a,   \
                                      const Polynomial<T,DIM>& b ); \
template Polynomial<T,DIM> operator+( const Polynomial<T,DIM>& a,   \
                                      const T b );                  \
template Polynomial<T,DIM> operator+( const T a,                    \
                                      const Polynomial<T,DIM>& b ); \
template Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a,   \
                                      const Polynomial<T,DIM>& b ); \
template  Polynomial<T,DIM> operator-( const Polynomial<T,DIM>& a,  \
                                       const T b );                 \
template Polynomial<T,DIM> operator-( const T a,                    \
                                      const Polynomial<T,DIM>& b ); \
template Polynomial<T,DIM> operator*( const Polynomial<T,DIM>& a,   \
                                      const Polynomial<T,DIM>& b ); \
template Polynomial<T,DIM> operator*( const Polynomial<T,DIM>& a,   \
                                      const T b );                  \
template Polynomial<T,DIM> operator*( const T a,                    \
                                      const Polynomial<T,DIM>& b ); \
template std::ostream& operator<<( std::ostream& out,               \
                                    const Polynomial<T,DIM> a );
INSTANTIATE_POLY(Double, 1)
INSTANTIATE_POLY(Double, 2)
INSTANTIATE_POLY(Double, 3)

#undef INSTANTIATE_POLY

// ========================================================================
//  TEST CASES
// ========================================================================
#ifdef TEST_POLYNOMIAL_CLASS
void PolynomialTest() {
  
  
  // define polynomial
  Polynomial<Double,2> x( 1.0, 0), y( 1.0, 1);
  Polynomial<Double,2> a( 1.0, 0), b( 1.0, 1); 
  Polynomial<Double,2> res;
  
  std::cerr << "x is " << x << std::endl << std::endl;
  std::cerr << "y is " << y << std::endl << std::endl;
  a = x;
    std::cerr << "a = x is " << a << std::endl << std::endl;
  a = -x;
  std::cerr << "a = -x is " << a << std::endl << std::endl;
  a = (x+y)+5.;
  std::cerr << "a = x+y+5 is " << a << std::endl << std::endl;

  a*=a;
  std::cerr << "a*a is " << a << "\n\n";
  
  a = x*5.0 + y*3.0 - x*y + 2.0;
  b = -1.0 + x*x*y;
  std::cerr << "x is: " << x << std::endl;
  std::cerr << "y is: " << y << std::endl;
  
  std::cerr << a.ToString();
  std::cerr << "a is: " << a << std::endl;
  std::cerr << "b is: " << b << std::endl;
  
  res = a*b;
  std::cerr << "a*b: " << res << std::endl;
  
  std::cerr << "=== CONST TEST === " << std::endl;
  a= 5.0;
  res = a*b;
  std::cerr << "a: " << a.ToString() << std::endl;
  std::cerr << "b: " << b.ToString() << std::endl;
  std::cerr << "a*b = " << res.ToString() << std::endl << std::endl;
  
  std::cerr << "=== ZERO TEST === " << std::endl;
  a = 2.0*x  + 1.*y;
  b = -2.0*x + 0.5*y*y;
  res = a+b;
  std::cerr << "a+b unnormalized: " << res.ToString(false) << std::endl;
  res.Cleanup();
  std::cerr << "a+b unnormalized: " << res.ToString(false) << std::endl;
  
  
  std::cerr << "maximum x order: " << res.GetMaxOrder(0) << std::endl;
  std::cerr << "maximum y order: " << res.GetMaxOrder(1) << std::endl;
  
  
  std::cerr << "=== LAMBDA TEST === " << std::endl;
  a= 0.25*(1.0-x)*(1.0-y);
  b = 0.25*(1.0+x)*(1.0-y);
  res = a+b;
  std::cerr << "a: " << a.ToString() << std::endl;
  std::cerr << "b: " << b.ToString() << std::endl;
  std::cerr << "a+b = " << res.ToString() << std::endl << std::endl;
  
  std::cerr << "=== SIGMA TEST === " << std::endl;
  a= 0.5*((1.0-x)+(1.0-y));
  b = 0.5*((1.0+x)+(1.0-y));
  res = b-a;
  std::cerr << "a: " << a.ToString() << std::endl;
  std::cerr << "b: " << b.ToString() << std::endl;
  std::cerr << "b-a = " << res.ToString() << std::endl << std::endl;
  
}
#endif

} // namespace CoupledField
