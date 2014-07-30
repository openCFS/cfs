#ifndef FILE_CFS_POLYNOMIALS_HH
#define FILE_CFS_POLYNOMIALS_HH

namespace CoupledField {

  //! \file Define different types of hierarchic Polynomials
  

// ========================================================================
//   L E G E N D R E    P O L Y N O M I A L    C L A S S E S
// ========================================================================
//
// Note: All recurrence formulas are according to the formulas of the
// PhD of Sabine Zaglmayr, 2006

// ------------------------
//  Legendre Polynomials
// ------------------------
//@{ \name Legendre Polynomials

//! Evaluate Legendre polynomials

//! This functions calculated the Legendre basis functions
//! from p=0 up to the given order. 
//!\param values (output) return values (length: order + 1)
//!\param order (input) order of the polynomials (>= 0)
//!\param loc (input) location of the polynomial to be evaluated [-1,1]
//!
//! \note The length of the values array is always (p+1) 
template <typename T_SCAL, class T_VEC>
inline void Legendre( T_VEC& values, UInt order, T_SCAL loc ) {
  T_SCAL p1 = 1.0;
  T_SCAL p2 =  0;
  T_SCAL p3 =  0;
  values.Resize(order+1);
  //if( order >= 0 ) {
    values[0] = 1.0;
  //}
  for (UInt j=1; j<=order; j++) {
    double inv = 1.0 / (double)j;
    p3=p2;
    p2=p1;
    p1=( (2*(double)j-1) * loc * p2 - ((double)j-1) * p3) *inv;
    values[j] = p1;
  }
}
//@}

// ---------------------------------
//  Integrated Legendre Polynomials
// ---------------------------------
//@{ \name Integrated Legendre Polynomials
//! Evaluate integrated Legendre polynomials

//! This functions calculated the integrated Legendre basis functions
//! from p=0 up to the given order. 
//! \param values (output) return values (length: order + 1)
//! \param order (input) order of the polynomials (>= 0)
//! \param loc (input) location of the polynomial to be evaluated [-1,1]
//!
//! \note The length of the values array is always (p+1) 
template <typename T_SCAL, class T_VEC>
inline void IntLegendre( T_VEC& values, UInt order, T_SCAL loc ) {
  T_SCAL p1 = -1.0;
  T_SCAL p2 =  0;
  T_SCAL p3 =  0;
  values.Resize(order+1);
  //  if( order >= 0 ) {
    values[0] = -1.0;
    //  }
  for (UInt j=1; j<=order; j++) {
    double inv = 1.0 / (double)j;
    p3=p2;
    p2=p1;
    p1=( (2*(double)j-3) * loc * p2 - ((double)j-3) * p3) *inv;
    values[j] = p1;
  }
}

//! Evaluate integrated Legendre polynomials (for p >= 2) 

//! This functions calculated the integrated Legendre basis functions
//! from p=2 up to the given \a order. 
//! \param values (output) return values (length: order-1)
//! \param order (input) order of the polynomials (>= 2)
//! \param loc (input) location of the polynomial to be evaluated [-1,1]
//!
//! \note The length of the values array is always (p-1), as the 0th order 
//!       and linear 1st order  shape functions are not included. For
//!       order =0,1, the return vector will be empty.   
template <typename T_SCAL, class T_VEC>
inline void IntLegendreP2( T_VEC& values, UInt order, T_SCAL loc ) {
  T_SCAL p1 =  loc;
  T_SCAL p2 = -1.0;
  T_SCAL p3 =  0.0;
  values.Resize(order-1);
  for (UInt j=2; j<=order; j++) {
    double inv = 1.0 / (double)j;
    p3=p2;
    p2=p1;
    p1=( (2*(double)j-3) * loc * p2 - ((double)j-3) * p3) *inv;
    values[j-2] = p1;
  }
}
//@}

// -----------------------------
//  Scaled Legendre Polynomials
// -----------------------------
//@{ \name Scaled Legendre Polynomials
//! Evaluate scaled Legendre polynomials

//! This functions calculates the scaled Legendre basis functions
//! from p=0 up to the given order. Those functions are used within
//! the triangular, tetrahedral, wedge and pyramidal elements
//! \param values (output) return values (length: order + 1)
//! \param order (input) order of the polynomials (>= 0)
//! \param scal (input) scaling parameter
//! \param loc (input) location of the polynomial to be evaluated [-1,1]
//!
//! \note The length of the values array is always (p+1) 

template <typename T_SCAL, class T_VEC>
inline void ScaledLegendre( T_VEC& values, UInt order, T_SCAL scal, T_SCAL loc ) {
  T_SCAL p1 = loc;
  T_SCAL p2 = 1.0;
  T_SCAL p3 = 0;
  T_SCAL tsquare = scal * scal;
  values.Resize(order+1);

  //  if( order >= 0 ) {
    values[0] = 1.0;
    //  }
  if( order >= 1 ) {
    values[1] = loc;
  }
  for (UInt j=2; j<=order; j++) {
    double inv = 1.0 / (double)j;
    p3=p2;
    p2=p1;
    p1=( (2*(double)j-1) * loc * p2 - ((double)j-1) * p3 * tsquare) *inv;
    values[j] = p1;
   }
}
//@}

// -------------------------------
//  Scaled Integrated Polynomials
// -------------------------------
//@{ \name Scaled Integrated Legendre Polynomials
//! Evaluate scaled integrated Legendre polynomials

//! This functions calculates the scaled integrated Legendre basis functions
//! from p=0 up to the given order. Those functions are used within
//! the triangular, tetrahedral, wedge and pyramidal elements
//! \param values (output) return values (length: order + 1)
//! \param order (input) order of the polynomials (>= 0)
//! \param scal (input) scaling parameter
//! \param loc (input) location of the polynomial to be evaluated [-1,1]
//!
//! \note The length of the values array is always (p+1) 
template <typename T_SCAL, class T_VEC>
inline void ScaledIntLegendre( T_VEC& values, UInt order, 
                               T_SCAL scal, T_SCAL loc ) {
  T_SCAL p1 = -1.0;
  T_SCAL p2 =  0;
  T_SCAL p3 =  0;
  T_SCAL tsquare = scal * scal;
  values.Resize(order+1);
  //  if( order >= 0 ) {
    values[0] = -1.0;
    //  }
  for (UInt j=1; j<=order; j++) {
    double inv = 1.0 / (double)j;
    p3=p2;
    p2=p1;
    p1=( (2*(double)j-3) * loc * p2 - ((double)j-3) * p3 * tsquare) *inv;
    values[j] = p1;
  }
}

//! Evaluate scaled integrated Legendre polynomials (for p >= 2) 

//! This functions calculates the scaled integrated Legendre basis functions
//! from p=2 up to the given order. Those functions are used within
//! the triangular, tetrahedral, wedge and pyramidal elements
//! \param values (output) return values (length: order + 1)
//! \param order (input) order of the polynomials (>= 0)
//! \param scal (input) scaling parameter
//! \param loc (input) location of the polynomial to be evaluated [-1,1]
//!
//! \note The length of the values array is always (p-1) 
template <typename T_SCAL, class T_VEC>
inline void ScaledIntLegendreP2( T_VEC& values, UInt order, 
                                 T_SCAL scal, T_SCAL loc ) {
  T_SCAL p1 =  loc;
  T_SCAL p2 = -1.0;
  T_SCAL p3 =  0;
  T_SCAL tsquare = scal * scal;
  values.Resize(order-1);
  for (UInt j=2; j<=order; j++) {
    double inv = 1.0 / (double)j;
    p3=p2;
    p2=p1;
    p1=( (2*(double)j-3) * loc * p2 - ((double)j-3) * p3 * tsquare) *inv;
    values[j-2] = p1;
  }
}
//@}


// =======================================================
//   T R I A N G U L A R   S H A P E   F U N C T I O N S
// =======================================================

template <typename T_SCAL, class T_VEC>
inline UInt TriaInnerLegendre( T_VEC& values,
                               const UInt& pos, UInt order,
                               const T_SCAL& lambda1, const T_SCAL& lambda2,
                               const T_SCAL& lambda3 ) {
  T_VEC f1, f2;
  UInt nfct = 0;
  UInt myPos = pos;
  ScaledIntLegendreP2(f1, order-1, lambda2+lambda1,lambda2-lambda1);
  Legendre(f2, order-3, lambda3*2.0 - T_SCAL(1));
  for( UInt i = 0; i <= order - 3; ++i ) {
    for( UInt j = 0; j <= order - 3 - i; ++j ) {
      values[myPos++] = f1[i] * f2[j] * lambda3;
      nfct++;
    }
  }
  return nfct;
}


template <typename T_SCAL, class T_VEC>
inline void TriaInnerLegendre2( T_VEC& values, UInt order,
                               T_SCAL y, T_SCAL x ) {
  
  // -----------
  //  Variant 2:
  // -----------
  // The genuine thing about the re.orientation is the following
  // the new x-axis is defined as: 
  //      x' = y-x  (135 deg  direction, i.e. along edge 2)
  // the new y-axis is defined as: 
  //      y' = 1 - x -y (45 deg direction)
  // so both axes are again orthogonal.
  //    
  //  ^
  //  | y-oringal
  //  +
  //  |\
  //  |  \            
  //  |    \      ... y' direction
  //  |     .\    \\\ x' direction (get scaled with increasing y')
  //  |  .     \         
  //  |.         \    
  //  +-----------+---> x-original
  //
  // The bubble function evaluates to 4*x*y*(1 - x - y), i.e. it
  // is zero on all edges, so we can use standard Legendre-based
  // integrals in the interior.
  //
  // Note: Although we could use the integrated Legendre polynomials,
  // the spatial variation using the non-integrated ones seems
  // to be more uniformly distributed. We should check, whether the 
  // condition number get influenced by this
  
  
  UInt nVals = (order-2)*(order-1)/2;
  values.Resize(nVals);
  UInt pos = 0;
  T_VEC xVals, yVals;
  
  ScaledLegendre( xVals, order-3, 1-y, x );
  Legendre( yVals, order-3, 2*y-1);
  
  T_SCAL blend = y * (1-x-y) * (1+x+y); 
  
  for( UInt i = 0; i <= order - 3; ++i ){
    for( UInt j = 0; j <= order - 3- i; ++j ) {
      values[pos++] = blend * xVals[i] * yVals[j];
    }
  }
}



#ifdef CFS_POLYNOMILAS_TEST
void TestPolys() {
  // ==============
  //  TEST SECTION FOR POLYNOMIALS
  // ==============
  StdVector<Double> Leg, IntLeg, IntLegP2, ScaledLeg, ScaledIntLeg; 
  StdVector<Double> ScaledIntLegP2, LegNetgen;
  UInt nPoints = 100;
  Double dx = 2.0/nPoints;
  UInt order = 5;
  Double scal = 1;
  std::ofstream f;
  f.open("./polys.txt");
  for( UInt i = 0; i < nPoints; ++i ) {
    Double x = -1.0 + i*dx;
    Legendre(Leg,order, x);
    IntLegendre(IntLeg,order, x);
    IntLegendreP2(IntLegP2,order, x);
    ScaledLegendre(ScaledLeg,order,scal, x);
    ScaledIntLegendre(ScaledIntLeg,order,scal, x);
    ScaledIntLegendreP2(ScaledIntLegP2,order,scal, x);
    f << x << "\t";
    f << Leg.Last() << "\t";
    f << IntLeg.Last() << "\t";
    f << IntLegP2.Last() << "\t";
    f << ScaledLeg.Last() << "\t";
    f << ScaledIntLeg.Last() << "\t";
    f << ScaledIntLegP2.Last() << "\n";
  }
  f.close();
}
#endif

} // namespace CoupledField


#endif
