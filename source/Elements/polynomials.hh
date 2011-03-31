// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
#ifndef FILE_CFS_POLYNOMIALS_HH
#define FILE_CFS_POLYNOMIALS_HH

namespace CoupledField {

  //! \file Define different types of hierarchic Polynomials
  

// ========================================================================
//  Integrated Legendre Polynomials
// ========================================================================

//! Integrated Legendre polynomials. Based on the shape functions of 
//! S. Zaglmayr.
/*!
  \param values (output) return values
  \param order (input) order of the polynomials
  \param loc (input) location of the polynomial to be evaluated [-1,1]
*/
//! \note If the length of the array is always (p-1), as the linear 1st order
//!       shape functions are not included 
//  template<typename T_SCAL, typename T_VEC>
//  void IntLegendrePoly( T_VEC& values, UInt order, T_SCAL loc ) {
//    
//    T_SCAL s1 = -1.0;
//    T_SCAL s2 = 0.0;
//    T_SCAL s3;
//
//    values.Resize(order);
//     
//    //values[0] = loc -1.0; 
//    for (int j=1; j<=order+1; j++){
//      Double inv = 1.0 / j;
//      s3=s2; 
//      s2=s1;
//      s1=( (2*(Double)j-3) * loc * s2 - ((Double)j-3) * s3) * inv;
//      values[j-1] = s1;
//    }
//  }
  
 
// this version corresponds to the class
// IntegratedLegendreMonomialExt
  template <typename T_SCAL, class T_VEC>
  void IntLegendrePoly( T_VEC& values, UInt order, T_SCAL loc ) {
     T_SCAL p3 = 0;
     T_SCAL p2 = -1;
     T_SCAL p1 = loc;
     values.Resize(order-1);
     for (int j=2; j<=order; j++) {
       double inv = 1.0 / (double)j;
       p3=p2;
       p2=p1;
       p1=( (2*(double)j-3) * loc * p2 - ((double)j-3) * p3) *inv;
       values[j-2] = p1;
     }
   }
   

} // namespace CoupledField

#endif
