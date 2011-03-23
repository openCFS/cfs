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
  template<typename SCAL, typename LOCTYPE>
  void IntLegendrePoly( SCAL& values, UInt order, LOCTYPE loc ) {
    
    SCAL s1 = -1.0;
    SCAL s2 = 0.0;
    SCAL s3;

    values.Resize(order+1);
     
    for (int j=1; j<=order+1; j++){
      s3=s2; 
      s2=s1;
      s1=( (2*j-3) * loc * s2 - (j-3) * s3) / j;
      values[j] = s1;
    }
  }






} // namespace CoupledField

#endif
