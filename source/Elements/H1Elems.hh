// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
#ifndef FILE_CFS_FEFUNCTION_HH
#define FILE_CFS_FEFUNCTION_HH

namespace CoupledField {

  //! Base class for H1-conforming reference elements (non hierarchical)
  class FeH1 : public BaseFE {
  public:
  
  protected:
  
    //! Evalute Lagrange polynomial 
    void EvaluatePolynomial( Vector<Double> & shape, Double coord );
    
    //! Evaluate derivative of Lagrange polynomial
    void EvaluateDerivPolynomial( Vector<Double> & deriv, Double coord );
     
  };
  
  
  // ========================================================================
  // Lagrange H1 elements of first order
  // ========================================================================
  
  class FeH1LoLine1 : public FeH1 {
  
  };
  
  class FeH1LoQuad1 : public FeH1 {
  
  };
  
  class FeH1LoHex1 : public FeH1 {
    
  };
  
  
  // ========================================================================
  // Lagrange H1 elements of second order
  // ========================================================================
  class FeH1LoLine2 : public FeH1 {

  };

  class FeH1LoQuad2 : public FeH1 {

  };

  class FeH1LoHex2 : public FeH1 {

  };


  // ========================================================================
  // Lagrange H1 elements of arbitrary order (evaluation using
  // general Lagrange polynomials)
  // ========================================================================
  class FeH1LoLineVar : public FeH1 {

  };

  class FeH1LoQuadVar : public FeH1 {

  };

  class FeH1LoHexVar : public FeH1 {

  };  
  
} // namespace CoupledField

#endif
