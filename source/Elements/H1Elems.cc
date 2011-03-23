// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "H1Elems.hh"

namespace CoupledField {


  // ========================================================================
  //  FeH1
  // ========================================================================
  
  FeH1::FeH1() {
  }
  
  FeH1::~FeH1() {
  }
  
  
  void FeH1::EvaluateLagrangePolynomial( Vector<Double> & shape, Double coord ){
    if(supportingPoints_.find(order_) == supportingPoints_.end() )
     supportingPoints_[order_] = CalcGaussLobattoPoints(order_);

    StdVector<Double> curSupPoints = supportingPoints_[order_];

    shape.Resize(order_+1);
    shape.Init();
    //get iutegration Pointes
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 122.
     for( UInt i=0; i<order_+1; i++)
     {
      shape[i] = 1.0;
      for( UInt p = 0;p< order_+1; p++)
      {
        if(p==i)
          continue;
        else
          shape[i] *= (coord - curSupPoints[p]) / (curSupPoints[i] - curSupPoints[p]);
      }
    }
  }

  void FeH1::EvaluateDerivLagrangePolynomial( Vector<Double> & deriv, Double coord ) {
    if(supportingPoints_.find(order_) == supportingPoints_.end() )
     supportingPoints_[order_] = CalcGaussLobattoPoints(order_);

    StdVector<Double> curSupPoints = supportingPoints_[order_];

    deriv.Resize(order_+1);
    deriv.Init();
    for ( UInt i = 0; i<order_+1  ; i++)
    {
      Double sum = 1.0;
      for ( UInt k=0;k<order_+1 ; k++ )
      {
        if ( k != i )
        {
          sum=1.0;
          for ( UInt l=0;l<order_+1 ; l++)
            if ( (l!=i) && (l != k))
              sum *= (coord - curSupPoints[l]);
          deriv[i] += sum;
        }
      }
    }
    for( UInt i=0; i< order_+1; i++)
    {
      for( UInt p = 0;p< order_+1; p++)
      {
        if(p==i)
          continue;
        else
          deriv[i] *= 1.0 / (curSupPoints[i] - curSupPoints[p]);
      }
    }

  }
  
  void FeH1::GetGlobDerivShFnc( Matrix<Double>& deriv, LocPointMapped& lp,
                                const Elem* elem, UInt comp ){
    
    // Get local derivative
    Matrix<Double> locDeriv;
    GetDerivShFnc( locDeriv, lp.lp, elem, comp );
    deriv = locDeriv * lp.jacInv;
    
    
  }
  
  void FeH1::CalcAllSupportingPoints(UInt maxOrder){
    if(maxOrder == 0){
      return;
    }else{
      for ( UInt i = 1; i <= maxOrder; i += 1 ){
        supportingPoints_[i] = CalcGaussLobattoPoints(i);
      }
    }
  }
  
  void FeH1::CalcSupportingPoints(UInt order){
     supportingPoints_[order] = CalcGaussLobattoPoints(order);
  }
  
  StdVector<Double> FeH1::CalcGaussLobattoPoints(UInt order){
    // Computes the Legendre-Gauss-Lobatto nodes by the LGL Vandermonde 
    // matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
    // integration and spectral methods. 
    Integer N1=order + 1;
    Integer i,j;
  
    StdVector<Double> x(N1);
    StdVector<Double> xold(N1,2);
    StdVector<Double> xUpdate(N1,2);
    Matrix<Double> P;
    P.Resize(N1,N1);
  
    for ( i = 0; i < N1; i += 1 ){
      x[i]=cos(M_PI * i/order);
    }
    // Compute P_N using the recursion relation
    // Compute its first and second derivatives and 
    // update x using the Newton-Raphson method.
  
    while (*std::max_element(xUpdate.Begin(), xUpdate.End() ) >  1e-18) {
  
      xold = x;
      
      for ( i = 0; i < N1; i += 1 ){
        P[i][0] = 1;
        P[i][1] = x[i];
      }
  
      for ( i = 0; i < N1; i += 1 ){
        for ( j = 1; j < (N1-1); j += 1 ){
          P[i][j+1] = ( (2*(j+1)-1)*x[i] * P[i][j] - j*P[i][j-1] )/(j+1);
        }
      }
  
      for ( i = 0; i < N1; i += 1 ){
        x[i] = xold[i] - (x[i]*P[i][N1-1] - P[i][N1-2])/(N1 * P[i][N1-1]);
        xUpdate[i] = abs(x[i] - xold[i]);
      }
  
    }
    return x;
  }
  } // namespace CoupledField
