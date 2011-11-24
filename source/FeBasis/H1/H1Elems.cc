#include "H1Elems.hh"

#include <algorithm>

namespace CoupledField {


  // ========================================================================
  //  FeH1
  // ========================================================================
  
  FeH1::FeH1() {
  }
  
  FeH1::~FeH1() {
  }
  
  
  void FeH1::EvaluateLagrangePolynomial( Vector<Double> & shape, Double coord, UInt order ){
    if(supportingPoints_.find(order) == supportingPoints_.end() )
     supportingPoints_[order] = CalcGaussLobattoPoints(order);

    StdVector<Double> curSupPoints = supportingPoints_[order];

    shape.Resize(order+1);
    shape.Init();
    //get iutegration Pointes
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 122.
     for( UInt i=0; i<order+1; i++)
     {
      shape[i] = 1.0;
      for( UInt p = 0;p< order+1; p++)
      {
        if(p==i)
          continue;
        else
          shape[i] *= (coord - curSupPoints[p]) / (curSupPoints[i] - curSupPoints[p]);
      }
    }
  }

  void FeH1::EvaluateDerivLagrangePolynomial( Vector<Double> & deriv, Double coord,
                                              UInt order) {
    if(supportingPoints_.find(order) == supportingPoints_.end() )
     supportingPoints_[order] = CalcGaussLobattoPoints(order);

    StdVector<Double> curSupPoints = supportingPoints_[order];

    deriv.Resize(order+1);
    deriv.Init();
    for ( UInt i = 0; i<order+1  ; i++)
    {
      Double sum = 1.0;
      for ( UInt k=0;k<order+1 ; k++ )
      {
        if ( k != i )
        {
          sum=1.0;
          for ( UInt l=0;l<order+1 ; l++)
            if ( (l!=i) && (l != k))
              sum *= (coord - curSupPoints[l]);
          deriv[i] += sum;
        }
      }
    }
    for( UInt i=0; i< order+1; i++)
    {
      for( UInt p = 0;p< order+1; p++)
      {
        if(p==i)
          continue;
        else
          deriv[i] *= 1.0 / (curSupPoints[i] - curSupPoints[p]);
      }
    }

  }
  
  void FeH1::GetShFnc( Vector<Double>& S, const LocPoint& lp,
                 const Elem* ptElem,  UInt comp ){
    
    //check if the shape function is already computed
    if(shapeFncsAtIp_.find(lp.number) == shapeFncsAtIp_.end() || comp !=1 ){
      CalcShFnc( S, lp.coord, ptElem, comp);
      shapeFncsAtIp_[lp.number] = S;
    }else{
      S = shapeFncsAtIp_[lp.number];
    }
  }
  
  void FeH1::GetGlobDerivShFnc( Matrix<Double>& deriv, const LocPointMapped& lpm,
                                const Elem* elem, UInt comp ){
    // Get local derivative
    Matrix<Double> locDeriv;
    
    //check if the shfunction is already computed
    if(shapeFncDerivsAtIp_.find(lpm.lp.number) == shapeFncDerivsAtIp_.end() || comp !=1 ){
      CalcLocDerivShFnc( locDeriv, lpm.lp.coord, elem, comp);
      //add them to the map
      shapeFncDerivsAtIp_[lpm.lp.number] = locDeriv;
    }else{
      locDeriv = shapeFncDerivsAtIp_[lpm.lp.number];
    }
    deriv = locDeriv * lpm.jacInv;
  }
  
  void FeH1::GetLocDerivShFnc( Matrix<Double>& deriv, const LocPoint& lp,
                               const Elem* elem, UInt comp  ) {
    if(shapeFncDerivsAtIp_.find(lp.number) == shapeFncDerivsAtIp_.end()){
      CalcLocDerivShFnc( deriv, lp.coord, elem, comp);
      //add them to the map
      shapeFncDerivsAtIp_[lp.number] = deriv;
    }else{
      deriv = shapeFncDerivsAtIp_[lp.number];
    }
  }

  void FeH1::SetFunctionsAtIp(const StdVector<LocPoint>& iPoints){
    
    
    // can only be performed for non-hierarchical elements
    //shapeFncsAtIp_.resize(iPoints.GetSize());
    //shapeFncDerivsAtIp_.resize(iPoints.GetSize());
    for(UInt aPoint = 0; aPoint < iPoints.GetSize();aPoint++){
     const LocPoint& lp = iPoints[aPoint];
      CalcShFnc( shapeFncsAtIp_[lp.number], lp.coord, NULL, 1);
      CalcLocDerivShFnc( shapeFncDerivsAtIp_[lp.number], lp.coord,
                         NULL, 1);
    }
  }
  
  void FeH1::SetFunctionsAtIp(const std::map<Integer,LocPoint>& iPoints){


    // can only be performed for non-hierarchical elements
    //shapeFncsAtIp_.resize(iPoints.GetSize());
    //shapeFncDerivsAtIp_.resize(iPoints.GetSize());
    std::map<Integer,LocPoint>::const_iterator pIt = iPoints.begin();
    while(pIt != iPoints.end()){
     const LocPoint& lp = pIt->second;
      CalcShFnc( shapeFncsAtIp_[lp.number], lp.coord, NULL, 1);
      CalcLocDerivShFnc( shapeFncDerivsAtIp_[lp.number], lp.coord,
                         NULL, 1);
      pIt++;
    }
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
    StdVector<Double> xold(N1);
    xold.Init(2);
    StdVector<Double> xUpdate(N1);
    xUpdate.Init(2);
    Matrix<Double> P;
    P.Resize(N1,N1);
  
    for ( i = 0; i < N1; i += 1 ){
      x[i]=cos(M_PI * i/order);
    }

    // Compute P_N using the recursion relation
    // Compute its first and second derivatives and 
    // update x using the Newton-Raphson method.
  
    while (*std::max_element(xUpdate.Begin(), xUpdate.End() ) >  1e-14) {
  
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
    //x is now a vector ranging from 1 to -1
    //in contradiction to our convention so we revert it
    std::reverse(x.Begin(),x.End());
    //std::cout << "ELEM" << std::endl;
    //std::cout << x << std::endl << std::endl;

    return x;
  }
  } // namespace CoupledField
