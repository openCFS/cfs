#include "FeNodal.hh"

namespace CoupledField {

FeNodal::FeNodal() {
  
}

FeNodal::FeNodal(const FeNodal & other){
  this->supportingPoints_ = other.supportingPoints_;
}

FeNodal::~FeNodal() {
  
}

void FeNodal::EvaluateLagrangePolynomial( Vector<Double> & shape, Double coord, UInt order ){
    if(supportingPoints_.find(order) == supportingPoints_.end() )
     supportingPoints_[order] = CalcGaussLobattoPoints(order);

    const StdVector<Double>& curSupPoints = supportingPoints_[order];

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

  void FeNodal::EvaluateDerivLagrangePolynomial( Vector<Double> & deriv, Double coord,
                                              UInt order) {
    if(supportingPoints_.find(order) == supportingPoints_.end() )
     supportingPoints_[order] = CalcGaussLobattoPoints(order);

    const StdVector<Double>& curSupPoints = supportingPoints_[order];

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
  
  void FeNodal::CalcAllSupportingPoints(UInt maxOrder){
      if(maxOrder == 0){
        return;
      }else{
        for ( UInt i = 1; i <= maxOrder; i += 1 ){
          supportingPoints_[i] = CalcGaussLobattoPoints(i);
        }
      }
    }
    
    void FeNodal::CalcSupportingPoints(UInt order){
       supportingPoints_[order] = CalcGaussLobattoPoints(order);
    }
    
    StdVector<Double> FeNodal::CalcGaussLobattoPoints(UInt order){
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
        x[i]=cos(4*atan(1.0) * i/order);
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
    
} // end of namespace
