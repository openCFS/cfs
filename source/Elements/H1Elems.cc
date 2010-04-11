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
  
  // ========================================================================
  //  FeH1LagrangeExpl 
  // ========================================================================
  
  FeH1LagrangeExpl::FeH1LagrangeExpl() {
    
  }
    
  FeH1LagrangeExpl::~FeH1LagrangeExpl() {
  }
  
  void FeH1LagrangeExpl::SetIntPoints( StdVector<LocPoint>& intPoints ) {
  }
  
  void FeH1LagrangeExpl::GetNumFncs( StdVector<UInt>& numFcns,
                                     EntityType fctEntityType,
                                     UInt dof ) {
    // Initialize explictily with number of nodes
    if( fctEntityType == NODE ) {
      numFcns.Resize( shape_.numNodes );
      numFcns.Init( 1 );
    }
  }

  void FeH1LagrangeExpl::GetFncPermutation( StdVector<UInt>& fncPermutation,
                                            const Elem* ptElem,
                                            EntityType fctEntityType,
                                            UInt entNumber){
    if( fctEntityType == NODE ) {
      fncPermutation.Resize(shape_.numNodes);
      for ( UInt i = 0; i < shape_.numNodes; i++ ){
        fncPermutation[i] = i;
      }
    }else{
      fncPermutation.Resize(0);
    }
    return;
  }
  UInt FeH1LagrangeExpl::GetNumFncsPerEntType( EntityType fctEntityType,
                                     UInt dof){
    UInt numFnc = 0;
    // Initialize explictily with number of nodes
    if( fctEntityType == NODE ) {
      numFnc = shape_.numNodes;
    }
    return numFnc;
  }

//  void FeH1LagrangeExpl::GetShFnc( Vector<Double> & shape, const LocPoint& lp,
//                                   const Elem* ptElem,  UInt comp ) {
//  
//    // In the future, we re-use the shape functions calculated at
//      // the integration points
//    //  if( lp.number != 0 ) {
//    //    //shape = shapeAtIp_[lp.number];
//    //  } else { 
//        CalcShFnc( shape, lp.coord);
//    //  }
//  }
  
//  void FeH1LagrangeExpl::GetDerivShFnc( Matrix<Double> & deriv, 
//                                        const LocPoint& lp,
//                                        const Elem * elem, 
//                                        UInt comp  ) {
//    
//    // In the future, we re-use the shape functions calculated at
//    // the integration points
//  //  if( lp.number != 0 ) {
//  //    //deriv = shapeDerivAtIp_[lp.number];
//  //  } else { 
//      CalcDerivShFnc( deriv, lp.coord);
//  //  }
//  }
    
  
  // ========================================================================
  //  Lagrangian Elements of 1st order
  // ========================================================================
  
  // --- Line 1st order ---
  
  FeH1LagrangeLine1::FeH1LagrangeLine1() {
    feType_ = Elem::ET_LINE2;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 2;
    order_ = 1;
  }
  FeH1LagrangeLine1::~FeH1LagrangeLine1() {
    
  }
  
  void FeH1LagrangeLine1::CalcShFnc( Vector<Double>& shape,
                                             const Vector<Double>& point ) {
     shape.Resize( 2 );
     shape[0] = 0.5 * ( 1.0 - point[0] );
     shape[1] = 0.5 * ( 1.0 + point[0] );
  }
  
  void FeH1LagrangeLine1::CalcDerivShFnc( Matrix<Double> & deriv, 
                                               const Vector<Double>& point ) {
   
      deriv.Resize(2, 1);
      deriv[0][0] = 0.5 * -1.0;
      deriv[1][0] = 0.5 *  1.0;
  }
  
  // --- Quad 1st order ---
   
  FeH1LagrangeQuad1::FeH1LagrangeQuad1() {
    feType_ = Elem::ET_QUAD4;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 4;
    order_ = 1; 
  }
    
  FeH1LagrangeQuad1::~FeH1LagrangeQuad1() {
    
  }
  
  void FeH1LagrangeQuad1::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point ) {
     shape.Resize( 4 );
     shape[0] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ); 
     shape[1] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] );
     shape[2] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] );
     shape[3] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] );
    
  }
  
  void FeH1LagrangeQuad1::CalcDerivShFnc( Matrix<Double> & deriv, 
                                               const Vector<Double>& point ) {
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    deriv.Resize( 4, 2 );
    for( UInt i = 0; i < 4; i++ ) {
      deriv[i][0] = 0.25 * coords[i][0] * (1 + coords[i][1] * point[1] );
      deriv[i][1] = 0.25 * (1 + coords[i][0] * point[0] ) * coords[i][1];
    }
  }
  
  
  // --- Hex 1st order ---
  FeH1LagrangeHex1::FeH1LagrangeHex1() {
    feType_ = Elem::ET_HEXA8;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 8;
  }
    
  FeH1LagrangeHex1::~FeH1LagrangeHex1() {
    
  }
  
  void FeH1LagrangeHex1::CalcShFnc( Vector<Double>& shape,
                                    const Vector<Double>& point ) {
    shape.Resize( 8 );
    shape[0] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]); 
    shape[1] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]);
    shape[2] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
    shape[3] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
    shape[4] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]); 
    shape[5] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]);
    shape[6] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
    shape[7] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
    
  }
  
  void FeH1LagrangeHex1::CalcDerivShFnc( Matrix<Double> & deriv, 
                                         const Vector<Double>& point ) {
    deriv.Resize( 8, 3 );
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    for( UInt i = 0; i < 8; i++ ) {
      deriv[i][0] = 0.125  * coords[i][0] 
                           * (1 + coords[i][1] * point[1] ) 
                           * (1 + coords[i][2] * point[2] );
      
      deriv[i][1] = 0.125  * (1 + coords[i][0] * point[0] )
                           * coords[i][1] 
                           * (1 + coords[i][2] * point[2] );
      
      deriv[i][2] = 0.125  * (1 + coords[i][0] * point[0] ) 
                           * (1 + coords[i][1] * point[1] )
                           * coords[i][2];
    }
    
  }
  
  // ========================================================================
  //  Lagrangian Elements of 2nd order
  // ========================================================================
  
  // --- Line 2nd order ---
  
  FeH1LagrangeLine2::FeH1LagrangeLine2() {
    feType_ = Elem::ET_LINE3;;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 3;
  }
  FeH1LagrangeLine2::~FeH1LagrangeLine2() {
    
  }
  
  void FeH1LagrangeLine2::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point ) {
    EXCEPTION("Implement me"); 
  //  shape.Resize( 2 );
  //   shape[0] = 0.5 * ( 1.0 - point[0] );
  //   shape[1] = 0.5 * ( 1.0 + point[0] );
  }
  
  void FeH1LagrangeLine2::CalcDerivShFnc( Matrix<Double> & deriv, 
                                          const Vector<Double>& point ) {
    EXCEPTION("Implement me");
  //    deriv.Resize(2, 1);
  //    deriv[0][0] = 0.5 * -1.0;
  //    deriv[1][0] = 0.5 *  1.0;
  }
  
  // --- Quad 2nd order ---
   
  FeH1LagrangeQuad2::FeH1LagrangeQuad2() {
    feType_ = Elem::ET_QUAD8;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 8;  
  }
    
  FeH1LagrangeQuad2::~FeH1LagrangeQuad2() {
    
  }
  
  void FeH1LagrangeQuad2::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point ) {
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    shape.Resize( 8 );
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 122.
    
    // corner nodesf
    for( UInt i = 0; i < 4; i++ ) {
      shape[i] = 0.25 * ( 1 + coords[i][0] * point[0] ) 
                      * ( 1 + coords[i][1] * point[1] )
                      * ( coords[i][0] * point[0] + coords[i][1] * point[1] - 1 );
    }
    
    // mid-side nodes
    for( UInt i = 4; i < 8; i = i + 2 ) {
      shape[i]   = 0.5 * ( 1 - point[0] * point[0] )
                       * ( 1 + coords[i][1] * point[1] );
      shape[i+1] = 0.5 * (1 - point[1] * point[1] )
                       * (1 + coords[i+1][0] * point[0] );
    }
  }
  
  void FeH1LagrangeQuad2::CalcDerivShFnc( Matrix<Double> & deriv, 
                                          const Vector<Double>& point ) {
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    deriv.Resize( 8, 2 );
    
    // corner nodes
    for( UInt i = 0; i < 4; i++ ) {
      deriv[i][0] = 0.25 * coords[i][0]  
                         * ( 1 + coords[i][1] * point[1] )
                         * ( 2 * coords[i][0] * point[0] 
                               + coords[i][1] * point[1] );
      deriv[i][1] = 0.25 * coords[i][1]
                         * ( 1 + coords[i][0] * point[0] )
                         * ( 2 * coords[i][1] * point[1] 
                               + coords[i][0] * point[0] );
    }
      
    // mid-side nodes
    for( UInt i = 4; i < 8; i = i + 2 ) {
      deriv[i][0] = - point[0] * ( 1 + coords[i][1] * point[1] );
      deriv[i][1] =  0.5 * coords[i][1] * ( 1 - point[0] * point[0] );
      
      deriv[i+1][0] = 0.5 * coords[i+1][0] * ( 1 - point[1] * point[1] );
      deriv[i+1][1] = - point[1] * ( 1 + coords[i+1][0] * point[0] );
    }
  }
  
  
  // --- Hex 2nd order ---
  FeH1LagrangeHex2::FeH1LagrangeHex2() {
    feType_ = Elem::ET_HEXA20;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 20;
  }
    
  FeH1LagrangeHex2::~FeH1LagrangeHex2() {
    
  }
  
  void FeH1LagrangeHex2::CalcShFnc( Vector<Double>& shape,
                                    const Vector<Double>& point ) {
  EXCEPTION("Implement me");
    //  shape.Resize( 8 );
  //  shape[0] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]); 
  //  shape[1] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]);
  //  shape[2] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
  //  shape[3] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
  //  shape[4] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]); 
  //  shape[5] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]);
  //  shape[6] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
  //  shape[7] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
    
  }
  
  void FeH1LagrangeHex2::CalcDerivShFnc( Matrix<Double> & deriv, 
                                         const Vector<Double>& point ) {
    EXCEPTION("Implement me");
    //  deriv.Resize( 8, 3 );
  //  StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
  //  for( UInt i = 0; i < 8; i++ ) {
  //    deriv[i][0] = 0.125  * coords[i][0] 
  //                         * (1 + coords[i][1] * point[1] ) 
  //                         * (1 + coords[i][2] * point[2] );
  //    
  //    deriv[i][1] = 0.125  * (1 + coords[i][0] * point[0] )
  //                         * coords[i][1] 
  //                         * (1 + coords[i][2] * point[2] );
  //    
  //    deriv[i][2] = 0.125  * (1 + coords[i][0] * point[0] ) 
  //                         * (1 + coords[i][1] * point[1] )
  //                         * coords[i][2];
  //  }
    
  }

} // namespace CoupledField
