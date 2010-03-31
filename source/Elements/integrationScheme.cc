#include "integrationScheme.hh"

namespace CoupledField {


IntegrationScheme::IntegrationScheme() {
 
  integMethod_ = ECONOMICAL;
  order_ = 2;
  // Very hard coded section at the moment
  // =====================================
 
  // ---- ET_QUAD4 ----
  
  // Gauss  quadrature  points  and  weights  on  the  reference  quadrilateral  order  p=2  3
  static Double a3[][3] = { 
                           {0.577350269189626,  0.577350269189626,  1.000000000000000},
                           {0.577350269189626,  -0.577350269189626,  1.000000000000000},
                           {-0.577350269189626,  0.577350269189626,  1.000000000000000},
                           {-0.577350269189626,  -0.577350269189626,  1.000000000000000}
  };

  StdVector<LocPoint> points (4);
  StdVector<Double> weights(4);

  for( UInt i = 0; i < 4; i++ ) {
    LocPoint lp;
    lp.coord.Resize(2);
    lp.coord[0] = a3[i][0];
    lp.coord[1] = a3[i][1];
    lp.number = numIntPts_[Elem::ET_QUAD4]++;
    points[i] = lp;
    weights[i] =  a3[i][2];
  }

  // store integration points / weights back to member
  intPoints_[ECONOMICAL][2][Elem::ET_QUAD4] = points;
  intWeights_[ECONOMICAL][2][Elem::ET_QUAD4] = weights;
  
  
  // ---- ET_HEXA8 ---
  static Double c1_2FE[][4] = { 
    { -0.57735026919,  -0.57735026919,  -0.57735026919,  1.000000000000000 },
    {  0.57735026919,  -0.57735026919,  -0.57735026919,  1.000000000000000 },
    {  0.57735026919,   0.57735026919,  -0.57735026919,  1.000000000000000 },                
    { -0.57735026919,   0.57735026919,  -0.57735026919,  1.000000000000000 },        
    { -0.57735026919,  -0.57735026919,   0.57735026919,  1.000000000000000 },        
    {  0.57735026919,  -0.57735026919,   0.57735026919,  1.000000000000000 },
    {  0.57735026919,   0.57735026919,   0.57735026919,  1.000000000000000 },
    { -0.57735026919,   0.57735026919,   0.57735026919,  1.000000000000000 },
  };
  points.Resize(8);
  weights.Resize(8);
  for( UInt i = 0; i < 8; i++ ) {
     LocPoint lp;
     lp.coord.Resize(3);
     lp.coord[0] = c1_2FE[i][0];
     lp.coord[1] = c1_2FE[i][1];
     lp.coord[2] = c1_2FE[i][2];
     lp.number = numIntPts_[Elem::ET_HEXA8]++;
     points[i] = lp;
     weights[i] =  c1_2FE[i][3];
   }
  // store integration points / weights back to member
  intPoints_[ECONOMICAL][2][Elem::ET_HEXA8] = points;
  intWeights_[ECONOMICAL][2][Elem::ET_HEXA8] = weights;
  FillGaussLobattoIntegPoints(11);
}


//================================================================
//Fill Gauss Lobatto Points and weights
//================================================================
void IntegrationScheme::FillGaussLobattoIntegPoints(UInt order){
  StdVector<LocPoint> points;
  StdVector<Double> weights;
  StdVector<Double> intPoints1D;
  StdVector<Double> weights1D;

  //loop over every order here 10
  for(UInt ord = 1; ord < order ; ord++){
    CalcGaussLobattoPointsWeights(ord,intPoints1D,weights1D);

    points.Resize(ord+1);
    weights.Resize(ord+1);
    //fill the lines
    for( UInt i = 0; i < ord+1; i++ ) {
      LocPoint lp;
      lp.coord.Resize(1);
      lp.coord[0] = intPoints1D[i];
      lp.number = numIntPts_[Elem::ET_LINE2]++;
      points[i] = lp;
      weights[i] =  weights1D[i];
    }
    intPoints_[LOBATTO][ord][Elem::ET_LINE2] = points;
    intWeights_[LOBATTO][ord][Elem::ET_LINE2] = weights;

    //fill rects
    UInt numPts = (ord+1)*(ord+1);
    points.Resize(numPts);
    weights.Resize(numPts);

    for (UInt i = 0; i < ord+1; i ++ ){
      for (UInt j = 0; j < ord+1; j ++ ){
        LocPoint lp;
        lp.coord.Resize(2);
        lp.coord[0] = intPoints1D[i];
        lp.coord[1] = intPoints1D[j];
        lp.number = numIntPts_[Elem::ET_QUAD4]++;
        points[(i*(ord+1) + j)] = lp;
        weights[(i*(ord+1) + j)] =  weights1D[i]*weights1D[j];
      }
    }
    intPoints_[LOBATTO][ord][Elem::ET_QUAD4] = points;
    intWeights_[LOBATTO][ord][Elem::ET_QUAD4] = weights;

    //fill hexas 
    numPts = (ord+1)*(ord+1)*(ord+1);
    points.Resize(numPts);
    weights.Resize(numPts);

    for (UInt i = 0; i < ord+1; i ++ ){
      for (UInt j = 0; j < ord+1; j ++ ){
        for (UInt k = 0; k < ord+1; k ++ ){
          LocPoint lp;
          lp.coord.Resize(3);
          lp.coord[0] = intPoints1D[i];
          lp.coord[1] = intPoints1D[j];
          lp.coord[2] = intPoints1D[k];
          lp.number = numIntPts_[Elem::ET_HEXA8]++;
          points[(i*((ord+1)*(ord+1)) + j*(ord+1) + k)] = lp;
          weights[(i*((ord+1)*(ord+1)) + j*(ord+1) + k)] =  weights1D[i]*weights1D[j]*weights1D[k];
        }
      }
    }
    intPoints_[LOBATTO][ord][Elem::ET_HEXA8] = points;
    intWeights_[LOBATTO][ord][Elem::ET_HEXA8] = weights;
  }

  
}


IntegrationScheme::~IntegrationScheme() {
  
}


void IntegrationScheme::SetOrder( std::string method, UInt order) {
  String2Enum(method,integMethod_);
  if(order == -1){
    order_ = 2;
  }else{
    order_ = order;
  }
}

void IntegrationScheme::SetOrder(IntegrationMethod method, UInt order) {
  integMethod_ = method;
  if(order == -1){
    order_ = 2;
  }else{
    order_ = order;
  }
}


void IntegrationScheme::GetIntPoints( Elem::FEType elemType,
                                      StdVector<LocPoint>& intPts, 
                                      StdVector<Double>& weights ) {

  if( intPoints_.find(integMethod_) == intPoints_.end() ) {
    std::string method;
    Enum2String(integMethod_,method);
    EXCEPTION( "No integration points defined for  '"
        << method 
        << "' integration!");
  }
  if( intPoints_[integMethod_].find(order_) == intPoints_[integMethod_].end() ) {
    EXCEPTION( "No integration points defined for  order '"
        << order_ 
        << "' !");
  }
  if( intPoints_[integMethod_][order_].find(elemType) == intPoints_[integMethod_][order_].end() ) {
    EXCEPTION( "No integration points defined for element of typed '"
        << Elem::feType.ToString( elemType ) 
        << "'!");
  }
  intPts = intPoints_[integMethod_][order_][elemType]; 
  weights = intWeights_[integMethod_][order_][elemType];
}


  void IntegrationScheme::CalcGaussLobattoPointsWeights(UInt order,StdVector<Double>& points, StdVector<Double>& weights){
    // Computes the Legendre-Gauss-Lobatto nodes by the LGL Vandermonde 
    // matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). 
    points.Resize(order+1);
    weights.Resize(order+1);

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
    points = x;
    Double divisor;
    for(i = 0 ; i<N1; i++){
      divisor = ((N1-1)*N1*P[i][N1-1]) * ((N1-1)*N1*P[i][N1-1]);
      weights[i] = 2 / divisor;
    }
    //now sort the arrays
    StdVector<Double> tmpPoints = points;
    StdVector<Double> tmpWeights = weights;
    //hier uu noch einen algorithmus einbauen um nachkommstellen abzuschneiden
    points[0] = tmpPoints[0];
    points[1] = tmpPoints[N1-1];
    weights[1] = tmpWeights[N1-1];
    for(i = 2; i < N1 ; i++){
      points[i] = tmpPoints[i-1];
      weights[i] = tmpWeights[i-1];
    }
  }


  void IntegrationScheme::GetAllIntegrationPoints(StdVector< LocPoint >& points,Elem::FEType type){ 

    points.Resize(numIntPts_[type]);
    points.Init();
    UInt counter = 0;
    //loop over all defined methods
    std::map<IntegrationMethod, std::map< UInt, IntegrationPoints > >::iterator methods;
    for(methods = intPoints_.begin();methods != intPoints_.end() ; methods++){
      //loop over every order available
      std::map< UInt, IntegrationPoints >::iterator orders;
      for(orders = methods->second.begin();orders != methods->second.end(); orders++){
        //get the vector for the given FeType and fill the vector
        if(orders->second[type].GetSize() > 0){
          StdVector<LocPoint> curPoints = orders->second[type];
          for(UInt curPoint = 0; curPoint < curPoints.GetSize();curPoint++){
            points[counter++] = orders->second[type][curPoint];
          }
        }
      }
    }
  }

} // namespace CoupledField
