#include "integrationScheme.hh"

namespace CoupledField {


IntegrationScheme::IntegrationScheme() {
 
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
    lp.number = 0;
    points[i] = lp;
    weights[i] =  a3[i][2];
  }

  // store integration points / weights back to member
  intPoints_[Elem::ET_QUAD4] = points;
  intWeights_[Elem::ET_QUAD4] = weights;
  
  
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
     lp.number = 0;
     points[i] = lp;
     weights[i] =  c1_2FE[i][3];
   }
  // store integration points / weights back to member
  intPoints_[Elem::ET_HEXA8] = points;
  intWeights_[Elem::ET_HEXA8] = weights;
  
}


IntegrationScheme::~IntegrationScheme() {
  
}


void IntegrationScheme::SetOrder( std::string method, UInt oder) {
  Warning("Implement me");
}


void IntegrationScheme::GetIntPoints( Elem::FEType elemType,
                                      StdVector<LocPoint>& intPts, 
                                      StdVector<Double>& weights ) {

  if( intPoints_.find(elemType) == intPoints_.end() ) {
    EXCEPTION( "No integration points defined for element of typed '"
        << Elem::feType.ToString( elemType ) 
        << "'!");
  }
  intPts = intPoints_[elemType]; 
  weights = intWeights_[elemType];
}




} // namespace CoupledField
