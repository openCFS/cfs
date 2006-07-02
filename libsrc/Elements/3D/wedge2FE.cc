#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "wedge2FE.hh"

namespace CoupledField
{

  Wedge2FE::Wedge2FE():WedgeFE()
  { 
    ENTER_FCN( "Wedge2FE::Wedge2FE" );
    Init();
  }


  Wedge2FE::~Wedge2FE()
  {
    ENTER_FCN( "Wedge2FE::~Wedge2FE" );
  }

  void Wedge2FE::Init()
  {
    ENTER_IFCN( "Wedge2FE::Init" );

    NumNodes_ = 15; 

    CommonInit();
  }


  void Wedge2FE::SetCornerCoords()
  {
    ENTER_IFCN( "Wedge2FE::SetCornerCoords" );

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =   0.0;
    LCornerCoords_[1][0] =   0.0;
    LCornerCoords_[2][0] =  -1.0;

    LCornerCoords_[0][1] =   1.0;
    LCornerCoords_[1][1] =   0.0;
    LCornerCoords_[2][1] =  -1.0;

    LCornerCoords_[0][2] =   0.0;
    LCornerCoords_[1][2] =   1.0;
    LCornerCoords_[2][2] =  -1.0;

    LCornerCoords_[0][3] =   0.0;
    LCornerCoords_[1][3] =   0.0;
    LCornerCoords_[2][3] =   1.0;

    LCornerCoords_[0][4] =   1.0;
    LCornerCoords_[1][4] =   0.0;
    LCornerCoords_[2][4] =   1.0;

    LCornerCoords_[0][5] =  0.0;
    LCornerCoords_[1][5] =  1.0;
    LCornerCoords_[2][5] =  1.0;

    LCornerCoords_[0][6] =  0.5;
    LCornerCoords_[1][6] =  0.0;
    LCornerCoords_[2][6] = -1.0;

    LCornerCoords_[0][7] =  0.5;
    LCornerCoords_[1][7] =  0.5;
    LCornerCoords_[2][7] = -1.0;

    LCornerCoords_[0][8] =  0.0;
    LCornerCoords_[1][8] =  0.5;
    LCornerCoords_[2][8] = -1.0;

    LCornerCoords_[0][9] =  0.5;
    LCornerCoords_[1][9] =  0.0;
    LCornerCoords_[2][9] =  1.0;

    LCornerCoords_[0][10] =  0.5;
    LCornerCoords_[1][10] =  0.5;
    LCornerCoords_[2][10] =  1.0;

    LCornerCoords_[0][11] =  0.0;
    LCornerCoords_[1][11] =  0.5;
    LCornerCoords_[2][11] =  1.0;

    LCornerCoords_[0][12] =  0.0;
    LCornerCoords_[1][12] =  0.0;
    LCornerCoords_[2][12] =  0.0;

    LCornerCoords_[0][13] =  0.5;
    LCornerCoords_[1][13] =  0.0;
    LCornerCoords_[2][13] =  0.0;

    LCornerCoords_[0][14] =  0.0;
    LCornerCoords_[1][14] =  0.5;
    LCornerCoords_[2][14] =  0.0;


  }


  void Wedge2FE :: CalcShapeFnc(Vector<Double> & Shape, 
                                const Vector<Double> & LCoord)
  {
    ENTER_IFCN( "Wedge2FE::CalcShapeFnc" );

    Shape.Resize(NumNodes_);

    Shape[0] =  0.5*LCoord[2] * (1 - LCoord[2]) * (1 - LCoord[0] - LCoord[1]) 
      * (2 * LCoord[0] + 2*LCoord[1] -1);
    Shape[1] =  0.5*LCoord[2] * (1 - LCoord[2]) *  LCoord[0] * (1 - 2*LCoord[0]);
    Shape[2] =  0.5*LCoord[2] * (1 - LCoord[2]) *  LCoord[1] * (1 - 2*LCoord[1]);

    Shape[3] =- 0.5*LCoord[2] * (1 + LCoord[2]) * (1 - LCoord[0] - LCoord[1]) 
      * (2*LCoord[0] + 2*LCoord[1] -1);
    Shape[4] = -0.5*LCoord[2] * (1 + LCoord[2]) * LCoord[0] * (1 - 2*LCoord[0]);
    Shape[5] = -0.5*LCoord[2] * (1 + LCoord[2]) * LCoord[1] * (1 - 2*LCoord[1]);

    Shape[6] = -2*LCoord[2] * (1 - LCoord[2]) * LCoord[0] * (1 - LCoord[0] - LCoord[1]);
    Shape[7] = -2*LCoord[2] * (1 - LCoord[2]) * LCoord[0] * LCoord[1];
    Shape[8] = -2*LCoord[2] * (1 - LCoord[2]) * LCoord[1] * (1 - LCoord[0] - LCoord[1]);

    Shape[9]  = 2*LCoord[2] * (1 + LCoord[2]) * LCoord[0] * (1 - LCoord[0] - LCoord[1]);
    Shape[10] = 2*LCoord[2] * (1 + LCoord[2]) * LCoord[0] * LCoord[1];
    Shape[11] = 2*LCoord[2] * (1 + LCoord[2]) * LCoord[1] * (1 - LCoord[0] - LCoord[1]);

    Shape[12] = (1 - LCoord[0] - LCoord[1]) * (1 - LCoord[2] * LCoord[2]);
    Shape[13] = LCoord[0] * (1 - LCoord[2] * LCoord[2]);
    Shape[14] = LCoord[1] * (1 - LCoord[2] * LCoord[2]);
  }


  void Wedge2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                          const Vector<Double> & LCoord)
  {
    ENTER_IFCN( "Wedge2FE::CalcLocalDerivShapeFnc" );

    LDeriv.Resize(NumNodes_,Dim_);

    LDeriv.Init();

    LDeriv[0][0] = +0.5*LCoord[2] * (1 - LCoord[2]) * (3 - 4*LCoord[0] - 4*LCoord[1]);
    LDeriv[0][1] = +0.5*LCoord[2] * (1 - LCoord[2]) * (3 - 4*LCoord[0] - 4*LCoord[1]);
    LDeriv[0][2] = +0.5 * (1 - 2*LCoord[2]) * (1 - LCoord[0] - LCoord[1]) 
      * (2*LCoord[0] + 2*LCoord[1] -1 );

    LDeriv[1][0] =  0.5*LCoord[2] * (1 - LCoord[2]) * (1 - 4*LCoord[0]);
    LDeriv[1][1] =  0.0;
    LDeriv[1][2] =  0.5 * (1 - 2*LCoord[2]) * LCoord[0] * (1 - 2*LCoord[0]);

    LDeriv[2][0] =  0.0;
    LDeriv[2][1] =  0.5*LCoord[2] * (1 - LCoord[2]) * (1 - 4*LCoord[1]);
    LDeriv[2][2] =  0.5 * (1 - 2*LCoord[2]) * LCoord[1] * (1 - 2*LCoord[1]);

    LDeriv[3][0] = -0.5*LCoord[2] * (1 + LCoord[2]) * (3 - 4*LCoord[0] - 4*LCoord[1]);
    LDeriv[3][1] = -0.5*LCoord[2] * (1 + LCoord[2]) * (3 - 4*LCoord[0] - 4*LCoord[1]);
    LDeriv[3][2] = -0.5 * (1 + 2*LCoord[2]) * (1 - LCoord[0] - LCoord[1]) 
      * (2*LCoord[0] + 2*LCoord[1] -1 );

    LDeriv[4][0] = -0.5*LCoord[2] * (1 + LCoord[2]) * (1 - 4*LCoord[0]);
    LDeriv[4][1] =  0.0;
    LDeriv[4][2] = -0.5 * (1 + 2*LCoord[2]) * LCoord[0] * (1 - 2*LCoord[0]);

    LDeriv[5][0] =  0.0;
    LDeriv[5][1] = -0.5*LCoord[2] * (1 + LCoord[2]) * (1 - 4*LCoord[1]);
    LDeriv[5][2] = -0.5 * (1 + 2*LCoord[2]) * LCoord[1] * (1 - 2*LCoord[1]);

    LDeriv[6][0] =  -2*LCoord[2] * (1 - LCoord[2]) * (1 - 2*LCoord[0] - LCoord[1]);
    LDeriv[6][1] =   2*LCoord[2] * (1 - LCoord[2]) * LCoord[0];
    LDeriv[6][2] =  -2*(1 - 2*LCoord[2]) * LCoord[0] * (1 - LCoord[0] - LCoord[1]); 

    LDeriv[7][0] =  -2*LCoord[2] * (1 - LCoord[2]) * LCoord[1];
    LDeriv[7][1] =  -2*LCoord[2] * (1 - LCoord[2]) * LCoord[0];
    LDeriv[7][2] =  -2*(1 - 2*LCoord[2]) * LCoord[0] * LCoord[1]; 

    LDeriv[8][0] =   2*LCoord[2] * (1 - LCoord[2]) * LCoord[1];
    LDeriv[8][1] =  -2*LCoord[2] * (1 - LCoord[2]) * (1 - LCoord[0] - 2*LCoord[1]);
    LDeriv[8][2] =  -2*(1 - 2*LCoord[2]) * LCoord[1] * (1 - LCoord[0] - LCoord[1]); 

    LDeriv[9][0] =   2*LCoord[2] * (1 + LCoord[2]) * (1 - 2*LCoord[0] - LCoord[1]);
    LDeriv[9][1] =  -2*LCoord[2] * (1 + LCoord[2]) * LCoord[0];
    LDeriv[9][2] =   2*(1 + 2*LCoord[2]) * LCoord[0] * (1 - LCoord[0] - LCoord[1]); 

    LDeriv[10][0] =  2*LCoord[2] * (1 + LCoord[2]) * LCoord[1];
    LDeriv[10][1] =  2*LCoord[2] * (1 + LCoord[2]) * LCoord[0];
    LDeriv[10][2] =  2*(1 + 2*LCoord[2]) * LCoord[0] * LCoord[1]; 

    LDeriv[11][0] =  -2*LCoord[2] * (1 + LCoord[2]) * LCoord[1];
    LDeriv[11][1] =   2*LCoord[2] * (1 + LCoord[2]) * (1 - LCoord[0] - 2*LCoord[1]);
    LDeriv[11][2] =   2*(1 + 2*LCoord[2]) * LCoord[1] * (1 - LCoord[0] - LCoord[1]); 

    LDeriv[12][0] =  -(1 - LCoord[2] * LCoord[2]);
    LDeriv[12][1] =  -(1 - LCoord[2] * LCoord[2]);
    LDeriv[12][2] =  -2*LCoord[2] * (1 - LCoord[0] - LCoord[1]); 

    LDeriv[13][0] =  (1 - LCoord[2] * LCoord[2]);
    LDeriv[13][1] =  0;
    LDeriv[13][2] =  -2*LCoord[2] * LCoord[0];

    LDeriv[14][0] =  0;
    LDeriv[14][1] =  (1 - LCoord[2] * LCoord[2]);
    LDeriv[14][2] = -2*LCoord[2] * LCoord[1];

  }
  



} // end of namespace

