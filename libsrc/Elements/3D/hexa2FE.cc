#include <iostream>
#include <fstream>

//#include <General/environment.hh>
#include "hexa2FE.hh"

namespace CoupledField
{

Hexa2FE::Hexa2FE():HexaFE()
{ 
  ENTER_FCN( "Hexa2FE::Hexa2FE" );

  Init();
}

Hexa2FE::~Hexa2FE()
{
  ENTER_FCN( "Hexa2FE::~Hexa2FE" );

}

void Hexa2FE::Init()
{
  ENTER_IFCN( "Hexa2FE::Init" );
  
  Dim_ = 3;
  NumNodes_ = 20;
  NumEdges_ = 12;
  //NumFaces_ = 6;

  // first set integration points and corner coords ...
  SetIntPoints();
  SetCornerCoords();


  // ... then calc shape function values at integration points
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();


}

void Hexa2FE::SetCornerCoords()
{
  ENTER_IFCN( "Hexa2FE::SetCornerCoords" );

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =  -1;
  LCornerCoords_[1][0] =  -1;
  LCornerCoords_[2][0] =  -1;

  LCornerCoords_[0][1] =   1;
  LCornerCoords_[1][1] =  -1;
  LCornerCoords_[2][1] =  -1;

  LCornerCoords_[0][2] =   1;
  LCornerCoords_[1][2] =   1;
  LCornerCoords_[2][2] =  -1;

  LCornerCoords_[0][3] =  -1;
  LCornerCoords_[1][3] =   1;
  LCornerCoords_[2][3] =  -1;

  LCornerCoords_[0][4] =  -1;
  LCornerCoords_[1][4] =  -1;
  LCornerCoords_[2][4] =   1;

  LCornerCoords_[0][5] =   1;
  LCornerCoords_[1][5] =  -1;
  LCornerCoords_[2][5] =   1;

  LCornerCoords_[0][6] =   1;
  LCornerCoords_[1][6] =   1;
  LCornerCoords_[2][6] =   1;

  LCornerCoords_[0][7] =  -1;
  LCornerCoords_[1][7] =   1;
  LCornerCoords_[2][7] =   1;

  LCornerCoords_[0][8] =   0;
  LCornerCoords_[1][8] =  -1;
  LCornerCoords_[2][8] =  -1;
 
  LCornerCoords_[0][9] =   1;
  LCornerCoords_[1][9] =   0;
  LCornerCoords_[2][9] =  -1;

 
  LCornerCoords_[0][10] =   0;
  LCornerCoords_[1][10] =   1;
  LCornerCoords_[2][10] =  -1;
 
  LCornerCoords_[0][11] =  -1;
  LCornerCoords_[1][11] =   0;
  LCornerCoords_[2][11] =  -1;
 
  LCornerCoords_[0][12] =   0;
  LCornerCoords_[1][12] =  -1;
  LCornerCoords_[2][12] =   1;
 
  LCornerCoords_[0][13] =   1;
  LCornerCoords_[1][13] =   0;
  LCornerCoords_[2][13] =   1;
 
  LCornerCoords_[0][14] =   0;
  LCornerCoords_[1][14] =   1;
  LCornerCoords_[2][14] =   1;
 
  LCornerCoords_[0][15] =  -1;
  LCornerCoords_[1][15] =   0;
  LCornerCoords_[2][15] =   1;
 
  LCornerCoords_[0][16] =  -1;
  LCornerCoords_[1][16] =  -1;
  LCornerCoords_[2][16] =   0;
 
  LCornerCoords_[0][17] =   1;
  LCornerCoords_[1][17] =  -1;
  LCornerCoords_[2][17] =   0;
 
  LCornerCoords_[0][18] =   1;
  LCornerCoords_[1][18] =   1;
  LCornerCoords_[2][18] =   0;
 
  LCornerCoords_[0][19] =  -1;
  LCornerCoords_[1][19] =   1;
  LCornerCoords_[2][19] =   0;

}



void Hexa2FE :: CalcShapeFnc(Vector<Double> & Shape, 
			     const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Hexa2FE::CalcShapeFnc" );

  Shape.Resize(NumNodes_);

 
  //nodes 1-8
  for( Integer i=0; i<8; i++)
    Shape[i] = 0.125 * (1 + LCornerCoords_[0][i] * LCoord[0])
                     * (1 + LCornerCoords_[1][i] * LCoord[1]) 
                     * (1 + LCornerCoords_[2][i] * LCoord[2])
                     * (-2 + LCornerCoords_[0][i] * LCoord[0]
                           + LCornerCoords_[1][i] * LCoord[1]
                           + LCornerCoords_[2][i] * LCoord[2]);
  
  //nodes 9,11,14,18
  Integer k =8;
  for( Integer i=0; i<4; i++) {
    Shape[k] = 0.25 * (1 - LCoord[0]*LCoord[0])
                    * (1 + LCornerCoords_[1][k] * LCoord[1]) 
                    * (1 + LCornerCoords_[2][k] * LCoord[2]);
    k +=2;
  }
  
  //nodes 10,12,14,16
  k = 9;
  for( Integer i=0; i<4; i++) {
    Shape[k] = 0.25 * (1 - LCoord[1]*LCoord[1])
                    * (1 + LCornerCoords_[0][k] * LCoord[0]) 
                    * (1 + LCornerCoords_[2][k] * LCoord[2]);
    k +=2;
  }

  //nodes 17, 18, 19, 20
  for( Integer i=16; i<20; i++) 
    Shape[i] = 0.25 * (1 - LCoord[2]*LCoord[2])
                    * (1 + LCornerCoords_[0][i] * LCoord[0]) 
                    * (1 + LCornerCoords_[1][i] * LCoord[1]);
}

void Hexa2FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				     const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Hexa2FE::CalcLocalDerivShapeFnc" );

  LDeriv.Resize(NumNodes_,Dim_);

  //nodes 1-8
  for( Integer i=0; i<8; i++)
    {
      LDeriv[i][0] = 0.125 * LCornerCoords_[0][i] 
	* (1 + LCornerCoords_[1][i] * LCoord[1] )* (1 + LCornerCoords_[2][i] * LCoord[2])
        * (-1 + 2 * LCornerCoords_[0][i] * LCoord[0] 
              + LCornerCoords_[1][i] * LCoord[1]
              + LCornerCoords_[2][i] * LCoord[2]);

      LDeriv[i][1] = 0.125 * LCornerCoords_[1][i] 
	* (1 + LCornerCoords_[0][i] * LCoord[0] )* (1 + LCornerCoords_[2][i] * LCoord[2])
        * (-1 + LCornerCoords_[0][i] * LCoord[0] 
              + 2 * LCornerCoords_[1][i] * LCoord[1]
              + LCornerCoords_[2][i] * LCoord[2]);

      LDeriv[i][2] = 0.125 * LCornerCoords_[2][i] 
	* (1 + LCornerCoords_[0][i] * LCoord[0] )* (1 + LCornerCoords_[1][i] * LCoord[1])
        * (-1 + LCornerCoords_[0][i] * LCoord[0] 
              + LCornerCoords_[1][i] * LCoord[1]
              + 2 * LCornerCoords_[2][i] * LCoord[2]);

    }

  //nodes 9,11,14,18
  Integer k =8;
  for( Integer i=0; i<4; i++) {
    LDeriv[k][0] = -0.5  * LCoord[0]
                         * (1 + LCornerCoords_[1][k] * LCoord[1]) 
                         * (1 + LCornerCoords_[2][k] * LCoord[2]);

    LDeriv[k][1] = 0.25 * LCornerCoords_[1][k]
                        * (1 - LCoord[0] * LCoord[0]) 
                        * (1 + LCornerCoords_[2][k] * LCoord[2]);

    LDeriv[k][2] = 0.25 * LCornerCoords_[2][k]
                       * (1 - LCoord[0] * LCoord[0]) 
                       * (1 + LCornerCoords_[1][k] * LCoord[1]);

    k +=2;
  }
  

  //nodes 10,12,14,16
  k = 9;
  for( Integer i=0; i<4; i++) {
    LDeriv[k][0] = 0.25 * LCornerCoords_[0][k]
                    * (1 - LCoord[1] * LCoord[1]) 
                    * (1 + LCornerCoords_[2][k] * LCoord[2]);

    LDeriv[k][1] = -0.5 * LCoord[1]
                    * (1 + LCornerCoords_[0][k] * LCoord[0]) 
                    * (1 + LCornerCoords_[2][k] * LCoord[2]);

    LDeriv[k][2] = 0.25 * LCornerCoords_[2][k]
                    * (1 + LCornerCoords_[0][k] * LCoord[0])
                    * (1 - LCoord[1] * LCoord[1]);

    k +=2;
  }


  //nodes 17, 18, 19, 20
  for( Integer i=16; i<20; i++) {
    LDeriv[i][0] = 0.25 * LCornerCoords_[0][i]
                    * (1 + LCornerCoords_[1][i] * LCoord[1]) 
                    * (1 - LCoord[2] * LCoord[2]);

    LDeriv[i][1] = 0.25 * LCornerCoords_[1][i]
                    * (1 + LCornerCoords_[0][i] * LCoord[0]) 
                    * (1 - LCoord[2] * LCoord[2]);

    LDeriv[i][2] = -0.5 * LCoord[2]
                    * (1 + LCornerCoords_[0][i] * LCoord[0]) 
                    * (1 + LCornerCoords_[1][i] * LCoord[1]);
  }
  
}

} // end of namespace

