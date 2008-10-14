// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "hexa27FE.hh"

namespace CoupledField
{

  Hexa27FE::Hexa27FE():HexaFE()
  { 

    Init();
  }

  Hexa27FE::~Hexa27FE()
  {

  }

  void Hexa27FE::Init()
  {
    NumNodes_ = 27;
    NumEdges_ = 12;

    CommonInit(); 
  }

  void Hexa27FE::SetCornerCoords()
  {

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



    LCornerCoords_[0][20] =   0;
    LCornerCoords_[1][20] =  -1;
    LCornerCoords_[2][20] =   0;
 
    LCornerCoords_[0][21] =   1;
    LCornerCoords_[1][21] =   0;
    LCornerCoords_[2][21] =   0;
 
    LCornerCoords_[0][22] =   0;
    LCornerCoords_[1][22] =   1;
    LCornerCoords_[2][22] =   0;
 
    LCornerCoords_[0][23] =  -1;
    LCornerCoords_[1][23] =   0;
    LCornerCoords_[2][23] =   0;
 


    LCornerCoords_[0][24] =   0;
    LCornerCoords_[1][24] =   0;
    LCornerCoords_[2][24] =  -1;
 
    LCornerCoords_[0][25] =   0;
    LCornerCoords_[1][25] =   0;
    LCornerCoords_[2][25] =   1;
 
    LCornerCoords_[0][26] =   0;
    LCornerCoords_[1][26] =   0;
    LCornerCoords_[2][26] =   0;
  }



  void Hexa27FE::CalcShapeFnc(Vector<Double> & Shape, 
                               const Vector<Double> & LCoord,
                               const Elem*, UInt dof,
                               AnsatzFct::FctEntityType fctEntityType )
  {
    static double shapeXi[3], shapeEta[3], shapeZeta[3];
    Shape.Resize(NumNodes_);

    shapeXi[0] = 0.5*LCoord[0]*(LCoord[0]-1);
    shapeXi[2] = 1.0 - LCoord[0]*LCoord[0];
    shapeXi[1] = 0.5*LCoord[0]*(LCoord[0]+1);

    shapeEta[0] = 0.5*LCoord[1]*(LCoord[1]-1);
    shapeEta[2] = 1.0 - LCoord[1]*LCoord[1];
    shapeEta[1] = 0.5*LCoord[1]*(LCoord[1]+1);

    shapeZeta[0] = 0.5*LCoord[2]*(LCoord[2]-1);
    shapeZeta[2] = 1.0 - LCoord[2]*LCoord[2];
    shapeZeta[1] = 0.5*LCoord[2]*(LCoord[2]+1);

    for(UInt xi=0; xi<3; xi++)
    {
      for(UInt eta=0; eta<3; eta++)
      {
        for(UInt zeta=0; zeta<3; zeta++)
        {
          Shape[xi*9+eta*3+zeta] = shapeXi[xi] * shapeEta[eta] * shapeZeta[zeta];
        }
      }
    }
  }

  void Hexa27FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                       const Vector<Double> & LCoord,
                                       const Elem*, UInt dof,
                                       AnsatzFct::FctEntityType fctEntityType )
  {
    static double shapeXi[3], shapeEta[3], shapeZeta[3];
    static double shapeDerivXi[3], shapeDerivEta[3], shapeDerivZeta[3];
    LDeriv.Resize(NumNodes_,Dim_);

    shapeXi[0] = 0.5*LCoord[0]*(LCoord[0]-1);
    shapeXi[2] = 1.0 - LCoord[0]*LCoord[0];
    shapeXi[1] = 0.5*LCoord[0]*(LCoord[0]+1);

    shapeEta[0] = 0.5*LCoord[1]*(LCoord[1]-1);
    shapeEta[2] = 1.0 - LCoord[1]*LCoord[1];
    shapeEta[1] = 0.5*LCoord[1]*(LCoord[1]+1);

    shapeZeta[0] = 0.5*LCoord[2]*(LCoord[2]-1);
    shapeZeta[2] = 1.0 - LCoord[2]*LCoord[2];
    shapeZeta[1] = 0.5*LCoord[2]*(LCoord[2]+1);


    shapeDerivXi[0] = 0.5*(2*LCoord[0] - 1);
    shapeDerivXi[2] = -2.0*LCoord[0];
    shapeDerivXi[1] = 0.5*(2*LCoord[0] + 1); 

    shapeDerivEta[0] = 0.5*(2*LCoord[1] - 1);
    shapeDerivEta[2] = -2.0*LCoord[1];
    shapeDerivEta[1] = 0.5*(2*LCoord[1] + 1); 

    shapeDerivZeta[0] = 0.5*(2*LCoord[2] - 1);
    shapeDerivZeta[2] = -2.0*LCoord[2];
    shapeDerivZeta[1] = 0.5*(2*LCoord[2] + 1); 

    for(UInt xi=0; xi<3; xi++)
    {
      for(UInt eta=0; eta<3; eta++)
      {
        for(UInt zeta=0; zeta<3; zeta++)
        {
          LDeriv[xi*9 + eta*3 + zeta][0] = shapeDerivXi[xi] * shapeEta[eta]      * shapeZeta[zeta];
          LDeriv[xi*9 + eta*3 + zeta][1] = shapeXi[xi]      * shapeDerivEta[eta] * shapeZeta[zeta];
          LDeriv[xi*9 + eta*3 + zeta][2] = shapeXi[xi]      * shapeEta[eta]      * shapeDerivZeta[zeta];
        }
      }
    }
  }


} // end of namespace

