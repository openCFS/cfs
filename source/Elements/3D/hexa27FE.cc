// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Elements/3D/hexaFE.hh"
#include "MatVec/matrix.hh"
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

    // Corners
    Shape[0]= shapeXi[0]     * shapeEta[0]      * shapeZeta[0];
    Shape[1]= shapeXi[1]     * shapeEta[0]      * shapeZeta[0];
    Shape[2]= shapeXi[1]     * shapeEta[1]      * shapeZeta[0];
    Shape[3]= shapeXi[0]     * shapeEta[1]      * shapeZeta[0];

    Shape[4]= shapeXi[0]     * shapeEta[0]      * shapeZeta[1];
    Shape[5]= shapeXi[1]     * shapeEta[0]      * shapeZeta[1];
    Shape[6]= shapeXi[1]     * shapeEta[1]      * shapeZeta[1];
    Shape[7]= shapeXi[0]     * shapeEta[1]      * shapeZeta[1];

    // Edges
    Shape[8]= shapeXi[2]     * shapeEta[0]      * shapeZeta[0];
    Shape[9]= shapeXi[1]     * shapeEta[2]      * shapeZeta[0];
    Shape[10]= shapeXi[2]     * shapeEta[1]      * shapeZeta[0];
    Shape[11]= shapeXi[0]     * shapeEta[2]      * shapeZeta[0];

    Shape[12]= shapeXi[2]     * shapeEta[0]      * shapeZeta[1];
    Shape[13]= shapeXi[1]     * shapeEta[2]      * shapeZeta[1];
    Shape[14]= shapeXi[2]     * shapeEta[1]      * shapeZeta[1];
    Shape[15]= shapeXi[0]     * shapeEta[2]      * shapeZeta[1];

    Shape[16]= shapeXi[0]     * shapeEta[0]      * shapeZeta[2];
    Shape[17]= shapeXi[1]     * shapeEta[0]      * shapeZeta[2];
    Shape[18]= shapeXi[1]     * shapeEta[1]      * shapeZeta[2];
    Shape[19]= shapeXi[0]     * shapeEta[1]      * shapeZeta[2];
    
    // Faces
    Shape[20]= shapeXi[2]     * shapeEta[0]      * shapeZeta[2];
    Shape[21]= shapeXi[1]     * shapeEta[2]      * shapeZeta[2];
    Shape[22]= shapeXi[2]     * shapeEta[1]      * shapeZeta[2];
    Shape[23]= shapeXi[0]     * shapeEta[2]      * shapeZeta[2];
    
    // Bottom
    Shape[24]= shapeXi[2]     * shapeEta[2]      * shapeZeta[0];
    
    // Top
    Shape[25]= shapeXi[2]     * shapeEta[2]      * shapeZeta[1];

    // Center
    Shape[26]= shapeXi[2]     * shapeEta[2]      * shapeZeta[2];
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

    // Corners bottom
    LDeriv[0][0]= shapeDerivXi[0]     * shapeEta[0]      * shapeZeta[0];
    LDeriv[0][1]= shapeXi[0]     * shapeDerivEta[0]      * shapeZeta[0];
    LDeriv[0][2]= shapeXi[0]     * shapeEta[0]      * shapeDerivZeta[0];
    
    LDeriv[1][0]= shapeDerivXi[1]     * shapeEta[0]      * shapeZeta[0];
    LDeriv[1][1]= shapeXi[1]     * shapeDerivEta[0]      * shapeZeta[0];
    LDeriv[1][2]= shapeXi[1]     * shapeEta[0]      * shapeDerivZeta[0];
    
    LDeriv[2][0]= shapeDerivXi[1]     * shapeEta[1]      * shapeZeta[0];
    LDeriv[2][1]= shapeXi[1]     * shapeDerivEta[1]      * shapeZeta[0];
    LDeriv[2][2]= shapeXi[1]     * shapeEta[1]      * shapeDerivZeta[0];

    LDeriv[3][0]= shapeDerivXi[0]     * shapeEta[1]      * shapeZeta[0];
    LDeriv[3][1]= shapeXi[0]     * shapeDerivEta[1]      * shapeZeta[0];
    LDeriv[3][2]= shapeXi[0]     * shapeEta[1]      * shapeDerivZeta[0];

    //Corners top
    LDeriv[4][0]= shapeDerivXi[0]     * shapeEta[0]      * shapeZeta[1];
    LDeriv[4][1]= shapeXi[0]     * shapeDerivEta[0]      * shapeZeta[1];
    LDeriv[4][2]= shapeXi[0]     * shapeEta[0]      * shapeDerivZeta[1];

    LDeriv[5][0]= shapeDerivXi[1]     * shapeEta[0]      * shapeZeta[1];
    LDeriv[5][1]= shapeXi[1]     * shapeDerivEta[0]      * shapeZeta[1];
    LDeriv[5][2]= shapeXi[1]     * shapeEta[0]      * shapeDerivZeta[1];

    LDeriv[6][0]= shapeDerivXi[1]     * shapeEta[1]      * shapeZeta[1];
    LDeriv[6][1]= shapeXi[1]     * shapeDerivEta[1]      * shapeZeta[1];
    LDeriv[6][2]= shapeXi[1]     * shapeEta[1]      * shapeDerivZeta[1];

    LDeriv[7][0]= shapeDerivXi[0]     * shapeEta[1]      * shapeZeta[1];
    LDeriv[7][1]= shapeXi[0]     * shapeDerivEta[1]      * shapeZeta[1];
    LDeriv[7][2]= shapeXi[0]     * shapeEta[1]      * shapeDerivZeta[1];

    // Edges
    LDeriv[8][0]= shapeDerivXi[2]     * shapeEta[0]      * shapeZeta[0];
    LDeriv[8][1]= shapeXi[2]     * shapeDerivEta[0]      * shapeZeta[0];
    LDeriv[8][2]= shapeXi[2]     * shapeEta[0]      * shapeDerivZeta[0];

    LDeriv[9][0]= shapeDerivXi[1]     * shapeEta[2]      * shapeZeta[0];
    LDeriv[9][1]= shapeXi[1]     * shapeDerivEta[2]      * shapeZeta[0];
    LDeriv[9][2]= shapeXi[1]     * shapeEta[2]      * shapeDerivZeta[0];

    LDeriv[10][0]= shapeDerivXi[2]     * shapeEta[1]      * shapeZeta[0];
    LDeriv[10][1]= shapeXi[2]     * shapeDerivEta[1]      * shapeZeta[0];
    LDeriv[10][2]= shapeXi[2]     * shapeEta[1]      * shapeDerivZeta[0];

    LDeriv[11][0]= shapeDerivXi[0]     * shapeEta[2]      * shapeZeta[0];
    LDeriv[11][1]= shapeXi[0]     * shapeDerivEta[2]      * shapeZeta[0];
    LDeriv[11][2]= shapeXi[0]     * shapeEta[2]      * shapeDerivZeta[0];

    
    LDeriv[12][0]= shapeDerivXi[2]     * shapeEta[0]      * shapeZeta[1];
    LDeriv[12][1]= shapeXi[2]     * shapeDerivEta[0]      * shapeZeta[1];
    LDeriv[12][2]= shapeXi[2]     * shapeEta[0]      * shapeDerivZeta[1];

    LDeriv[13][0]= shapeDerivXi[1]     * shapeEta[2]      * shapeZeta[1];
    LDeriv[13][1]= shapeXi[1]     * shapeDerivEta[2]      * shapeZeta[1];
    LDeriv[13][2]= shapeXi[1]     * shapeEta[2]      * shapeDerivZeta[1];

    LDeriv[14][0]= shapeDerivXi[2]     * shapeEta[1]      * shapeZeta[1];
    LDeriv[14][1]= shapeXi[2]     * shapeDerivEta[1]      * shapeZeta[1];
    LDeriv[14][2]= shapeXi[2]     * shapeEta[1]      * shapeDerivZeta[1];

    LDeriv[15][0]= shapeDerivXi[0]     * shapeEta[2]      * shapeZeta[1];
    LDeriv[15][1]= shapeXi[0]     * shapeDerivEta[2]      * shapeZeta[1];
    LDeriv[15][2]= shapeXi[0]     * shapeEta[2]      * shapeDerivZeta[1];

    
    LDeriv[16][0]= shapeDerivXi[0]     * shapeEta[0]      * shapeZeta[2];
    LDeriv[16][1]= shapeXi[0]     * shapeDerivEta[0]      * shapeZeta[2];
    LDeriv[16][2]= shapeXi[0]     * shapeEta[0]      * shapeDerivZeta[2];

    LDeriv[17][0]= shapeDerivXi[1]     * shapeEta[0]      * shapeZeta[2];
    LDeriv[17][1]= shapeXi[1]     * shapeDerivEta[0]      * shapeZeta[2];
    LDeriv[17][2]= shapeXi[1]     * shapeEta[0]      * shapeDerivZeta[2];

    LDeriv[18][0]= shapeDerivXi[1]     * shapeEta[1]      * shapeZeta[2];
    LDeriv[18][1]= shapeXi[1]     * shapeDerivEta[1]      * shapeZeta[2];
    LDeriv[18][2]= shapeXi[1]     * shapeEta[1]      * shapeDerivZeta[2];

    LDeriv[19][0]= shapeDerivXi[0]     * shapeEta[1]      * shapeZeta[2];
    LDeriv[19][1]= shapeXi[0]     * shapeDerivEta[1]      * shapeZeta[2];
    LDeriv[19][2]= shapeXi[0]     * shapeEta[1]      * shapeDerivZeta[2];
    
    // Faces
    LDeriv[20][0]= shapeDerivXi[2]     * shapeEta[0]      * shapeZeta[2];
    LDeriv[20][1]= shapeXi[2]     * shapeDerivEta[0]      * shapeZeta[2];
    LDeriv[20][2]= shapeXi[2]     * shapeEta[0]      * shapeDerivZeta[2];

    LDeriv[21][0]= shapeDerivXi[1]     * shapeEta[2]      * shapeZeta[2];
    LDeriv[21][1]= shapeXi[1]     * shapeDerivEta[2]      * shapeZeta[2];
    LDeriv[21][2]= shapeXi[1]     * shapeEta[2]      * shapeDerivZeta[2];

    LDeriv[22][0]= shapeDerivXi[2]     * shapeEta[1]      * shapeZeta[2];
    LDeriv[22][1]= shapeXi[2]     * shapeDerivEta[1]      * shapeZeta[2];
    LDeriv[22][2]= shapeXi[2]     * shapeEta[1]      * shapeDerivZeta[2];

    LDeriv[23][0]= shapeDerivXi[0]     * shapeEta[2]      * shapeZeta[2];
    LDeriv[23][1]= shapeXi[0]     * shapeDerivEta[2]      * shapeZeta[2];
    LDeriv[23][2]= shapeXi[0]     * shapeEta[2]      * shapeDerivZeta[2];
    
    // Bottom
    LDeriv[24][0]= shapeDerivXi[2]     * shapeEta[2]      * shapeZeta[0];
    LDeriv[24][1]= shapeXi[2]     * shapeDerivEta[2]      * shapeZeta[0];
    LDeriv[24][2]= shapeXi[2]     * shapeEta[2]      * shapeDerivZeta[0];
    
    // Top
    LDeriv[25][0]= shapeDerivXi[2]     * shapeEta[2]      * shapeZeta[1];
    LDeriv[25][1]= shapeXi[2]     * shapeDerivEta[2]      * shapeZeta[1];
    LDeriv[25][2]= shapeXi[2]     * shapeEta[2]      * shapeDerivZeta[1];

    // Center
    LDeriv[26][0]= shapeDerivXi[2]     * shapeEta[2]      * shapeZeta[2];
    LDeriv[26][1]= shapeXi[2]     * shapeDerivEta[2]      * shapeZeta[2];
    LDeriv[26][2]= shapeXi[2]     * shapeEta[2]      * shapeDerivZeta[2];

  }


} // end of namespace

