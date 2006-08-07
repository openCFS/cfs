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
  
    NumNodes_ = 20;
    NumEdges_ = 12;

    CommonInit(); 
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

    //helping variables
    Double  rgauss,sgauss,tgauss;
    Double  onehalf,one,two;
    Double  oneplusr,oneminusr,oneminusr2;
    Double  onepluss,oneminuss,oneminuss2;
    Double  oneplust,oneminust,oneminust2;

    Shape.Resize(NumNodes_);

    //integration points
    rgauss = LCoord[0];
    sgauss = LCoord[1];
    tgauss = LCoord[2];

    onehalf = 0.5;
    one     = 1.0;
    two     = 2.0;

    oneplusr   = onehalf * (one + rgauss);
    oneminusr  = onehalf * (one - rgauss);
    onepluss   = onehalf * (one + sgauss);
    oneminuss  = onehalf * (one - sgauss);
    oneplust   = onehalf * (one + tgauss);
    oneminust  = onehalf * (one - tgauss);
    oneminusr2 = (one + rgauss) * (one - rgauss);
    oneminuss2 = (one + sgauss) * (one - sgauss);
    oneminust2 = (one + tgauss) * (one - tgauss);


    //compute the shape functions
    Shape[0]  = oneminusr  * oneminuss * oneminust;
    Shape[1]  = oneplusr   * oneminuss * oneminust;
    Shape[2]  = oneplusr   * onepluss  * oneminust;
    Shape[3]  = oneminusr  * onepluss  * oneminust;
    Shape[4]  = oneminusr  * oneminuss * oneplust;
    Shape[5]  = oneplusr   * oneminuss * oneplust;
    Shape[6]  = oneplusr   * onepluss  * oneplust;
    Shape[7]  = oneminusr  * onepluss  * oneplust;
    Shape[8]  = oneminusr2 * oneminuss * oneminust;
    Shape[9]  = oneminuss2 * oneplusr  * oneminust;
    Shape[10] = oneminusr2 * onepluss  * oneminust;
    Shape[11] = oneminuss2 * oneminusr * oneminust;
    Shape[12] = oneminusr2 * oneminuss * oneplust;
    Shape[13] = oneminuss2 * oneplusr  * oneplust;
    Shape[14] = oneminusr2 * onepluss  * oneplust;
    Shape[15] = oneminuss2 * oneminusr * oneplust;
    Shape[16] = oneminusr  * oneminuss * oneminust2;
    Shape[17] = oneminuss  * oneplusr  * oneminust2;
    Shape[18] = oneplusr   * onepluss  * oneminust2;
    Shape[19] = onepluss   * oneminusr * oneminust2;


    Shape[0] -= onehalf * (Shape[11] + Shape[8] + Shape[16]);
    Shape[1] -= onehalf * (Shape[8]  + Shape[9] + Shape[17]);
    Shape[2] -= onehalf * (Shape[9]  + Shape[10] + Shape[18]);
    Shape[3] -= onehalf * (Shape[13] + Shape[10] + Shape[19]);
    Shape[4] -= onehalf * (Shape[15] + Shape[12] + Shape[16]);
    Shape[5] -= onehalf * (Shape[12] + Shape[13] + Shape[19]);
    Shape[6] -= onehalf * (Shape[13] + Shape[14] + Shape[16]);
    Shape[7] -= onehalf * (Shape[14] + Shape[15] + Shape[19]);
  }

  void Hexa2FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                       const Vector<Double> & LCoord)
  {
    ENTER_IFCN( "Hexa2FE::CalcLocalDerivShapeFnc" );

    //helper variables
    Double  rgauss,sgauss,tgauss;
    Double  onehalf,one,two;
    Double  oneplusr,oneminusr,oneminusr2,two_r;
    Double  onepluss,oneminuss,oneminuss2,two_s;
    Double  oneplust,oneminust,oneminust2,two_t;


    LDeriv.Resize(NumNodes_,Dim_);
    LDeriv.Init(0);

    //integration points
    rgauss = LCoord[0];
    sgauss = LCoord[1];
    tgauss = LCoord[2];

    onehalf = 0.5;
    one     = 1.0;
    two     = 2.0;
    two_r   = -two * rgauss;
    two_s   = -two * sgauss;
    two_t   = -two * tgauss;

    oneplusr   = onehalf * (one + rgauss);
    oneminusr  = onehalf * (one - rgauss);
    onepluss   = onehalf * (one + sgauss);
    oneminuss  = onehalf * (one - sgauss);
    oneplust   = onehalf * (one + tgauss);
    oneminust  = onehalf * (one - tgauss);
    oneminusr2 = (one + rgauss) * (one - rgauss);
    oneminuss2 = (one + sgauss) * (one - sgauss);
    oneminust2 = (one + tgauss) * (one - tgauss);


    //deriavtives of shape functions
    LDeriv[6][0] =  onehalf * onepluss  * oneplust;
    LDeriv[5][0] =  onehalf * oneminuss * oneplust;
    LDeriv[2][0] =  onehalf * onepluss  * oneminust;
    LDeriv[1][0] =  onehalf * oneminuss * oneminust;
    LDeriv[7][0] = -LDeriv[6][0];
    LDeriv[4][0] = -LDeriv[5][0];
    LDeriv[3][0] = -LDeriv[2][0];
    LDeriv[0][0] = -LDeriv[1][0];
  
    LDeriv[6][1] = onehalf * oneplusr  * oneplust;
    LDeriv[7][1] = onehalf * oneminusr * oneplust;
    LDeriv[2][1] = onehalf * oneplusr  * oneminust;
    LDeriv[3][1] = onehalf * oneminusr * oneminust;
    LDeriv[4][1] = -LDeriv[7][1];
    LDeriv[5][1] = -LDeriv[6][1];
    LDeriv[0][1] = -LDeriv[3][1];
    LDeriv[1][1] = -LDeriv[2][1];
  
    LDeriv[6][2] = onehalf * oneplusr  * onepluss;
    LDeriv[7][2] = onehalf * oneminusr * onepluss;
    LDeriv[4][2] = onehalf * oneminusr * oneminuss;
    LDeriv[5][2] = onehalf * oneplusr  * oneminuss;
    LDeriv[2][2] = -LDeriv[6][2];
    LDeriv[3][2] = -LDeriv[7][2];
    LDeriv[0][2] = -LDeriv[4][2];
    LDeriv[1][2] = -LDeriv[5][2];


    LDeriv[14][0] =  two_r   * onepluss   * oneplust;
    LDeriv[15][0] = -onehalf * oneminuss2 * oneplust;
    LDeriv[12][0] =  two_r   * oneminuss  * oneplust;
    LDeriv[13][0] = -LDeriv[15][0];
    LDeriv[10][0] =  two_r   * onepluss   * oneminust;
    LDeriv[11][0] = -onehalf * oneminuss2 * oneminust;
    LDeriv[8][0]  =  two_r   * oneminuss  * oneminust;
    LDeriv[9][0]  = -LDeriv[11][0];
    LDeriv[18][0] =  onehalf * onepluss   * oneminust2;
    LDeriv[19][0] = -LDeriv[18][0];
    LDeriv[16][0] = -onehalf * oneminuss  * oneminust2;
    LDeriv[17][0] = -LDeriv[16][0];
  
    LDeriv[14][1] =  onehalf * oneminusr2 * oneplust;
    LDeriv[15][1] =  two_s   * oneminusr  * oneplust;
    LDeriv[12][1] = -LDeriv[14][1];
    LDeriv[13][1] =  two_s   * oneplusr   * oneplust;
    LDeriv[10][1] =  onehalf * oneminusr2 * oneminust;
    LDeriv[11][1] =  two_s   * oneminusr  * oneminust;
    LDeriv[8][1]  = -LDeriv[10][1];
    LDeriv[9][1]  =  two_s   * oneplusr   * oneminust;
    LDeriv[18][1] =  onehalf * oneplusr   * oneminust2;
    LDeriv[19][1] =  onehalf * oneminusr  * oneminust2;
    LDeriv[16][1] = -LDeriv[19][1];
    LDeriv[17][1] = -LDeriv[18][1];

    LDeriv[14][2] =  oneminusr2 * onepluss  * onehalf;
    LDeriv[15][2] =  oneminuss2 * oneminusr * onehalf;
    LDeriv[12][2] =  oneminusr2 * oneminuss * onehalf;
    LDeriv[13][2] =  oneminuss2 * oneplusr  * onehalf;
    LDeriv[10][2] = -LDeriv[14][2];
    LDeriv[11][2] = -LDeriv[15][2];
    LDeriv[8][2]  = -LDeriv[12][2];
    LDeriv[9][2]  = -LDeriv[13][2];
    LDeriv[18][2] =  oneplusr  * onepluss  * two_t;
    LDeriv[19][2] =  onepluss  * oneminusr * two_t;
    LDeriv[16][2] =  oneminusr * oneminuss * two_t;
    LDeriv[17][2] =  oneminuss * oneplusr  * two_t;


    for (UInt i=0; i<3; i++) {
      LDeriv[0][i] -= onehalf * (LDeriv[11][i] + LDeriv[8][i]  + LDeriv[16][i]);
      LDeriv[1][i] -= onehalf * (LDeriv[8][i]  + LDeriv[9][i]  + LDeriv[17][i]);
      LDeriv[2][i] -= onehalf * (LDeriv[9][i]  + LDeriv[10][i] + LDeriv[18][i]);
      LDeriv[3][i] -= onehalf * (LDeriv[10][i] + LDeriv[11][i] + LDeriv[19][i]);
      LDeriv[4][i] -= onehalf * (LDeriv[15][i] + LDeriv[12][i] + LDeriv[16][i]);
      LDeriv[5][i] -= onehalf * (LDeriv[12][i] + LDeriv[13][i] + LDeriv[17][i]);
      LDeriv[6][i] -= onehalf * (LDeriv[13][i] + LDeriv[14][i] + LDeriv[18][i]);
      LDeriv[7][i] -= onehalf * (LDeriv[14][i] + LDeriv[15][i] + LDeriv[19][i]);
    }

  }

} // end of namespace

