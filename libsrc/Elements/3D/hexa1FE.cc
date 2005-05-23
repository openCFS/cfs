#include <iostream>
#include <fstream>

//#include <General/environment.hh>
#include "hexa1FE.hh"

namespace CoupledField
{

  Hexa1FE::Hexa1FE():HexaFE()
  { 
    ENTER_FCN( "Hexa1FE::Hexa1FE" );

    Init();
  }

  Hexa1FE::~Hexa1FE()
  {
    ENTER_FCN( "Hexa1FE::~Hexa1FE" );

  }

  void Hexa1FE::Init()
  {
    ENTER_IFCN( "Hexa1FE::Init" );
  
    Dim_ = 3;
    NumNodes_ = 8;
    NumEdges_ = 12;
    //NumFaces_ = 6;

    // first set integration points and corner coords ...
    SetIntPoints();
    SetCornerCoords();


    // ... then calc shape function values at integration points
    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();


  }

  void Hexa1FE::SetCornerCoords()
  {
    ENTER_IFCN( "Hexa1FE::SetCornerCoords" );

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =  -1;
    LCornerCoords_[1][0] =  -1;
    LCornerCoords_[2][0] =  -1;

    LCornerCoords_[0][1] =  1;
    LCornerCoords_[1][1] =  -1;
    LCornerCoords_[2][1] =  -1;

    LCornerCoords_[0][2] =  1;
    LCornerCoords_[1][2] =  1;
    LCornerCoords_[2][2] =  -1;

    LCornerCoords_[0][3] =  -1;
    LCornerCoords_[1][3] =  1;
    LCornerCoords_[2][3] =  -1;

    LCornerCoords_[0][4] =  -1;
    LCornerCoords_[1][4] =  -1;
    LCornerCoords_[2][4] =  1;

    LCornerCoords_[0][5] =  1;
    LCornerCoords_[1][5] =  -1;
    LCornerCoords_[2][5] =  1;

    LCornerCoords_[0][6] =  1;
    LCornerCoords_[1][6] =  1;
    LCornerCoords_[2][6] =  1;

    LCornerCoords_[0][7] =  -1;
    LCornerCoords_[1][7] =  1;
    LCornerCoords_[2][7] =  1;

  }



  void Hexa1FE :: CalcShapeFnc(Vector<Double> & Shape, 
                               const Vector<Double> & LCoord)
  {
    ENTER_IFCN( "Hexa1FE::CalcShapeFnc" );

    Shape.Resize(NumNodes_);

 
    for( Integer i=0; i<NumNodes_; i++)
      Shape[i] = 0.125 * (1 + LCornerCoords_[0][i] * LCoord[0])
        * (1 + LCornerCoords_[1][i] * LCoord[1]) * (1 + LCornerCoords_[2][i] * LCoord[2]);

  }

  void Hexa1FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                       const Vector<Double> & LCoord)
  {
    ENTER_IFCN( "Hexa1FE::CalcLocalDerivShapeFnc" );

    LDeriv.Resize(NumNodes_,Dim_);

    for( Integer i=0; i<NumNodes_; i++)
      {
        LDeriv[i][0] = 0.125 * LCornerCoords_[0][i] 
          * (1 + LCornerCoords_[1][i] * LCoord[1] )* (1 + LCornerCoords_[2][i] * LCoord[2]);

        LDeriv[i][1] = 0.125 * (1 + LCornerCoords_[0][i] * LCoord[0] )
          * LCornerCoords_[1][i] * (1 + LCornerCoords_[2][i] * LCoord[2]);

        LDeriv[i][2] = 0.125 * (1 + LCornerCoords_[0][i] * LCoord[0])
          * (1 + LCornerCoords_[1][i] * LCoord[1]) * LCornerCoords_[2][i];
      }
  
  }

} // end of namespace

