#include <iostream>
#include <fstream>

#include "line1fe.hh"

namespace CoupledField
{


  Line1FE :: Line1FE() : LineFE()
  {
    ENTER_FCN( "Line1FE::Line1FE" );

    Init();
  }
  
  Line1FE :: ~Line1FE()
  {
    ENTER_FCN( "Line1FE::~Line1FE" );
  }

  void Line1FE :: Init()
  {
    ENTER_IFCN( "Line1FE::Init" );

    NumNodes_ = 2;
    SetIntPoints();
    SetCornerCoords();
    SetShapeFncAtIp();
    SetShapeFncDerivAtIp();  
  }

  void Line1FE :: SetCornerCoords()
  {
    ENTER_IFCN( "Line1FE::SetCornerCoords" );

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =   1;
    LCornerCoords_[0][1] =  -1;

  }

  void Line1FE :: CalcShapeFnc(Vector<Double> & Shape, 
                               const Vector<Double> & LCoord)
  {
    ENTER_IFCN( "Line1FE::CalcShapeFnc" );

    Shape.Resize(NumNodes_);
  
    for( UInt i=0; i<NumNodes_; i++)
      Shape[i] = 0.5*(1+LCornerCoords_[0][i] * LCoord[0]);

  }

  void Line1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                         const Vector<Double> & LCoord)
  {
    ENTER_IFCN( "Line1FE::CalcLocalDerivShapeFnc" );

    LDeriv.Resize(NumNodes_,Dim_);

    for( UInt i=0; i<NumNodes_; i++)
      {
        LDeriv[i][0] = 0.5*LCornerCoords_[0][i];

      }
  
  }


} // end of namespace
