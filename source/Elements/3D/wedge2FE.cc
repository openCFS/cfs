// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Elements/3D/wedgeFE.hh"
#include "MatVec/matrix.hh"
#include "wedge2FE.hh"

namespace CoupledField
{

  Wedge2FE::Wedge2FE():WedgeFE()
  {
    Init();
  }


  Wedge2FE::~Wedge2FE()
  {
  }

  void Wedge2FE::Init()
  {

    NumNodes_ = 15;

    CommonInit();
  }


  void Wedge2FE::SetCornerCoords()
  {

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

    LCornerCoords_[0][13] =  1.0;
    LCornerCoords_[1][13] =  0.0;
    LCornerCoords_[2][13] =  0.0;

    LCornerCoords_[0][14] =  0.0;
    LCornerCoords_[1][14] =  1.0;
    LCornerCoords_[2][14] =  0.0;


  }


  void Wedge2FE :: CalcShapeFnc(Vector<Double> & shape, 
                                const Vector<Double> & LCoord,
                                const Elem*, UInt dof,
                                AnsatzFct::FctEntityType fctEntityType )
  {

    const Double xi = LCoord[0];
    const Double eta = LCoord[1];
    const Double zeta = LCoord[2];

    shape.Resize( 15 );
    shape[0] = (-0.5*xi - 0.5*eta + 0.5)*(zeta*zeta + (zeta - 1)*(2*xi + 2*eta - 1) - 1);
    shape[1] = xi*(-1.0*xi*zeta + 1.0*xi + 0.5*zeta*zeta + 0.5*zeta - 1.0);
    shape[2] = eta*(-1.0*eta*zeta + 1.0*eta + 0.5*zeta*zeta + 0.5*zeta - 1.0);
    shape[3] = (xi + eta - 1)*(-0.5*zeta*zeta + 0.5*(zeta + 1)*(2*xi + 2*eta - 1) + 0.5);
    shape[4] = xi*(1.0*xi*zeta + 1.0*xi + 0.5*zeta*zeta - 0.5*zeta - 1.0);
    shape[5] = eta*(1.0*eta*zeta + 1.0*eta + 0.5*zeta*zeta - 0.5*zeta - 1.0);
    shape[6] = 2*xi*(zeta - 1)*(xi + eta - 1);
    shape[7] = 2*xi*eta*(-zeta + 1);
    shape[8] = 2*eta*(zeta - 1)*(xi + eta - 1);
    shape[9] = -2*xi*(zeta + 1)*(xi + eta - 1);
    shape[10] = 2*xi*eta*(zeta + 1);
    shape[11] = -2*eta*(zeta + 1)*(xi + eta - 1);
    shape[12] = (zeta*zeta - 1)*(xi + eta - 1);
    shape[13] = xi*(-zeta*zeta + 1);
    shape[14] = eta*(-zeta*zeta + 1);
  }


  void Wedge2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & deriv, 
                                          const Vector<Double> & LCoord,
                                          const Elem*, UInt dof,
                                          AnsatzFct::FctEntityType fctEntityType )
  {

    const Double xi = LCoord[0];
    const Double eta = LCoord[1];
    const Double zeta = LCoord[2];

    deriv.Resize(NumNodes_,Dim_);
    deriv.Init();

    // symbolically solved using sympy

    deriv[0][0] = -2.0*xi*zeta + 2.0*xi - 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta + 1.5*zeta - 1.0;
    deriv[0][1] = -2.0*xi*zeta + 2.0*xi - 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta + 1.5*zeta - 1.0;
    deriv[0][2] = (xi + eta - 1)*(-1.0*xi - 1.0*eta - 1.0*zeta + 0.5);

    deriv[1][0] = -2.0*xi*zeta + 2.0*xi + 0.5*zeta*zeta + 0.5*zeta - 1.0;
    deriv[1][1] = 0;
    deriv[1][2] = xi*(-1.0*xi + 1.0*zeta + 0.5);

    deriv[2][0] = 0;
    deriv[2][1] = -2.0*eta*zeta + 2.0*eta + 0.5*zeta*zeta + 0.5*zeta - 1.0;
    deriv[2][2] = eta*(-1.0*eta + 1.0*zeta + 0.5);

    deriv[3][0] = 2.0*xi*zeta + 2.0*xi + 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta - 1.5*zeta - 1.0;
    deriv[3][1] = 2.0*xi*zeta + 2.0*xi + 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta - 1.5*zeta - 1.0;
    deriv[3][2] = (xi + eta - 1)*(1.0*xi + 1.0*eta - 1.0*zeta - 0.5);

    deriv[4][0] = 2.0*xi*zeta + 2.0*xi + 0.5*zeta*zeta - 0.5*zeta - 1.0;
    deriv[4][1] = 0;
    deriv[4][2] = xi*(1.0*xi + 1.0*zeta - 0.5);

    deriv[5][0] = 0;
    deriv[5][1] = 2.0*eta*zeta + 2.0*eta + 0.5*zeta*zeta - 0.5*zeta - 1.0;
    deriv[5][2] = eta*(1.0*eta + 1.0*zeta - 0.5);

    deriv[6][0] = 4*xi*zeta - 4*xi + 2*eta*zeta - 2*eta - 2*zeta + 2;
    deriv[6][1] = 2*xi*(zeta - 1);
    deriv[6][2] = 2*xi*(xi + eta - 1);

    deriv[7][0] = 2*eta*(-zeta + 1);
    deriv[7][1] = 2*xi*(-zeta + 1);
    deriv[7][2] = -2*xi*eta;

    deriv[8][0] = 2*eta*(zeta - 1);
    deriv[8][1] = 2*xi*zeta - 2*xi + 4*eta*zeta - 4*eta - 2*zeta + 2;
    deriv[8][2] = 2*eta*(xi + eta - 1);

    deriv[9][0] = -4*xi*zeta - 4*xi - 2*eta*zeta - 2*eta + 2*zeta + 2;
    deriv[9][1] = 2*xi*(-zeta - 1);
    deriv[9][2] = 2*xi*(-xi - eta + 1);

    deriv[10][0] = 2*eta*(zeta + 1);
    deriv[10][1] = 2*xi*(zeta + 1);
    deriv[10][2] = 2*xi*eta;

    deriv[11][0] = 2*eta*(-zeta - 1);
    deriv[11][1] = -2*xi*zeta - 2*xi - 4*eta*zeta - 4*eta + 2*zeta + 2;
    deriv[11][2] = 2*eta*(-xi - eta + 1);

    deriv[12][0] = zeta*zeta - 1;
    deriv[12][1] = zeta*zeta - 1;
    deriv[12][2] = 2*zeta*(xi + eta - 1);

    deriv[13][0] = -zeta*zeta + 1;  
    deriv[13][1] = 0;
    deriv[13][2] = -2*xi*zeta;

    deriv[14][0] = 0;
    deriv[14][1] = -zeta*zeta + 1;
    deriv[14][2] = -2*eta*zeta;

  }



} // end of namespace

