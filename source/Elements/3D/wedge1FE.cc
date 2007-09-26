// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "wedge1FE.hh"

namespace CoupledField
{

  Wedge1FE::Wedge1FE():WedgeFE()
  { 
    Init();
  }


  Wedge1FE::~Wedge1FE()
  {
  }

  void Wedge1FE::Init()
  {

    Dim_ = 3;
    NumNodes_ = 6;
    NumEdges_ = 9;

    CommonInit();
  }


  void Wedge1FE::SetCornerCoords()
  {

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =  0.0;
    LCornerCoords_[1][0] =  0.0;
    LCornerCoords_[2][0] = -1.0;

    LCornerCoords_[0][1] =  1.0;
    LCornerCoords_[1][1] =  0.0;
    LCornerCoords_[2][1] = -1.0;

    LCornerCoords_[0][2] =  0.0;
    LCornerCoords_[1][2] =  1.0;
    LCornerCoords_[2][2] = -1.0;

    LCornerCoords_[0][3] =  0.0;
    LCornerCoords_[1][3] =  0.0;
    LCornerCoords_[2][3] =  1.0;

    LCornerCoords_[0][4] =  1.0;
    LCornerCoords_[1][4] =  0.0;
    LCornerCoords_[2][4] =  1.0;

    LCornerCoords_[0][5] =  0.0;
    LCornerCoords_[1][5] =  1.0;
    LCornerCoords_[2][5] =  1.0;

  }


  void Wedge1FE :: CalcShapeFnc(Vector<Double> & Shape, 
                                const Vector<Double> & LCoord,
                                const Elem*, UInt dof,
                                AnsatzFct::FctEntityType fctEntityType )
  {

    Shape.Resize(NumNodes_);

    //"Wedge Elements"
    // from "Dhatt, G.: The Finite Element Method Displayed, p. 120"
    // corner nodes

    Shape[0] = 0.5 * (1 - LCoord[2]) * (1 - LCoord[0] - LCoord[1]);
    Shape[1] = 0.5 * (1 - LCoord[2]) * LCoord[0];
    Shape[2] = 0.5 * (1 - LCoord[2]) * LCoord[1];
    Shape[3] = 0.5 * (1 + LCoord[2]) * (1 - LCoord[0] - LCoord[1]); 
    Shape[4] = 0.5 * (1 + LCoord[2]) * LCoord[0];
    Shape[5] = 0.5 * (1 + LCoord[2]) * LCoord[1];

  }


  void Wedge1FE :: 
  CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                         const Vector<Double> & LCoord,
                         const Elem*, UInt dof,
                         AnsatzFct::FctEntityType fctEntityType )
  {

    LDeriv.Resize(NumNodes_,Dim_);

    LDeriv.Init();

    LDeriv[0][0] = -0.5 * (1 - LCoord[2]);
    LDeriv[0][1] = -0.5 * (1 - LCoord[2]);
    LDeriv[0][2] = -0.5 * (1 - LCoord[0] - LCoord[1]);
  
    LDeriv[1][0] =  0.5 * (1 - LCoord[2]);
    LDeriv[1][1] =  0.0;
    LDeriv[1][2] = -0.5 * LCoord[0];
  
    LDeriv[2][0] =  0.0;
    LDeriv[2][1] =  0.5 * (1 - LCoord[2]);
    LDeriv[2][2] = -0.5 * LCoord[1];
  
    LDeriv[3][0] = -0.5 * (1+ LCoord[2]);
    LDeriv[3][1] = -0.5 * (1+ LCoord[2]);
    LDeriv[3][2] =  0.5 * (1- LCoord[0] - LCoord[1]);
  
    LDeriv[4][0] =  0.5 * (1+ LCoord[2]);
    LDeriv[4][1] =  0.0;
    LDeriv[4][2] =  0.5 * LCoord[0];
  
    LDeriv[5][0] =  0.0;
    LDeriv[5][1] =  0.5 * (1 + LCoord[2]);
    LDeriv[5][2] =  0.5 * LCoord[1];
  }
  



} // end of namespace

