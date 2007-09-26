// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "line2fe.hh"

namespace CoupledField
{


  Line2FE :: Line2FE() : LineFE()
  {

    Init();
  }
  
  Line2FE :: ~Line2FE()
  {
  }

  void Line2FE :: Init()
  {
     NumNodes_ = 3;
     
     CommonInit();
  }

  void Line2FE :: SetCornerCoords()
  {

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =   -1;
    LCornerCoords_[0][1] =   1;
    LCornerCoords_[0][2] =  0;

  }

  void Line2FE :: CalcShapeFnc(Vector<Double> & Shape, 
                               const Vector<Double> & LCoord,
                               const Elem*, UInt dof,
                               AnsatzFct::FctEntityType )
  {

    Shape.Resize(NumNodes_);

    Shape[0] = 0.5*LCoord[0]*(LCoord[0]-1);

    Shape[2] = 1.0 - LCoord[0]*LCoord[0];
  
    Shape[1] = 0.5*LCoord[0]*(LCoord[0]+1);

  }


  void Line2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                         const Vector<Double> & LCoord,
                                         const Elem*, UInt dof,
                                         AnsatzFct::FctEntityType )
  {

    LDeriv.Resize(NumNodes_,Dim_);

    LDeriv[0][0] = 0.5*(2*LCoord[0] - 1);
    LDeriv[2][0] = -2.0*LCoord[0];
    LDeriv[1][0] = 0.5*(2*LCoord[0] + 1); 

  }


} // end of namespace
