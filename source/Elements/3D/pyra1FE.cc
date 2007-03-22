// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "pyra1FE.hh"

namespace CoupledField
{

  Pyra1FE::Pyra1FE():PyraFE()
  { 
    ENTER_FCN( "Pyra1FE::Pyra1FE" );
  
    Init();
  }


  Pyra1FE::~Pyra1FE()
  {
    ENTER_FCN( "Pyra1FE::~Pyra1FE" );
  }

  void Pyra1FE::Init()
  {
    ENTER_IFCN( "Pyra1FE::Init" );
  
    NumNodes_ = 5;
    NumEdges_ = 8;

    CommonInit();
  }


  void Pyra1FE::SetCornerCoords()
  {
    ENTER_IFCN( "Pyra1FE::SetCornerCoords" );

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =  1;
    LCornerCoords_[1][0] =  1;
    LCornerCoords_[2][0] =  0;

    LCornerCoords_[0][1] =  -1;
    LCornerCoords_[1][1] =  1;
    LCornerCoords_[2][1] =  0;

    LCornerCoords_[0][2] =  -1;
    LCornerCoords_[1][2] =  -1;
    LCornerCoords_[2][2] =  0;

    LCornerCoords_[0][3] =  1;
    LCornerCoords_[1][3] =  -1;
    LCornerCoords_[2][3] =  0;

    LCornerCoords_[0][4] =  0;
    LCornerCoords_[1][4] =  0;
    LCornerCoords_[2][4] =  1;

  }




  // std::ostream& operator<< (std::ostream & outStr, Vector<Double> xOut)
  // {
  //   for (UInt i=0; i<xOut.size(); i++)
  //     outStr <<  " " << xOut[i];
  //   return outStr;
  // }


  void Pyra1FE :: CalcShapeFnc(Vector<Double> & Shape, 
                               const Vector<Double> & LCoord,
                               const Elem*, UInt dof,
                               AnsatzFct::FctEntityType fctEntityType )
  {
    ENTER_IFCN( "Pyra1FE::CalcShapeFnc" );

    Shape.Resize(NumNodes_);

    //"Pyramidal Elements"
    // F. Zgainski, J.L. Coulomb, Y. Marechal. IEEE Transactions on Magnetics, Vol. 32, No. 3, May 1996


    Shape[4] = LCoord[2];


    if (LCoord[2]==1)
      {
        Shape[0] = 0.25*((1+LCoord[0])*(1+LCoord[1])-LCoord[2]);
        Shape[1] = 0.25*((1-LCoord[0])*(1+LCoord[1])-LCoord[2]);
        Shape[2] = 0.25*((1-LCoord[0])*(1-LCoord[1])-LCoord[2]);
        Shape[3] = 0.25*((1+LCoord[0])*(1-LCoord[1])-LCoord[2]);

      }
    else
      {
        Shape[0] = 0.25*((1+LCoord[0])*(1+LCoord[1])-LCoord[2]+(LCoord[0]*LCoord[1]*LCoord[2])/(1-LCoord[2]));
        Shape[1] = 0.25*((1-LCoord[0])*(1+LCoord[1])-LCoord[2]-(LCoord[0]*LCoord[1]*LCoord[2])/(1-LCoord[2]));
        Shape[2] = 0.25*((1-LCoord[0])*(1-LCoord[1])-LCoord[2]+(LCoord[0]*LCoord[1]*LCoord[2])/(1-LCoord[2]));
        Shape[3] = 0.25*((1+LCoord[0])*(1-LCoord[1])-LCoord[2]-(LCoord[0]*LCoord[1]*LCoord[2])/(1-LCoord[2]));

      }


    if (Shape[4] < 0)
      std::cerr << "There would be 'Local coordinates are not inside pyramidal element!' in Pyra1FE :: CalcShapeFnc() ?? - Fabian\n";
      // Killme - check this!! Fabian 14.06.06
      // Error("Local coordinates are not inside pyramidal element!",__FILE__,__LINE__);

  }


  void Pyra1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                         const Vector<Double> & LCoord,
                                         const Elem*, UInt dof,
                                         AnsatzFct::FctEntityType fctEntityType )
  {
    ENTER_IFCN( "Pyra1FE::CalcLocalDerivShapeFnc" );

    LDeriv.Resize(NumNodes_,Dim_);

    LDeriv.Init();
 
    LDeriv[4][0] = 0;
    LDeriv[4][1] = 0;
    LDeriv[4][2] = 1; 

    if (LCoord[2]==1)
      {
        LDeriv[0][0] = 0.25 * (1 + LCoord[1]);
        LDeriv[0][1] = 0.25 * (1 + LCoord[0]);
        LDeriv[0][2] = -0.25;

        LDeriv[1][0] = 0.25 * (-1 - LCoord[1]);
        LDeriv[1][1] = 0.25 * (1 - LCoord[0]);
        LDeriv[1][2] = -0.25;

        LDeriv[2][0] = 0.25 * (-1 + LCoord[1]);
        LDeriv[2][1] = 0.25 * (-1 + LCoord[0]);
        LDeriv[2][2] = -0.25;
      
        LDeriv[3][0] = 0.25 * (1 - LCoord[1]);
        LDeriv[3][1] = 0.25 * (-1 - LCoord[0]);
        LDeriv[3][2] = -0.25;    

      }
    else
      {
        LDeriv[0][0] = 0.25 * (1 + LCoord[1] + (LCoord[1]*LCoord[2])/(1-LCoord[2]));
        LDeriv[0][1] = 0.25 * (1 + LCoord[0] + (LCoord[0]*LCoord[2])/(1-LCoord[2]));
        LDeriv[0][2] = 0.25 * (-1 + (LCoord[0]*LCoord[1])/((1-LCoord[2])*(1-LCoord[2])));

        LDeriv[1][0] = 0.25 * (-1 - LCoord[1] - (LCoord[1]*LCoord[2])/(1-LCoord[2]));
        LDeriv[1][1] = 0.25 * (1 - LCoord[0] - (LCoord[0]*LCoord[2])/(1-LCoord[2]));
        LDeriv[1][2] = 0.25 * (-1 - (LCoord[0]*LCoord[1])/((1-LCoord[2])*(1-LCoord[2])));

        LDeriv[2][0] = 0.25 * (-1 + LCoord[1] + (LCoord[1]*LCoord[2])/(1-LCoord[2]));
        LDeriv[2][1] = 0.25 * (-1 + LCoord[0] + (LCoord[0]*LCoord[2])/(1-LCoord[2]));
        LDeriv[2][2] = 0.25 * (-1 + (LCoord[0]*LCoord[1])/((1-LCoord[2])*(1-LCoord[2])));
      
        LDeriv[3][0] = 0.25 * (1 - LCoord[1] - (LCoord[1]*LCoord[2])/(1-LCoord[2]));
        LDeriv[3][1] = 0.25 * (-1 - LCoord[0] - (LCoord[0]*LCoord[2])/(1-LCoord[2]));
        LDeriv[3][2] = 0.25 * (-1 - (LCoord[0]*LCoord[1])/((1-LCoord[2])*(1-LCoord[2])));

      }
      

  }
  



} // end of namespace

