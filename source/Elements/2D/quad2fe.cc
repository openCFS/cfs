// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "quad2fe.hh"

namespace CoupledField
{


  Quad2FE :: Quad2FE() : RectangleFE()
  {
    ENTER_FCN( "Quad2FE::Quad2FE" );

    Init();
  }
  
  Quad2FE :: ~Quad2FE()
  {
    ENTER_FCN( "Quad2FE::~Quad2FE" );
  }

  void Quad2FE :: Init()
  {
    ENTER_IFCN( "Quad2FE::Init" );
    NumNodes_ = 8;

    CommonInit();   
  }
  // Should be called SetNodalCoords!!
  void Quad2FE :: SetCornerCoords()
  {
    ENTER_IFCN( "Quad2FE::SetCornerCoords" );

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] = -1;
    LCornerCoords_[1][0] = -1;
    LCornerCoords_[0][1] =  1;
    LCornerCoords_[1][1] = -1;
    LCornerCoords_[0][2] =  1;
    LCornerCoords_[1][2] =  1;
    LCornerCoords_[0][3] = -1;
    LCornerCoords_[1][3] =  1;

    LCornerCoords_[0][4] =  0;
    LCornerCoords_[1][4] = -1;
    LCornerCoords_[0][5] =  1;
    LCornerCoords_[1][5] =  0;
    LCornerCoords_[0][6] =  0;
    LCornerCoords_[1][6] =  1;
    LCornerCoords_[0][7] = -1;
    LCornerCoords_[1][7] =  0;

  }

  void Quad2FE :: CalcShapeFnc(Vector<Double> & Shape, 
                               const Vector<Double> & LCoord,
                               const Elem*, UInt dof,
                               AnsatzFct::FctEntityType )
  {
    ENTER_IFCN( "Quad2FE::CalcShapeFnc" );

    Shape.Resize(NumNodes_);
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 122.
    // corner nodes
    for( UInt i=0; i<(NumNodes_/2); i++)
      Shape[i] = 0.25 * (1 + LCornerCoords_[0][i] * LCoord[0])
        * (1 + LCornerCoords_[1][i] * LCoord[1])
        * (LCornerCoords_[0][i] * LCoord[0] +
           LCornerCoords_[1][i] * LCoord[1] - 1);
    // midside node
    for( UInt i=(NumNodes_/2); i<NumNodes_; i=i+2)
      {
        Shape[i]   = 0.5 * (1 - LCoord[0]*LCoord[0])*
          (1 + LCornerCoords_[1][i] * LCoord[1]);
        Shape[i+1] = 0.5 * (1 - LCoord[1]*LCoord[1])*
          (1 + LCornerCoords_[0][i+1] * LCoord[0]);
      }
  }


  void Quad2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                         const Vector<Double> & LCoord,
                                         const Elem*, UInt dof,
                                         AnsatzFct::FctEntityType )
  {
    ENTER_IFCN( "Quad2FE::CalcLocalDerivShapeFnc" );

    LDeriv.Resize(NumNodes_,Dim_);

    // corner nodes
    for( UInt i=0; i<(NumNodes_/2); i++)
      {
        LDeriv[i][0] = 0.25 * LCornerCoords_[0][i] * 
          (1 + LCornerCoords_[1][i] * LCoord[1]) * 
          (2 * LCornerCoords_[0][i] * LCoord[0] +
           LCornerCoords_[1][i] * LCoord[1]);

        LDeriv[i][1] = 0.25 * LCornerCoords_[1][i] * 
          (1 + LCornerCoords_[0][i] * LCoord[0]) * 
          (2 * LCornerCoords_[1][i] * LCoord[1] +
           LCornerCoords_[0][i] * LCoord[0]);
      }
  
    // midside node
    for( UInt i=(NumNodes_/2); i<NumNodes_; i=i+2)
      {
        LDeriv[i][0]   = - LCoord[0]*
          (1 + LCornerCoords_[1][i] * LCoord[1]);

        LDeriv[i][1]   = 0.5 * LCornerCoords_[1][i] * 
          (1 - LCoord[0]*LCoord[0]);

        LDeriv[i+1][0] = 0.5 *  LCornerCoords_[0][i+1] * 
          (1 - LCoord[1]*LCoord[1]);

        LDeriv[i+1][1] = -LCoord[1]*
          (1 + LCornerCoords_[0][i+1] * LCoord[0]);
      }

  }

  Double Quad2FE::CalcMeanStrain(Matrix<Double> &cornerCoords, Matrix<Double> &displacements)
  {
    ENTER_IFCN( "Quad2FE::CalcDistortion" );

    (*error) << "Quad2FE::CalcDistortion. This function has not yet been implemented in " 
             << " quadratic quadrilaterals" << std::endl;
    Error( __FILE__, __LINE__ );
    return -1;

    //   Double factor;
    //   Double eps1, eps2, eps4, eps5, eps11, eps12, eps21, eps22, eps41, eps42, eps51, eps52;
    //   Double length1, length2, length11, length12, length21, length22;


    //   // calculate inital size of element
    //   length11 = abs(cornerCoords[0][0] - cornerCoords[0][1]);
    //   length12 = abs(cornerCoords[0][3] - cornerCoords[0][2]);
    //   length1 = (length11+length12) * 0.5;
  
    //   length21 = abs(cornerCoords[1][0] - cornerCoords[1][3]);
    //   length22 = abs(cornerCoords[1][1] - cornerCoords[1][2]);
    //   length2 = (length21+length22) * 0.5;

    //    // calculate absolute change of size
    //   eps11 = displacements[0][0] - displacements[0][1];
    //   eps12 = displacements[0][3] - displacements[0][2];
    //   eps1 = (eps11+eps12) * 0.5;

    //   eps21 = displacements[1][0] - displacements[1][3];
    //   eps22 = displacements[1][1] - displacements[1][2];
    //   eps2 = (eps21+eps22) * 0.5;
  
    //   eps41 = displacements[0][2] - displacements[0][1];
    //   eps42 = displacements[0][3] - displacements[0][0]; 
    //   eps4 = (eps41+eps42)*0.5;
    
    //   eps51 = displacements[1][1] - displacements[1][0];
    //   eps52 = displacements[1][3] - displacements[1][2];
    //   eps5= (eps51+eps52)*0.5;  

    //   factor =  0.2 * ((eps1*eps1/(length1*length1))
    //               + (eps2*eps2/(length2*length2)) 
    //               + (eps5*eps5/(length1*length2))
    //               + (eps4*eps4/(length1*length1)));

    //   return factor;


  }


} // end of namespace


