// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "quad9fe.hh"

namespace CoupledField
{


  Quad9FE::Quad9FE() : RectangleFE()
  {

    Init();
  }

  Quad9FE::~Quad9FE()
  {
  }

  void Quad9FE::Init()
  {
    NumNodes_ = 9;

    CommonInit();
  }
  // Should be called SetNodalCoords!!
  void Quad9FE::SetCornerCoords()
  {

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

    LCornerCoords_[0][8] =  0;
    LCornerCoords_[1][8] =  0;

  }

  void Quad9FE::CalcShapeFnc(Vector<Double> & Shape,
                               const Vector<Double> & LCoord,
                               const Elem*, UInt dof,
                               AnsatzFct::FctEntityType )
  {
    static double shapeXi[3], shapeEta[3];
    Shape.Resize(NumNodes_);

    shapeXi[0] = 0.5*LCoord[0]*(LCoord[0]-1);
    shapeXi[2] = 1.0 - LCoord[0]*LCoord[0];
    shapeXi[1] = 0.5*LCoord[0]*(LCoord[0]+1);

    shapeEta[0] = 0.5*LCoord[1]*(LCoord[1]-1);
    shapeEta[2] = 1.0 - LCoord[1]*LCoord[1];
    shapeEta[1] = 0.5*LCoord[1]*(LCoord[1]+1);

    for(UInt xi=0; xi<3; xi++)
    {
      for(UInt eta=0; eta<3; eta++)
      {
        Shape[xi*3+eta] = shapeXi[xi] * shapeEta[eta];
      }
    }
  }


  void Quad9FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
                                         const Vector<Double> & LCoord,
                                         const Elem*, UInt dof,
                                         AnsatzFct::FctEntityType )
  {
    static double shapeXi[3], shapeEta[3];
    static double shapeDerivXi[3], shapeDerivEta[3];
    LDeriv.Resize(NumNodes_,Dim_);

    shapeXi[0] = 0.5*LCoord[0]*(LCoord[0]-1);
    shapeXi[2] = 1.0 - LCoord[0]*LCoord[0];
    shapeXi[1] = 0.5*LCoord[0]*(LCoord[0]+1);

    shapeEta[0] = 0.5*LCoord[1]*(LCoord[1]-1);
    shapeEta[2] = 1.0 - LCoord[1]*LCoord[1];
    shapeEta[1] = 0.5*LCoord[1]*(LCoord[1]+1);

    shapeDerivXi[0] = 0.5*(2*LCoord[0] - 1);
    shapeDerivXi[2] = -2.0*LCoord[0];
    shapeDerivXi[1] = 0.5*(2*LCoord[0] + 1); 

    shapeDerivEta[0] = 0.5*(2*LCoord[1] - 1);
    shapeDerivEta[2] = -2.0*LCoord[1];
    shapeDerivEta[1] = 0.5*(2*LCoord[1] + 1); 


    for(UInt xi=0; xi<3; xi++)
    {
      for(UInt eta=0; eta<3; eta++)
      {
        LDeriv[xi*3+eta][0] = shapeDerivXi[xi] * shapeEta[eta];
        LDeriv[xi*3+eta][1] = shapeXi[xi] * shapeDerivEta[eta];
      }
    }

  }

  Double Quad9FE::CalcMeanStrain(Matrix<Double> &cornerCoords, Matrix<Double> &displacements)
  {

    EXCEPTION("Quad9FE::CalcDistortion. This function has not yet been"
              << "implemented in quadratic (9-node) quadrilaterals");
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


