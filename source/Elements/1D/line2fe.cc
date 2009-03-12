// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "DataInOut/Logging/cfslog.hh"

#include "line2fe.hh"

namespace CoupledField
{

  DECLARE_LOG(line2fe)
  DEFINE_LOG(line2fe, "line2fe")

  Line2FE::Line2FE() : LineFE()
  {

    Init();
  }

  Line2FE::~Line2FE()
  {
  }

  void Line2FE::Init()
  {
     NumNodes_ = 3;

     CommonInit();
  }

  void Line2FE::SetCornerCoords()
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


  void Line2FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
                                         const Vector<Double> & LCoord,
                                         const Elem*, UInt dof,
                                         AnsatzFct::FctEntityType )
  {

    LDeriv.Resize(NumNodes_,Dim_);

    LDeriv[0][0] = 0.5*(2*LCoord[0] - 1);
    LDeriv[2][0] = -2.0*LCoord[0];
    LDeriv[1][0] = 0.5*(2*LCoord[0] + 1);

  }

  void Line2FE :: Global2LocalCoords(Matrix<Double> & localCoords,
                                     const Matrix<Double> & globalCoords,
                                     const Matrix<Double> & coordMat)
  {
    Vector<Double> c0, c1; // endpoint-coordinates
    Vector<Double> diff0, diff1;
    Vector<Double> dummy;
    Double s;
    Double dist, dist1, fac;
    UInt globDim = globalCoords.GetNumRows();

    LOG_DBG(line2fe) << "WARNING: Line2FE just uses linear mapping "
                     << "for Global2LocalCoords!";
    // Get coordinates of the endpoints
    c0.Resize(globDim);
    //    c1.Resize(Dim_);
    //    dummy.Resize(Dim_);
    diff0.Resize(globDim);
    diff1.Resize(globDim);

    //    std::cout << "SIMON> Line1FE->Global2Local(): coordMat " << coordMat << std::endl;
    //    std::cout << "SIMON> Line1FE->Global2Local(): globalCoords " << globalCoords << std::endl;

    for(UInt i = 0; i < globDim; i++)
    {
        c0[i] = coordMat[i][0];
        //        c1[i] = coordMat[i][1];
        diff1[i] = coordMat[i][1] - c0[i];
    }


    //    diff1 = c1 - c0;
    dist = diff1.NormL2();
    fac = 1.0 / dist;
    diff1 *= fac;

    localCoords.Resize(Dim_, globalCoords.GetNumCols());

    for(UInt i=0; i < globalCoords.GetNumCols(); i++)
    {
        for(UInt j = 0; j < globDim; j++)
        {
            diff0[j] = globalCoords[j][i] - c0[j];
        }
        //        diff0 = globalCoords[i] - c0;
        dist1 = diff0.NormL2();
        diff0.Inner(diff1, s);
        s *= fac;

        //        std::cout << "SIMON> Line1FE->Global2Local(): s " << s << std::endl;

        for(UInt j = 0; j < Dim_; j++)
        {
            localCoords[j][i] = LCornerCoords_[j][0] * s + LCornerCoords_[j][1] * (1-s);
            //            std::cout << "SIMON> Line1FE->Global2Local(): localcoord " <<localCoords[j][i]<< std::endl;
        }
    }
  }

} // end of namespace
