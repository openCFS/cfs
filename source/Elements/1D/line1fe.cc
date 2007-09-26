// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "line1fe.hh"

namespace CoupledField
{


  Line1FE :: Line1FE(IntegrationMethod method, int order) : LineFE()
  {

    Init(method, order);
  }
  
  Line1FE :: ~Line1FE()
  {
  }

  void Line1FE :: Init(IntegrationMethod method, int order)
  {
    NumNodes_ = 2;

    CommonInit(method, order);
    SetEdgeIndices();
  }

  void Line1FE :: SetCornerCoords()
  {

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =   1;
    LCornerCoords_[0][1] =  -1;

  }

  void Line1FE :: SetEdgeIndices () {

    edgeIndices_ = new StdVector<UInt>[NumEdges_];
    edgeIndices_[0].Resize(2);

    // edge 1
    edgeIndices_[0][0] = 1;
    edgeIndices_[0][1] = 2;

  }
    

  void Line1FE :: CalcShapeFnc(Vector<Double> & Shape, 
                               const Vector<Double> & LCoord,
                               const Elem* elem, UInt dof,
                               AnsatzFct::FctEntityType )
  {

    Shape.Resize(NumNodes_);
  
    for( UInt i=0; i<NumNodes_; i++)
      Shape[i] = 0.5*(1+LCornerCoords_[0][i] * LCoord[0]);

  }

  void Line1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
                                         const Vector<Double> & LCoord,
                                         const Elem* elem, UInt dof,
                                         AnsatzFct::FctEntityType )
  {

    LDeriv.Resize(NumNodes_,Dim_);

    for( UInt i=0; i<NumNodes_; i++)
      {
        LDeriv[i][0] = 0.5*LCornerCoords_[0][i];

      }
  
  }

  void Line1FE :: Global2LocalCoords(Matrix<Double> & localCoords,
                                     const Matrix<Double> & globalCoords,
                                     const Matrix<Double> & coordMat)
  {
    Vector<Double> c0, c1; // endpoint-coordinates
    Vector<Double> diff0, diff1;
    Vector<Double> dummy;
    Double s;
    Double dist, dist1, fac;
    UInt globDim = globalCoords.GetSizeRow();
    
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
    
    localCoords.Resize(Dim_, globalCoords.GetSizeCol());
    
    for(UInt i=0; i < globalCoords.GetSizeCol(); i++)
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
