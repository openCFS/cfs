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


  void Line2FE::Global2LocalCoords(Matrix<Double> & localCoords,
                                   const Matrix<Double> & globalCoords,
                                   const Matrix<Double> & coordMat)
  {
    Vector<Double> c0; // endpoint-coordinates
    Vector<Double> diffVecWholeElem;
    Vector<Double> diffVecElemMidPoint;
    Vector<Double> diffVecToPoint;
    Double s;
    Double lengthOfWholeElem, lengthToPoint, fac;
    UInt globDim = globalCoords.GetNumRows();

    // x-----------x-----o------x
    // c0                point  coordMat[i][1]
    //
    //  ------------------------>
    //  diffVecWholeElem
    //
    //  ---------->
    //  diffVecElemMidPoint
    //
    //  ---------------->
    //  diffVecToPoint
    //
    //  lengthOfWholeElem = length(diffVecWholeElem)
    //  lengthToPoint = length(diffVecToPoint)
    //  fac = 1 / lengthOfWholeElem
    //  s = lengthToPoint / lengthWholeElem

    // Get coordinates of the endpoints
    c0.Resize(globDim);
    diffVecToPoint.Resize(globDim);
    diffVecWholeElem.Resize(globDim);
    diffVecElemMidPoint.Resize(globDim);

    for(UInt i = 0; i < globDim; i++)
    {
        c0[i] = coordMat[i][0];
        diffVecWholeElem[i] = coordMat[i][1] - c0[i];
        diffVecElemMidPoint[i] = coordMat[i][2] - c0[i];
    }

    lengthOfWholeElem = diffVecWholeElem.NormL2();
    fac = 1.0 / lengthOfWholeElem;

    // Check if element is curved and issue a warning
    diffVecElemMidPoint.Inner(diffVecWholeElem, s);
    s *= fac / diffVecElemMidPoint.NormL2();
    if(fabs(s - 1) > 1e-6) {
      WARN("Line2FE just uses linear mapping " \
           "for Global2LocalCoords and element is curved!");
    }

    localCoords.Resize(Dim_, globalCoords.GetNumCols());

    for(UInt i=0; i < globalCoords.GetNumCols(); i++)
    {
        for(UInt j = 0; j < globDim; j++)
        {
            diffVecToPoint[j] = globalCoords[j][i] - c0[j];
        }

        lengthToPoint = diffVecToPoint.NormL2();
        diffVecToPoint.Inner(diffVecWholeElem, s);
        s = fac * lengthToPoint;

        localCoords[0][i] = LCornerCoords_[0][0] + 
                            (LCornerCoords_[0][1] - LCornerCoords_[0][0]) * s;
    }
  }


} // end of namespace
