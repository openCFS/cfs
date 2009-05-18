// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <iostream>
#include <fstream>
#include "Domain/elem.hh"
#include <General/environment.hh>
#include "wedge1FE.hh"

namespace CoupledField
{

  Wedge1FE::Wedge1FE():WedgeFE()
  {
    ShFnc2AtIp_       = NULL;
    ShFnc2DerivAtIp_  = NULL; 
    
    ShFnc2AtIp_       = NULL;
    ShFnc2DerivAtIp_  = NULL; 
    
    Init();
  }


  Wedge1FE::~Wedge1FE()
  {
    if(ShFnc2DerivAtIp_) delete[] ShFnc2DerivAtIp_;
    if(ShFnc2AtIp_)      delete[] ShFnc2AtIp_;
  }

  void Wedge1FE::Init()
  {

    Dim_ = 3;
    NumNodes_ = 6;
    NumEdges_ = 9;

    CommonInit();
    SetEdgeIndices();
    SetShapeFnc2AtIp();
    SetShapeFnc2DerivAtIp();
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
  
  void Wedge1FE::SetEdgeIndices() {
    edgeIndices_ = new StdVector<UInt>[NumEdges_];
    for (UInt i=0; i<NumEdges_; i++) {
      edgeIndices_[i].Resize(2);
    }
    
    // edge 1
    edgeIndices_[0][0] = 1;
    edgeIndices_[0][1] = 2;
    
    // edge 2
    edgeIndices_[1][0] = 1;
    edgeIndices_[1][1] = 3;
    
    // edge 3
    edgeIndices_[2][0] = 2;
    edgeIndices_[2][1] = 3;
    
    // edge 4
    edgeIndices_[3][0] = 4;
    edgeIndices_[3][1] = 5;
    
    // edge 5
    edgeIndices_[4][0] = 4;
    edgeIndices_[4][1] = 6;
    
    // edge 6
    edgeIndices_[5][0] = 5;
    edgeIndices_[5][1] = 6;
    
    // edge 7
    edgeIndices_[6][0] = 1;
    edgeIndices_[6][1] = 4;
    
    // edge 8
    edgeIndices_[7][0] = 2;
    edgeIndices_[7][1] = 5;
    
    // edge 9
    edgeIndices_[8][0] = 3;
    edgeIndices_[8][1] = 6;
  }
  
  void Wedge1FE::CalcShapeFnc(Vector<Double> & Shape, 
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


  void Wedge1FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
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

  // see Ph.D. Sabine Zagelmayr
  void Wedge1FE::CalcEdgeShapeFnc(Matrix<Double> & edgeShape, 
                                    const Vector<Double> & LCoord, 
                                    const Matrix<Double> & cornernodes,
                                    const Elem* elem)
  {

    edgeShape.Resize(NumEdges_, Dim_);
    edgeShape.Init();

    // nodal shape functions of a tet
    Vector<Double> nodeShape;
    CalcShapeFnc2(nodeShape, LCoord);

    // local derivates for special shape functions
    Matrix<Double> xDxi;  
    GetGlobDerivShFnc2(xDxi, LCoord, cornernodes, NULL, 1);    

    Double factor;

    //edge shape functions: 1-6
    for (UInt actEdge=0; actEdge<6; actEdge++)
    {
      UInt node1 = edgeIndices_[actEdge][0] - 1;
      UInt node2 = edgeIndices_[actEdge][1] - 1;

      factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;
      if ( actEdge < 3) 
        factor *= nodeShape[3];
      else 
        factor *= nodeShape[4];

      for (UInt actDim=0; actDim<Dim_; actDim++)
        edgeShape[actEdge][actDim] = 
          ( nodeShape[node2%3] * xDxi[node1%3][actDim] - 
              nodeShape[node1%3] * xDxi[node2%3][actDim] ) * factor;
    }  

    //edge shape functions: 6-9
    for (UInt actEdge=6; actEdge<NumEdges_; actEdge++)
    {
      UInt node1 = edgeIndices_[actEdge][0] - 1;
      //        UInt node2 = edgeIndices_[actEdge][1] - 1;

      factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;
      for (UInt actDim=0; actDim<Dim_; actDim++)
        edgeShape[actEdge][actDim] = 
          nodeShape[node1] * xDxi[3][actDim] * factor;
    }  
  }


  // calculated the Nedelec shape function in an arbitrary point
  void Wedge1FE::GetEdgeGlobalDerivShapeFnc(StdVector< Matrix<Double> > & shapeDeriv, 
                                              const Vector<Double> & lCoord,
                                              const Matrix<Double> & cornerCoords,
                                              const Elem* elem)
  {

    shapeDeriv.Resize(NumEdges_);

    Vector<Double> nodeShape;
    CalcShapeFnc2(nodeShape, lCoord);

    Matrix<Double> xDxi;  
    GetGlobDerivShFnc2(xDxi, lCoord, cornerCoords, NULL, 1);

    // derivative of edge functions: 1-6  
    for (UInt actEdge=0; actEdge<6; actEdge++) {
      shapeDeriv[actEdge].Resize(Dim_,Dim_);

      UInt node1 = edgeIndices_[actEdge][0] -1;
      UInt node2 = edgeIndices_[actEdge][1] -1;

      Double fnc;
      UInt idx;

      if ( actEdge < 3) {
        fnc = nodeShape[3];
        idx = 3;
      }
      else {
        fnc = nodeShape[4];
        idx = 4;
      }

      Double factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;  

      for (UInt dim1=0; dim1<Dim_; dim1++)
        for (UInt dim2=0; dim2<Dim_; dim2++) {
          (shapeDeriv[actEdge])[dim2][dim1] = 
            ( ( xDxi[node2%3][dim1] * xDxi[node1%3][dim2] -
                xDxi[node2%3][dim2] * xDxi[node1%3][dim1] ) * fnc
                + ( nodeShape[node2%3] * xDxi[node1%3][dim2] - 
                    nodeShape[node1%3] * xDxi[node2%3][dim2] ) * xDxi[idx][dim1] 
            ) * factor;
        }
    }


    // derivative of edge functions: 7-9  
    for (UInt actEdge=6; actEdge<NumEdges_; actEdge++) {
      shapeDeriv[actEdge].Resize(Dim_,Dim_);

      UInt node1 = edgeIndices_[actEdge][0] -1;
      //      UInt node2 = edgeIndices_[actEdge][1] -1;

      Double factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;  

      for (UInt dim1=0; dim1<Dim_; dim1++)
        for (UInt dim2=0; dim2<Dim_; dim2++) {
          (shapeDeriv[actEdge])[dim2][dim1] = 
            ( xDxi[node1][dim1] * xDxi[3][dim2] ) * factor;
        }
    }

  }


  void Wedge1FE::GetGlobDerivShFnc2(Matrix<Double> & Deriv, 
                                      const Vector<Double> & LCoord,
                                      const Matrix<Double> & CornerCoords,
                                      const Elem * elem, 
                                      UInt dof )
  {

    CalcInvJacobian(JInv, LCoord, CornerCoords, elem );

    CalcLocalDerivShapeFnc2(LDeriv, LCoord, elem, dof); 

    Deriv = LDeriv * JInv;
  }


  void Wedge1FE::SetShapeFnc2AtIp()
  {

    if (!ShFnc2AtIp_) {
      ShFnc2AtIp_ = new Vector<Double>[NumIntPoints_];
    } 
    else { 
      delete[] ShFnc2AtIp_ ;
      ShFnc2AtIp_ = new Vector<Double>[NumIntPoints_];
    }

    for( UInt i=0; i<NumIntPoints_; i++ ) {
      CalcShapeFnc2( ShFnc2AtIp_[i], IntPoints_[i]);
    }

  }

  void Wedge1FE::CalcShapeFnc2(Vector<Double> & Shape, 
                                 const Vector<Double> & LCoord)
  {

    Shape.Resize(5);

    //"Wedge Elements"
    // from Ph.D. thesis of Sabine Zagelmayer

    Shape[0] = 1 - LCoord[0] - LCoord[1];
    Shape[1] = LCoord[0];
    Shape[2] = LCoord[1];
    Shape[3] = 0.5 * (1 - LCoord[2]);
    Shape[4] = 0.5 * (1 + LCoord[2]);
  }


  void Wedge1FE::SetShapeFnc2DerivAtIp()
  {

    if( !ShFnc2DerivAtIp_) {
      ShFnc2DerivAtIp_ = new Matrix<Double>[NumIntPoints_]; 
    } 
    else{ 
      delete[] ShFnc2DerivAtIp_ ;
      ShFnc2DerivAtIp_ = new Matrix<Double>[NumIntPoints_];
    }

    for( UInt i=0; i<NumIntPoints_; i++ )
      CalcLocalDerivShapeFnc2( ShFnc2DerivAtIp_[i], IntPoints_[i], 
                               NULL, 1);
  }

  void Wedge1FE::CalcLocalDerivShapeFnc2( Matrix<Double> & LDeriv, 
                                            const Vector<Double> & LCoord,
                                            const Elem*, UInt dof)
  {

    LDeriv.Resize(5,Dim_);

    LDeriv.Init();

    LDeriv[0][0] = -1.0;
    LDeriv[0][1] = -1.0;
    LDeriv[1][0] =  1.0;
    LDeriv[2][1] =  1.0;
    LDeriv[3][2] = -0.5;
    LDeriv[4][2] =  0.5;
  }

} // end of namespace

