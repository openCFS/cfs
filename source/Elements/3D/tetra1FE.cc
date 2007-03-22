// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "tetra1FE.hh"

namespace CoupledField
{

  Tetra1FE::Tetra1FE():TetraFE()
  { 
    ENTER_FCN( "Tetra1FE::Tetra1FE" );

    Init();
  }


  Tetra1FE::~Tetra1FE()
  {
    ENTER_FCN( "Tetra1FE::~Tetra1FE" );
  
  }

  void Tetra1FE::Init()
  {
    ENTER_IFCN( "Tetra1FE::Init" );
  
    NumNodes_ = 4;
    NumEdges_ = 6;
    
    CommonInit();
  }


  void Tetra1FE::SetCornerCoords()
  {
    ENTER_IFCN( "Tetra1FE::SetCornerCoords" );

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] =  0;
    LCornerCoords_[1][0] =  0;
    LCornerCoords_[2][0] =  0;

    LCornerCoords_[0][1] =  1;
    LCornerCoords_[1][1] =  0;
    LCornerCoords_[2][1] =  0;

    LCornerCoords_[0][2] =  0;
    LCornerCoords_[1][2] =  1;
    LCornerCoords_[2][2] =  0;

    LCornerCoords_[0][3] =  0;
    LCornerCoords_[1][3] =  0;
    LCornerCoords_[2][3] =  1;

    //   LCornerCoords_[0][0] =  1;
    //   LCornerCoords_[1][0] =  0;
    //   LCornerCoords_[2][0] =  0;

    //   LCornerCoords_[0][1] =  0;
    //   LCornerCoords_[1][1] =  1;
    //   LCornerCoords_[2][1] =  0;

    //   LCornerCoords_[0][2] =  0;
    //   LCornerCoords_[1][2] =  0;
    //   LCornerCoords_[2][2] =  1;

    //   LCornerCoords_[0][3] =  0;
    //   LCornerCoords_[1][3] =  0;
    //   LCornerCoords_[2][3] =  0;
  }

  /// defines the connection between nodes with "their" edge 
  void Tetra1FE :: SetEdgeVertices()
  {
    ENTER_IFCN( "SetEdgeVertices" );
    const UInt nrNodesPerEdge = 2;
  
    edgeVertices_.Resize(NumEdges_, nrNodesPerEdge);

    edgeVertices_[0][0] = 0;
    edgeVertices_[0][1] = 1;

    edgeVertices_[1][0] = 0;
    edgeVertices_[1][1] = 2;

    edgeVertices_[2][0] = 0;
    edgeVertices_[2][1] = 3;

    edgeVertices_[3][0] = 1;
    edgeVertices_[3][1] = 2;

    edgeVertices_[4][0] = 1;
    edgeVertices_[4][1] = 3;

    edgeVertices_[5][0] = 2;
    edgeVertices_[5][1] = 3;

    //   edgeVertices_[0][0] = 3;
    //   edgeVertices_[0][1] = 0;

    //   edgeVertices_[1][0] = 3;
    //   edgeVertices_[1][1] = 1;

    //   edgeVertices_[2][0] = 3;
    //   edgeVertices_[2][1] = 2;

    //   edgeVertices_[3][0] = 0;
    //   edgeVertices_[3][1] = 1;

    //   edgeVertices_[4][0] = 0;
    //   edgeVertices_[4][1] = 2;

    //   edgeVertices_[5][0] = 1;
    //   edgeVertices_[5][1] = 2;
  }



  // std::ostream& operator<< (std::ostream & outStr, Vector<Double> xOut)
  // {
  //   for (UInt i=0; i<xOut.size(); i++)
  //     outStr <<  " " << xOut[i];
  //   return outStr;
  // }


  void Tetra1FE :: CalcShapeFnc( Vector<Double> & Shape, 
                                 const Vector<Double> & LCoord,
                                 const Elem*, UInt dof,
                                 AnsatzFct::FctEntityType fctEntityType )
  {
    ENTER_IFCN( "Tetra1FE::CalcShapeFnc" );


    Shape.Resize(NumNodes_);

    // see Hughes p. 170
  
    Shape[0] = 1.0 - LCoord[0] - LCoord[1] - LCoord[2];

    if (Shape[0] < 0)
      Error("Local coordinates are not inside tetrahedral element!",__FILE__,__LINE__);

    for( UInt i=1; i<NumNodes_; i++)
      Shape[i] = LCoord[i-1];

  }


  void Tetra1FE :: CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
                                           const Vector<Double> & LCoord,
                                           const Elem*, UInt dof,
                                           AnsatzFct::FctEntityType fctEntityType )
  {
    ENTER_IFCN( "Tetra1FE::CalcLocalDerivShapeFnc" );

    LDeriv.Resize(NumNodes_,Dim_);

    LDeriv.Init();
  
    for( UInt i=0; i<Dim_; i++)
      LDeriv[0][i] = -1.0;

    for( UInt i=1; i < NumNodes_; i++)
      LDeriv[i][i-1] = 1.0;
  }

  // see Kaltenbacher: "Numerical Sim. of Mechatronic Sensors and Actuators" p. 25
  // calculates the edge shape function of a tetrahedral of first order.
  void Tetra1FE :: CalcEdgeShapeFnc(Matrix<Double> & edgeShape, 
                                    const Vector<Double> & LCoord, 
                                    const Matrix<Double> & cornernodes)
  {
    ENTER_IFCN( "Tetra1FE::CalcShapeFnc" );

    edgeShape.Resize(NumEdges_, Dim_);
    edgeShape.Init();

    // nodal shape functions of a tet
    Vector<Double> nodeShape;
    CalcShapeFnc(nodeShape, LCoord, NULL, 1, AnsatzFct::NODE);


    // local derivates of nodal tet, dimension: nrNodes x Dim_
    Matrix<Double> xDxi;  
    //  CalcLocalDerivShapeFnc(xDxi, localcoord);
    GetGlobDerivShFnc(xDxi, LCoord, cornernodes, NULL, 1);    

    for (UInt actEdge=0; actEdge<NumEdges_; actEdge++)
      {
        UInt node1 = edgeVertices_[actEdge][0];
        UInt node2 = edgeVertices_[actEdge][1];
      
        for (UInt actDim=0; actDim<Dim_; actDim++)
          edgeShape[actEdge][actDim] = 
            nodeShape[node1] * xDxi[node2][actDim] - 
            nodeShape[node2] * xDxi[node1][actDim];
      }  
  }



  // calculated the Nedelec shape function in an arbitrary point
  void Tetra1FE :: GetEdgeGlobalDerivShapeFnc(Vector< Matrix<Double>* > & shapeDeriv, 
                                              const Vector<Double> & lCoord,
                                              const Matrix<Double> & cornerCoords)
  {
    ENTER_IFCN( "Tetra1FE::GetEdgeGlobalDerivShapeFnc" );

    shapeDeriv.Resize(NumEdges_);


    // local derivates of nodal tet, dimension: nrNodes x Dim_
    Matrix<Double> xDxi;  
    GetGlobDerivShFnc(xDxi, lCoord, cornerCoords, NULL);
  
  
    for (UInt actEdge=0; actEdge<NumEdges_; actEdge++)
      {
        shapeDeriv[actEdge]->Resize(Dim_,Dim_);

        UInt node1 = edgeVertices_[actEdge][0];
        UInt node2 = edgeVertices_[actEdge][1];

        for (UInt dim1=0; dim1<Dim_; dim1++)
          for (UInt dim2=0; dim2<Dim_; dim2++)
            (*shapeDeriv[actEdge]) [dim2][dim1] = 
              xDxi[node1][dim1] * xDxi[node2][dim2] -
              xDxi[node1][dim2] * xDxi[node2][dim1];
      }
  }


} // end of namespace

