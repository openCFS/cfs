// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "Domain/elem.hh"
#include "tetra1FE.hh"

namespace CoupledField
{

  Tetra1FE::Tetra1FE():TetraFE()
  { 

    Init();
  }


  Tetra1FE::~Tetra1FE()
  {
  
  }

  void Tetra1FE::Init()
  {
  
    NumNodes_ = 4;
    NumEdges_ = 6;
    
    CommonInit();
    SetEdgeIndices();
  }


  void Tetra1FE::SetCornerCoords()
  {

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

  void Tetra1FE :: SetEdgeIndices() {
    
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
    edgeIndices_[2][0] = 1;
    edgeIndices_[2][1] = 4;

    // edge 4
    edgeIndices_[3][0] = 2;
    edgeIndices_[3][1] = 3;

    // edge 5
    edgeIndices_[4][0] = 2;
    edgeIndices_[4][1] = 4;

    // edge 6
    edgeIndices_[5][0] = 3;
    edgeIndices_[5][1] = 4;
        
  }
  
  void Tetra1FE::GetNumFncs(Vector<UInt>& numFcns,  
                            const shared_ptr<AnsatzFct>& fcnType, 
                            AnsatzFct::FctEntityType fctEntityType, 
                            UInt dof) {


    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      numFcns.Resize( NumNodes_ );
      numFcns.Init(1);
    } else if( fcnType->GetType() == AnsatzFct::NEDELEC ) {
      numFcns.Resize( NumEdges_ );
      numFcns.Init(1);
    }
    else {
      EXCEPTION("In base class only implemented for Lagrange functions!");
    }
  }


  void Tetra1FE::CalcShapeFnc( Vector<Double> & Shape, 
                                 const Vector<Double> & LCoord,
                                 const Elem*, UInt dof,
                                 AnsatzFct::FctEntityType fctEntityType )
  {


    Shape.Resize(NumNodes_);

    // see Hughes p. 170
  
    Shape[0] = 1.0 - LCoord[0] - LCoord[1] - LCoord[2];

    if (Shape[0] < 0)
      Error("Local coordinates are not inside tetrahedral element!",__FILE__,__LINE__);

    for( UInt i=1; i<NumNodes_; i++)
      Shape[i] = LCoord[i-1];

  }


  void Tetra1FE::CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
                                           const Vector<Double> & LCoord,
                                           const Elem*, UInt dof,
                                           AnsatzFct::FctEntityType fctEntityType )
  {

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
                                    const Matrix<Double> & cornernodes,
                                    const Elem* elem)
  {

    edgeShape.Resize(NumEdges_, Dim_);
    edgeShape.Init();

    // nodal shape functions of a tet
    Vector<Double> nodeShape;
    CalcShapeFnc(nodeShape, LCoord, NULL, 1, AnsatzFct::NODE);


    // local derivates of nodal tet, dimension: nrNodes x Dim_
    Matrix<Double> xDxi;  
    //  CalcLocalDerivShapeFnc(xDxi, localcoord);
    GetGlobDerivShFnc(xDxi, LCoord, cornernodes, NULL, 1);    

    Double factor;
    for (UInt actEdge=0; actEdge<NumEdges_; actEdge++)
    {
      UInt node1 = edgeIndices_[actEdge][0] - 1;
      UInt node2 = edgeIndices_[actEdge][1] - 1;

      factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;  
      for (UInt actDim=0; actDim<Dim_; actDim++)
        edgeShape[actEdge][actDim] = 
          ( nodeShape[node1] * xDxi[node2][actDim] - 
              nodeShape[node2] * xDxi[node1][actDim] ) * factor;
    }  
  }


  // calculated the Nedelec shape function in an arbitrary point
  void Tetra1FE :: GetEdgeGlobalDerivShapeFnc(StdVector< Matrix<Double> > & shapeDeriv, 
                                              const Vector<Double> & lCoord,
                                              const Matrix<Double> & cornerCoords,
                                              const Elem* elem)
  {

    shapeDeriv.Resize(NumEdges_);


    // local derivates of nodal tet, dimension: nrNodes x Dim_
    Matrix<Double> xDxi;  
    GetGlobDerivShFnc(xDxi, lCoord, cornerCoords, NULL);

    for (UInt actEdge=0; actEdge<NumEdges_; actEdge++)
    {
      shapeDeriv[actEdge].Resize(Dim_,Dim_);

      UInt node1 = edgeIndices_[actEdge][0] -1;
      UInt node2 = edgeIndices_[actEdge][1] -1;

      Double factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;  
      for (UInt dim1=0; dim1<Dim_; dim1++)
        for (UInt dim2=0; dim2<Dim_; dim2++)
          (shapeDeriv[actEdge])[dim2][dim1] = 
            ( xDxi[node1][dim1] * xDxi[node2][dim2] -
                xDxi[node1][dim2] * xDxi[node2][dim1] ) * factor;
    }
  }



} // end of namespace

