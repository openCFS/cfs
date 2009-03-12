// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "Domain/elem.hh"
#include "pyra1FE.hh"

namespace CoupledField
{

  Pyra1FE::Pyra1FE():PyraFE()
  {
    Init();
  }


  Pyra1FE::~Pyra1FE()
  {
  }

  void Pyra1FE::Init()
  {

    NumNodes_ = 5;
    NumEdges_ = 8;

    CommonInit();
    SetEdgeIndices();
  }


  void Pyra1FE::SetCornerCoords()
  {

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

  void Pyra1FE :: SetEdgeIndices() {
     
     edgeIndices_ = new StdVector<UInt>[NumEdges_];
     for (UInt i=0; i<NumEdges_; i++) {
       edgeIndices_[i].Resize(2);
     }

     // edge 1
     edgeIndices_[0][0] = 1;
     edgeIndices_[0][1] = 2;

     // edge 2
     edgeIndices_[1][0] = 2;
     edgeIndices_[1][1] = 3;

     // edge 3
     edgeIndices_[2][0] = 3;
     edgeIndices_[2][1] = 4;

     // edge 4
     edgeIndices_[3][0] = 4;
     edgeIndices_[3][1] = 1;

     // edge 5
     edgeIndices_[4][0] = 1;
     edgeIndices_[4][1] = 5;

     // edge 6
     edgeIndices_[5][0] = 2;
     edgeIndices_[5][1] = 5;

     // edge 7
     edgeIndices_[6][0] = 3;
     edgeIndices_[6][1] = 5;

     // edge 8
     edgeIndices_[7][0] = 4;
     edgeIndices_[7][1] = 5;

   }
   
  void Pyra1FE::GetNumFncs(StdVector<UInt>& numFcns,  
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
    
  void Pyra1FE::CalcShapeFnc(Vector<Double> & Shape, 
                             const Vector<Double> & LCoord,
                             const Elem*, UInt dof,
                             AnsatzFct::FctEntityType fctEntityType )
  {

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
      std::cerr << "There would be 'Local coordinates are not inside pyramidal element!' in Pyra1FE::CalcShapeFnc() ?? - Fabian\n";
      // Killme - check this!! Fabian 14.06.06
      // Error("Local coordinates are not inside pyramidal element!",__FILE__,__LINE__);

  }


  void Pyra1FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
                                         const Vector<Double> & LCoord,
                                         const Elem*, UInt dof,
                                         AnsatzFct::FctEntityType fctEntityType )
  {

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

  //"Pyramidal Edge Element"
   // J.L. Coulomb, F. Zgainski, Y. Marechal. IEEE Transactions on Magnetics, Vol. 33, No. 2, March 1997
   void Pyra1FE :: CalcEdgeShapeFnc(Matrix<Double> & edgeShape, 
                                     const Vector<Double> & LCoord, 
                                     const Matrix<Double> & cornernodes,
                                     const Elem* elem)
   {

     edgeShape.Resize(NumEdges_, Dim_);
     edgeShape.Init();

     // nodal shape functions of a tet
     Vector<Double> nodeShape;
     CalcShapeFnc(nodeShape, LCoord, NULL, 1, AnsatzFct::NODE);

     // local derivates for special shape functions
     Matrix<Double> xDxi;  
     GetGlobDerivShFnc(xDxi, LCoord, cornernodes, NULL, 1);    

     Double factor;

     //edge shape functions: 1-4
     for (UInt actEdge=0; actEdge<5; actEdge++) {
       UInt node1 = edgeIndices_[actEdge][0] - 1;
       UInt node2 = edgeIndices_[actEdge][1] - 1;
       
       factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;      
       for (UInt actDim=0; actDim<Dim_; actDim++)
         edgeShape[actEdge][actDim] = 
           (   nodeShape[node1] * ( xDxi[node2][actDim] +  xDxi[(node2+1)%4][actDim] )
             - nodeShape[node2] * ( xDxi[node1][actDim] +  xDxi[(node1+3)%4][actDim] )
           ) * factor;
     }  
     
     //edge shape functions: 5-8
     for (UInt actEdge=5; actEdge<NumEdges_; actEdge++) {
       UInt node1 = edgeIndices_[actEdge][0] - 1;
       UInt node2 = edgeIndices_[actEdge][1] - 1;
       
       factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;
       for (UInt actDim=0; actDim<Dim_; actDim++)
         edgeShape[actEdge][actDim] = 
           ( nodeShape[node1] * xDxi[node2][actDim] 
             - nodeShape[node2] * xDxi[node1][actDim] ) * factor;
     }  
   }
   
   
   // calculated the Nedelec shape function in an arbitrary point
   void Pyra1FE :: GetEdgeGlobalDerivShapeFnc(StdVector< Matrix<Double> > & shapeDeriv, 
                                              const Vector<Double> & LCoord,
                                              const Matrix<Double> & cornerCoords,
                                              const Elem* elem)
   {

     shapeDeriv.Resize(NumEdges_);

     Vector<Double> nodeShape;
     CalcShapeFnc(nodeShape, LCoord, NULL, 1, AnsatzFct::NODE);

     Matrix<Double> xDxi;  
     GetGlobDerivShFnc(xDxi, LCoord, cornerCoords, NULL, 1);    

     //edge shape functions: 1-4
     for (UInt actEdge=0; actEdge<5; actEdge++) {
       shapeDeriv[actEdge].Resize(Dim_,Dim_);

       UInt node1 = edgeIndices_[actEdge][0] - 1;
       UInt node2 = edgeIndices_[actEdge][1] - 1;

       Double factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;    
       for (UInt dim1=0; dim1<Dim_; dim1++)
         for (UInt dim2=0; dim2<Dim_; dim2++)
           (shapeDeriv[actEdge])[dim2][dim1] = 
             (  xDxi[node1][dim1] * ( xDxi[node2][dim2] +  xDxi[(node2+1)%4][dim2] )
                 - xDxi[node2][dim1] * ( xDxi[node1][dim2] +  xDxi[(node1+3)%4][dim2] )
             ) * factor;
     }

     //edge shape functions: 5-8
     for (UInt actEdge=5; actEdge<NumEdges_; actEdge++) {
       shapeDeriv[actEdge].Resize(Dim_,Dim_);

       UInt node1 = edgeIndices_[actEdge][0] - 1;
       UInt node2 = edgeIndices_[actEdge][1] - 1;

       Double factor = elem->edges[actEdge] < 0 ? -1.0 : 1.0;
       for (UInt dim1=0; dim1<Dim_; dim1++)
         for (UInt dim2=0; dim2<Dim_; dim2++)
           (shapeDeriv[actEdge])[dim2][dim1] = 
             ( xDxi[node1][dim1] *  xDxi[node2][dim2]
                                                - xDxi[node2][dim1] *  xDxi[node1][dim2] ) * factor;
     }  
   }
} // end of namespace

