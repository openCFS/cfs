// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "Utils/mathfunctions.hh"
#include "Domain/elem.hh"
#include "triangle1fe.hh"

namespace CoupledField
{

  Triangle1FE :: Triangle1FE(IntegrationMethod method, int order) : TriangleFE()
  {
    ENTER_FCN( "Triangle1FE::Triangle1FE" );

    Init(method, order);
  }
  
  Triangle1FE :: ~Triangle1FE()
  {
    ENTER_FCN( "Triangle1FE::~Triangle1FE" );
  }

  void Triangle1FE :: Init(IntegrationMethod method, int order)
  {
    ENTER_FCN( "Triangle1FE::Init" );
    NumNodes_ = 3;

    CommonInit(method, order);
  }

  void Triangle1FE :: SetCornerCoords()
  {
    ENTER_IFCN( "Triangle1FE::SetCornerCoords" );

    LCornerCoords_.Resize(Dim_,NumNodes_);
  
    LCornerCoords_[0][0] = 0;
    LCornerCoords_[1][0] = 0;
    LCornerCoords_[0][1] = 1;
    LCornerCoords_[1][1] = 0;
    LCornerCoords_[0][2] = 0;
    LCornerCoords_[1][2] = 1;

  }

  void Triangle1FE :: SetEdgeIndices() {
    ENTER_IFCN( "Triangle1FE::SetEdgeIndices" );
    
    edgeIndices_ = new StdVector<UInt>[NumEdges_];
    for (UInt i=0; i<NumEdges_; i++) {
      edgeIndices_[i].Resize(2);
    }
    
    edgeIndices_[0][0] = 1;
    edgeIndices_[0][1] = 2;
    edgeIndices_[1][0] = 2;
    edgeIndices_[1][1] = 3;
    edgeIndices_[2][0] = 3;
    edgeIndices_[2][1] = 1;
  }


  void Triangle1FE :: CalcShapeFnc( Vector<Double> & Shape, 
                                    const Vector<Double> & LCoord,
                                    const Elem* el, UInt dof,
                                    AnsatzFct::FctEntityType )
  {
    ENTER_IFCN( "Triangle1FE::CalcShapeFnc" );

    Shape.Resize(NumNodes_);

    Shape[0] = 1.0 - LCoord[0] - LCoord[1];

    if (Shape[0] < 0)
        EXCEPTION("Local coordinates are not inside triangular element!");
    
    for( UInt i=1; i<NumNodes_; i++)
      Shape[i] = LCoord[i-1];

#ifdef DEBUG
    //  (*debug) << "LCoord \n " << LCoord << std::endl;
    //  (*debug) << "Shape \n " << Shape << std::endl;
#endif

  }


  void Triangle1FE :: CalcLocalDerivShapeFnc( Matrix<Double> & LDeriv, 
                                              const Vector<Double> & LCoord,
                                              const Elem*, UInt dof,
                                              AnsatzFct::FctEntityType )
  {
    ENTER_IFCN( "Triangle1FE::CalcLocalDerivShapeFnc" );

    LDeriv.Resize(NumNodes_,Dim_);

    LDeriv[0][0] = -1;
    LDeriv[0][1] = -1;
    LDeriv[1][0] =  1;
    LDeriv[1][1] =  0;
    LDeriv[2][0] =  0;
    LDeriv[2][1] =  1; 
  }

  void Triangle1FE :: Global2LocalCoords(Matrix<Double> & localCoords,
                                         const Matrix<Double> & globalCoords,
                                         const Matrix<Double> & coordMat)
  {
      Vector<Double> c0, c1, c2; // endpoint-coordinates
      Vector<Double> dummy;
      Vector<Double> bCoords;
      UInt globDim = globalCoords.GetSizeRow();

    // Get coordinates of the endpoints
    c0.Resize(3);
    c1.Resize(3);
    c2.Resize(3);
    dummy.Resize(3);
    bCoords.Resize(3);
    
    //    std::cout << "SIMON> Line1FE->Global2Local(): coordMat " << coordMat << std::endl;
    //    std::cout << "SIMON> Line1FE->Global2Local(): globalCoords " << globalCoords << std::endl;

    for(UInt i = 0; i < globDim; i++)
    {
        c0[i] = coordMat[i][0];
        c1[i] = coordMat[i][1];
        c2[i] = coordMat[i][2];
    }


    localCoords.Resize(Dim_, globalCoords.GetSizeCol());
    
    for(UInt i=0; i < globalCoords.GetSizeCol(); i++)
    {
        for(UInt j = 0; j < globDim; j++)
        {
            dummy[j] = globalCoords[j][i];
        }

        GetBarycentricCoords(c0, c1, c2, dummy, bCoords);
        

        localCoords[0][i] = bCoords[0] * LCornerCoords_[0][0] +
                            bCoords[1] * LCornerCoords_[0][1] +
                            bCoords[2] * LCornerCoords_[0][2];
        
                            
        localCoords[1][i] = bCoords[0] * LCornerCoords_[1][0] +
                            bCoords[1] * LCornerCoords_[1][1] +
                            bCoords[2] * LCornerCoords_[1][2];

        //        std::cout << "SIMON> Triangle1FE->Global2Local(): localCoord[" << i << "]: (" <<localCoords[0][i]<< ", " << localCoords[1][i] << ")" << std::endl;
        /*        
        std::cout << "SIMON> Line1FE->Global2Local(): s " << s << std::endl;

        for(UInt j = 0; j < Dim_; j++)
        {
            localCoords[j][i] = LCornerCoords_[j][0] * s + LCornerCoords_[j][1] * (1-s);
            std::cout << "SIMON> Line1FE->Global2Local(): localcoord " <<localCoords[j][i]<< std::endl;
        }
        */
    }
    
  }


} // end of namespace


