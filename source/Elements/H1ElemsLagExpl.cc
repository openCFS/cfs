// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "H1ElemsLagExpl.hh"

namespace CoupledField {

  // ========================================================================
  //  FeH1LagrangeExpl 
  // ========================================================================
  
  FeH1LagrangeExpl::FeH1LagrangeExpl() {
   order_ = 0; 
  }
    
  FeH1LagrangeExpl::~FeH1LagrangeExpl() {
  }
  
  void FeH1LagrangeExpl::SetIntPoints( StdVector<LocPoint>& intPoints ) {
  }
  
  void FeH1LagrangeExpl::GetNumFncs( StdVector<UInt>& numFcns,
                                     EntityType fctEntityType,
                                     UInt dof ) {
    // Initialize explicitly with number of nodes
    if( fctEntityType == VERTEX  ) {
      numFcns.Resize( shape_.numVertices );
      numFcns.Init( 1 );
    }else if( fctEntityType == EDGE ) {
      numFcns.Resize((order_-1)*shape_.numEdges);
      numFcns.Init(1);
    }else if( fctEntityType == ALL){
      numFcns.Resize(shape_.numNodes);
      numFcns.Init(1);
    }else{
      numFcns.Resize(0);
    }
  }

  void FeH1LagrangeExpl::GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                             const Elem* ptElem,
                                             EntityType fctEntityType,
                                             UInt entNumber){
    if( fctEntityType == VERTEX ) {
      fncPermutation.Resize(1);
      fncPermutation.Init(0);
    }else if( fctEntityType == EDGE ) {
      fncPermutation.Resize(order_-1);
      for ( UInt i = 0; i < order_-1 ; i++ ) {
        fncPermutation[i] = i;
      }
    }else if( fctEntityType == FACE && ptElem->faces.GetSize() > 0) {
      fncPermutation.Resize((order_-1) * (order_-1));
      for(UInt i = 0; i< order_-1 ; i++){
        for(UInt j = 0; j< order_-1 ; j++){
          fncPermutation[(i*(order_-1)) + j] = i*(order_-1) + j;
        }
      }
    }else{
      fncPermutation.Resize(0);
    }
  }


  UInt FeH1LagrangeExpl::GetNumFncsPerEntType( EntityType fctEntityType,
                                     UInt dof){
    UInt numFnc = 0;
    // Initialize explictily with number of nodes
    if( fctEntityType == VERTEX ) {
      numFnc = shape_.numVertices;
    }else if( fctEntityType == EDGE ) {
      numFnc = shape_.numEdges*(order_-1);
    }
    return numFnc;
  }
  
  void FeH1LagrangeExpl::GetMaxOrderLocDir(StdVector<UInt>& order ) {
    
    // resize vector to dimenson of local element
    order.Resize( Elem::shapes[feType_].dim);
    order.Init(order_);
  }
  
  // ========================================================================
  //  Lagrangian Elements of 1st order
  // ========================================================================
  
  // --- Line 1st order ---
  
  FeH1LagrangeLine1::FeH1LagrangeLine1() {
    feType_ = Elem::ET_LINE2;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 2;
    order_ = 1;
    elemDim_ = 1;
  }
  FeH1LagrangeLine1::~FeH1LagrangeLine1() {
    
  }
  
  void FeH1LagrangeLine1::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
     shape.Resize( 2 );
     shape[0] = 0.5 * ( 1.0 - point[0] );
     shape[1] = 0.5 * ( 1.0 + point[0] );
  }
  
  void FeH1LagrangeLine1::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {

      deriv.Resize(2, 1);
      deriv[0][0] = 0.5 * -1.0;
      deriv[1][0] = 0.5 *  1.0;
  }
  
  
  bool FeH1LagrangeLine1::CoordIsInsideElem( const Vector<Double>& point,
                                             Double tolerance )  {
   const Double & xi = point[0];
   return (xi >= (-1.0 - tolerance) &&
           xi <= (1.0 + tolerance));
  }
  
  void FeH1LagrangeLine1::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                                    const StdVector<UInt> & volConnect,
                                                    const LocPoint & surfIntPoint,
                                                    LocPoint & volIntPoint,
                                                    Vector<Double>& locNormal ) {
    EXCEPTION("Not implemented");
  }
  
  // --- Tria 1st order ---

  FeH1LagrangeTria1::FeH1LagrangeTria1() {
    feType_ = Elem::ET_TRIA3;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 3;
    order_ = 1; 
    elemDim_ = 2;
  }

  FeH1LagrangeTria1::~FeH1LagrangeTria1() {

  }

  void FeH1LagrangeTria1::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
    shape.Resize( 3 );
    shape[0] = 1.0 - point[0] - point[1]; 
    shape[1] = point[0];
    shape[2] = point[1];
  }

  void FeH1LagrangeTria1::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
    
    deriv.Resize( 3, 2 );
    deriv[0][0] = -1;
    deriv[0][1] = -1;
    
    deriv[1][0] =  1;
    deriv[1][1] =  0;
    
    deriv[2][0] =  0;
    deriv[2][1] =  1;
  }


  bool FeH1LagrangeTria1::CoordIsInsideElem( const Vector<Double>& point,
                                             Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[1];
    return ( xi >= (0 - tolerance)) &&
           ( eta >= (0 - tolerance)) &&
           ((xi + eta) <= (1 + tolerance));
  }

  void FeH1LagrangeTria1::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                            const StdVector<UInt> & volConnect,
                            const LocPoint & surfIntPoint,
                            LocPoint & volIntPoint,
                            Vector<Double>& locNormal ) {

    // Try to find out, which vertices are in common with
    // the surface element. Then calculate the product of both
    // and compare them
    //
    //            
    //            
    // 3 +        
    //   |\         eta
    //   | \       ^        REFERENCE VOLUME ELEMENT
    //   |  \      | 
    // 1 +---+ 2   +--> xi



    // NOTE: Since the line element is defined in the range [-1;+1]
    // we have to calculate (1+surfCoord)/2 in order to get the right
    // position on the triangular element

    StdVector<UInt> commonIndex(2);
    UInt found = 0;
    UInt indexProduct = 0;
    std::string errMsg;

    volIntPoint.coord.Resize(2);
    locNormal.Resize(2);

    // loop over surface connect
    for (UInt iSurf=0; iSurf<2; iSurf++)
      // loop over volume connect
      for (UInt iVol=0; iVol<3; iVol++)
        if (surfConnect[iSurf] == volConnect[iVol])
        {
          commonIndex[found++] = iVol+1;
        }

    indexProduct= commonIndex[0] * commonIndex[1];
    switch(indexProduct)
    {
      case 2:
        // Edge[1,2] is common
        volIntPoint[0] = 0.5 + (surfIntPoint[0] / 2.0);
        volIntPoint[1] = 0.0;
        locNormal[0] = 0;
        locNormal[1] = -1.0;
        break;

      case 3:
        // Edge[1,3] is common
        volIntPoint[0] = 0.0;
        volIntPoint[1] = 0.5 + (surfIntPoint[0] / 2.0);
        locNormal[0] = -1.0;
        locNormal[1] =  0.0;
        break;

      case 6:
        // Edge[2,3] is common
        volIntPoint[0] = 0.5 - (surfIntPoint[0] / 2.0);
        volIntPoint[1] = 0.5 + (surfIntPoint[0] / 2.0);
        locNormal[0] = sqrt(.5);
        locNormal[1] = sqrt(.5);
        break;

      default:
        EXCEPTION( "TriangleFE::GetLocalIntPoints4Surface: surface and volume element "
            << "have not two nodes in common. Check your .mesh-file.");
    }
  }

  
  // --- Quad 1st order ---
   
  FeH1LagrangeQuad1::FeH1LagrangeQuad1() {
    feType_ = Elem::ET_QUAD4;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 4;
    order_ = 1; 
    elemDim_ = 2;
  }
    
  FeH1LagrangeQuad1::~FeH1LagrangeQuad1() {
    
  }
  
  void FeH1LagrangeQuad1::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
     shape.Resize( 4 );
     shape[0] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ); 
     shape[1] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] );
     shape[2] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] );
     shape[3] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] );
    
  }
  
  void FeH1LagrangeQuad1::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    deriv.Resize( 4, 2 );
    for( UInt i = 0; i < 4; i++ ) {
      deriv[i][0] = 0.25 * coords[i][0] * (1 + coords[i][1] * point[1] );
      deriv[i][1] = 0.25 * (1 + coords[i][0] * point[0] ) * coords[i][1];
    }
  }
  
  
  bool FeH1LagrangeQuad1::CoordIsInsideElem( const Vector<Double>& point,
                                             Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[1];
    return  ( xi >= (-1.0 - tolerance)) &&
            (eta >= (-1.0 - tolerance)) &&
            ( xi <= (1.0 + tolerance)) &&
            (eta <= (1.0 + tolerance));  
  }
  
  void FeH1LagrangeQuad1::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                            const StdVector<UInt> & volConnect,
                            const LocPoint & surfIntPoint,
                            LocPoint & volIntPoint,
                            Vector<Double>& locNormal ) {
    // Try to find out, which vertices are in common with
      // the surface element. Then calculate the product of both
      // and compare them
      //
      //      eta
      //       ^
      // 4 +---|---+ 3    
      //   |   |   |      
      //   |   0---|-> xi     REFERENCE VOLUME ELEMENT
      //   |       |
      // 1 +-------+ 2



      StdVector<UInt> commonIndex(2);
      UInt found = 0;
      UInt indexProduct = 0;
      std::string errMsg;
    
      volIntPoint.coord.Resize(2);
      locNormal.Resize(2);
    
      // loop over surface connect
      for (UInt iSurf=0; iSurf<2; iSurf++)
        // loop over volume connect
        for (UInt iVol=0; iVol<4; iVol++)
          if (surfConnect[iSurf] == volConnect[iVol])
            {
              commonIndex[found++] = iVol+1;
            }

      indexProduct= commonIndex[0] * commonIndex[1];
      switch(indexProduct)
        {
        case 2:
          // Edge[1,2] is common
          volIntPoint[0] = surfIntPoint[0];
          volIntPoint[1] = -1.0;
          locNormal[0] =  0.0;
          locNormal[1] = -1.0;
          break;

        case 12:
          // Edge[4,3] is common
          volIntPoint[0] = surfIntPoint[0];
          volIntPoint[1] = 1.0;
          locNormal[0] =  0.0;
          locNormal[1] =  1.0;
          break;

        case 4:
          // Edge[1,4] is common
          volIntPoint[0] = -1.0;
          volIntPoint[1] = surfIntPoint[0];
          locNormal[0] = -1.0;
          locNormal[1] =  0.0;
          break;

        case 6:
          // Edge[2,3] is common
          volIntPoint[0] = 1.0;
          volIntPoint[1] = surfIntPoint[0];
          locNormal[0] =  1.0;
          locNormal[1] =  0.0;
          break;

        default:
          EXCEPTION( "RectangleFE::GetLocalIntPoints4Surface: surface and volume element "
                     <<  "have not two nodes in common. Check your .mesh-file.");
        }
  }
  
  // --- Hex 1st order ---
  FeH1LagrangeHex1::FeH1LagrangeHex1() {
    feType_ = Elem::ET_HEXA8;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 8;
    order_ = 1; 
    elemDim_ = 3;
  }
    
  FeH1LagrangeHex1::~FeH1LagrangeHex1() {
    
  }
  
  void FeH1LagrangeHex1::CalcShFnc( Vector<Double>& shape,
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp ) {
    shape.Resize( 8 );
    shape[0] = 0.125 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]); 
    shape[1] = 0.125 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]);
    shape[2] = 0.125 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
    shape[3] = 0.125 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
    shape[4] = 0.125 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]); 
    shape[5] = 0.125 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]);
    shape[6] = 0.125 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
    shape[7] = 0.125 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
    
  }
  
  void FeH1LagrangeHex1::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                            const Vector<Double>& point,
                                            const Elem* ptElem,
                                            UInt comp ) {
    deriv.Resize( 8, 3 );
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    for( UInt i = 0; i < 8; i++ ) {
      deriv[i][0] = 0.125  * coords[i][0] 
                           * (1 + coords[i][1] * point[1] ) 
                           * (1 + coords[i][2] * point[2] );
      
      deriv[i][1] = 0.125  * (1 + coords[i][0] * point[0] )
                           * coords[i][1] 
                           * (1 + coords[i][2] * point[2] );
      
      deriv[i][2] = 0.125  * (1 + coords[i][0] * point[0] ) 
                           * (1 + coords[i][1] * point[1] )
                           * coords[i][2];
    }
    
  }
  
  bool FeH1LagrangeHex1::CoordIsInsideElem( const Vector<Double>& point,
                                            Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[1];
    const Double & zeta = point[2];
    return (  xi >= (-1.0 - tolerance)) &&
           ( eta >= (-1.0 - tolerance)) &&
           (zeta >= (-1.0 - tolerance)) &&
           (  xi <= (1.0 + tolerance)) &&
           ( eta <= (1.0 + tolerance)) &&
           (zeta <= (1.0 + tolerance));  
  }
  
  void FeH1LagrangeHex1::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                            const StdVector<UInt> & volConnect,
                            const LocPoint & surfIntPoint,
                            LocPoint & volIntPoint,
                            Vector<Double>& locNormal ) {

    // Try to find out, which vertices are in common with
    // the surface element. Then calculate the product of all four
    // and compare them
    //
    //                    zeta 
    //                     ^ eta 
    //    8 +-------+ 7    |/
    //     /|      /|      0--> xi 
   //    / |     / |
    // 5 +--+----+6 |   
    //   |  +-- -|- + 3    
    //   | / 4   | /    REFERENCE VOLUME ELEMENT
    //   |/      |/
    // 1 +-------+ 2



    StdVector<UInt> commonIndex(4);
    UInt found = 0;
    UInt indexProduct = 0;
    std::string errMsg;
  
    volIntPoint.coord.Resize(3);
    locNormal.Resize(3);
  
    // loop over surface connect
    for (UInt iSurf=0; iSurf<4; iSurf++)
      // loop over volume connect
      for (UInt iVol=0; iVol<8; iVol++)
        if (surfConnect[iSurf] == volConnect[iVol])
          {
            commonIndex[found++] = iVol+1;
          }

    // std::cerr << std::endl << std::endl;
    //std::cerr << "commonIndex = " << std::endl << commonIndex << std::endl << std::endl;
    indexProduct =  commonIndex[0] * commonIndex[1];
    indexProduct *= commonIndex[2] * commonIndex[3];

    //std::cerr << "indexProduct = " << indexProduct << std::endl;
    switch(indexProduct)
      {
      case 24:
        // Surface[1,2,3,4] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = surfIntPoint[1];
        volIntPoint[2] = -1.0;
        locNormal[0] =  0.0;
        locNormal[1] =  0.0;
        locNormal[2] = -1.0;
        break;

      case 1680:
        // Surface[5,6,7,8] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = surfIntPoint[1];
        volIntPoint[2] = 1.0;
        locNormal[0] =  0.0;
        locNormal[1] =  0.0;
        locNormal[2] =  1.0;
        break;

      case 252:
        // Surface[2,3,7,6] is common
        volIntPoint[0] = 1.0;
        volIntPoint[1] = surfIntPoint[0];
        volIntPoint[2] = surfIntPoint[1];
        locNormal[0] =  1.0;
        locNormal[1] =  0.0;
        locNormal[2] =  0.0;
        break;
      
      case 160:
        // Surface[1,5,8,4] is common
        volIntPoint[0] = -1.0;
        volIntPoint[1] = surfIntPoint[0];
        volIntPoint[2] = surfIntPoint[1];
        locNormal[0] = -1.0;
        locNormal[1] =  0.0;
        locNormal[2] =  0.0;
        break;
      
      case 672:
        // Surface[4,3,7,8] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = 1.0;
        volIntPoint[2] = surfIntPoint[1];
        locNormal[0] =  0.0;
        locNormal[1] =  1.0;
        locNormal[2] =  0.0;
        break;
      
      case 60:
        // Surface[1,2,6,5] is common
        volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = -1.0;
        volIntPoint[2] = surfIntPoint[1];
        locNormal[0] =  0.0;
        locNormal[1] = -1.0;
        locNormal[2] =  0.0;
        break;
      
      default:
        EXCEPTION("HexaFE::GetLocalIntPoints4Surface: surface and volume element "
                  << "have not four nodes in common. Check your .mesh-file.");
      }
   }
  
  
  // --- Wedge 1st order ---
  FeH1LagrangeWedge1::FeH1LagrangeWedge1() {
    feType_ = Elem::ET_WEDGE6;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 6;
    order_ = 1; 
  }
    
  FeH1LagrangeWedge1::~FeH1LagrangeWedge1() {
    
  }
  
  void FeH1LagrangeWedge1::CalcShFnc( Vector<Double>& shape,
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp ) {
    shape.Resize( 6 );
    //"Wedge Elements"
    // from "Dhatt, G.: The Finite Element Method Displayed, p. 120"
    // corner nodes

    shape[0] = 0.5 * (1 - point[2]) * (1 - point[0] - point[1]);
    shape[1] = 0.5 * (1 - point[2]) * point[0];
    shape[2] = 0.5 * (1 - point[2]) * point[1];
    shape[3] = 0.5 * (1 + point[2]) * (1 - point[0] - point[1]);
    shape[4] = 0.5 * (1 + point[2]) * point[0];
    shape[5] = 0.5 * (1 + point[2]) * point[1];
  }
  
  void FeH1LagrangeWedge1::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                            const Vector<Double>& point,
                                            const Elem* ptElem,
                                            UInt comp ) {
    deriv.Resize( 6, 3 );

    deriv[0][0] = -0.5 * (1 - point[2]);
    deriv[0][1] = -0.5 * (1 - point[2]);
    deriv[0][2] = -0.5 * (1 - point[0] - point[1]);

    deriv[1][0] =  0.5 * (1 - point[2]);
    deriv[1][1] =  0.0;
    deriv[1][2] = -0.5 * point[0];

    deriv[2][0] =  0.0;
    deriv[2][1] =  0.5 * (1 - point[2]);
    deriv[2][2] = -0.5 * point[1];

    deriv[3][0] = -0.5 * (1+ point[2]);
    deriv[3][1] = -0.5 * (1+ point[2]);
    deriv[3][2] =  0.5 * (1- point[0] - point[1]);

    deriv[4][0] =  0.5 * (1+ point[2]);
    deriv[4][1] =  0.0;
    deriv[4][2] =  0.5 * point[0];

    deriv[5][0] =  0.0;
    deriv[5][1] =  0.5 * (1 + point[2]);
    deriv[5][2] =  0.5 * point[1];
  }

  bool FeH1LagrangeWedge1::CoordIsInsideElem( const Vector<Double>& point,
                                              Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[1];
    const Double & zeta = point[2];
    bool isInside = 
        (        xi >= ( 0 - tolerance)) &&
        (       eta >= ( 0 - tolerance)) &&
        ((xi + eta) <= ( 1 + tolerance)) &&
        (      zeta >= (-1 - tolerance)) &&
        (      zeta <= ( 1 + tolerance));
    return isInside;
  }
  
  void FeH1LagrangeWedge1::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                            const StdVector<UInt> & volConnect,
                            const LocPoint & surfIntPoint,
                            LocPoint & volIntPoint,
                            Vector<Double>& locNormal ) {

    // Try to find out, which vertices are in common with
    // the surface element. Then calculate the product of all four
    // and compare them
    //      + 6    
    //     /|\   
    //    / |  \           zeta
    // 4 +----- + 5         ^  eta
    //   |  + 3 |           |/ 
    //   | / \  |           0--> xi
    //   |/    \|   
    // 1 +------+ 2

    // Check if surface element is triangle 
    // or quadrilateral
    if (surfConnect.GetSize() == 3 ||
        surfConnect.GetSize() == 6) 
    {
      // ---- Triangle Surface ---
      StdVector<Integer> commonIndex(3);
      Integer found = 0;
      Integer indexSum = 0;
      std::string errMsg;

      volIntPoint.coord.Resize(3);
      locNormal.Resize(3);

      // loop over surface connect
      for (Integer iSurf=0; iSurf<3; iSurf++)
        // loop over volume connect
        for (Integer iVol=0; iVol<6; iVol++)
          if (surfConnect[iSurf] == volConnect[iVol])
          {
            commonIndex[found++] = iVol+1;
          }


      indexSum =  commonIndex[0] + commonIndex[1] + commonIndex[2];

      switch(indexSum)
      {
        case 6:
          // Surface[1,2,3] is common
          volIntPoint[0] = surfIntPoint[0];
          volIntPoint[1] = surfIntPoint[1];
          volIntPoint[2] = -1.0;
          locNormal[0] =  0.0;
          locNormal[1] =  0.0;
          locNormal[2] = -1.0;
          break;

        case 15:
          // Surface[4,5,6] is common
          volIntPoint[0] = surfIntPoint[0];
          volIntPoint[1] = surfIntPoint[1];
          volIntPoint[2] = 1.0;
          locNormal[0] =  0.0;
          locNormal[1] =  0.0;
          locNormal[2] =  1.0;
          break;

        default:
          EXCEPTION( "WedgeFE::GetLocalIntPoints4Surface: surface and volume element "
              << "have not three nodes in common. Check your .mesh-file.");
      }

    } else {
      // ---- Quadrilateral Surface ---
      StdVector<Integer> commonIndex(4);
      Integer found = 0;
      Integer indexSum = 0;
      std::string errMsg;

      volIntPoint.coord.Resize(3);
      locNormal.Resize(3);

      // loop over surface connect
      for (Integer iSurf=0; iSurf<4; iSurf++)
        // loop over volume connect
        for (Integer iVol=0; iVol<6; iVol++)
          if (surfConnect[iSurf] == volConnect[iVol])
          {
            commonIndex[found++] = iVol+1;
          }


      indexSum = 
          commonIndex[0] + commonIndex[1] 
          + commonIndex[2] + commonIndex[3];

      // NOTE: Since the line quad-element is defined in the range [-1;+1]
      // we have to calculate (1+surfCoord)/2 in order to get the right
      // position on the wedge element

      switch(indexSum)
      {
        case 16:
          // Surface[2,3,5,6] is common
          volIntPoint[0] = 0.5 - (surfIntPoint[0] / 2.0);
          volIntPoint[1] = 0.5 + (surfIntPoint[0] / 2.0);
          volIntPoint[2] = surfIntPoint[1];
          locNormal[0] = sqrt(.5);
          locNormal[1] = sqrt(.5);
          locNormal[2] = 0.0;
          break;

        case 14:
          // Surface[1,3,4,6] is common
          volIntPoint[0] = 0.0;
          volIntPoint[1] = 0.5 * (1 + surfIntPoint[0]);
          volIntPoint[2] = 0.5 * (1 + surfIntPoint[1]);
          locNormal[0] = -1.0;
          locNormal[1] =  0.0;
          locNormal[2] =  0.0;
          break;

        case 12:
          // Surface[1,2,4,5] is common
          volIntPoint[0] = 0.5 * (1 + surfIntPoint[1]);
          volIntPoint[1] = 0.0;
          volIntPoint[2] = 0.5 * (1 + surfIntPoint[0]);
          locNormal[0] =  0.0;
          locNormal[1] = -1.0;
          locNormal[2] =  0.0;
          break;

        default:
          EXCEPTION("WedgeFE::GetLocalIntPoints4Surface: surface and volume element "
              << "have not four nodes in common. Check your .mesh-file.");
      } // switch
    } // if
  }
  
  // ========================================================================
  //  Lagrangian Elements of 2nd order
  // ========================================================================
  
  // --- Line 2nd order ---
  
  FeH1LagrangeLine2::FeH1LagrangeLine2() {
    feType_ = Elem::ET_LINE3;;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 3;
    order_ = 2; 
    elemDim_ = 1;
  }
  FeH1LagrangeLine2::~FeH1LagrangeLine2() {
    
  }
  
  void FeH1LagrangeLine2::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
    //WARN("FeH1LagrangeLine2::CalcShFnc: Implement me");
  //  shape.Resize( 2 );
  //   shape[0] = 0.5 * ( 1.0 - point[0] );
  //   shape[1] = 0.5 * ( 1.0 + point[0] );
  }
  
  void FeH1LagrangeLine2::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
    //WARN("FeH1LagrangeLine2::CalcLocDerivShFnc: Implement me");
  //    deriv.Resize(2, 1);
  //    deriv[0][0] = 0.5 * -1.0;
  //    deriv[1][0] = 0.5 *  1.0;
  }
  
  bool FeH1LagrangeLine2::CoordIsInsideElem( const Vector<Double>& point,
                                             Double tolerance )  {
    const Double & xi = point[0];
    return (xi >= (-1.0 - tolerance) &&
            xi <= (1.0 + tolerance));
  }
  
  void FeH1LagrangeLine2::
   GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                                     const StdVector<UInt> & volConnect,
                                                     const LocPoint & surfIntPoint,
                                                     LocPoint & volIntPoint,
                                                     Vector<Double>& locNormal ) {
     EXCEPTION("Not implemented");
   }
  
  // --- Quad 2nd order ---
   
  FeH1LagrangeQuad2::FeH1LagrangeQuad2() {
    feType_ = Elem::ET_QUAD8;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 8;
    order_ = 2;
  }
    
  FeH1LagrangeQuad2::~FeH1LagrangeQuad2() {
    
  }
  
  void FeH1LagrangeQuad2::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    shape.Resize( 8 );
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 122.
    
    // corner nodesf
    for( UInt i = 0; i < 4; i++ ) {
      shape[i] = 0.25 * ( 1 + coords[i][0] * point[0] ) 
                      * ( 1 + coords[i][1] * point[1] )
                      * ( coords[i][0] * point[0] + coords[i][1] * point[1] - 1 );
    }
    
    // mid-side nodes
    for( UInt i = 4; i < 8; i = i + 2 ) {
      shape[i]   = 0.5 * ( 1 - point[0] * point[0] )
                       * ( 1 + coords[i][1] * point[1] );
      shape[i+1] = 0.5 * (1 - point[1] * point[1] )
                       * (1 + coords[i+1][0] * point[0] );
    }
  }
  
  void FeH1LagrangeQuad2::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    deriv.Resize( 8, 2 );
    
    // corner nodes
    for( UInt i = 0; i < 4; i++ ) {
      deriv[i][0] = 0.25 * coords[i][0]  
                         * ( 1 + coords[i][1] * point[1] )
                         * ( 2 * coords[i][0] * point[0] 
                               + coords[i][1] * point[1] );
      deriv[i][1] = 0.25 * coords[i][1]
                         * ( 1 + coords[i][0] * point[0] )
                         * ( 2 * coords[i][1] * point[1] 
                               + coords[i][0] * point[0] );
    }
      
    // mid-side nodes
    for( UInt i = 4; i < 8; i = i + 2 ) {
      deriv[i][0] = - point[0] * ( 1 + coords[i][1] * point[1] );
      deriv[i][1] =  0.5 * coords[i][1] * ( 1 - point[0] * point[0] );
      
      deriv[i+1][0] = 0.5 * coords[i+1][0] * ( 1 - point[1] * point[1] );
      deriv[i+1][1] = - point[1] * ( 1 + coords[i+1][0] * point[0] );
    }
  }
  
  bool FeH1LagrangeQuad2::CoordIsInsideElem( const Vector<Double>& point,
                                             Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[0];
    return  ( xi >= (-1.0 - tolerance)) &&
            (eta >= (-1.0 - tolerance)) &&
            ( xi <= (1.0 + tolerance)) &&
            (eta <= (1.0 + tolerance));  
    }
  
  void FeH1LagrangeQuad2::
   GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                                                     const StdVector<UInt> & volConnect,
                                                     const LocPoint & surfIntPoint,
                                                     LocPoint & volIntPoint,
                                                     Vector<Double>& locNormal ) {
     EXCEPTION("Not implemented");
   }
  // --- Hex 2nd order ---
  FeH1LagrangeHex2::FeH1LagrangeHex2() {
    feType_ = Elem::ET_HEXA20;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 20;
    order_ = 2;
    elemDim_ = 3;
  }
    
  FeH1LagrangeHex2::~FeH1LagrangeHex2() {
    

  }
  
  void FeH1LagrangeHex2::CalcShFnc( Vector<Double>& shape,
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp ) {
    //WARN("FeH1LagrangeHex2::CalcShFnc: Implement me");
    //  shape.Resize( 8 );
  //  shape[0] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]); 
  //  shape[1] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 - point[2]);
  //  shape[2] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
  //  shape[3] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 - point[2]);
  //  shape[4] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]); 
  //  shape[5] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] ) * (1.0 + point[2]);
  //  shape[6] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
  //  shape[7] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] ) * (1.0 + point[2]);
    
  }
  
  void FeH1LagrangeHex2::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                            const Vector<Double>& point,
                                            const Elem* ptElem,
                                            UInt comp ) {
    //WARN("FeH1LagrangeHex2::CalcLocDerivShFnc: Implement me");
    //  deriv.Resize( 8, 3 );
  //  StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
  //  for( UInt i = 0; i < 8; i++ ) {
  //    deriv[i][0] = 0.125  * coords[i][0] 
  //                         * (1 + coords[i][1] * point[1] ) 
  //                         * (1 + coords[i][2] * point[2] );
  //    
  //    deriv[i][1] = 0.125  * (1 + coords[i][0] * point[0] )
  //                         * coords[i][1] 
  //                         * (1 + coords[i][2] * point[2] );
  //    
  //    deriv[i][2] = 0.125  * (1 + coords[i][0] * point[0] ) 
  //                         * (1 + coords[i][1] * point[1] )
  //                         * coords[i][2];
  //  }
    
  }

  bool FeH1LagrangeHex2::CoordIsInsideElem( const Vector<Double>& point,
                                            Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[1];
    const Double & zeta = point[2];
    return (  xi >= (-1.0 - tolerance)) &&
        ( eta >= (-1.0 - tolerance)) &&
        (zeta >= (-1.0 - tolerance)) &&
        (  xi <= (1.0 + tolerance)) &&
        ( eta <= (1.0 + tolerance)) &&
        (zeta <= (1.0 + tolerance));  
  }

  void FeH1LagrangeHex2::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                            const StdVector<UInt> & volConnect,
                            const LocPoint & surfIntPoint,
                            LocPoint & volIntPoint,
                            Vector<Double>& locNormal ) {
    EXCEPTION("Not implemented");
  }

} // namespace CoupledField
