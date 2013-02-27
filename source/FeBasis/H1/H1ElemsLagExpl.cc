#include "H1ElemsLagExpl.hh"

namespace CoupledField {

  // ========================================================================
  //  FeH1LagrangeExpl 
  // ========================================================================
  
  FeH1LagrangeExpl::FeH1LagrangeExpl()
  : FeH1(), FeNodal() {
   order_ = 0;
   preComputShFnc_ = true;
  }
    
  FeH1LagrangeExpl::~FeH1LagrangeExpl() {
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


  void FeH1LagrangeExpl::SetFunctionsAtIp(const StdVector<LocPoint>& iPoints){
    
    //! Precompute shape functions only, if we are allowed to, i.e.
    //! if we have no higher order shape functions, depending on 
    //! the global orientation of the element
    if( !preComputShFnc_ ) 
      return;
    //shapeFncsAtIp_.resize(iPoints.GetSize());
    //shapeFncDerivsAtIp_.resize(iPoints.GetSize());
    for(UInt aPoint = 0; aPoint < iPoints.GetSize();aPoint++){
     const LocPoint& lp = iPoints[aPoint];
      CalcShFnc( shapeFncsAtIp_[lp.number], lp.coord, NULL, 1);
      CalcLocDerivShFnc( shapeFncDerivsAtIp_[lp.number], lp.coord,
                         NULL, 1);
    }
  }
  
  void FeH1LagrangeExpl::SetFunctionsAtIp(const std::map<Integer,LocPoint>& iPoints){


    // can only be performed for non-hierarchical elements
    //shapeFncsAtIp_.resize(iPoints.GetSize());
    //shapeFncDerivsAtIp_.resize(iPoints.GetSize());
    std::map<Integer,LocPoint>::const_iterator pIt = iPoints.begin();
    while(pIt != iPoints.end()){
     const LocPoint& lp = pIt->second;
      CalcShFnc( shapeFncsAtIp_[lp.number], lp.coord, NULL, 1);
      CalcLocDerivShFnc( shapeFncDerivsAtIp_[lp.number], lp.coord,
                         NULL, 1);
      pIt++;
    }
  }
  
  
  // ========================================================================
  //  Lagrangian Elements of 1st order
  // ========================================================================
  
  // --- Line 1st order ---
  
  FeH1LagrangeLine1::FeH1LagrangeLine1() : FeH1LagrangeLine() {
    feType_ = Elem::ET_LINE2;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 2;
    order_ = 1;
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
  
  
  bool FeH1LagrangeLine::CoordIsInsideElem( const Vector<Double>& point,
                                            Double tolerance )  {
   const Double & xi = point[0];
   return (xi >= (-1.0 - tolerance) &&
           xi <= (1.0 + tolerance));
  }
  
  void FeH1LagrangeLine::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                            const StdVector<UInt> & volConnect,
                            const LocPoint & surfIntPoint,
                            LocPoint & volIntPoint,
                            Vector<Double>& locNormal ) {
    EXCEPTION("Not implemented");
  }
  
  // --- Tria 1st order ---

  FeH1LagrangeTria1::FeH1LagrangeTria1()  : FeH1LagrangeTria(){
    feType_ = Elem::ET_TRIA3;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 3;
    order_ = 1; 
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


  bool FeH1LagrangeTria::CoordIsInsideElem( const Vector<Double>& point,
                                            Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[1];
    return ( xi >= (0 - tolerance)) &&
           ( eta >= (0 - tolerance)) &&
           ((xi + eta) <= (1 + tolerance));
  }

  void FeH1LagrangeTria::
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
   
  FeH1LagrangeQuad1::FeH1LagrangeQuad1() : FeH1LagrangeQuad() {
    feType_ = Elem::ET_QUAD4;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 4;
    order_ = 1; 
    // This element supports incompatible modes
    hasICModes_ = true;
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
    StdVector<Vector<Double> >& coords = shape_.nodeCoords;
    deriv.Resize( 4, 2 );
    for( UInt i = 0; i < 4; i++ ) {
      deriv[i][0] = 0.25 * coords[i][0] * (1 + coords[i][1] * point[1] );
      deriv[i][1] = 0.25 * (1 + coords[i][0] * point[0] ) * coords[i][1];
    }
  }
  
  
   void FeH1LagrangeQuad1::CalcShFncICModes( Vector<Double>& shape,
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp  ) {
     shape.Resize( 2 );
     shape[0] = 1.0 - point[0] * point[0];
     shape[1] = 1.0 - point[1] * point[1];
   }

   void FeH1LagrangeQuad1::CalcLocDerivShFncICModes( Matrix<Double> & deriv, 
                                                     const Vector<Double>& point,
                                                     const Elem* ptElem,
                                                     UInt comp ) {
     deriv.Resize( 2, 2);
     deriv.Init();
     deriv[0][0] = -2.0 * point[0];
     deriv[1][1] = -2.0 * point[1];
   }



  bool FeH1LagrangeQuad::CoordIsInsideElem( const Vector<Double>& point,
                                            Double tolerance )  {
    const Double & xi = point[0];
    const Double & eta = point[1];
    return  ( xi >= (-1.0 - tolerance)) &&
            (eta >= (-1.0 - tolerance)) &&
            ( xi <= (1.0 + tolerance)) &&
            (eta <= (1.0 + tolerance));  
  }
  
  void FeH1LagrangeQuad::
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
          if(commonIndex[0] == 2){
            volIntPoint[0] = surfIntPoint[0]*-1.0;
            locNormal[0] =  0.0;
            locNormal[1] = 1.0;
          }else{
            volIntPoint[0] = surfIntPoint[0];
            locNormal[0] =  0.0;
            locNormal[1] = -1.0;
          }
          volIntPoint[1] = -1.0;

          break;

        case 12:
          // Edge[4,3] is common
          if(commonIndex[0] == 3){
            volIntPoint[0] = surfIntPoint[0]*-1.0;
            locNormal[0] =  0.0;
            locNormal[1] = -1.0;
          }else{
            volIntPoint[0] = surfIntPoint[0];
            locNormal[0] =  0.0;
            locNormal[1] =  1.0;
          }
          volIntPoint[1] = 1.0;
          break;

        case 4:
          // Edge[1,4] is common
          if(commonIndex[0] == 4){
            volIntPoint[1] = surfIntPoint[0]*-1.0;
            locNormal[0] =  1.0;
            locNormal[1] =  0.0;
          }else{
            volIntPoint[1] = surfIntPoint[0];
            locNormal[0] = -1.0;
            locNormal[1] =  0.0;
          }
          volIntPoint[0] = -1.0;
          break;

        case 6:
          // Edge[2,3] is common
          if(commonIndex[0] == 3){
            volIntPoint[1] = surfIntPoint[0]*-1.0;
            locNormal[0] = -1.0;
            locNormal[1] =  0.0;
          }else{
            volIntPoint[1] = surfIntPoint[0];
            locNormal[0] =  1.0;
            locNormal[1] =  0.0;
          }
          volIntPoint[0] = 1.0;
          break;

        default:
          EXCEPTION( "RectangleFE::GetLocalIntPoints4Surface: surface and volume element "
                     <<  "have not two nodes in common. Check your .mesh-file.");
        }
  }
  
  // --- Hex 1st order ---
  FeH1LagrangeHex1::FeH1LagrangeHex1() : FeH1LagrangeHex() {
    feType_ = Elem::ET_HEXA8;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 8;
    order_ = 1; 
    
    // This element supports incompatible modes
    hasICModes_ = true;
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
    StdVector<Vector<Double> >& coords = shape_.nodeCoords;
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
  
  void FeH1LagrangeHex1::CalcShFncICModes( Vector<Double>& shape,
                                           const Vector<Double>& point,
                                           const Elem* ptElem,
                                           UInt comp  ) {
    shape.Resize( 3 );
    shape[0] = 1.0 - point[0] * point[0];
    shape[1] = 1.0 - point[1] * point[1];
    shape[2] = 1.0 - point[2] * point[2];

  }

  void FeH1LagrangeHex1::CalcLocDerivShFncICModes( Matrix<Double> & deriv, 
                                                   const Vector<Double>& point,
                                                   const Elem* ptElem,
                                                   UInt comp ) {
    deriv.Resize( 3, 3);
    deriv.Init();
    deriv[0][0] = -2.0 * point[0];
    deriv[1][1] = -2.0 * point[1];
    deriv[2][2] = -2.0 * point[2];
  }

  bool FeH1LagrangeHex::CoordIsInsideElem( const Vector<Double>& point,
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
  
  void FeH1LagrangeHex::
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
  FeH1LagrangeWedge1::FeH1LagrangeWedge1() : FeH1LagrangeWedge() {
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

  bool FeH1LagrangeWedge::CoordIsInsideElem( const Vector<Double>& point,
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
  
  void FeH1LagrangeWedge::
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
  
  FeH1LagrangeLine2::FeH1LagrangeLine2() : FeH1LagrangeLine() {
    feType_ = Elem::ET_LINE3;;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 3;
    order_ = 2; 
  }
  FeH1LagrangeLine2::~FeH1LagrangeLine2() {
    
  }
  
  void FeH1LagrangeLine2::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
    shape.Resize(actNumFncs_);
    shape[0] = 0.5*point[0]*(point[0]-1);
    shape[2] = 1.0 - point[0]*point[0];
    shape[1] = 0.5*point[0]*(point[0]+1);
  }
  
  void FeH1LagrangeLine2::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
    deriv.Resize(actNumFncs_,1);

    deriv[0][0] = 0.5*(2*point[0] - 1);
    deriv[2][0] = -2.0*point[0];
    deriv[1][0] = 0.5*(2*point[0] + 1);

  }
  
  // --- Tria 2nd order ---

  FeH1LagrangeTria2::FeH1LagrangeTria2()  : FeH1LagrangeTria() {
    feType_ = Elem::ET_TRIA6;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 6;
    order_ = 2;
  }

  FeH1LagrangeTria2::~FeH1LagrangeTria2() {

  }

  void FeH1LagrangeTria2::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
    // From Zienkiewicz, The Finite Element Method. Vol 1, page 128.
    // corner nodes
    shape.Resize(6);

    // Define the third component of the triangular coord.
    Double t = 1.0 - point[0] - point[1];

    shape[0] = t * (2*t - 1);
    shape[1] = point[0]*(2*point[0] - 1);
    shape[2] = point[1]*(2*point[1] - 1);
    shape[3] = 4 * point[0] * t;
    shape[4] = 4 * point[0] * point[1];
    shape[5] = 4 * point[1] * t;
  }

  void FeH1LagrangeTria2::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {

    deriv.Resize( 6, 2 );

    deriv[0][0] =  4*point[0] + 4*point[1] - 3.0;
    deriv[0][1] =  4*point[0] + 4*point[1] - 3.0;

    deriv[1][0] =  4*point[0]-1;
    deriv[1][1] =  0;

    deriv[2][0] =  0;
    deriv[2][1] =  4*point[1]-1;

    deriv[3][0] =  4*(1 - 2*point[0] - point[1]);
    deriv[3][1] = -4*point[0];

    deriv[4][0] =  4* point[1];
    deriv[4][1] =  4* point[0];

    deriv[5][0] = -4*point[1];
    deriv[5][1] =  4*(1 - 2*point[1] - point[0]);
  }

  // --- Quad 2nd order ---
  FeH1LagrangeQuad2::FeH1LagrangeQuad2() : FeH1LagrangeQuad() {
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
    StdVector<Vector<Double> >& coords = shape_.nodeCoords;
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
    StdVector<Vector<Double> >& coords = shape_.nodeCoords;
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
  
  // --- Quad 2nd order tensor product ---
  FeH1LagrangeQuad9::FeH1LagrangeQuad9() : FeH1LagrangeQuad(){
    feType_ = Elem::ET_QUAD9;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 9;
    order_ = 2;
  }

  FeH1LagrangeQuad9::~FeH1LagrangeQuad9() {

  }

  void FeH1LagrangeQuad9::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
    static Double shapeXi[3], shapeEta[3];
    shape.Resize(9);

    shapeXi[0] = 0.5*point[0]*(point[0]-1);
    shapeXi[2] = 1.0 - point[0]*point[0];
    shapeXi[1] = 0.5*point[0]*(point[0]+1);

    shapeEta[0] = 0.5*point[1]*(point[1]-1);
    shapeEta[2] = 1.0 - point[1]*point[1];
    shapeEta[1] = 0.5*point[1]*(point[1]+1);

    shape[0]= shapeXi[0]     * shapeEta[0];
    shape[1]= shapeXi[1]     * shapeEta[0];
    shape[2]= shapeXi[1]     * shapeEta[1];
    shape[3]= shapeXi[0]     * shapeEta[1];

    shape[4]= shapeXi[2]     * shapeEta[0];
    shape[5]= shapeXi[1]     * shapeEta[2];
    shape[6]= shapeXi[2]     * shapeEta[1];
    shape[7]= shapeXi[0]     * shapeEta[2];

    shape[8]= shapeXi[2]     * shapeEta[2];
  }

  void FeH1LagrangeQuad9::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
	    static double shapeXi[3], shapeEta[3];
	    static double shapeDerivXi[3], shapeDerivEta[3];
	    deriv.Resize(9,2);

	    shapeXi[0] = 0.5*point[0]*(point[0]-1);
	    shapeXi[2] = 1.0 - point[0]*point[0];
	    shapeXi[1] = 0.5*point[0]*(point[0]+1);

	    shapeEta[0] = 0.5*point[1]*(point[1]-1);
	    shapeEta[2] = 1.0 - point[1]*point[1];
	    shapeEta[1] = 0.5*point[1]*(point[1]+1);

	    shapeDerivXi[0] = 0.5*(2*point[0] - 1);
	    shapeDerivXi[2] = -2.0*point[0];
	    shapeDerivXi[1] = 0.5*(2*point[0] + 1);

	    shapeDerivEta[0] = 0.5*(2*point[1] - 1);
	    shapeDerivEta[2] = -2.0*point[1];
	    shapeDerivEta[1] = 0.5*(2*point[1] + 1);

	    deriv[0][0]= shapeDerivXi[0]     * shapeEta[0];
	    deriv[0][1]= shapeXi[0]     * shapeDerivEta[0];

	    deriv[1][0]= shapeDerivXi[1]     * shapeEta[0];
	    deriv[1][1]= shapeXi[1]     * shapeDerivEta[0];

	    deriv[2][0]= shapeDerivXi[1]     * shapeEta[1];
	    deriv[2][1]= shapeXi[1]     * shapeDerivEta[1];

	    deriv[3][0]= shapeDerivXi[0]     * shapeEta[1];
	    deriv[3][1]= shapeXi[0]     * shapeDerivEta[1];


	    deriv[4][0]= shapeDerivXi[2]     * shapeEta[0];
	    deriv[4][1]= shapeXi[2]     * shapeDerivEta[0];

	    deriv[5][0]= shapeDerivXi[1]     * shapeEta[2];
	    deriv[5][1]= shapeXi[1]     * shapeDerivEta[2];

	    deriv[6][0]= shapeDerivXi[2]     * shapeEta[1];
	    deriv[6][1]= shapeXi[2]     * shapeDerivEta[1];

	    deriv[7][0]= shapeDerivXi[0]     * shapeEta[2];
	    deriv[7][1]= shapeXi[0]     * shapeDerivEta[2];


	    deriv[8][0]= shapeDerivXi[2]     * shapeEta[2];
	    deriv[8][1]= shapeXi[2]     * shapeDerivEta[2];  }

  // --- Hex 2nd order ---
  FeH1LagrangeHex2::FeH1LagrangeHex2() : FeH1LagrangeHex() {
    feType_ = Elem::ET_HEXA20;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 20;
    order_ = 2;
  }
    
  FeH1LagrangeHex2::~FeH1LagrangeHex2() {
  }
  
  void FeH1LagrangeHex2::CalcShFnc( Vector<Double>& shape,
                                    const Vector<Double>& point,
                                    const Elem* ptElem,
                                    UInt comp ) {
    Double  xi, eta, zeta;
    UInt i;
    shape.Resize(20);

    //integration points
    xi   = point[0];
    eta  = point[1];
    zeta = point[2];

    StdVector<Vector<Double> >& coords = shape_.nodeCoords;
    
    //Corner coordinates:
    // Ni
    for (i=0;i<8; i++) {
      shape[i] =  0.125
          *(1.0+xi  *coords[i][0])
          *(1.0+eta *coords[i][1])
          *(1.0+zeta*coords[i][2])
          *(xi*coords[i][0]+eta*coords[i][1]+zeta*coords[i][2]-2.0);
    }

    //Midside nodes xi_i=0:
    // Ni
    for (i=8;i<=14; i+=2) {
      shape[i]  = 0.25*(1.0-xi*xi)*
                       (1.0+eta*coords[i][1])*
                       (1.0+zeta*coords[i][2]);
    }

    //Midside nodes eta_i=0:
    // Ni
    for (i=9;i<=15; i+=2) {
      shape[i]  = 0.25*(1.0+xi*coords[i][0])*
                       (1.0-eta*eta)*
                       (1.0+zeta*coords[i][2]);
    }

    //Midside nodes zeta_i=0:
    // Ni
    for (i=16;i<20; i++) {
      shape[i] = 0.25*(1.0+xi*coords[i][0])*
                      (1.0+eta*coords[i][1])*(1.0-zeta*zeta);
    }
  }
  
  void FeH1LagrangeHex2::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                            const Vector<Double>& point,
                                            const Elem* ptElem,
                                            UInt comp ) {
	// shape function according to Zienkiewicz 2000 (Volume 1 - The Basis)
	// "The Finite Element Method" S. 185 Kap. 8.11
    Double  xi, eta, zeta;
    UInt i;
    deriv.Resize(20, 3);
    deriv.Init();
    xi   = point[0];
    eta  = point[1];
    zeta = point[2];

    StdVector<Vector<Double> >& coords = shape_.nodeCoords;
    for (i=0;i<8; i++) {
      //Corner coordinates: Ni,x
      deriv[i][0] =  0.25 * xi * coords[i][0] * coords[i][0]
        - 0.125 * coords[i][0]
        + 0.25 * zeta * xi * coords[i][2] * coords[i][0] * coords[i][0]
        + 0.125 * zeta * zeta * coords[i][0] * coords[i][2] * coords[i][2]
        + 0.25 * eta * xi * coords[i][1] * coords[i][0] * coords[i][0]
        + 0.125 * eta * eta * coords[i][0] * coords[i][1] * coords[i][1]
        + 0.25 * xi * eta * zeta * coords[i][1] * coords[i][2] * coords[i][0] * coords[i][0]
        + 0.125 * eta * eta * zeta * coords[i][0] * coords[i][1] * coords[i][1] * coords[i][2]
        + 0.125 * eta * zeta * zeta * coords[i][0] * coords[i][1] * coords[i][2] * coords[i][2]
        + 0.125 * eta * zeta * coords[i][0] * coords[i][1] * coords[i][2];
      
      //Corner coordinates: Ni,y
      deriv[i][1] =  0.25 * eta * coords[i][1] * coords[i][1]
        - 0.125 * coords[i][1]
        + 0.25 * eta * zeta * coords[i][1] * coords[i][1] * coords[i][2]
        + 0.125 * zeta * zeta * coords[i][1] * coords[i][2] * coords[i][2]
        + 0.125 * xi * xi * coords[i][0] * coords[i][0] * coords[i][1]
        + 0.25 * eta * xi * coords[i][0] * coords[i][1] * coords[i][1]
        + 0.125 * xi * xi * zeta * coords[i][0] * coords[i][0] * coords[i][1] * coords[i][2]
        + 0.25 * xi * eta * zeta * coords[i][0] * coords[i][1] * coords[i][1] * coords[i][2]
        + 0.125 * xi * zeta * zeta * coords[i][0] * coords[i][1] * coords[i][2] * coords[i][2]
        + 0.125 * xi * zeta * coords[i][0] * coords[i][1] * coords[i][2];

      //Corner coordinates: Ni,z
      deriv[i][2] =  0.25 * zeta * coords[i][2] * coords[i][2]
        - 0.125 * coords[i][2]
        + 0.125 * eta * eta * coords[i][1] * coords[i][1] * coords[i][2]
        + 0.25 * eta * zeta * coords[i][1] * coords[i][2] * coords[i][2]
        + 0.125 * xi * xi * coords[i][0] * coords[i][0] * coords[i][2]
        + 0.25 * xi * zeta * coords[i][0] * coords[i][2] * coords[i][2]
        + 0.125 * xi * xi * eta * coords[i][0] * coords[i][0] * coords[i][1] * coords[i][2]
        + 0.125 * xi * eta * eta * coords[i][0] * coords[i][1] * coords[i][1] * coords[i][2]
        + 0.25 * xi * eta * zeta * coords[i][0] * coords[i][1] * coords[i][2] * coords[i][2]
        + 0.125 * xi * eta * coords[i][0] * coords[i][1] * coords[i][2];
    }

    for (i=8;i<=14; i+=2) {
      //Midside nodes xi_i=0: Ni,x
      deriv[i][0]  = -0.5 * xi * (1.0+eta * coords[i][1]) * (1.0+zeta * coords[i][2]);
      //Midside nodes xi_i=0: Ni,y
      deriv[i][1]  = -0.25 * (xi * xi - 1.0) * coords[i][1] * (1.0+zeta * coords[i][2]);
      //Midside nodes xi_i=0: Ni,z
      deriv[i][2]  = -0.25 * (xi * xi - 1.0) * (1.0+eta * coords[i][1]) * coords[i][2];
    }

    for (i=9;i<=15; i+=2) {
      //Midside nodes eta_i=0: Ni,x
      deriv[i][0]  = -0.25 * coords[i][0]  * (eta * eta-1.0) * (1.0+zeta * coords[i][2]);
      //Midside nodes eta_i=0: Ni,y
      deriv[i][1]  = -0.5 * (1.0+xi * coords[i][0]) * eta * (1.0+zeta * coords[i][2]);
      //Midside nodes eta_i=0: Ni,z
      deriv[i][2]  = -0.25 * (1.0+xi * coords[i][0]) * (eta * eta-1.0) * coords[i][2];
    }
    for (i=16;i<20; i++) {
      //Midside nodes zeta_i=0: Ni,x
      deriv[i][0] = -0.25 * coords[i][0] * (1+eta * coords[i][1]) * (zeta * zeta-1.0);
      //Midside nodes zeta_i=0: Ni,y
      deriv[i][1] = -0.25 * (1.0+xi * coords[i][0]) * coords[i][1] * (zeta * zeta-1.0);
      //Midside nodes zeta_i=0:  Ni,z
      deriv[i][2] = -0.5 * (1.0+xi * coords[i][0]) * (1.0+eta * coords[i][1]) * zeta;
    }
  }

  // --- Hex 2nd order ---
  FeH1LagrangeHex27::FeH1LagrangeHex27() : FeH1LagrangeHex() {
    feType_ = Elem::ET_HEXA27;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 27;
    order_ = 2;
  }

  FeH1LagrangeHex27::~FeH1LagrangeHex27() {
  }

  void FeH1LagrangeHex27::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
    static double shapeXi[3], shapeEta[3], shapeZeta[3];
    shape.Resize(actNumFncs_);

    shapeXi[0] = 0.5*point[0]*(point[0]-1);
    shapeXi[2] = 1.0 - point[0]*point[0];
    shapeXi[1] = 0.5*point[0]*(point[0]+1);

    shapeEta[0] = 0.5*point[1]*(point[1]-1);
    shapeEta[2] = 1.0 - point[1]*point[1];
    shapeEta[1] = 0.5*point[1]*(point[1]+1);

    shapeZeta[0] = 0.5*point[2]*(point[2]-1);
    shapeZeta[2] = 1.0 - point[2]*point[2];
    shapeZeta[1] = 0.5*point[2]*(point[2]+1);

    // Corners
    shape[0]= shapeXi[0]     * shapeEta[0]      * shapeZeta[0];
    shape[1]= shapeXi[1]     * shapeEta[0]      * shapeZeta[0];
    shape[2]= shapeXi[1]     * shapeEta[1]      * shapeZeta[0];
    shape[3]= shapeXi[0]     * shapeEta[1]      * shapeZeta[0];

    shape[4]= shapeXi[0]     * shapeEta[0]      * shapeZeta[1];
    shape[5]= shapeXi[1]     * shapeEta[0]      * shapeZeta[1];
    shape[6]= shapeXi[1]     * shapeEta[1]      * shapeZeta[1];
    shape[7]= shapeXi[0]     * shapeEta[1]      * shapeZeta[1];

    // Edges
    shape[8]= shapeXi[2]     * shapeEta[0]      * shapeZeta[0];
    shape[9]= shapeXi[1]     * shapeEta[2]      * shapeZeta[0];
    shape[10]= shapeXi[2]     * shapeEta[1]      * shapeZeta[0];
    shape[11]= shapeXi[0]     * shapeEta[2]      * shapeZeta[0];

    shape[12]= shapeXi[2]     * shapeEta[0]      * shapeZeta[1];
    shape[13]= shapeXi[1]     * shapeEta[2]      * shapeZeta[1];
    shape[14]= shapeXi[2]     * shapeEta[1]      * shapeZeta[1];
    shape[15]= shapeXi[0]     * shapeEta[2]      * shapeZeta[1];

    shape[16]= shapeXi[0]     * shapeEta[0]      * shapeZeta[2];
    shape[17]= shapeXi[1]     * shapeEta[0]      * shapeZeta[2];
    shape[18]= shapeXi[1]     * shapeEta[1]      * shapeZeta[2];
    shape[19]= shapeXi[0]     * shapeEta[1]      * shapeZeta[2];

    // Faces
    shape[20]= shapeXi[2]     * shapeEta[0]      * shapeZeta[2];
    shape[21]= shapeXi[1]     * shapeEta[2]      * shapeZeta[2];
    shape[22]= shapeXi[2]     * shapeEta[1]      * shapeZeta[2];
    shape[23]= shapeXi[0]     * shapeEta[2]      * shapeZeta[2];

    // Bottom
    shape[24]= shapeXi[2]     * shapeEta[2]      * shapeZeta[0];

    // Top
    shape[25]= shapeXi[2]     * shapeEta[2]      * shapeZeta[1];

    // Center
    shape[26]= shapeXi[2]     * shapeEta[2]      * shapeZeta[2];
  }

  void FeH1LagrangeHex27::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
    static double shapeXi[3], shapeEta[3], shapeZeta[3];
    static double shapeDerivXi[3], shapeDerivEta[3], shapeDerivZeta[3];
    deriv.Resize(actNumFncs_,3);

    shapeXi[0] = 0.5*point[0]*(point[0]-1);
    shapeXi[2] = 1.0 - point[0]*point[0];
    shapeXi[1] = 0.5*point[0]*(point[0]+1);

    shapeEta[0] = 0.5*point[1]*(point[1]-1);
    shapeEta[2] = 1.0 - point[1]*point[1];
    shapeEta[1] = 0.5*point[1]*(point[1]+1);

    shapeZeta[0] = 0.5*point[2]*(point[2]-1);
    shapeZeta[2] = 1.0 - point[2]*point[2];
    shapeZeta[1] = 0.5*point[2]*(point[2]+1);

    shapeDerivXi[0] = 0.5*(2*point[0] - 1);
    shapeDerivXi[2] = -2.0*point[0];
    shapeDerivXi[1] = 0.5*(2*point[0] + 1);

    shapeDerivEta[0] = 0.5*(2*point[1] - 1);
    shapeDerivEta[2] = -2.0*point[1];
    shapeDerivEta[1] = 0.5*(2*point[1] + 1);

    shapeDerivZeta[0] = 0.5*(2*point[2] - 1);
    shapeDerivZeta[2] = -2.0*point[2];
    shapeDerivZeta[1] = 0.5*(2*point[2] + 1);

    // Corners bottom
    deriv[0][0]= shapeDerivXi[0]     * shapeEta[0]      * shapeZeta[0];
    deriv[0][1]= shapeXi[0]     * shapeDerivEta[0]      * shapeZeta[0];
    deriv[0][2]= shapeXi[0]     * shapeEta[0]      * shapeDerivZeta[0];

    deriv[1][0]= shapeDerivXi[1]     * shapeEta[0]      * shapeZeta[0];
    deriv[1][1]= shapeXi[1]     * shapeDerivEta[0]      * shapeZeta[0];
    deriv[1][2]= shapeXi[1]     * shapeEta[0]      * shapeDerivZeta[0];

    deriv[2][0]= shapeDerivXi[1]     * shapeEta[1]      * shapeZeta[0];
    deriv[2][1]= shapeXi[1]     * shapeDerivEta[1]      * shapeZeta[0];
    deriv[2][2]= shapeXi[1]     * shapeEta[1]      * shapeDerivZeta[0];

    deriv[3][0]= shapeDerivXi[0]     * shapeEta[1]      * shapeZeta[0];
    deriv[3][1]= shapeXi[0]     * shapeDerivEta[1]      * shapeZeta[0];
    deriv[3][2]= shapeXi[0]     * shapeEta[1]      * shapeDerivZeta[0];

    //Corners top
    deriv[4][0]= shapeDerivXi[0]     * shapeEta[0]      * shapeZeta[1];
    deriv[4][1]= shapeXi[0]     * shapeDerivEta[0]      * shapeZeta[1];
    deriv[4][2]= shapeXi[0]     * shapeEta[0]      * shapeDerivZeta[1];

    deriv[5][0]= shapeDerivXi[1]     * shapeEta[0]      * shapeZeta[1];
    deriv[5][1]= shapeXi[1]     * shapeDerivEta[0]      * shapeZeta[1];
    deriv[5][2]= shapeXi[1]     * shapeEta[0]      * shapeDerivZeta[1];

    deriv[6][0]= shapeDerivXi[1]     * shapeEta[1]      * shapeZeta[1];
    deriv[6][1]= shapeXi[1]     * shapeDerivEta[1]      * shapeZeta[1];
    deriv[6][2]= shapeXi[1]     * shapeEta[1]      * shapeDerivZeta[1];

    deriv[7][0]= shapeDerivXi[0]     * shapeEta[1]      * shapeZeta[1];
    deriv[7][1]= shapeXi[0]     * shapeDerivEta[1]      * shapeZeta[1];
    deriv[7][2]= shapeXi[0]     * shapeEta[1]      * shapeDerivZeta[1];

    // Edges
    deriv[8][0]= shapeDerivXi[2]     * shapeEta[0]      * shapeZeta[0];
    deriv[8][1]= shapeXi[2]     * shapeDerivEta[0]      * shapeZeta[0];
    deriv[8][2]= shapeXi[2]     * shapeEta[0]      * shapeDerivZeta[0];

    deriv[9][0]= shapeDerivXi[1]     * shapeEta[2]      * shapeZeta[0];
    deriv[9][1]= shapeXi[1]     * shapeDerivEta[2]      * shapeZeta[0];
    deriv[9][2]= shapeXi[1]     * shapeEta[2]      * shapeDerivZeta[0];

    deriv[10][0]= shapeDerivXi[2]     * shapeEta[1]      * shapeZeta[0];
    deriv[10][1]= shapeXi[2]     * shapeDerivEta[1]      * shapeZeta[0];
    deriv[10][2]= shapeXi[2]     * shapeEta[1]      * shapeDerivZeta[0];

    deriv[11][0]= shapeDerivXi[0]     * shapeEta[2]      * shapeZeta[0];
    deriv[11][1]= shapeXi[0]     * shapeDerivEta[2]      * shapeZeta[0];
    deriv[11][2]= shapeXi[0]     * shapeEta[2]      * shapeDerivZeta[0];


    deriv[12][0]= shapeDerivXi[2]     * shapeEta[0]      * shapeZeta[1];
    deriv[12][1]= shapeXi[2]     * shapeDerivEta[0]      * shapeZeta[1];
    deriv[12][2]= shapeXi[2]     * shapeEta[0]      * shapeDerivZeta[1];

    deriv[13][0]= shapeDerivXi[1]     * shapeEta[2]      * shapeZeta[1];
    deriv[13][1]= shapeXi[1]     * shapeDerivEta[2]      * shapeZeta[1];
    deriv[13][2]= shapeXi[1]     * shapeEta[2]      * shapeDerivZeta[1];

    deriv[14][0]= shapeDerivXi[2]     * shapeEta[1]      * shapeZeta[1];
    deriv[14][1]= shapeXi[2]     * shapeDerivEta[1]      * shapeZeta[1];
    deriv[14][2]= shapeXi[2]     * shapeEta[1]      * shapeDerivZeta[1];

    deriv[15][0]= shapeDerivXi[0]     * shapeEta[2]      * shapeZeta[1];
    deriv[15][1]= shapeXi[0]     * shapeDerivEta[2]      * shapeZeta[1];
    deriv[15][2]= shapeXi[0]     * shapeEta[2]      * shapeDerivZeta[1];


    deriv[16][0]= shapeDerivXi[0]     * shapeEta[0]      * shapeZeta[2];
    deriv[16][1]= shapeXi[0]     * shapeDerivEta[0]      * shapeZeta[2];
    deriv[16][2]= shapeXi[0]     * shapeEta[0]      * shapeDerivZeta[2];

    deriv[17][0]= shapeDerivXi[1]     * shapeEta[0]      * shapeZeta[2];
    deriv[17][1]= shapeXi[1]     * shapeDerivEta[0]      * shapeZeta[2];
    deriv[17][2]= shapeXi[1]     * shapeEta[0]      * shapeDerivZeta[2];

    deriv[18][0]= shapeDerivXi[1]     * shapeEta[1]      * shapeZeta[2];
    deriv[18][1]= shapeXi[1]     * shapeDerivEta[1]      * shapeZeta[2];
    deriv[18][2]= shapeXi[1]     * shapeEta[1]      * shapeDerivZeta[2];

    deriv[19][0]= shapeDerivXi[0]     * shapeEta[1]      * shapeZeta[2];
    deriv[19][1]= shapeXi[0]     * shapeDerivEta[1]      * shapeZeta[2];
    deriv[19][2]= shapeXi[0]     * shapeEta[1]      * shapeDerivZeta[2];

    // Faces
    deriv[20][0]= shapeDerivXi[2]     * shapeEta[0]      * shapeZeta[2];
    deriv[20][1]= shapeXi[2]     * shapeDerivEta[0]      * shapeZeta[2];
    deriv[20][2]= shapeXi[2]     * shapeEta[0]      * shapeDerivZeta[2];

    deriv[21][0]= shapeDerivXi[1]     * shapeEta[2]      * shapeZeta[2];
    deriv[21][1]= shapeXi[1]     * shapeDerivEta[2]      * shapeZeta[2];
    deriv[21][2]= shapeXi[1]     * shapeEta[2]      * shapeDerivZeta[2];

    deriv[22][0]= shapeDerivXi[2]     * shapeEta[1]      * shapeZeta[2];
    deriv[22][1]= shapeXi[2]     * shapeDerivEta[1]      * shapeZeta[2];
    deriv[22][2]= shapeXi[2]     * shapeEta[1]      * shapeDerivZeta[2];

    deriv[23][0]= shapeDerivXi[0]     * shapeEta[2]      * shapeZeta[2];
    deriv[23][1]= shapeXi[0]     * shapeDerivEta[2]      * shapeZeta[2];
    deriv[23][2]= shapeXi[0]     * shapeEta[2]      * shapeDerivZeta[2];

    // Bottom
    deriv[24][0]= shapeDerivXi[2]     * shapeEta[2]      * shapeZeta[0];
    deriv[24][1]= shapeXi[2]     * shapeDerivEta[2]      * shapeZeta[0];
    deriv[24][2]= shapeXi[2]     * shapeEta[2]      * shapeDerivZeta[0];

    // Top
    deriv[25][0]= shapeDerivXi[2]     * shapeEta[2]      * shapeZeta[1];
    deriv[25][1]= shapeXi[2]     * shapeDerivEta[2]      * shapeZeta[1];
    deriv[25][2]= shapeXi[2]     * shapeEta[2]      * shapeDerivZeta[1];

    // Center
    deriv[26][0]= shapeDerivXi[2]     * shapeEta[2]      * shapeZeta[2];
    deriv[26][1]= shapeXi[2]     * shapeDerivEta[2]      * shapeZeta[2];
    deriv[26][2]= shapeXi[2]     * shapeEta[2]      * shapeDerivZeta[2];
  }

  // --- Wedge 2nd order ---
   FeH1LagrangeWedge2::FeH1LagrangeWedge2() : FeH1LagrangeWedge() {
     feType_ = Elem::ET_WEDGE15;
     shape_ = Elem::shapes[feType_];
     actNumFncs_ = 15;
     order_ = 2; 
   }
     
   FeH1LagrangeWedge2::~FeH1LagrangeWedge2() {
     
   }
   
   void FeH1LagrangeWedge2::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
     
     shape.Resize( 15 );

     shape[0] =  0.5 * point[2] * (1 - point[2]) * (1 - point[0] - point[1])
                      * (2 * point[0] + 2*point[1] -1);
     shape[1] =  0.5 * point[2] * (1 - point[2]) *  point[0] * (1 - 2 * point[0]);
     shape[2] =  0.5 * point[2] * (1 - point[2]) *  point[1] * (1 - 2 * point[1]);

     shape[3] =- 0.5 * point[2] * (1 + point[2]) * (1 - point[0] - point[1])
                     * (2*point[0] + 2 * point[1] -1);
     shape[4] = -0.5*point[2] * (1 + point[2]) * point[0] * (1 - 2 * point[0]);
     shape[5] = -0.5*point[2] * (1 + point[2]) * point[1] * (1 - 2 * point[1]);

     shape[6] = -2 * point[2] * (1 - point[2]) * point[0] * (1 - point[0] - point[1]);
     shape[7] = -2 * point[2] * (1 - point[2]) * point[0] * point[1];
     shape[8] = -2 * point[2] * (1 - point[2]) * point[1] * (1 - point[0] - point[1]);

     shape[9]  = 2 * point[2] * (1 + point[2]) * point[0] * (1 - point[0] - point[1]);
     shape[10] = 2 * point[2] * (1 + point[2]) * point[0] * point[1];
     shape[11] = 2 * point[2] * (1 + point[2]) * point[1] * (1 - point[0] - point[1]);

     shape[12] = (1 - point[0] - point[1]) * (1 - point[2] * point[2]);
     shape[13] = point[0] * (1 - point[2] * point[2]);
     shape[14] = point[1] * (1 - point[2] * point[2]);
   }
   
   void FeH1LagrangeWedge2::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
     deriv.Resize(15,3);
     deriv.Init();

     deriv[0][0] = +0.5*point[2] * (1 - point[2]) * (3 - 4*point[0] - 4*point[1]);
     deriv[0][1] = +0.5*point[2] * (1 - point[2]) * (3 - 4*point[0] - 4*point[1]);
     deriv[0][2] = +0.5 * (1 - 2*point[2]) * (1 - point[0] - point[1])
               * (2*point[0] + 2*point[1] -1 );

     deriv[1][0] =  0.5*point[2] * (1 - point[2]) * (1 - 4*point[0]);
     deriv[1][1] =  0.0;
     deriv[1][2] =  0.5 * (1 - 2*point[2]) * point[0] * (1 - 2*point[0]);

     deriv[2][0] =  0.0;
     deriv[2][1] =  0.5*point[2] * (1 - point[2]) * (1 - 4*point[1]);
     deriv[2][2] =  0.5 * (1 - 2*point[2]) * point[1] * (1 - 2*point[1]);

     deriv[3][0] = -0.5*point[2] * (1 + point[2]) * (3 - 4*point[0] - 4*point[1]);
     deriv[3][1] = -0.5*point[2] * (1 + point[2]) * (3 - 4*point[0] - 4*point[1]);
     deriv[3][2] = -0.5 * (1 + 2*point[2]) * (1 - point[0] - point[1])
               * (2*point[0] + 2*point[1] -1 );

     deriv[4][0] = -0.5*point[2] * (1 + point[2]) * (1 - 4*point[0]);
     deriv[4][1] =  0.0;
     deriv[4][2] = -0.5 * (1 + 2*point[2]) * point[0] * (1 - 2*point[0]);

     deriv[5][0] =  0.0;
     deriv[5][1] = -0.5*point[2] * (1 + point[2]) * (1 - 4*point[1]);
     deriv[5][2] = -0.5 * (1 + 2*point[2]) * point[1] * (1 - 2*point[1]);

     deriv[6][0] =  -2*point[2] * (1 - point[2]) * (1 - 2*point[0] - point[1]);
     deriv[6][1] =   2*point[2] * (1 - point[2]) * point[0];
     deriv[6][2] =  -2*(1 - 2*point[2]) * point[0] * (1 - point[0] - point[1]);

     deriv[7][0] =  -2*point[2] * (1 - point[2]) * point[1];
     deriv[7][1] =  -2*point[2] * (1 - point[2]) * point[0];
     deriv[7][2] =  -2*(1 - 2*point[2]) * point[0] * point[1];

     deriv[8][0] =   2*point[2] * (1 - point[2]) * point[1];
     deriv[8][1] =  -2*point[2] * (1 - point[2]) * (1 - point[0] - 2*point[1]);
     deriv[8][2] =  -2*(1 - 2*point[2]) * point[1] * (1 - point[0] - point[1]);

     deriv[9][0] =   2*point[2] * (1 + point[2]) * (1 - 2*point[0] - point[1]);
     deriv[9][1] =  -2*point[2] * (1 + point[2]) * point[0];
     deriv[9][2] =   2*(1 + 2*point[2]) * point[0] * (1 - point[0] - point[1]);

     deriv[10][0] =  2*point[2] * (1 + point[2]) * point[1];
     deriv[10][1] =  2*point[2] * (1 + point[2]) * point[0];
     deriv[10][2] =  2*(1 + 2*point[2]) * point[0] * point[1];

     deriv[11][0] =  -2*point[2] * (1 + point[2]) * point[1];
     deriv[11][1] =   2*point[2] * (1 + point[2]) * (1 - point[0] - 2*point[1]);
     deriv[11][2] =   2*(1 + 2*point[2]) * point[1] * (1 - point[0] - point[1]);

     deriv[12][0] =  -(1 - point[2] * point[2]);
     deriv[12][1] =  -(1 - point[2] * point[2]);
     deriv[12][2] =  -2*point[2] * (1 - point[0] - point[1]);

     deriv[13][0] =  (1 - point[2] * point[2]);
     deriv[13][1] =  0;
     deriv[13][2] =  -2*point[2] * point[0];

     deriv[14][0] =  0;
     deriv[14][1] =  (1 - point[2] * point[2]);
     deriv[14][2] = -2*point[2] * point[1];

   }

   // --- Wedge 2nd order ---
    FeH1LagrangeWedge18::FeH1LagrangeWedge18() : FeH1LagrangeWedge() {
      feType_ = Elem::ET_WEDGE18;
      shape_ = Elem::shapes[feType_];
      actNumFncs_ = 18;
      order_ = 2;
    }

    FeH1LagrangeWedge18::~FeH1LagrangeWedge18() {

    }

    void FeH1LagrangeWedge18::CalcShFnc( Vector<Double>& shape,
                                         const Vector<Double>& point,
                                         const Elem* ptElem,
                                         UInt comp ) {
      static bool issueWarning = true;
      if(issueWarning) {
        WARN("CalcShFnc for ET_WEDGE18 implemented but not tested!");
        issueWarning = false;
      }

      shape.Resize( actNumFncs_ );
      Double x = point[0];
      Double y = point[1];
      Double z = point[2];

      // corners
      shape[0] =-0.25 * (x + y) * (x + y + 1) * z * (1 - z);
      shape[1] =-0.25 *  x      * (x + 1)     * z * (1 - z);
      shape[2] =-0.25 *      y  * (1 + y)     * z * (1 - z);
      shape[3] = 0.25 * (x + y) * (x + y + 1) * z * (1 + z);
      shape[4] = 0.25 *  x      * (x + 1)     * z * (1 + z);
      shape[5] = 0.25 *      y  * (1 + y)     * z * (1 + z);

      // midsides of quadratic triangles
      shape[6] =  (x + 1)*(x + y) *  0.5 * z * (1 - z);
      shape[7] = -(x + 1)*(y + 1) *  0.5 * z * (1 - z);
      shape[8] =  (y + 1)*(x + y) *  0.5 * z * (1 - z);
      shape[9] = -(x + 1)*(x + y) *  0.5 * z * (1 + z);
      shape[10]=  (x + 1)*(y + 1) *  0.5 * z * (1 + z);
      shape[11]= -(y + 1)*(x + y) *  0.5 * z * (1 + z);

      // midsides of edges between the two triangles
      shape[12] = 0.5 * (x + y) * (x + y + 1) * (1 + z)*(1 - z);
      shape[13] = 0.5 *  x      * (x + 1)     * (1 + z)*(1 - z);
      shape[14] = 0.5 *      y  * (1 + y)     * (1 + z)*(1 - z);

      //Centerpoints of the biquadratic quads
      shape[15] = -(x + 1)*(x + y) * (1 + z)*(1 - z);
      shape[16] =  (x + 1)*(y + 1) * (1 + z)*(1 - z);
      shape[17] = -(y + 1)*(x + y) * (1 + z)*(1 - z);
    }

    void FeH1LagrangeWedge18::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                                 const Vector<Double>& point,
                                                 const Elem* ptElem,
                                                 UInt comp ) {
      static bool issueWarning = true;
      if(issueWarning) {
        WARN("CalcLocDerivShFnc for ET_WEDGE18 implemented but not tested!");
        issueWarning = false;
      }
      deriv.Resize(actNumFncs_,3);
      deriv.Init();
      Double x = point[0];
      Double y = point[1];
      Double z = point[2];

      //Derivatives in x-direction
      // corners
      deriv[0][0] = -0.25 * (2 * x + 2 * y + 1) * z * (1 - z);
      deriv[1][0] = -0.25 * (2 * x + 1)         * z * (1 - z);
      deriv[2][0] =  0;
      deriv[3][0] =  0.25 * (2 * x + 2 * y + 1) * z * (1 + z);
      deriv[4][0] =  0.25 * (2 * x + 1)         * z * (1 + z);
      deriv[5][0] =  0;
      // midsides of quadratic triangles
      deriv[6][0] =  (2 * x + y + 1) *  0.5 * z * (1 - z);
      deriv[7][0] = -(y + 1)         *  0.5 * z * (1 - z);
      deriv[8][0] =  (y + 1)         *  0.5 * z * (1 - z);
      deriv[9][0] = -(2 * x + y + 1) *  0.5 * z * (1 + z);
      deriv[10][0] = (y + 1)         *  0.5 * z * (1 + z) ;
      deriv[11][0] =-(y + 1)         *  0.5 * z * (1 + z) ;
      // midsides of edges between the two triangles
      deriv[12][0] = 0.5 * (2 * x + 2 * y + 1) * (1 + z)*(1 - z);
      deriv[13][0] = 0.5 * (2 * x + 1)         * (1 + z)*(1 - z);
      deriv[14][0] = 0;
      //Centerpoints of the biquadratic quads
      deriv[15][0] = -(2 * x + y + 1) * (1 + z)*(1 - z);
      deriv[16][0] =  (y + 1)         * (1 + z)*(1 - z);
      deriv[17][0] = -(y + 1)         * (1 + z)*(1 - z);

      //Derivatives in y-direction
      // corners
      deriv[0][1] = -0.25 * (2 * y + 2 * x + 1)   * z * (1 - z);
      deriv[1][1] =  0;
      deriv[2][1] = -0.25 * (2 * y + 1)           * z * (1 - z);
      deriv[3][1] =  0.25 * (2 * y + 2 * x + 1)   * z * (1 + z);
      deriv[4][1] =  0;
      deriv[5][1] =  0.25 * (2 * y + 1)           * z * (1 + z);
      // midsides of quadratic triangles
      deriv[6][1] =  (x + 1)         *  0.5 * z * (1 - z);
      deriv[7][1] = -(x + 1)         *  0.5 * z * (1 - z);
      deriv[8][1] =  (2 * y + x + 1) *  0.5 * z * (1 - z);
      deriv[9][1] = -(x + 1)         *  0.5 * z * (1 + z);
      deriv[10][1] =  (x + 1)         *  0.5 * z * (1 + z);
      deriv[11][1] = -(2 * y + x + 1) *  0.5 * z * (1 + z);
      // midsides of edges between the two triangles
      deriv[12][1] = 0.5 * (2 * y + 2 * x + 1) * (1 + z)*(1 - z);
      deriv[13][1] = 0;
      deriv[14][1] = 0.5 * (2 * y + 1)         * (1 + z)*(1 - z);
      //Centerpoints of the biquadratic quads
      deriv[15][1] = -(x + 1)         * (1 + z)*(1 - z);
      deriv[16][1] =  (x + 1)         * (1 + z)*(1 - z);
      deriv[17][1] = -(2 * y + x + 1) * (1 + z)*(1 - z);

      //Derivatives in z-direction
      // corners
      deriv[0][2] = -0.25 * (x + y) * (x + y + 1) * (1 - 2 * z);
      deriv[1][2] = -0.25 *  x      * (x + 1)     * (1 - 2 * z);
      deriv[2][2] = -0.25 *      y  * (1 + y)     * (1 - 2 * z);
      deriv[3][2] =  0.25 * (x + y) * (x + y + 1) * (1 + 2 * z);
      deriv[4][2] =  0.25 *  x      * (x + 1)     * (1 + 2 * z);
      deriv[5][2] =  0.25 *      y  * (1 + y)     * (1 + 2 * z);
      // midsides of quadratic triangles
      deriv[6][2] =  (x + 1)*(x + y) *  0.5 * (1 - 2 * z);
      deriv[7][2] = -(x + 1)*(y + 1) *  0.5 * (1 - 2 * z);
      deriv[8][2] =  (y + 1)*(x + y) *  0.5 * (1 - 2 * z);
      deriv[9][2] = -(x + 1)*(x + y) *  0.5 * (1 + 2 * z);
      deriv[10][2] =  (x + 1)*(y + 1) *  0.5 * (1 + 2 * z);
      deriv[11][2] = -(y + 1)*(x + y) *  0.5 * (1 + 2 * z);
      // midsides of edges between the two triangles
      deriv[12][2] = 0.5 * (x + y) * (x + y + 1) * (-2 * z);
      deriv[13][2] = 0.5 *  x      * (x + 1)     * (-2 * z);
      deriv[14][2] = 0.5 *      y  * (1 + y)     * (-2 * z);
      //Centerpoints of the biquadratic quads
      deriv[15][2] = -(x + 1)*(x + y) * (-2 * z);
      deriv[16][2] =  (x + 1)*(y + 1) * (-2 * z);
      deriv[17][2] = -(y + 1)*(x + y) * (-2 * z);
    }

   // --- Tetra 1st order ---
    FeH1LagrangeTet1::FeH1LagrangeTet1() : FeH1LagrangeTet() {
      feType_ = Elem::ET_TET4;
      shape_ = Elem::shapes[feType_];
      actNumFncs_ = 4;
      order_ = 1;
    }

    FeH1LagrangeTet1::~FeH1LagrangeTet1() {

    }

    void FeH1LagrangeTet1::CalcShFnc( Vector<Double>& shape,
                                      const Vector<Double>& point,
                                      const Elem* ptElem,
                                      UInt comp ) {

    	shape.Resize(4);

        // see Hughes p. 170

        shape[0] = 1.0 - point[0] - point[1] - point[2];

        for( UInt i=1; i<4; i++ )
          shape[i] = point[i-1];
    }

    void FeH1LagrangeTet1::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                              const Vector<Double>& point,
                                              const Elem* ptElem,
                                              UInt comp ) {
      deriv.Resize(4,3);
      deriv.Init();

      for( UInt i=0; i<3; i++)
        deriv[0][i] = -1.0;

      for( UInt i=1; i < 4; i++)
        deriv[i][i-1] = 1.0;
    }

   // --- Tetra 2nd order ---
    FeH1LagrangeTet2::FeH1LagrangeTet2() : FeH1LagrangeTet() {
      feType_ = Elem::ET_TET10;
      shape_ = Elem::shapes[feType_];
      actNumFncs_ = 10;
      order_ = 2;
    }

    FeH1LagrangeTet2::~FeH1LagrangeTet2() {

    }

    void FeH1LagrangeTet2::CalcShFnc( Vector<Double>& shape,
                                      const Vector<Double>& point,
                                      const Elem* ptElem,
                                      UInt comp ) {

    	shape.Resize(10);

    	//See Finite Element Procedures :Klaus Juergen Bathe, Prentice Hall,
    	//page 375 Sec. 5.3.
    	//Definition of the shape functions from 5 to 10

    	shape[4]=4*point[0]*(1 - point[0] - point[1] - point[2]);
    	shape[5]=4*point[0]*point[1];
    	shape[6]=4*point[1]*(1 - point[0] - point[1] - point[2]);
    	shape[7]=4*point[2]*(1 - point[0] - point[1] - point[2]);
    	shape[8]=4*point[0]*point[2];
    	shape[9]=4*point[1]*point[2];


    	//definition of the shape functions from 1 to 4

    	shape[0]= 1. - point[0] - point[1] - point[2] - 0.5*shape[4] -
    			0.5*shape[6] - 0.5*shape[7];
    	shape[1]= point[0] - 0.5*shape[4] - 0.5*shape[5] - 0.5*shape[8];
    	shape[2]= point[1] - 0.5*shape[5] - 0.5*shape[6] - 0.5*shape[9];
    	shape[3]= point[2] - 0.5*shape[8] - 0.5*shape[7] - 0.5*shape[9];
    }

    void FeH1LagrangeTet2::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                              const Vector<Double>& point,
                                              const Elem* ptElem,
                                              UInt comp ) {
      deriv.Resize(10,3);
      deriv.Init();

      //Calculation of the local derivatives from 5 to 10

      deriv[4][0]=4. - 8*point[0] - 4.*point[1] - 4.*point[2];
      deriv[4][1]=-4.*point[0];
      deriv[4][2]=-4.*point[0];

      deriv[5][0]=4.*point[1];
      deriv[5][1]=4.*point[0];
      deriv[5][2]=0.0;

      deriv[6][0]=-4.*point[1];
      deriv[6][1]=4. - 8*point[1] - 4.*point[0] - 4.*point[2];
      deriv[6][2]=-4.*point[1];

      deriv[8][0]=4.*point[2];
      deriv[8][1]=0.0;
      deriv[8][2]=4.*point[0];

      deriv[9][0]=0.0;
      deriv[9][1]=4.*point[2];
      deriv[9][2]=4.*point[1];

      deriv[7][0]=-4.*point[2];
      deriv[7][1]=-4.*point[2];
      deriv[7][2]=4. - 8.*point[2] - 4.*point[0] - 4.*point[1];

      //Calculation of the local derivatives from 1 to 4

      deriv[0][0]= -1.-0.5*deriv[4][0]-0.5*deriv[6][0]-0.5*deriv[7][0];
      deriv[0][1]= -1.-0.5*deriv[4][1]-0.5*deriv[6][1]-0.5*deriv[7][1];
      deriv[0][2]= -1.-0.5*deriv[4][2]-0.5*deriv[6][2]-0.5*deriv[7][2];

      deriv[1][0]= 1.-0.5*deriv[4][0]-0.5*deriv[5][0]-0.5*deriv[8][0];
      deriv[1][1]= -0.5*deriv[4][1]-0.5*deriv[5][1]-0.5*deriv[8][1];
      deriv[1][2]= -0.5*deriv[4][2]-0.5*deriv[5][2]-0.5*deriv[8][2];

      deriv[2][0]= -0.5*deriv[5][0]-0.5*deriv[6][0]-0.5*deriv[9][0];
      deriv[2][1]= 1.-0.5*deriv[5][1]-0.5*deriv[6][1]-0.5*deriv[9][1];
      deriv[2][2]= -0.5*deriv[5][2]-0.5*deriv[6][2]-0.5*deriv[9][2];

      deriv[3][0]= -0.5*deriv[9][0]-0.5*deriv[7][0]-0.5*deriv[8][0];
      deriv[3][1]= -0.5*deriv[9][1]-0.5*deriv[7][1]-0.5*deriv[8][1];
      deriv[3][2]= 1.-0.5*deriv[9][2]-0.5*deriv[7][2]-0.5*deriv[8][2];
    }

    bool FeH1LagrangeTet::CoordIsInsideElem( const Vector<Double>& point,
                                             Double tolerance )  {
      const Double & xi = point[0];
      const Double & eta = point[1];
      const Double & zeta = point[2];
      bool isInside =
    		  (               xi >= (0 - tolerance)) &&
    		  (              eta >= (0 - tolerance)) &&
    		  (             zeta >= (0 - tolerance)) &&
    		  ((xi + eta + zeta) <= (1 + tolerance));
      return isInside;
    }

    void FeH1LagrangeTet::
    GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
    		const StdVector<UInt> & volConnect,
    		const LocPoint & surfIntPoint,
    		LocPoint & volIntPoint,
    		Vector<Double>& locNormal ) {

    	// Try to find out, which vertices are in common with
    	// the surface element. Then calculate the product of all four
    	// and compare them
    	//
    	//
    	// 4+\
    	//  |\ \           zeta
    	//  | \  \ 	      ^ eta
    	//  |  \  +3	      |/
    	//  |   \ |	      0--> xi
    	//  |    \ \
    	//  |     \|     REFERENCE TETRAHEDRAL ELEMENT
    	//  +------+
    	//  1      2

    	volIntPoint.coord.Resize(3);
    	locNormal.Resize(3);

    	StdVector<UInt> commonIndex(3);
    	UInt found = 0;
    	UInt indexProduct = 0;
    	std::string errMsg;

    	// loop over surface connect
    	for (UInt iSurf=0; iSurf<3; iSurf++)
    		// loop over volume connect
    		for (UInt iVol=0; iVol<4; iVol++)
    			if (surfConnect[iSurf] == volConnect[iVol])
    			{
    				commonIndex[found++] = iVol+1;
    			}

    	indexProduct =  commonIndex[0] * commonIndex[1] * commonIndex[2];

    	//std::cerr << "indexProduct = " << indexProduct << std::endl;
    	switch(indexProduct)
    	{
    	case 8:
    		// Surface[1,2,4] is common
    		volIntPoint[0] = surfIntPoint[0];
    		volIntPoint[1] = 0.0;
    		volIntPoint[2] = surfIntPoint[1];

    		locNormal[0] =  0.0;
    		locNormal[1] = -1.0;
    		locNormal[2] =  0.0;
    		break;

    	case 24:
    		// Surface[2,3,4] is common
    		volIntPoint[0] = surfIntPoint[0];
    		volIntPoint[1] = surfIntPoint[1];
    		volIntPoint[2] = 1.0 - surfIntPoint[0] - surfIntPoint[1];

    		locNormal[0] =  1.0;
    		locNormal[1] =  1.0;
    		locNormal[2] =  1.0;
    		break;

    	case 12:
    		// Surface[1,3,4] is common
    		volIntPoint[0] = 0.0;
    		volIntPoint[1] = surfIntPoint[0];
    		volIntPoint[2] = surfIntPoint[1];

    		locNormal[0] = -1.0;
    		locNormal[1] =  0.0;
    		locNormal[2] =  0.0;
    		break;

    	case 6:
    		// Surface[1,2,3] is common
    		volIntPoint[0] = surfIntPoint[0];
    		volIntPoint[1] = surfIntPoint[1];
    		volIntPoint[2] = 0.0;

    		locNormal[0] =  0.0;
    		locNormal[1] =  0.0;
    		locNormal[2] = -1.0;
    		break;
    	default:
    		EXCEPTION("FeH1LagrangeTet::GetLocalIntPoints4Surface: surface "
    				<< "and volume element have not three nodes in common. "
    				<< "Check your .mesh-file.");
    		break;
    	}
    }

    // --- Pyramid 1st order ---
     FeH1LagrangePyra1::FeH1LagrangePyra1() : FeH1LagrangePyra() {
       feType_ = Elem::ET_PYRA5;
       shape_ = Elem::shapes[feType_];
       actNumFncs_ = 5;
       order_ = 1;
     }

     FeH1LagrangePyra1::~FeH1LagrangePyra1() {

     }

     void FeH1LagrangePyra1::CalcShFnc( Vector<Double>& shape,
                                       const Vector<Double>& point,
                                       const Elem* ptElem,
                                       UInt comp ) {

    	 shape.Resize(5);

    	 // "A New Family of Finite Elements: The Pyramidal Elements"
    	 // F. Zgainski, J.L. Coulomb, Y. Marechal.
    	 // IEEE Transactions on Magnetics, Vol. 32, No. 3, May 1996, p. 1394

    	 shape[4] = point[2];

    	 shape[0] = 0.25*((1+point[0])*(1+point[1])-point[2]);
    	 shape[1] = 0.25*((1-point[0])*(1+point[1])-point[2]);
    	 shape[2] = 0.25*((1-point[0])*(1-point[1])-point[2]);
    	 shape[3] = 0.25*((1+point[0])*(1-point[1])-point[2]);

    	 if (point[2] != 1.0) {
    		 shape[0] += 0.25*(+(point[0]*point[1]*point[2])/(1-point[2]));
    		 shape[1] += 0.25*(-(point[0]*point[1]*point[2])/(1-point[2]));
    		 shape[2] += 0.25*(+(point[0]*point[1]*point[2])/(1-point[2]));
    		 shape[3] += 0.25*(-(point[0]*point[1]*point[2])/(1-point[2]));
    	 }
     }

     void FeH1LagrangePyra1::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                               const Vector<Double>& point,
                                               const Elem* ptElem,
                                               UInt comp ) {
    	 deriv.Resize(5,3);

    	 deriv.Init();

    	 deriv[4][0] = 0;
    	 deriv[4][1] = 0;
    	 deriv[4][2] = 1;

    	 deriv[0][0] =  0.25 * ( 1 + point[1]);
    	 deriv[0][1] =  0.25 * ( 1 + point[0]);
    	 deriv[0][2] = -0.25;

    	 deriv[1][0] =  0.25 * (-1 - point[1]);
    	 deriv[1][1] =  0.25 * ( 1 - point[0]);
    	 deriv[1][2] = -0.25;

    	 deriv[2][0] =  0.25 * (-1 + point[1]);
    	 deriv[2][1] =  0.25 * (-1 + point[0]);
    	 deriv[2][2] = -0.25;

    	 deriv[3][0] =  0.25 * ( 1 - point[1]);
    	 deriv[3][1] =  0.25 * (-1 - point[0]);
    	 deriv[3][2] = -0.25;

    	 if (point[2] != 1.0)
    	 {
    		 deriv[0][0] += 0.25 * (+(point[1]*point[2])/(1-point[2]));
    		 deriv[0][1] += 0.25 * (+(point[0]*point[2])/(1-point[2]));
    		 deriv[0][2] += 0.25 * (+(point[0]*point[1])/((1-point[2])*(1-point[2])));

    		 deriv[1][0] += 0.25 * (-(point[1]*point[2])/(1-point[2]));
    		 deriv[1][1] += 0.25 * (-(point[0]*point[2])/(1-point[2]));
    		 deriv[1][2] += 0.25 * (-(point[0]*point[1])/((1-point[2])*(1-point[2])));

    		 deriv[2][0] += 0.25 * (+(point[1]*point[2])/(1-point[2]));
    		 deriv[2][1] += 0.25 * (+(point[0]*point[2])/(1-point[2]));
    		 deriv[2][2] += 0.25 * (+(point[0]*point[1])/((1-point[2])*(1-point[2])));

    		 deriv[3][0] += 0.25 * (-(point[1]*point[2])/(1-point[2]));
    		 deriv[3][1] += 0.25 * (-(point[0]*point[2])/(1-point[2]));
    		 deriv[3][2] += 0.25 * (-(point[0]*point[1])/((1-point[2])*(1-point[2])));
    	 }
     }

    // --- Pyra 2nd order ---
     FeH1LagrangePyra2::FeH1LagrangePyra2() : FeH1LagrangePyra() {
       feType_ = Elem::ET_PYRA13;
       shape_ = Elem::shapes[feType_];
       actNumFncs_ = 13;
       order_ = 2;
     }

     FeH1LagrangePyra2::~FeH1LagrangePyra2() {

     }

     void FeH1LagrangePyra2::CalcShFnc( Vector<Double>& shape,
                                       const Vector<Double>& point,
                                       const Elem* ptElem,
                                       UInt comp ) {

    	 shape.Resize(13);

    	 // Order is the same as in source paper, just index ist changed to follow
    	 // our elements standard numbering
    	 shape[1]  = 0.25 * ( point[0]+point[1]-1) *
    			 ( (1+point[0]) * (1+point[1]) - point[2] );
    	 shape[0]  = 0.25 * ( point[0]-point[1]-1) *
    			 ( (1+point[0]) * (1-point[1]) - point[2] );
    	 shape[3]  = 0.25 * (-point[0]-point[1]-1) *
    			 ( (1-point[0]) * (1-point[1]) - point[2] );
    	 shape[2]  = 0.25 * (-point[0]+point[1]-1) *
    			 ( (1-point[0]) * (1+point[1]) - point[2] );
    	 shape[4]  = point[2] * ( 2 * point[2] - 1);
    	 shape[5]  = 0.0;
    	 shape[8]  = 0.0;
    	 shape[7]  = 0.0;
    	 shape[6]  = 0.0;
    	 shape[10] = 0.0;
    	 shape[9]  = 0.0;
    	 shape[12] = 0.0;
    	 shape[11] = 0.0;

    	 if (point[2] != 1.0)
    	 {
    		 Double fac = (point[0] * point[1] * point[2]) / (1-point[2]);
    		 shape[1]  += 0.25 * ( point[0]+point[1]-1) * (+fac);
    		 shape[0]  += 0.25 * ( point[0]-point[1]-1) * (-fac);
    		 shape[3]  += 0.25 * (-point[0]-point[1]-1) * (+fac);
    		 shape[2]  += 0.25 * (-point[0]+point[1]-1) * (-fac);
    		 shape[4]  += 0.0;
    		 shape[5]  += (1+point[0]-point[2]) * (1-point[1]-point[2])*
    				 (1+point[1]-point[2])/(2*(1-point[2]));
    		 shape[8]  += (1+point[0]-point[2])*(1-point[0]-point[2])*
    				 (1-point[1]-point[2])/(2*(1-point[2]));
    		 shape[7]  += (1-point[0]-point[2])*(1-point[1]-point[2])*
    				 (1+point[1]-point[2])/(2*(1-point[2]));
    		 shape[6]  += (1+point[0]-point[2])*(1-point[0]-point[2])*
    				 (1+point[1]-point[2])/(2*(1-point[2]));
    		 shape[10] += point[2]*(1+point[0]-point[2])*
    				 (1+point[1]-point[2])/(1-point[2]);
    		 shape[9]  += point[2]*(1+point[0]-point[2])*
    				 (1-point[1]-point[2])/(1-point[2]);
    		 shape[12] += point[2]*(1-point[0]-point[2])*
    				 (1-point[1]-point[2])/(1-point[2]);
    		 shape[11] += point[2]*(1-point[0]-point[2])*
    				 (1+point[1]-point[2])/(1-point[2]);
    	 }
     }

     void FeH1LagrangePyra2::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                               const Vector<Double>& point,
                                               const Elem* ptElem,
                                               UInt comp ) {
    	 deriv.Resize(13,3);

    	 deriv.Init();

    	 // Derivatives for the quadratic case.
    	 // Calculated symbolically with Maple.

    	 deriv[4][0] = 0;
    	 deriv[4][1] = 0;
    	 deriv[4][2] = 4*point[2]-1;

    	 if (point[2]==1)
    	 {
    		 deriv[1][0] = 0.25*( (1+point[0])*(1+point[1])-point[2]+
    				 (point[0]+point[1]-1)*(1+point[1]));
    		 deriv[1][1] = 0.25*( (1+point[0])*(1+point[1])-point[2]+
    				 (point[0]+point[1]-1)*(1+point[0]));
    		 deriv[1][2] = 0.25*( -point[0]-point[1]+1);
    		 deriv[0][0] = 0.25*( (1+point[0])*(1-point[1])-point[2]+
    				 (point[0]-point[1]-1)*(1-point[1]));
    		 deriv[0][1] = 0.25*( -(1+point[0])*(1-point[1])+point[2]+
    				 (point[0]-point[1]-1)*(-1-point[0]));
    		 deriv[0][2] = 0.25*( -point[0]+point[1]+1);
    		 deriv[3][0] = 0.25*( -(1-point[0])*(1-point[1])+point[2]+
    				 (-point[0]-point[1]-1)*(-1+point[1]));
    		 deriv[3][1] = 0.25*( -(1-point[0])*(1-point[1])+point[2]+
    				 (-point[0]-point[1]-1)*(-1+point[0]));
    		 deriv[3][2] = 0.25*( point[0]+point[1]+1);
    		 deriv[2][0] = 0.25*( -(1-point[0])*(1+point[1])+point[2]+
    				 (-point[0]+point[1]-1)*(-1-point[1]));
    		 deriv[2][1] = 0.25*( (1-point[0])*(1+point[1])-point[2]+
    				 (-point[0]+point[1]-1)*(1-point[0]));
    		 deriv[2][2] = 0.25*( point[0]-point[1]+1);

    		 deriv[6][0] = 0.0;
    		 deriv[6][1] = 0.0;
    		 deriv[6][2] = 0.0;

    		 deriv[7][0] = 0.0;
    		 deriv[7][1] = 0.0;
    		 deriv[7][2] = 0.0;

    		 deriv[8][0] = 0.0;
    		 deriv[8][1] = 0.0;
    		 deriv[8][2] = 0.0;

    		 deriv[5][0] = 0.0;
    		 deriv[5][1] = 0.0;
    		 deriv[5][2] = 0.0;

    		 deriv[10][0] = 0.0;
    		 deriv[10][1] = 0.0;
    		 deriv[10][2] = 0.0;

    		 deriv[11][0] = 0.0;
    		 deriv[11][1] = 0.0;
    		 deriv[11][2] = 0.0;

    		 deriv[12][0] = 0.0;
    		 deriv[12][1] = 0.0;
    		 deriv[12][2] = 0.0;

    		 deriv[9][0] = 0.0;
    		 deriv[9][1] = 0.0;
    		 deriv[9][2] = 0.0;
    	 }
    	 else
    	 {
    		 deriv[1][0] = 0.25*((1+point[0])*(1+point[1])-point[2]+
    				 point[0]*point[1]*point[2]/(1-point[2])+
    				 (point[0]+point[1]-1)*(1+point[1]+
    						 point[1]*point[2]/(1-point[2])));
    		 deriv[1][1] = 0.25*((1+point[0])*(1+point[1])-point[2]+
    				 point[0]*point[1]*point[2]/(1-point[2])+
    				 (point[0]+point[1]-1)*(1+point[0]+point[0]*
    						 point[2]/(1-point[2])));
    		 deriv[1][2] = 0.25*((point[0]+point[1]-1)*(-1+point[0]*point[1]/
    				 (1-point[2])+point[0]*point[1]*point[2]/
    				 ((1-point[2])*(1-point[2]))));

    		 deriv[0][0] = 0.25*((1+point[0])*(1-point[1])-point[2]-point[0]*
    				 point[1]*point[2]/(1-point[2])+
    				 (point[0]-point[1]-1)*(1-point[1]-
    						 point[1]*point[2]/(1-point[2])));
    		 deriv[0][1] = 0.25*(-(1+point[0])*(1-point[1])+point[2]+point[0]*
    				 point[1]*point[2]/(1-point[2])+
    				 (point[0]-point[1]-1)*(-1-point[0]-
    						 point[0]*point[2]/(1-point[2])));
    		 deriv[0][2] = 0.25*((point[0]-point[1]-1)*(-1-point[0]*point[1]/
    				 (1-point[2])-point[0]*point[1]*point[2]/
    				 ((1-point[2])*(1-point[2]))));

    		 deriv[3][0] = 0.25*(-(1-point[0])*(1-point[1])+point[2]-point[0]*
    				 point[1]*point[2]/(1-point[2])+(-point[0]-
    						 point[1]-1)*(-1+point[1]+point[1]*point[2]/
    								 (1-point[2])));
    		 deriv[3][1] = 0.25*(-(1-point[0])*(1-point[1])+point[2]-point[0]*
    				 point[1]*point[2]/(1-point[2])+
    				 (-point[0]-point[1]-1)*(-1+point[0]+point[0]*
    						 point[2]/(1-point[2])));
    		 deriv[3][2] = 0.25*((-point[0]-point[1]-1)*(-1+point[0]*point[1]/
    				 (1-point[2])+point[0]*point[1]*point[2]/
    				 ((1-point[2])*(1-point[2]))));

    		 deriv[2][0] = 0.25*(-(1-point[0])*(1+point[1])+point[2]+point[0]*
    				 point[1]*point[2]/(1-point[2])+(-point[0]+
    						 point[1]-1)*(-1-point[1]-point[1]*
    								 point[2]/(1-point[2])));
    		 deriv[2][1] = 0.25*((1-point[0])*(1+point[1])-point[2]-point[0]*
    				 point[1]*point[2]/(1-point[2])+
    				 (-point[0]+point[1]-1)*(1-point[0]-point[0]*
    						 point[2]/(1-point[2])));
    		 deriv[2][2] = 0.25*((-point[0]+point[1]-1)*(-1-point[0]*point[1]/
    				 (1-point[2])-point[0]*point[1]*point[2]/
    				 ((1-point[2])*(1-point[2]))));

    		 deriv[5][0] = .5*(1-point[1]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2]);
    		 deriv[5][1] = -.5*(1+point[0]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])+.5*(1+point[0]-point[2])*
    				 (1-point[1]-point[2])/(1-point[2]);
    		 deriv[5][2] = -.5*(1-point[1]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])-.5*(1+point[0]-point[2])*
    				 (1+point[1]-point[2])/(1-point[2])-.5*
    				 (1+point[0]-point[2])*(1-point[1]-point[2])/
    				 (1-point[2])+.5*(1+point[0]-point[2])*
    				 (1-point[1]-point[2])*(1+point[1]-point[2])/
    				 ((1-point[2])*(1-point[2]));

    		 deriv[8][0] = .5*(1-point[0]-point[2])*(1-point[1]-point[2])/
    				 (1-point[2])-.5*(1+point[0]-point[2])*
    				 (1-point[1]-point[2])/(1-point[2]);
    		 deriv[8][1] = -.5*(1+point[0]-point[2])*(1-point[0]-point[2])/
    				 (1-point[2]);
    		 deriv[8][2] = -.5*(1-point[0]-point[2])*(1-point[1]-point[2])/
    				 (1-point[2])-.5*(1+point[0]-point[2])*
    				 (1-point[1]-point[2])/(1-point[2])-.5*
    				 (1+point[0]-point[2])*(1-point[0]-point[2])/
    				 (1-point[2])+.5*(1+point[0]-point[2])*
    				 (1-point[0]-point[2])*(1-point[1]-point[2])/
    				 ((1-point[2])*(1-point[2]));

    		 deriv[7][0] = -.5*(1-point[1]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2]);
    		 deriv[7][1] = -.5*(1-point[0]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])+.5*(1-point[0]-point[2])*
    				 (1-point[1]-point[2])/(1-point[2]);
    		 deriv[7][2] = -.5*(1-point[1]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])-.5*(1-point[0]-point[2])*
    				 (1+point[1]-point[2])/(1-point[2])-.5*(1-point[0]-
    						 point[2])*(1-point[1]-point[2])/(1-point[2])+
    						 .5*(1-point[0]-point[2])*(1-point[1]-point[2])*
    						 (1+point[1]-point[2])/((1-point[2])*(1-point[2]));

    		 deriv[6][0] = .5*(1-point[0]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])-.5*(1+point[0]-point[2])*
    				 (1+point[1]-point[2])/(1-point[2]);
    		 deriv[6][1] = .5*(1+point[0]-point[2])*(1-point[0]-point[2])/
    				 (1-point[2]);
    		 deriv[6][2] = -.5*(1-point[0]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])-.5*(1+point[0]-point[2])*
    				 (1+point[1]-point[2])/(1-point[2])-.5*(1+point[0]-
    						 point[2])*(1-point[0]-point[2])/
    						 (1-point[2])+.5*(1+point[0]-point[2])*
    						 (1-point[0]-point[2])*(1+point[1]-point[2])/
    						 ((1-point[2])*(1-point[2]));

    		 deriv[10][0] = point[2]*(1+point[1]-point[2])/(1-point[2]);
    		 deriv[10][1] = point[2]*(1+point[0]-point[2])/(1-point[2]);
    		 deriv[10][2] = (1+point[0]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])-point[2]*(1+point[1]-point[2])/
    				 (1-point[2])-point[2]*(1+point[0]-point[2])/
    				 (1-point[2])+point[2]*(1+point[0]-point[2])*
    				 (1+point[1]-point[2])/((1-point[2])*(1-point[2]));

    		 deriv[9][0] = point[2]*(1-point[1]-point[2])/(1-point[2]);
    		 deriv[9][1] = -point[2]*(1+point[0]-point[2])/(1-point[2]);
    		 deriv[9][2] = (1+point[0]-point[2])*(1-point[1]-point[2])/
    				 (1-point[2])-point[2]*(1-point[1]-point[2])/
    				 (1-point[2])-point[2]*(1+point[0]-point[2])/
    				 (1-point[2])+point[2]*(1+point[0]-point[2])*
    				 (1-point[1]-point[2])/((1-point[2])*(1-point[2]));

    		 deriv[12][0] = -point[2]*(1-point[1]-point[2])/(1-point[2]);
    		 deriv[12][1] = -point[2]*(1-point[0]-point[2])/(1-point[2]);
    		 deriv[12][2] = (1-point[0]-point[2])*(1-point[1]-point[2])/
    				 (1-point[2])-point[2]*(1-point[1]-point[2])/
    				 (1-point[2])-point[2]*(1-point[0]-point[2])/
    				 (1-point[2])+point[2]*(1-point[0]-point[2])*
    				 (1-point[1]-point[2])/((1-point[2])*(1-point[2]));

    		 deriv[11][0] = -point[2]*(1+point[1]-point[2])/(1-point[2]);
    		 deriv[11][1] = point[2]*(1-point[0]-point[2])/(1-point[2]);
    		 deriv[11][2] = (1-point[0]-point[2])*(1+point[1]-point[2])/
    				 (1-point[2])-point[2]*(1+point[1]-point[2])/
    				 (1-point[2])-point[2]*(1-point[0]-point[2])/
    				 (1-point[2])+point[2]*(1-point[0]-point[2])*
    				 (1+point[1]-point[2])/((1-point[2])*(1-point[2]));
    	 }
     }

     // --- Pyra 2nd order ---
      FeH1LagrangePyra14::FeH1LagrangePyra14() : FeH1LagrangePyra() {
        feType_ = Elem::ET_PYRA14;
        shape_ = Elem::shapes[feType_];
        actNumFncs_ = 14;
        order_ = 2;
      }

      FeH1LagrangePyra14::~FeH1LagrangePyra14() {
      }

      void FeH1LagrangePyra14::CalcShFnc( Vector<Double>& shape,
                                          const Vector<Double>& point,
                                          const Elem* ptElem,
                                          UInt comp ) {
        static bool issueWarning = true;
        if(issueWarning) {
          WARN("CalcShFnc for ET_PYRA14 implemented but not tested!");
          issueWarning = false;
        }

        // Shape functions for 14 node pyramid taken from
        // http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch19.d/AFEM.Ch19.pdf
        shape.Resize(actNumFncs_);

        // Transform coordinates to unit cube [-1,1]^3.
        Double x,y;
        x = 1.0 / (1.0-(point[2]-1e-10)) * point[0];
        //x = point[0] + (point[0] > 0.0 ? 1.0 : -1.0)*point[2];
        x *= -1.0;


        y = 1.0 / (1.0-(point[2]-1e-10)) * point[1];
        //y = point[1] + (point[1] > 0.0 ? 1.0 : -1.0)*point[2];
        y *= -1.0;

        Double z = 2*point[2]-1;

        Double t1 = 3*x;
        Double t2 = 3*y;
        Double t3 = 2*x*y;
        Double t4 = 2*z;
        Double t5 = x*z;
        Double t6 = y*z;
        Double t7 = 2*x*y*z;


        shape[4] =   1/2.0 * z * (1+z);
        if(z<1){
          shape[13] = 1/2.0 * (1 - x*x) * (1 - y*y) * (1 - z);

          shape[2] = -1/16.0 * (1-x) * (1-y) * (1-z) *
                    (4 + t1 + t2 + t3 + t4 + t5 + t6 + t7) + 1/4.0 * shape[13];
          shape[3] = -1/16.0 * (1+x) * (1-y) * (1-z) *
                    (4 - t1 + t2 - t3 + t4 - t5 + t6 - t7) + 1/4.0 * shape[13];
          shape[0] = -1/16.0 * (1+x) * (1+y) * (1-z) *
                    (4 - t1 - t2 + t3 + t4 - t5 - t6 + t7) + 1/4.0 * shape[13];
          shape[1] = -1/16.0 * (1-x) * (1+y) * (1-z) *
                    (4 + t1 - t2 - t3 + t4 + t5 - t6 - t7) + 1/4.0 * shape[13];



          shape[7] =   1/8.0 * (1 - x*x) * (1 - y) * (1 - z) *
                    (2 + y + y*z) - 1/2.0 * shape[13];
          shape[8] =   1/8.0 * (1 + x ) * (1 - y*y) * (1 - z) *
                    (2 - x - x * z) - 1/2.0 * shape[13];
          shape[5] =   1/8.0 * (1 - x*x) * (1 + y) * (1 - z) *
                    (2 - y - y*z) - 1/2.0 * shape[13];
          shape[6] =   1/8.0 * (1 - x ) * (1 - y*y) * (1 - z) *
                    (2 + x + x*z) - 1/2.0 * shape[13];

          shape[11] =  1/4.0 * (1 - x ) * (1 - y) * (1 - z*z);
          shape[12] =  1/4.0 * (1 + x ) * (1 - y) * (1 - z*z);
          shape[9] =  1/4.0 * (1 + x ) * (1 + y) * (1 - z*z);
          shape[10] =  1/4.0 * (1 - x ) * (1 + y) * (1 - z*z);
        }

#if 0
        // An alternative approach is described in:
        // On Higher Order Pyramidal Finite Elements
        // Liping Liu, Kevin B. Davies, Michal Krizek and Li Guan
        // DOI: 10.4208/aamm.09-m0989
        // cf. http://www.dune-project.org/flyspray/index.php?getfile=369

        Double x = 2*point[0] + point[2]-1.0;
        Double y = 2*point[1] + point[2]-1.0;
        Double z = zeta;

        if (x > y)
        {
          // vertices
          shape[0] = 0.25*(x + z)*(x + z - 1)*(y - z - 1)*(y - z);
          shape[1] = -0.25*(x + z)*(y - z)*((x + z + 1)*(-y + z + 1) - 4*z) - z*(x - y);
          shape[2] = 0.25*(x + z)*(y - z)*(y - z + 1)*(x + z - 1);
          shape[3] = 0.25*(y - z)*(x + z)*(y - z + 1)*(x + z + 1);
          shape[4] = z*(2*z - 1);

          // lower edges
          shape[5] = -0.5*(y - z + 1)*(x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1));
          shape[6] = -0.5*(y - z + 1)*(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1));
          shape[7] = -0.5*(x + z - 1)*(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1));
          shape[8] = -0.5*(y - z + 1)*(x + z - 1)*(x + 1)*y;

          // upper edges
          shape[9] = z*(x + z - 1)*(y - z - 1);
          shape[10] = -z*((x + z + 1)*(y - z - 1) + 4*z);
          shape[11] = -z*(y - z + 1)*(x + z - 1);
          shape[12] = z*(y - z + 1)*(x + z + 1);

          // base face
          shape[13] = (y - z + 1)*(x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1));
        }
        else
        {
          // vertices
          shape[0] = 0.25*(y + z)*(y + z - 1)*(x - z - 1)*(x - z);
          shape[1] = -0.25*(x - z)*(y + z)*(x - z + 1)*(-y - z + 1);
          shape[2] = 0.25*(x - z)*(y + z)*((x - z - 1)*(y + z + 1) + 4*z) + z*(x - y);
          shape[3] = 0.25*(y + z)*(x - z)*(x - z + 1)*(y + z + 1);
          shape[4] = z*(2*z - 1);

          // lower edges
          shape[5] = -0.5*(y + z - 1)*(((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1));
          shape[6] = -0.5*(x - z + 1)*(y + z - 1)*(y + 1)*x;
          shape[7] = -0.5*(x - z + 1)*(y + z - 1)*(x - 1)*y;
          shape[8] = -0.5*(x - z + 1)*(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1));

          // upper edges
          shape[9] = z*(y + z - 1)*(x - z - 1);
          shape[10] = -z*(x - z + 1)*(y + z - 1);
          shape[11] = -z*((y + z + 1)*(x - z - 1) + 4*z);
          shape[12] = z*(x - z + 1)*(y + z + 1);

          // base face
          shape[13] = (x - z + 1)*(y + z - 1)*((y + 1)*(x - 1) - z*(x - y - z - 1));
        }
#endif
      }

      void FeH1LagrangePyra14::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                                  const Vector<Double>& point,
                                                  const Elem* ptElem,
                                                  UInt comp ) {
        static bool issueWarning = true;
        if(issueWarning) {
          WARN("CalcLocDerivShFnc for ET_PYRA14 implemented but not tested!");
          issueWarning = false;
        }
        deriv.Resize(actNumFncs_, 3);
        // Transform coordinates to unit cube [-1,1]^3.
        Double x,y;

        x = 1.0 / (1.0-(point[2]-1e-10)) * point[0];
        //x = point[0] + (point[0] > 0.0 ? 1.0 : -1.0)*point[2];
        x *= -1.0;

        y = 1.0 / (1.0-(point[2]-1e-10)) * point[1];
        //y = point[1] + (point[1] > 0.0 ? 1.0 : -1.0)*point[2];
        y *= -1.0;

        Double z = 2*point[2]-1;
        deriv.Init();

        // http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch19.d/AFEM.Ch19.pdf

        // func 5
        deriv[4][0] =          0;
        deriv[4][1] =          0;
        deriv[4][2] = 1.0*z + 0.5;

        if(z<1){
          Double N_cx = -x*(1-y*y)*(1-z);
          Double N_cy = -y*(1-x*x)*(1-z);
          Double N_cz = -0.5*(1-x*x)*(1-y*y);
          // func 1
            deriv[2][0] = (1.0/16.0) * (1-y)*(1-z)*(1+6*x+y+4*x*y+z+2*x*z-y*z+4*x*z*y)+0.25*N_cx;
            deriv[2][1] = (1.0/16.0) * (1-x)*(1-z)*(1+x+6*y+4*x*y+z-x*z+2*y*z+4*x*z*y)+0.25*N_cy;
            deriv[2][2] = (1.0/8.0)  * (1-x)*(1-y)*(1+x+y+2*z+x*z+y*z+2*x*z*y)+0.25*N_cz;

            // func 2
            deriv[3][0] = (-1.0/16.0) * (1-y)*(1-z)*(1-6*x+y-4*x*y+z-2*x*z-y*z-4*x*z*y)+0.25*N_cx;
            deriv[3][1] = (1.0/16.0) * (1+x)*(1-z)*(1-x+6*y-4*x*y+z+x*z+2*y*z-4*x*z*y)+0.25*N_cy;
            deriv[3][2] = (1.0/8.0)  * (1+x)*(1-y)*(1-x+y+2*z-x*z+y*z-2*x*z*y)+0.25*N_cz;

            // func 3
            deriv[0][0] = (-1.0/16.0) * (1+y)*(1-z)*(1-6*x-y+4*x*y+z-2*x*z+y*z+4*x*z*y)+0.25*N_cx;
            deriv[0][1] = (-1.0/16.0) * (1+x)*(1-z)*(1-x-6*y+4*x*y+z+x*z-2*y*z+4*x*z*y)+0.25*N_cy;
            deriv[0][2] = (1.0/8.0)  * (1+x)*(1+y)*(1-x-y+2*z-x*z-y*z+2*x*z*y)+0.25*N_cz;

            // func 4
            deriv[1][0] = (1.0/16.0) * (1+y)*(1-z)*(1+6*x-y-4*x*y+z+2*x*z+y*z-4*x*z*y)+0.25*N_cx;
            deriv[1][1] = (-1.0/16.0) * (1-x)*(1-z)*(1+x-6*y-4*x*y+z-x*z-2*y*z-4*x*z*y)+0.25*N_cy;
            deriv[1][2] = (1.0/8.0)  * (1-x)*(1+y)*(1+x-y+2*z+x*z-y*z-2*x*z*y)+0.25*N_cz;

            // func 6
            deriv[7][0] = -0.25*x*(1-y)*(1-z)*(2+y+y*z)-0.5*N_cx;
            deriv[7][1] = -0.125*(1-x*x)*(1-z)*(1+2*y-z+2*y*z)-0.5*N_cy;
            deriv[7][2] = -0.25*(1-x*x)*(1-y)*(1+y*z)-0.5*N_cz;

            // func 7
            deriv[8][0] =  0.125*(1-y*y)*(1-z)*(1-2*x-z-2*x*z)-0.5*N_cx;
            deriv[8][1] = -0.25*y*(1+x)*(1-z)*(2-x-x*z)-0.5*N_cy;
            deriv[8][2] = -0.25*(1+x)*(1-y*y)*(1-x*z)-0.5*N_cz;

            // func 8
            deriv[5][0] = -0.25*x*(1+y)*(1-z)*(2-y-y*z)-0.5*N_cx;
            deriv[5][1] =  0.125*(1-x*x)*(1-z)*(1-2*y-z-2*y*z)-0.5*N_cy;
            deriv[5][2] = -0.25*(1+x*x)*(1+y)*(1-y*z)-0.5*N_cz;

            // func 9
            deriv[6][0] = -0.125*(1-y*y)*(1-z)*(1+2*x-z+2*x*z)-0.5*N_cx;
            deriv[6][1] = -0.25*y*(1-x)*(1-z)*(2+x+x*z)-0.5*N_cy;
            deriv[6][2] = -0.25*(1-x)*(1-y*y)*(1+x*z)-0.5*N_cz;

            // func 10
            deriv[11][0] = -0.25*(1-y)*(1-z*z);
            deriv[11][1] = -0.25*(1-x)*(1-z*z);
            deriv[11][2] = -0.5*(1-x)*(1-y)*z;

            // func 11
            deriv[12][0] =  0.25*(1-y)*(1-z*z);
            deriv[12][1] = -0.25*(1+x)*(1-z*z);
            deriv[12][2] = -0.5*(1+x)*(1-y)*z;

            // func 12
            deriv[9][0] = 0.25*(1+y)*(1-z*z);
            deriv[9][1] = 0.25*(1+x)*(1-z*z);
            deriv[9][2] = -0.5*(1+x)*(1+y)*z;

            // func 13
            deriv[10][0] = -0.25*(1+y)*(1-z*z);
            deriv[10][1] =  0.25*(1-y)*(1-z*z);
            deriv[10][2] = -0.5*(1-x)*(1+y)*z;

            // func 14
            deriv[13][0] = N_cx;
            deriv[13][1] = N_cy;
            deriv[13][2] = N_cz;


///++++++++++++++++++++++++++++++++++++OLD IMPLEMENTATION
          //
          //// func 1
          //deriv[2][0] =         -0.25*x*(-y*y + 1)*(-z + 1) + (0.0625*x - 0.0625)*(-y + 1)*(-z + 1)*(2*y*z + 2*y + z + 3) + 0.0625*(-y + 1)*(-z + 1)*(2*x*y*z + 2*x*y + x*z + 3*x + y*z + 3*y + 2*z + 4);
          //deriv[2][1] =-0.5*y*(-0.5*x*x + 0.5)*(-z + 1) + (0.0625*x - 0.0625)*(-y + 1)*(-z + 1)*(2*x*z + 2*x + z + 3) - (0.0625*x - 0.0625)*(-z + 1)*(2*x*y*z + 2*x*y + x*z + 3*x + y*z + 3*y + 2*z + 4);
          //deriv[2][2] = (0.0625*x - 0.0625)*(-y + 1)*(-z + 1)*(2*x*y + x + y + 2) - (0.0625*x - 0.0625)*(-y + 1)*(2*x*y*z + 2*x*y + x*z + 3*x + y*z + 3*y + 2*z + 4) - 0.25*(-0.5*x*x + 0.5)*(-y*y + 1);
          //
          //// func 2
          //deriv[3][0] =          -0.25*x*(-y*y + 1)*(-z + 1) + (-0.0625*x - 0.0625)*(-y + 1)*(-z + 1)*(-2*y*z - 2*y - z - 3) - 0.0625*(-y + 1)*(-z + 1)*(-2*x*y*z - 2*x*y - x*z - 3*x + y*z + 3*y + 2*z + 4);
          //deriv[3][1] =-0.5*y*(-0.5*x*x + 0.5)*(-z + 1) + (-0.0625*x - 0.0625)*(-y + 1)*(-z + 1)*(-2*x*z - 2*x + z + 3) - (-0.0625*x - 0.0625)*(-z + 1)*(-2*x*y*z - 2*x*y - x*z - 3*x + y*z + 3*y + 2*z + 4);
          //deriv[3][2] = (-0.0625*x - 0.0625)*(-y + 1)*(-z + 1)*(-2*x*y - x + y + 2) - (-0.0625*x - 0.0625)*(-y + 1)*(-2*x*y*z - 2*x*y - x*z - 3*x + y*z + 3*y + 2*z + 4) - 0.25*(-0.5*x*x + 0.5)*(-y*y + 1);
          //
          //// func 3
          //deriv[0][0] =           -0.25*x*(-y*y + 1)*(-z + 1) + (-0.0625*x - 0.0625)*(y + 1)*(-z + 1)*(2*y*z + 2*y - z - 3) - 0.0625*(y + 1)*(-z + 1)*(2*x*y*z + 2*x*y - x*z - 3*x - y*z - 3*y + 2*z + 4);
          //deriv[0][1] =-0.5*y*(-0.5*x*x + 0.5)*(-z + 1) + (-0.0625*x - 0.0625)*(y + 1)*(-z + 1)*(2*x*z + 2*x - z - 3) + (-0.0625*x - 0.0625)*(-z + 1)*(2*x*y*z + 2*x*y - x*z - 3*x - y*z - 3*y + 2*z + 4);
          //deriv[0][2] =  (-0.0625*x - 0.0625)*(y + 1)*(-z + 1)*(2*x*y - x - y + 2) - (-0.0625*x - 0.0625)*(y + 1)*(2*x*y*z + 2*x*y - x*z - 3*x - y*z - 3*y + 2*z + 4) - 0.25*(-0.5*x*x + 0.5)*(-y*y + 1);
          //
          //// func 4
          //deriv[1][0] =          -0.25*x*(-y*y + 1)*(-z + 1) + (0.0625*x - 0.0625)*(y + 1)*(-z + 1)*(-2*y*z - 2*y + z + 3) + 0.0625*(y + 1)*(-z + 1)*(-2*x*y*z - 2*x*y + x*z + 3*x - y*z - 3*y + 2*z + 4);
          //deriv[1][1] =-0.5*y*(-0.5*x*x + 0.5)*(-z + 1) + (0.0625*x - 0.0625)*(y + 1)*(-z + 1)*(-2*x*z - 2*x - z - 3) + (0.0625*x - 0.0625)*(-z + 1)*(-2*x*y*z - 2*x*y + x*z + 3*x - y*z - 3*y + 2*z + 4);
          //deriv[1][2] =  (0.0625*x - 0.0625)*(y + 1)*(-z + 1)*(-2*x*y + x - y + 2) - (0.0625*x - 0.0625)*(y + 1)*(-2*x*y*z - 2*x*y + x*z + 3*x - y*z - 3*y + 2*z + 4) - 0.25*(-0.5*x*x + 0.5)*(-y*y + 1);
          //
          //// func 6
          //deriv[7][0] =                                                             -0.25*x*(-y + 1)*(-z + 1)*(y*z + y + 2) + 0.5*x*(-y*y + 1)*(-z + 1);
          //deriv[7][1] =1.0*y*(-0.5*x*x + 0.5)*(-z + 1) + (-0.125*x*x + 0.125)*(-y + 1)*(-z + 1)*(z + 1) - (-0.125*x*x + 0.125)*(-z + 1)*(y*z + y + 2);
          //deriv[7][2] =     y*(-0.125*x*x + 0.125)*(-y + 1)*(-z + 1) + 0.5*(-0.5*x*x + 0.5)*(-y*y + 1) - (-0.125*x*x + 0.125)*(-y + 1)*(y*z + y + 2);
          //
          //// func 7
          //deriv[8][0] =    0.5*x*(-y*y + 1)*(-z + 1) + (0.125*x + 0.125)*(-y*y + 1)*(-z - 1)*(-z + 1) + 0.125*(-y*y + 1)*(-z + 1)*(-x*z - x + 2);
          //deriv[8][1] =                                           -2*y*(0.125*x + 0.125)*(-z + 1)*(-x*z - x + 2) + 1.0*y*(-0.5*x*x + 0.5)*(-z + 1);
          //deriv[8][2] = -x*(0.125*x + 0.125)*(-y*y + 1)*(-z + 1) - (0.125*x + 0.125)*(-y*y + 1)*(-x*z - x + 2) + 0.5*(-0.5*x*x + 0.5)*(-y*y + 1);
          //
          //// func 8
          //deriv[5][0] =                                                              -0.25*x*(y + 1)*(-z + 1)*(-y*z - y + 2) + 0.5*x*(-y*y + 1)*(-z + 1);
          //deriv[5][1] =1.0*y*(-0.5*x*x + 0.5)*(-z + 1) + (-0.125*x*x + 0.125)*(y + 1)*(-z - 1)*(-z + 1) + (-0.125*x*x + 0.125)*(-z + 1)*(-y*z - y + 2);
          //deriv[5][2] =      -y*(-0.125*x*x + 0.125)*(y + 1)*(-z + 1) + 0.5*(-0.5*x*x + 0.5)*(-y*y + 1) - (-0.125*x*x + 0.125)*(y + 1)*(-y*z - y + 2);
          //
          //// func 9
          //deriv[6][0] =     0.5*x*(-y*y + 1)*(-z + 1) + (-0.125*x + 0.125)*(-y*y + 1)*(-z + 1)*(z + 1) - 0.125*(-y*y + 1)*(-z + 1)*(x*z + x + 2);
          //deriv[6][1] =                                           -2*y*(-0.125*x + 0.125)*(-z + 1)*(x*z + x + 2) + 1.0*y*(-0.5*x*x + 0.5)*(-z + 1);
          //deriv[6][2] =x*(-0.125*x + 0.125)*(-y*y + 1)*(-z + 1) - (-0.125*x + 0.125)*(-y*y + 1)*(x*z + x + 2) + 0.5*(-0.5*x*x + 0.5)*(-y*y + 1);
          //
          //// func 10
          //deriv[11][0] =    -0.25*(-y + 1)*(-z*z + 1);
          //deriv[11][1] = -(-0.25*x + 0.25)*(-z*z + 1);
          //deriv[11][2] =-2*z*(-0.25*x + 0.25)*(-y + 1);
          //
          //// func 11
          //deriv[12][0] =    0.25*(-y + 1)*(-z*z + 1);
          //deriv[12][1] = -(0.25*x + 0.25)*(-z*z + 1);
          //deriv[12][2] =-2*z*(0.25*x + 0.25)*(-y + 1);
          //
          //// func 12
          //deriv[9][0] =    0.25*(y + 1)*(-z*z + 1);
          //deriv[9][1] = (0.25*x + 0.25)*(-z*z + 1);
          //deriv[9][2] =-2*z*(0.25*x + 0.25)*(y + 1);
          //
          //// func 13
          //deriv[10][0] =    -0.25*(y + 1)*(-z*z + 1);
          //deriv[10][1] = (-0.25*x + 0.25)*(-z*z + 1);
          //deriv[10][2] =-2*z*(-0.25*x + 0.25)*(y + 1);
          //
          //// func 14
          //deriv[13][0] =    -1.0*x*(-y*y + 1)*(-z + 1);
          //deriv[13][1] =-2*y*(-0.5*x*x + 0.5)*(-z + 1);
          //deriv[13][2] = -(-0.5*x*x + 0.5)*(-y*y + 1);
////++++++++++++++++++++++++++++++ENDOLD
        }
#if 0
        // An alternative approach is described in:
        // On Higher Order Pyramidal Finite Elements
        // Liping Liu, Kevin B. Davies, Michal Krizek and Li Guan
        // DOI: 10.4208/aamm.09-m0989
        // cf. http://www.dune-project.org/flyspray/index.php?getfile=369

        Double x = 2*point[0] + point[2]-1.0;
        Double y = 2*point[1] + point[2]-1.0;
        Double z = zeta;

            if (x > y)
              {
                // vertices
                deriv[0][0] = 0.5*(y - z - 1)*(y - z)*(2*x + 2*z - 1);
                deriv[0][1] = 0.5*(x + z)*(x + z - 1)*(2*y - 2*z - 1);
                deriv[0][2] = 0.5*(deriv[0][0] + deriv[0][1])
                  + 0.25*((2*x + 2*z - 1)*(y - z - 1)*(y - z)
                          + (x + z)*(x + z - 1)*(-2*y + 2*z + 1));

                deriv[1][0] = 2*(-0.25*((y - z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                      + (x + z)*(y - z)*(-y + z + 1)) - z);
                deriv[1][1] = 2*(-0.25*((x + z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                      + (x + z)*(y - z)*(-(x + z + 1))) + z);
                deriv[1][2] = 0.5*(deriv[1][0] + deriv[1][1])
                  - 0.25*((y - z)*((x + z + 1)*(-y + z + 1) - 4*z)
                          - (x + z)*((x + z + 1)*(-y + z + 1) - 4*z)
                          + (x + z)*(y - z)*(x - y + 2*z - 2))
                  - (x - y);

                deriv[2][0] = 0.5*(y - z)*(y - z + 1)*(2*x + 2*z - 1);
                deriv[2][1] = 0.5*(x + z)*(2*y - 2*z + 1)*(x + z - 1);
                deriv[2][2] = 0.5*(deriv[2][0] + deriv[2][1])
                  + 0.25*((y - x - 2*z)*(y - z + 1)*(x + z - 1)
                          + (x + z)*(y - z)*(y - x - 2*z + 2));

                deriv[3][0] = 0.5*(y - z)*(2*x + 2*z + 1)*(y - z + 1);
                deriv[3][1] = 0.5*(2*y - 2*z + 1)*(x + z)*(x + z + 1);
                deriv[3][2] = 0.5*(deriv[3][0] + deriv[3][1])
                  + 0.25*((y - x - 2*z)*(y - z + 1)*(x + z + 1)
                          + (y - z)*(x + z)*(y - x - 2*z));

                deriv[4][0] = 0;
                deriv[4][1] = 0;
                deriv[4][2] = 4*z - 1;

                // lower edges
                deriv[5][0] = -((y - z + 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                              + (y - z + 1)*(x + z - 1)*((y - 1) + z));
                deriv[5][1] = -((x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                              + (y - z + 1)*(x + z - 1)*((x + 1) - z));
                deriv[5][2] = 0.5*(deriv[5][0] + deriv[5][1])
                  - 0.5*((-x + y - 2*z + 2)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                         + (y - z + 1)*(x + z - 1)*(x - y + 2*z + 1));

                deriv[6][0] = -(y - z + 1)*(2*x + z + 1)*(y - 1);
                deriv[6][1] = -(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1)
                              + (y - z + 1)*((x + z + 1)*x + 2*z));
                deriv[6][2] = 0.5*(deriv[6][0] + deriv[6][1])
                  - 0.5*(-(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1))
                         + (y - z + 1)*(((y - 1)*x - 1) + 2*y + 1));

                deriv[7][0] = -(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1)
                              + (x + z - 1)*((y - z - 1)*y + 2*z));
                deriv[7][1] = -(x + z - 1)*(2*y - z - 1)*(x + 1);
                deriv[7][2] = 0.5*(deriv[7][0] + deriv[7][1])
                  - 0.5*(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1)
                         + (x + z - 1)*((-(x + 1)*y - 1) + 2*x + 1));

                deriv[8][0] = -(y - z + 1)*(2*x + z)*y;
                deriv[8][1] = -(2*y - z + 1)*(x + z - 1)*(x + 1);
                deriv[8][2] = 0.5*(deriv[8][0] + deriv[8][1])
                  - 0.5*(-x + y - 2*z + 2)*(x + 1)*y;

                // upper edges
                deriv[9][0] = 2*z*(y - z - 1);
                deriv[9][1] = 2*z*(x + z - 1);
                deriv[9][2] = 0.5*(deriv[9][0] + deriv[9][1])
                  + (x + z - 1)*(y - z - 1) + z*(-x + y - 2*z);

                deriv[10][0] = -2*z*(y - z - 1);
                deriv[10][1] = -2*z*(x + z + 1);
                deriv[10][2] = 0.5*(deriv[10][0] + deriv[10][1])
                  - ((x + z + 1)*(y - z - 1) + 4*z)
                  - z*(-x + y - 2*z + 2);

                deriv[11][0] = -2*z*(y - z + 1);
                deriv[11][1] = -2*z*(x + z - 1);
                deriv[11][2] = 0.5*(deriv[11][0] + deriv[11][1])
                  - (y - z + 1)*(x + z - 1) - z*(-x + y - 2*z + 2);

                deriv[12][0] = 2*z*(y - z + 1);
                deriv[12][1] = 2*z*(x + z + 1);
                deriv[12][2] = 0.5*(deriv[12][0] + deriv[12][1])
                  + (y - z + 1)*(x + z + 1) + z*(-x + y - 2*z);

                // base face
                deriv[13][0] = 2*((y - z + 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                                + (y - z + 1)*(x + z - 1)*(y - 1 + z));
                deriv[13][1] = 2*((x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                                + (y - z + 1)*(x + z - 1)*(x + 1 - z));
                deriv[13][2] = 0.5*(deriv[13][0] + deriv[13][1])
                  + ((-x + y - 2*z + 2)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                     + (y - z + 1)*(x + z - 1)*(x - y + 2*z + 1));
              }
            else
              {
                // vertices
                deriv[0][0] = 0.5*(y + z)*(y + z - 1)*(2*x - 2*z - 1);
                deriv[0][1] = 0.5*(2*y + 2*z - 1)*(x - z - 1)*(x - z);
                deriv[0][2] = 0.5*(deriv[0][0] + deriv[0][1])
                  + 0.25*((2*y + 2*z - 1)*(x - z - 1)*(x - z)
                          + (y + z)*(y + z - 1)*(-2*x + 2*z + 1));

                deriv[1][0] = -0.5*(y + z)*(2*x - 2*z + 1)*(-y - z + 1);
                deriv[1][1] = -0.5*(x - z)*(x - z + 1)*(-2*y - 2*z + 1);
                deriv[1][2] = 0.5*(deriv[1][0] + deriv[1][1])
                  - 0.25*((x - y - 2*z)*(x - z + 1)*(-y - z + 1)
                          + (x - z)*(y + z)*(-x + y + 2*z - 2));

                deriv[2][0] = 0.5*((y + z)*((x - z - 1)*(y + z + 1) + 4*z)
                                 + (x - z)*(y + z)*(y + z + 1) + 4*z);
                deriv[2][1] = 0.5*((x - z)*((x - z - 1)*(y + z + 1) + 4*z)
                                 + (x - z)*(y + z)*(x - z - 1) - 4*z);
                deriv[2][2] = 0.5*(deriv[2][0] + deriv[2][1])
                  + 0.25*((x - y - 2*z)*((x - z - 1)*(y + z + 1) + 4*z)
                          + (x - z)*(y + z)*(x - y - 2*z + 2) + 4*(x - y));

                deriv[3][0] = 0.5*(y + z)*(2*x - 2*z + 1)*(y + z + 1);
                deriv[3][1] = 0.5*(x - z)*(x - z + 1)*(2*y + 2*z + 1);
                deriv[3][2] = 0.5*(deriv[3][0] + deriv[3][1])
                  + 0.25*((x - y - 2*z)*(x - z + 1)*(y + z + 1)
                          + (y + z)*(x - z)*(x - y - 2*z));

                deriv[4][0] = 0;
                deriv[4][1] = 0;
                deriv[4][2] = 4*z - 1;

                // lower edges
                deriv[5][0] = -(y + z - 1)*(2*x - z - 1)*(y + 1);
                deriv[5][1] = -(((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1)
                              + (y + z - 1)*((x - z - 1)*x + 2*z));
                deriv[5][2] = 0.5*(deriv[5][0] + deriv[5][1])
                  - 0.5*((((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1))
                         + (y + z - 1)*((-(y + 1)*x - 1) + 2*y + 1));

                deriv[6][0] = -(2*x - z + 1)*(y + z - 1)*(y + 1);
                deriv[6][1] = -(x - z + 1)*(2*y + z)*x;
                deriv[6][2] = 0.5*(deriv[6][0] + deriv[6][1])
                  - 0.5*(x - y - 2*z + 2)*(y + 1)*x;

                deriv[7][0] = -(2*x - z)*(y + z - 1)*y;
                deriv[7][1] = -(x - z + 1)*(2*y + z - 1)*(x - 1);
                deriv[7][2] = 0.5*(deriv[7][0] + deriv[7][1])
                  - 0.5*(x - y - 2*z + 2)*(x - 1)*y;

                deriv[8][0] = -(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1)
                              + (x - z + 1)*((y + z + 1)*y + 2*z));
                deriv[8][1] = -(x - z + 1)*(2*y + z + 1)*(x - 1);
                deriv[8][2] = 0.5*(deriv[8][0] + deriv[8][1])
                  - 0.5*(-(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1))
                         + (x - z + 1)*(((x - 1)*y - 1) + 2*x + 1));

                // upper edges
                deriv[9][0] = 2*z*(y + z - 1);
                deriv[9][1] = 2*z*(x - z - 1);
                deriv[9][2] = 0.5*(deriv[9][0] + deriv[9][1])
                    + (y + z - 1)*(x - z - 1) + z*(x - y - 2*z);

                deriv[10][0] = -2*z*(y + z - 1);
                deriv[10][1] = -2*z*(x - z + 1);
                deriv[10][2] = 0.5*(deriv[10][0] + deriv[10][1])
                  - (x - z + 1)*(y + z - 1) - z*(x - y - 2*z + 2);

                deriv[11][0] = -2*z*(y + z + 1);
                deriv[11][1] = -2*z*(x - z - 1);
                deriv[11][2] = 0.5*(deriv[11][0] + deriv[11][1])
                  - ((y + z + 1)*(x - z - 1) + 4*z) - z*(x - y - 2*z + 2);

                deriv[12][0] = 2*z*(y + z + 1);
                deriv[12][1] = 2*z*(x - z + 1);
                deriv[12][2] = 0.5*(deriv[12][0] + deriv[12][1])
                  + (x - z + 1)*(y + z + 1) + z*(x - y - 2*z);

                // base face
                deriv[13][0] = 2*((y + z - 1)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                                + (x - z + 1)*(y + z - 1)*(y + 1 - z));
                deriv[13][1] = 2*((x - z + 1)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                                + (x - z + 1)*(y + z - 1)*(x - 1 + z));
                deriv[13][2] = 0.5*(deriv[13][0] + deriv[13][1])
                  + (x - y - 2*z + 2)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                  + (x - z + 1)*(y + z - 1)*(-(x - y - 2*z - 1));
              }
#endif
     }

     bool FeH1LagrangePyra::CoordIsInsideElem( const Vector<Double>& point,
                                               Double tolerance )  {
         const Double & xi = point[0];
         const Double & eta = point[1];
         const Double & zeta = point[2];
         double threshold = 1 - zeta;

         bool isInside =
      		   (zeta >= (0 - tolerance)) &&
      		   (zeta <= (1.0 + tolerance)) &&
      		   (  xi >= (-threshold - tolerance)) &&
      		   ( eta >= (-threshold - tolerance)) &&
      		   (  xi <= (threshold + tolerance)) &&
      		   ( eta <= (threshold + tolerance));
         return isInside;
     }

     void FeH1LagrangePyra::
     GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
     		const StdVector<UInt> & volConnect,
     		const LocPoint & surfIntPoint,
     		LocPoint & volIntPoint,
     		Vector<Double>& locNormal ) {
    	 // Try to find out, which vertices are in common with
    	 // the surface element. Then calculate the product of all four
    	 // and compare them
    	 //                                   zeta
    	 //             5                     ^  eta
    	 //             +                     |/
    	 //           // \                    0--> xi
    	 //          // \ \
    	 //         / / \  \
    	 //        / /   \  \
    	 //     2 +-/----\ ---+ 1
    	 //      / /     \   /
    	 //     / /      \  /     REFERENCE VOLUME ELEMENT
    	 //    //        \ /
    	 //  3+-----------+ 4
    	 //

         volIntPoint.coord.Resize(3);
	     locNormal.Resize(3);

    	 // Check if surface element is triangle
    	 // or quadrilateral
    	 if (surfConnect.GetSize() == 3 ||
    			 surfConnect.GetSize() == 6)
    	 {
    		 // ---- Triangle Surface ---
    		 StdVector<Integer> commonIndex(3);
    		 Integer found = 0;
    		 Integer indexProduct = 0;
    		 std::string errMsg;

    		 // loop over surface connect
    		 for (Integer iSurf=0; iSurf<3; iSurf++)
    			 // loop over volume connect
    			 for (Integer iVol=0; iVol<5; iVol++)
    				 if (surfConnect[iSurf] == volConnect[iVol]) {
    					 commonIndex[found++] = iVol+1;
    				 }
    		 indexProduct =  commonIndex[0] * commonIndex[1] * commonIndex[2];

    		 // Now we have to consider the following:
    		 // - The extension of the triangular element is from [0..1] in both
    		 //   local directions xi and eta
    		 // - The side length of the base-rectangular side is in each direction
    		 //   [-1..+1], so we need a mapping in
    		 // - The

    		 // General rule:


    		 switch( indexProduct ) {
    		 case 10:
    			 // Surface[1,2,5] is common
    			 volIntPoint[0] = - 2.0 * (surfIntPoint[0] - 0.5);
    			 volIntPoint[1] =   1.0 - surfIntPoint[1];
    			 volIntPoint[2] =   surfIntPoint[1];

			 locNormal[0] =  0;
                         locNormal[1] = -1;
                         locNormal[2] =  1;
    			 break;
    		 case 30:
    			 // Surface[2,3,5] is common
    			 volIntPoint[0] =   surfIntPoint[1] - 1.0;
    			 volIntPoint[1] = - 2.0 * (surfIntPoint[0] - 0.5);
    			 volIntPoint[2] =   surfIntPoint[1];
			 
			 locNormal[0] =  1;
                         locNormal[1] =  0;
                         locNormal[2] =  1;
			 break;
    		 case 60:
    			 // Surface[3,4,5] is common
    			 volIntPoint[0] =   2.0 * (surfIntPoint[0] - 0.5);
    			 volIntPoint[1] =   surfIntPoint[1] - 1.0;
    			 volIntPoint[2] =   surfIntPoint[1];
			 
			 locNormal[0] =  0;
                         locNormal[1] =  1;
                         locNormal[2] =  1;

    			 break;
    		 case 20:
    			 // Surface[4,1,5] is common
    			 volIntPoint[0] =  1.0 - surfIntPoint[1];
    			 volIntPoint[1] =  2.0 * (surfIntPoint[0] - 0.5);
    			 volIntPoint[2] =  surfIntPoint[1];
			 
			 locNormal[0] = -1;
                         locNormal[1] =  0;
                         locNormal[2] =  1;
    			 break;
    		 default:
    			 EXCEPTION("FeH1LagrangePyra::GetLocalIntPoints4Surface: surface and volume element "
    					 << "have not three nodes in common. Check your mesh.");
    			 break;
    		 }
    	 } else {
    		 // ---- Quadrilateral Surface ---
    		 StdVector<Integer> commonIndex(4);
    		 Integer found = 0;
    		 Integer indexSum = 0;
    		 std::string errMsg;

    		 // loop over surface connect
    		 for (Integer iSurf=0; iSurf<4; iSurf++)
    			 // loop over volume connect
    			 for (Integer iVol=0; iVol<5; iVol++)
    				 if (surfConnect[iSurf] == volConnect[iVol])
    				 {
    					 commonIndex[found++] = iVol+1;
    				 }
    		 indexSum =  commonIndex[0] + commonIndex[1]
    		                                          + commonIndex[2] + commonIndex[3];

    		 // Safety check: Check, that the surface element is
    		 // really located on the bottom of the pyramid.
    		 if( indexSum != 10 ) {
    			 EXCEPTION("FeH1LagrangePyra::GetLocalIntPoints4Surface: surface and volume element "
    					 << "have not four nodes in common. Check your mesh.");
    		 }

    		 volIntPoint[0] = surfIntPoint[0];
    		 volIntPoint[1] = surfIntPoint[1];
    		 volIntPoint[2] = 0.0; // always on bottom

                 locNormal[0] =  0;
                 locNormal[1] =  0;
                 locNormal[2] = -1;
    	 }
     }
} // namespace CoupledField
