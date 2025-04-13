#include "H1ElemsLagExpl.hh"

#include "Utils/AutoDiff.hh"


namespace CoupledField {

  // ========================================================================
  //  FeH1LagrangeExpl 
  // ========================================================================
  
  FeH1LagrangeExpl::FeH1LagrangeExpl()
  : FeH1(), FeNodal() {
   order_ = 0;
   preComputShFnc_ = true;
  }

  FeH1LagrangeExpl::FeH1LagrangeExpl(const FeH1LagrangeExpl & other)
                   : FeH1(other), FeNodal(other){
    this->order_ = other.order_;
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
      numFcns.Resize(shape_.numEdges);
      numFcns.Init((order_-1));
    }else if( fctEntityType == FACE 
        && completeType_ == TENSOR_TYPE) {
      // we only have face nodes, if we have the full tensorial
      // elements. Here only 4-sided faces have interior nodes.
      // collect all faces with 4 nodes
      numFcns.Resize(shape_.numFaces);
      numFcns.Init(0);
      for( UInt i = 0; i < shape_.numFaces; ++i ) {
        if( shape_.faceVertices[i].GetSize() == 4 ) {
          numFcns[i] = (order_-1)*(order_-1);
        }
      }
    }else if( fctEntityType == INTERIOR
            && feType_ == Elem::ET_HEXA27 ){
      numFcns.Resize(1);
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
    }else if( fctEntityType == FACE && ptElem->extended->faces.GetSize() > 0) {
      if(completeType_ == TENSOR_TYPE){
        //WARNING: for order > 2 we would need to check for orientation. See H1ElemsLagVar.cc
        fncPermutation.Resize((order_-1) * (order_-1));
        for(UInt i = 0; i< order_-1 ; i++){
          for(UInt j = 0; j< order_-1 ; j++){
            fncPermutation[(i*(order_-1)) + j] = i*(order_-1) + j;
          }
        }
      }else{
        //in case of serendipity elements, we apply an offset to face and interior
        //equation numbers as they have edge and face DOFs only for order > 3
        //as those are not implemented, we just throw an error
        if(order_>3){
          Exception("Function FeH1LagrangeExpl::GetNodalPermutation needs to be extended for higher order serendipity elements!");
        }
        fncPermutation.Resize(0);
      }
    }else if( fctEntityType == INTERIOR && ptElem->extended->faces.GetSize() > 3){
      if(completeType_ == SERENDIPITY_TYPE){
        if(order_>3){
          Exception("Function FeH1LagrangeExpl::GetNodalPermutation needs to be extended for higher order serendipity elements!");
        }else{
          fncPermutation.Resize(0);
        }
      }else{
        UInt numIFncs = (order_-1)*(order_-1)*(order_-1);
        fncPermutation.Resize(numIFncs);
        for(UInt k = 0; k< numIFncs ; k++){
          fncPermutation[k] = k;
        }
      }
    }else{
      fncPermutation.Resize(0);
    }
  }
  
  bool FeH1LagrangeExpl::operator==( const FeH1LagrangeExpl& comp) const {
    bool ret = true;
    ret &= this->feType_ == comp.feType_;
    ret &= this->order_ == comp.order_;
    return ret;
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
  
  
  void FeH1LagrangeExpl::GetLocalDOFCoordinates(Matrix<Double> & coordMat) {
    coordMat.Resize(actNumFncs_, shape_.dim);
    coordMat.Init();
    
    for( UInt iNode = 0; iNode < shape_.numNodes; ++iNode ) {
      for( UInt iDim = 0; iDim < shape_.dim; ++iDim ) {
        coordMat[iNode][iDim] = shape_.nodeCoords[iNode][iDim];  
      }
    }
  }

  void FeH1LagrangeExpl::MapQuadSurfOrientation(const std::map<UInt, UInt>& commonIndexMap, 
                                               const StdVector<UInt>& commonIndex, 
                                               const LocPoint& surfIntPoint, 
                                               StdVector<Double>& volIntPoint) {
    // consider orientation of the surface
    switch (commonIndexMap.at(commonIndex[0])) {
      // first index is 1
      case 1: {
        switch (commonIndexMap.at(commonIndex[1])) {
          case 2: {
            // Orientation [1,2,3,4]:
            volIntPoint[0] = surfIntPoint[0];
            volIntPoint[1] = surfIntPoint[1];
            break;
          }
          case 4: {
            // Orientation [1,4,3,2]:
            // here the surface needs to be flipped along the 1-3 line
            volIntPoint[0] = surfIntPoint[1];
            volIntPoint[1] = surfIntPoint[0];
            break;
          }
          default:
            EXCEPTION("FeH1LagrangeExpl::MapQuadSurfOrientation: Unexpected orientation of surface element!");
        }
        break;
      }
      // first index is 2
      case 2: {
        switch (commonIndexMap.at(commonIndex[1])) {
          case 3: {
            // Orientation [2,3,4,1]:
            // here the surface needs to be rotated by 90° counter-clockwise
            volIntPoint[0] = -surfIntPoint[1];
            volIntPoint[1] = surfIntPoint[0];
            break;
          }
          case 1: {
            // Orientation [2,1,4,3]:
            // here the surface needs to be flipped left-to-right
            volIntPoint[0] = -surfIntPoint[0];
            volIntPoint[1] = surfIntPoint[1];
            break;
          }
          default:
            EXCEPTION("FeH1LagrangeExpl::MapQuadSurfOrientation: Unexpected orientation of surface element!");
        }
        break;
      }
      // first index is 3
      case 3: {
        switch (commonIndexMap.at(commonIndex[1])) {
          case 4: {
            // Orientation [3,4,1,2]:
            // here the surface needs to be rotated by 180°
            volIntPoint[0] = -surfIntPoint[0];
            volIntPoint[1] = -surfIntPoint[1];
            break;
          }
          case 2: {
            // Orientation [3,2,1,4]:
            // here the surface needs to be flipped along the 2-4 line
            volIntPoint[0] = -surfIntPoint[1];
            volIntPoint[1] = -surfIntPoint[0];
            break;
          }
          default:
            EXCEPTION("FeH1LagrangeExpl::MapQuadSurfOrientation: Unexpected orientation of surface element!");
        }
        break;
      }
      // first index is 4
      case 4: {
        switch (commonIndexMap.at(commonIndex[1])) {
          case 1: {
            // Orientation [4,1,2,3]:
            // here the surface needs to be rotated by 90° clockwise
            volIntPoint[0] = surfIntPoint[1];
            volIntPoint[1] = -surfIntPoint[0];
            break;
          }
          case 3: {
            // Orientation [4,3,2,1]:
            // here the surface needs to be flipped up-to-down
            volIntPoint[0] = surfIntPoint[0];
            volIntPoint[1] = -surfIntPoint[1];
            break;
          }
          default:
            EXCEPTION("FeH1LagrangeExpl::MapQuadSurfOrientation: Unexpected orientation of surface element!");
        }
        break;
      }
      default:
        EXCEPTION("FeH1LagrangeExpl::MapQuadSurfOrientation: Unexpected orientation of surface element!");
    }
  }

  void FeH1LagrangeExpl::MapTriaSurfOrientation(const std::map<UInt, UInt>& commonIndexMap, 
                                                const StdVector<UInt>& commonIndex, 
                                                const LocPoint& surfIntPoint, 
                                                StdVector<Double>& volIntPoint,
                                                bool isEquilateral) {
    if (!isEquilateral) {
      // consider orientation of the rectangular triangles
      switch (commonIndexMap.at(commonIndex[0])) {
        // first index is 1
        case 1: {
          switch (commonIndexMap.at(commonIndex[1])) {
            case 2: {
              // Orientation [1,2,3]:
              volIntPoint[0] = surfIntPoint[0];
              volIntPoint[1] = surfIntPoint[1];
              break;
            }
            case 3: {
              // Orientation [1,3,2]:
              // here 2 and 3 are flipped
              volIntPoint[0] = surfIntPoint[1];
              volIntPoint[1] = surfIntPoint[0];
              break;
            }
            default:
              EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
          }
          break;
        }
        // first index is 2
        case 2: {
          switch (commonIndexMap.at(commonIndex[1])) {
            case 1: {
              // Orientation [2,1,3]:
              // here we flip left and right and shear by (-eta', 1) afterwards
              volIntPoint[0] = 1 - surfIntPoint[0] - surfIntPoint[1];
              volIntPoint[1] = surfIntPoint[1];
              break;
            }
            case 3: {
              // Orientation [2,3,1]:
              // here we shear +(1,xi), rotate 90° counter-clockwise, and add (0,1)
              volIntPoint[0] = - surfIntPoint[0] - surfIntPoint[1] + 1.0;
              volIntPoint[1] = surfIntPoint[0];
              break;
            }
            default:
              EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
          }
          break;
        }
        // first index is 3
        case 3: {
          switch (commonIndexMap.at(commonIndex[1])) {
            case 1: {
              // Orientation [3,1,2]:
              // here we rotate by -90°, add (0,1), and shear by (1, -xi') afterwards
              volIntPoint[0] = surfIntPoint[1];
              volIntPoint[1] = 1 - surfIntPoint[0] - surfIntPoint[1];
              break;
            }
            case 2: {
              // Orientation [3,2,1]:
              // here we flip up-down and shear by (1, -xi') afterwards
              volIntPoint[0] = surfIntPoint[0];
              volIntPoint[1] = 1 - surfIntPoint[0] - surfIntPoint[1];
              break;
            }
            default:
              EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
          }
          break;
        }
        default:
          EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
      }
    } else { //!isEquilateral
      // here we need to map the standard triangle to the equilateral triangle with coords 1: (-1,0), 2: (1,0), 3: (0,sqrt(3))
      // consider orientation of the equilateral triangles
      switch (commonIndexMap.at(commonIndex[0])) {
        // first index is 1
        case 1: {
          switch (commonIndexMap.at(commonIndex[1])) {
            case 2: {
              // Orientation [1,2,3]:
              volIntPoint[0] = 2.0 * surfIntPoint[0] + surfIntPoint[1] - 1.0;
              volIntPoint[1] = sqrt(3.0) * surfIntPoint[1];
              break;
            }
            case 3: {
              // Orientation [1,3,2]:
              volIntPoint[0] = 2.0 * surfIntPoint[1] + surfIntPoint[0] - 1.0;
              volIntPoint[1] = sqrt(3.0) * surfIntPoint[0];
              break;
            }
            default:
              EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
          }
          break;
        }
        // first index is 2
        case 2: {
          switch (commonIndexMap.at(commonIndex[1])) {
            case 1: {
              // Orientation [2,1,3]:
              volIntPoint[0] = -2.0 * surfIntPoint[0] - surfIntPoint[1] + 1.0;
              volIntPoint[1] = sqrt(3.0) * surfIntPoint[1];
              break;
            }
            case 3: {
              // Orientation [2,3,1]:
              volIntPoint[0] = -2.0 * surfIntPoint[1] - surfIntPoint[0] + 1.0;
              volIntPoint[1] = sqrt(3.0) * surfIntPoint[0];
              break;
            }
            default:
              EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
          }
          break;
        }
        // first index is 3
        case 3: {
          switch (commonIndexMap.at(commonIndex[1])) {
            case 1: {
              // Orientation [3,1,2]:
              volIntPoint[0] = surfIntPoint[1] - surfIntPoint[0];
              volIntPoint[1] = sqrt(3.0) * (1.0 - surfIntPoint[0] - surfIntPoint[1]);
              break;
            }
            case 2: {
              // Orientation [3,2,1]:
              volIntPoint[0] = surfIntPoint[0] - surfIntPoint[1];
              volIntPoint[1] = sqrt(3.0) * (1.0 - surfIntPoint[0] - surfIntPoint[1]);
              break;
            }
            default:
              EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
          }
          break;
        }
        default:
          EXCEPTION("FeH1LagrangeExpl::MapTriaSurfOrientation: Unexpected orientation of surface element!");
      } //switch (commonIndexMap.at(commonIndex[0]))
    } //if(!isEquilateral)
  }

  // ========================================================================
  //  Lagrangian Elements of 0th order
  // ========================================================================
  
  // --- Point 0th order ---
  
  FeH1LagrangePoint::FeH1LagrangePoint() {
    feType_ = Elem::ET_POINT;
    completeType_ = TENSOR_TYPE;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 1;
    order_ = 0;
  }
  FeH1LagrangePoint::~FeH1LagrangePoint() {
    
  }
  
  void FeH1LagrangePoint::CalcShFnc( Vector<Double>& shape,
                                     const Vector<Double>& point,
                                     const Elem* ptElem,
                                     UInt comp ) {
     shape.Resize( 1 );
     shape[0] = 1;
  }
  
  void FeH1LagrangePoint::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {

    deriv.Resize(1,1);
    deriv[0][0] = 1;
  }
  
  
  bool FeH1LagrangePoint::CoordIsInsideElem( const Vector<Double>& point,
                                            Double tolerance )  {
   const Double & xi = point[0];
   return (xi >= (- tolerance) &&
           xi <= (tolerance));
  }
  
  void FeH1LagrangePoint::
  GetLocalIntPoints4Surface(const StdVector<UInt> & surfConnect,
                            const StdVector<UInt> & volConnect,
                            const LocPoint & surfIntPoint,
                            LocPoint & volIntPoint,
                            Vector<Double>& locNormal ) {
    EXCEPTION("Not implemented");
  }

  
  // ========================================================================
  //  Lagrangian Elements of 1st order
  // ========================================================================
  
  // --- Line 1st order ---
  
  FeH1LagrangeLine1::FeH1LagrangeLine1() : FeH1LagrangeLine() {
    feType_ = Elem::ET_LINE2;
    completeType_ = TENSOR_TYPE;
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
    completeType_ = TENSOR_TYPE;
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

    // Node coords: 1: (0,0), 2: (1,0), 3: (0,1)
    // NOTE: Since the line element is defined in the range [-1;+1]
    // we have to calculate (1+surfCoord)/2 in order to get the right
    // position on the triangular element

    StdVector<UInt> commonIndex(2);
    UInt found = 0;
    UInt indexProduct = 0;

    volIntPoint.coord.Resize(2);
    locNormal.Resize(2);

    // loop over surface connect
    for (UInt iSurf=0; iSurf<2; iSurf++) {
      // loop over volume connect
      for (UInt iVol=0; iVol<3; iVol++) {
        if (surfConnect[iSurf] == volConnect[iVol])
          commonIndex[found++] = iVol+1;
      }
    }

    indexProduct= commonIndex[0] * commonIndex[1];
    switch(indexProduct)
    {
      case 2:
        // Edge[1,2] is common
        if (commonIndex[0] == 1) {
          volIntPoint[0] = 0.5 + (surfIntPoint[0] / 2.0);
        } else {
          volIntPoint[0] = 0.5 - (surfIntPoint[0] / 2.0);
        }
        volIntPoint[1] = 0.0;
        locNormal[0] = 0;
        locNormal[1] = -1.0;
        break;

      case 3:
        // Edge[1,3] is common
        volIntPoint[0] = 0.0;
        if (commonIndex[0] == 1) {
          volIntPoint[1] = 0.5 + (surfIntPoint[0] / 2.0);
        } else {
          volIntPoint[1] = 0.5 - (surfIntPoint[0] / 2.0);
        }
        locNormal[0] = -1.0;
        locNormal[1] =  0.0;
        break;

      case 6:
        // Edge[2,3] is common
        if (commonIndex[0] == 2) {
          volIntPoint[0] = 0.5 - (surfIntPoint[0] / 2.0);
          volIntPoint[1] = 0.5 + (surfIntPoint[0] / 2.0);
        } else {
          volIntPoint[0] = 0.5 + (surfIntPoint[0] / 2.0);
          volIntPoint[1] = 0.5 - (surfIntPoint[0] / 2.0);
        }
        locNormal[0] = sqrt(.5);
        locNormal[1] = sqrt(.5);
        break;

      default:
        EXCEPTION( "TriangleFE::GetLocalIntPoints4Surface: surface and volume element "
            << "have not two nodes in common. Check your .mesh-file.");
    }
  }


  void FeH1LagrangeTria::
  ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C){

  }

  
  void FeH1LagrangeTria1::Triangulate(StdVector< StdVector<UInt> > & triConnect){

    //trivial
    triConnect.Resize(1);
    triConnect.Init(StdVector<UInt>(3));
    triConnect[0][0] = 0;
    triConnect[0][1] = 1;
    triConnect[0][2] = 2;
  }

  // --- Quad 1st order ---
   
  FeH1LagrangeQuad1::FeH1LagrangeQuad1() : FeH1LagrangeQuad() {
    feType_ = Elem::ET_QUAD4;
    completeType_ = TENSOR_TYPE;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 4;
    order_ = 1; 
    // This element supports incompatible modes
    hasICModes_ = true;
  }
    
  FeH1LagrangeQuad1::~FeH1LagrangeQuad1() {
    
  }
  
  void FeH1LagrangeQuad1::ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C){
    P.Resize(4,3);
    P.Init();
    C.Resize(4,4);
    C.Init();
    C[0][0] = 0.25;
    C[0][1] = -0.25;
    C[0][2] = 0.25;
    C[0][3] = -0.25;

    C[1][0] = 0.25;
    C[1][1] = 0.25;
    C[1][2] = -0.25;
    C[1][3] = -0.25;

    C[2][0] = 0.25;
    C[2][1] = 0.25;
    C[2][2] = 0.25;
    C[2][3] = 0.25;

    C[3][0] = 0.25;
    C[3][1] = -0.25;
    C[3][2] = -0.25;
    C[3][3] = 0.25;

    P[1][0] = 1;
    P[2][1] = 1;
    P[3][0] = 1;
    P[3][1] = 1;
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
    //
    // Node coords: 1: (-1,-1), 2: (1,-1), 3: (1,1), 4: (-1,1)

    StdVector<UInt> commonIndex(2);
    UInt found = 0;
    UInt indexProduct = 0;
  
    volIntPoint.coord.Resize(2);
    locNormal.Resize(2);
  
    // loop over surface connect
    for (UInt iSurf=0; iSurf<2; iSurf++) {
      // loop over volume connect
      for (UInt iVol=0; iVol<4; iVol++) {
        if (surfConnect[iSurf] == volConnect[iVol])
          commonIndex[found++] = iVol+1;
      }
    }

    indexProduct= commonIndex[0] * commonIndex[1];
    switch(indexProduct) {
      case 2: {
        // Edge[1,2] is common
        if(commonIndex[0] == 2)
          volIntPoint[0] = -surfIntPoint[0];
        else
          volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = -1.0;
        locNormal[0] =  0.0;
        locNormal[1] = -1.0;
        break;
      }
      case 12: {
        // Edge[4,3] is common
        if(commonIndex[0] == 3)
          volIntPoint[0] = -surfIntPoint[0];
        else
          volIntPoint[0] = surfIntPoint[0];
        volIntPoint[1] = 1.0;
        locNormal[0] =  0.0;
        locNormal[1] =  1.0;
        break;
      }
      case 4: {
        // Edge[1,4] is common
        if(commonIndex[0] == 4)
          volIntPoint[1] = -surfIntPoint[0];
        else
          volIntPoint[1] = surfIntPoint[0];
        volIntPoint[0] = -1.0;
        locNormal[0] = -1.0;
        locNormal[1] =  0.0;
        break;
      }
      case 6: {
        // Edge[2,3] is common
        if(commonIndex[0] == 3)
          volIntPoint[1] = -surfIntPoint[0];
        else
          volIntPoint[1] = surfIntPoint[0];
        volIntPoint[0] = 1.0;
        locNormal[0] =  1.0;
        locNormal[1] =  0.0;
        break;
      }
      default:
        EXCEPTION( "RectangleFE::GetLocalIntPoints4Surface: surface and volume element "
                    <<  "have not two nodes in common. Check your mesh-file.");
    }
  }
  
  void FeH1LagrangeQuad1::Triangulate(StdVector< StdVector<UInt> > & triConnect){
    triConnect.Resize(2);
    triConnect.Init(StdVector<UInt>(3));

    //create two triangles in counterclockwise orientation
    //diagonal is in both cases third edge
    triConnect[0][0] = 3;
    triConnect[0][1] = 0;
    triConnect[0][2] = 1;

    triConnect[1][0] = 1;
    triConnect[1][1] = 2;
    triConnect[1][2] = 3;
  }

  // --- Hex 1st order ---
  FeH1LagrangeHex1::FeH1LagrangeHex1() : FeH1LagrangeHex() {
    feType_ = Elem::ET_HEXA8;
    completeType_ = TENSOR_TYPE;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 8;
    order_ = 1; 
    
    // This element supports incompatible modes
    hasICModes_ = true;
  }
    
  FeH1LagrangeHex1::~FeH1LagrangeHex1() {
    
  }
  
  void FeH1LagrangeHex1::ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C){
          P.Resize(actNumFncs_,3);
          C.Resize(actNumFncs_,actNumFncs_);

          P[0][0] = 0; P[0][1] = 0; P[0][2] = 0;
          P[1][0] = 0; P[1][1] = 0; P[1][2] = 1;
          P[2][0] = 0; P[2][1] = 1; P[2][2] = 0;
          P[3][0] = 0; P[3][1] = 1; P[3][2] = 1;
          P[4][0] = 1; P[4][1] = 0; P[4][2] = 0;
          P[5][0] = 1; P[5][1] = 0; P[5][2] = 1;
          P[6][0] = 1; P[6][1] = 1; P[6][2] = 0;
          P[7][0] = 1; P[7][1] = 1; P[7][2] = 1;

          C[0][0] =  1.2500000e-01; C[0][1] = -1.2500000e-01; C[0][2] = -1.2500000e-01;
          C[0][3] =  1.2500000e-01; C[0][4] = -1.2500000e-01; C[0][5] =  1.2500000e-01;
          C[0][6] =  1.2500000e-01; C[0][7] = -1.2500000e-01;

          C[1][0] =  1.2500000e-01; C[1][1] = -1.2500000e-01; C[1][2] = -1.2500000e-01;
          C[1][3] =  1.2500000e-01; C[1][4] =  1.2500000e-01; C[1][5] = -1.2500000e-01;
          C[1][6] = -1.2500000e-01; C[1][7] =  1.2500000e-01;

          C[2][0] =  1.2500000e-01; C[2][1] = -1.2500000e-01; C[2][2] =  1.2500000e-01;
          C[2][3] = -1.2500000e-01; C[2][4] =  1.2500000e-01; C[2][5] = -1.2500000e-01;
          C[2][6] =  1.2500000e-01; C[2][7] = -1.2500000e-01;

          C[3][0] =  1.2500000e-01; C[3][1] = -1.2500000e-01; C[3][2] =  1.2500000e-01;
          C[3][3] = -1.2500000e-01; C[3][4] = -1.2500000e-01; C[3][5] =  1.2500000e-01;
          C[3][6] = -1.2500000e-01; C[3][7] =  1.2500000e-01;

          C[4][0] =  1.2500000e-01; C[4][1] =  1.2500000e-01; C[4][2] = -1.2500000e-01;
          C[4][3] = -1.2500000e-01; C[4][4] = -1.2500000e-01; C[4][5] = -1.2500000e-01;
          C[4][6] =  1.2500000e-01; C[4][7] =  1.2500000e-01;

          C[5][0] =  1.2500000e-01; C[5][1] =  1.2500000e-01; C[5][2] = -1.2500000e-01;
          C[5][3] = -1.2500000e-01; C[5][4] =  1.2500000e-01; C[5][5] =  1.2500000e-01;
          C[5][6] = -1.2500000e-01; C[5][7] = -1.2500000e-01;

          C[6][0] =  1.2500000e-01; C[6][1] =  1.2500000e-01; C[6][2] =  1.2500000e-01;
          C[6][3] =  1.2500000e-01; C[6][4] =  1.2500000e-01; C[6][5] =  1.2500000e-01;
          C[6][6] =  1.2500000e-01; C[6][7] =  1.2500000e-01;

          C[7][0] =  1.2500000e-01; C[7][1] =  1.2500000e-01; C[7][2] =  1.2500000e-01;
          C[7][3] =  1.2500000e-01; C[7][4] = -1.2500000e-01; C[7][5] = -1.2500000e-01;
          C[7][6] = -1.2500000e-01; C[7][7] = -1.2500000e-01;
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
    //
    // Node coords: 1: (-1,-1,-1), 2: (1,-1,-1), 3: (1,1,-1), 4: (-1,1,-1),
    //              5: (-1,-1,1), 6: (1,-1,1), 7: (1,1,1), 8: (-1,1,1)

    StdVector<UInt> commonIndex(4);
    UInt found = 0;
    UInt indexProduct = 0;

    volIntPoint.coord.Resize(3);
    locNormal.Resize(3);

    // loop over surface connect
    for (UInt iSurf=0; iSurf<4; iSurf++) {
      // loop over volume connect
      for (UInt iVol=0; iVol<8; iVol++) {
        if (surfConnect[iSurf] == volConnect[iVol])
          commonIndex[found++] = iVol+1;
      }
    }

    indexProduct =  commonIndex[0] * commonIndex[1];
    indexProduct *= commonIndex[2] * commonIndex[3];

    switch (indexProduct) {
      // Surface [1,2,3,4] is common: downside surface, normal in -zeta direction
      case 24: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {2, 2}, {3, 3}, {4, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = tmpVolIntPoint[1];
        volIntPoint[2] = -1.0;
        locNormal[0] =  0.0;
        locNormal[1] =  0.0;
        locNormal[2] = -1.0;
        break;
      }
      // Surface [5,6,7,8] is common: upside surface normal in +zeta direction
      case 1680: {
        const std::map<UInt, UInt> commonIndexMap{{5, 1}, {6, 2}, {7, 3}, {8, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = tmpVolIntPoint[1];
        volIntPoint[2] = 1.0;
        locNormal[0] =  0.0;
        locNormal[1] =  0.0;
        locNormal[2] =  1.0;
        break;
      }
      // Surface [2,3,7,6] is common: right-side surface, normal in +xi direction
      case 252: {
        const std::map<UInt, UInt> commonIndexMap{{2, 1}, {3, 2}, {7, 3}, {6, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = 1.0;
        volIntPoint[1] = tmpVolIntPoint[0];
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] =  1.0;
        locNormal[1] =  0.0;
        locNormal[2] =  0.0;
        break;
      }
      // Surface [1,4,8,5] is common: left-side surface, normal in -xi direction
      case 160: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {4, 2}, {8, 3}, {5, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = -1.0;
        volIntPoint[1] = tmpVolIntPoint[0];
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] = -1.0;
        locNormal[1] =  0.0;
        locNormal[2] =  0.0;
        break;
      }
      // Surface[4,3,7,8] is common: back-side surface, normal in +eta direction
      case 672: {
        const std::map<UInt, UInt> commonIndexMap{{4, 1}, {3, 2}, {7, 3}, {8, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = 1.0;
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] =  0.0;
        locNormal[1] =  1.0;
        locNormal[2] =  0.0;
        break;
      }
      // Surface[1,2,6,5] is common: front-side surface, normal in -eta direction
      case 60: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {2, 2}, {6, 3}, {5, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = -1.0;
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] =  0.0;
        locNormal[1] = -1.0;
        locNormal[2] =  0.0;
        break;
      }
      default:
        EXCEPTION("HexaFE::GetLocalIntPoints4Surface: surface and volume element "
                  << "have not four nodes in common. Check your .mesh-file.");
    }
  }

  void FeH1LagrangeHex1::Triangulate(StdVector< StdVector<UInt> > & triConnect){
    //TODO: check orientation!
    triConnect.Resize(6);
    triConnect.Init(StdVector<UInt>(4));
    triConnect[0][0] = 0;
    triConnect[0][1] = 4;
    triConnect[0][2] = 1;
    triConnect[0][3] = 2;

    triConnect[1][0] = 1;
    triConnect[1][1] = 5;
    triConnect[1][2] = 2;
    triConnect[1][3] = 4;

    triConnect[2][0] = 0;
    triConnect[2][1] = 4;
    triConnect[2][2] = 2;
    triConnect[2][3] = 3;

    triConnect[3][0] = 2;
    triConnect[3][1] = 6;
    triConnect[3][2] = 3;
    triConnect[3][3] = 4;

    triConnect[4][0] = 3;
    triConnect[4][1] = 7;
    triConnect[4][2] = 4;
    triConnect[4][3] = 6;

    triConnect[5][0] = 2;
    triConnect[5][1] = 6;
    triConnect[5][2] = 4;
    triConnect[5][3] = 5;
  }

  
  // --- Wedge 1st order ---
  FeH1LagrangeWedge1::FeH1LagrangeWedge1() : FeH1LagrangeWedge() {
    feType_ = Elem::ET_WEDGE6;
    completeType_ = TENSOR_TYPE;
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
    // the surface element. Then calculate the sum of all four
    // and compare them
    //      + 6    
    //     /|\
    //    / |  \           zeta
    // 4 +----- + 5         ^  eta
    //   |  + 3 |           |/ 
    //   | / \  |           0--> xi
    //   |/    \|   
    // 1 +------+ 2
    //
    // Node coords: 1: (0,0,-1), 2: (1,0,-1), 3: (0,1,-1), 4: (0,0,1), 5: (1,0,1), 6: (0,1,1)

    // Check if surface element is triangle or quadrilateral
    UInt numNodes;
    if (surfConnect.GetSize() == 3 || surfConnect.GetSize() == 6)
      numNodes = 3;
    else
      numNodes = 4;

    StdVector<UInt> commonIndex(numNodes);
    UInt found = 0;
    UInt indexSum = 0;
    volIntPoint.coord.Resize(3);
    locNormal.Resize(3);

    // loop over surface connect
    for (UInt iSurf=0; iSurf<numNodes; iSurf++) {
      // loop over volume connect
      for (UInt iVol=0; iVol<6; iVol++) {
        if (surfConnect[iSurf] == volConnect[iVol]) {
          commonIndex[found] = iVol+1;
          indexSum +=  commonIndex[found++];
        }
      }
    }

    switch(indexSum) {
      // Surface[1,2,3] is common: downside triangular surface, normal in -zeta direction
      case 6: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {2, 2}, {3, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = tmpVolIntPoint[1];
        volIntPoint[2] = -1.0;
        locNormal[0] = 0.0;
        locNormal[1] = 0.0;
        locNormal[2] = -1.0;
        break;
      }
      // Surface[4,5,6] is common: upside triangular surfa, normal in +zeta direction
      case 15: {
        const std::map<UInt, UInt> commonIndexMap{{4, 1}, {5, 2}, {6, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = tmpVolIntPoint[1];
        volIntPoint[2] = 1.0;
        locNormal[0] = 0.0;
        locNormal[1] = 0.0;
        locNormal[2] = 1.0;
        break;
      }
      // Surface[2,3,5,6] is common: right-side quadrilateral surface, normal in +xi+eta direction
      case 16: {
        const std::map<UInt, UInt> commonIndexMap{{3, 1}, {2, 2}, {5, 3}, {6, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = 0.5 + tmpVolIntPoint[0] / 2.0;
        volIntPoint[1] = 0.5 - tmpVolIntPoint[0] / 2.0;
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] = sqrt(0.5);
        locNormal[1] = sqrt(0.5);
        locNormal[2] = 0.0;
        break;
      }
      // Surface[1,3,6,4] is common: left-side quadrilateral surface, normal in -xi direction
      case 14: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {3, 2}, {6, 3}, {4, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = 0.0;
        volIntPoint[1] = 0.5 + tmpVolIntPoint[0] / 2.0;
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] = -1.0;
        locNormal[1] = 0.0;
        locNormal[2] = 0.0;
        break;
      }
      // Surface[1,2,5,4] is common: back-side quadrilateral surface, normal in -eta direction
      case 12: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {2, 2}, {5, 3}, {4, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = 0.5 + tmpVolIntPoint[0] / 2.0;
        volIntPoint[1] = 0.0;
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] = 0.0;
        locNormal[1] = -1.0;
        locNormal[2] = 0.0;
        break;
      }
      default:
        EXCEPTION("FeH1LagrangeWedge::GetLocalIntPoints4Surface: surface "
            << "and volume element have not the expected nodes in common. "
            << "Check your mesh file.");
    }
  }
  
  void FeH1LagrangeWedge1::Triangulate(StdVector< StdVector<UInt> > & triConnect){

    triConnect.Resize(3);
    triConnect.Init(StdVector<UInt>(4));
    triConnect[0][0] = 0;
    triConnect[0][1] = 2;
    triConnect[0][2] = 1;
    triConnect[0][3] = 3;

    triConnect[1][0] = 1;
    triConnect[1][1] = 3;
    triConnect[1][2] = 5;
    triConnect[1][3] = 4;

    triConnect[2][0] = 1;
    triConnect[2][1] = 2;
    triConnect[2][2] = 5;
    triConnect[2][3] = 3;

  }

  // ========================================================================
  //  Lagrangian Elements of 2nd order
  // ========================================================================
  
  // --- Line 2nd order ---
  
  FeH1LagrangeLine2::FeH1LagrangeLine2() : FeH1LagrangeLine() {
    feType_ = Elem::ET_LINE3;
    completeType_ = TENSOR_TYPE;
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
    completeType_ = SERENDIPITY_TYPE;
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
    // From:
    // Zienkiewicz, The Finite Element Method.-  Vol 1, 5th ed., page 182
    
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

    deriv[0][0] =  4.0*point[0] + 4.0*point[1] - 3.0;
    deriv[0][1] =  4.0*point[0] + 4.0*point[1] - 3.0;

    deriv[1][0] =  4.0 * point[0] - 1.0;
    deriv[1][1] =  0;

    deriv[2][0] =  0;
    deriv[2][1] =  4.0 * point[1] - 1.0;

    deriv[3][0] =  4.0 * (1.0 - 2.0*point[0] - point[1]);
    deriv[3][1] = -4.0 * point[0];

    deriv[4][0] =  4.0 * point[1];
    deriv[4][1] =  4.0 * point[0];

    deriv[5][0] = -4.0 * point[1];
    deriv[5][1] =  4.0 * (1.0 - 2.0*point[1] - point[0]);
  }

  // --- Quad 2nd order ---
  FeH1LagrangeQuad2::FeH1LagrangeQuad2() : FeH1LagrangeQuad() {
    feType_ = Elem::ET_QUAD8;
    completeType_ = SERENDIPITY_TYPE;
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
    shape.Init();
    //intermediate Storage for important terms
    Double xProd = 0.0;
    Double yProd = 0.0;
    // From:
    // Zienkiewicz, The Finite Element Method.-  Vol 1, 5th ed., page 174
    
    // corner nodes
    for( UInt i = 0; i < 4; i++ ) {
      xProd = coords[i][0] * point[0];
      yProd = coords[i][1] * point[1];
      shape[i] = 0.25 * ( 1.0 + xProd ) * ( 1.0 + yProd ) * ( xProd + yProd - 1.0 );
    }
    
    // mid-side nodes
    for( UInt i = 4; i < 8; i = i + 2 ) {
      shape[i]   = 0.5 * ( 1.0 - point[0] * point[0] )
                       * ( 1.0 + coords[i][1] * point[1] );
      shape[i+1] = 0.5 * (1.0 - point[1] * point[1] )
                       * (1.0 + coords[i+1][0] * point[0] );
    }
  }
  
  void FeH1LagrangeQuad2::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                             const Vector<Double>& point,
                                             const Elem* ptElem,
                                             UInt comp ) {
    StdVector<Vector<Double> >& coords = shape_.nodeCoords;
    deriv.Resize( 8, 2 );
    deriv.Init();
    Double xProd = 0.0;
    Double yProd = 0.0;
    
    // corner nodes
    for( UInt i = 0; i < 4; i++ ) {
      xProd = coords[i][0] * point[0];
      yProd = coords[i][1] * point[1];
      deriv[i][0] = 0.25 * coords[i][0] * ( 1.0 + yProd ) * ( 2.0 * xProd + yProd );
      deriv[i][1] = 0.25 * coords[i][1] * ( 1.0 + xProd ) * ( 2.0 * yProd + xProd );
    }
      
    // mid-side nodes
    for( UInt i = 4; i < 8; i = i + 2 ) {
      deriv[i][0] = - point[0] * ( 1.0 + coords[i][1] * point[1] );
      deriv[i][1] =  0.5 * coords[i][1] * ( 1.0 - point[0] * point[0] );
      
      deriv[i+1][0] = 0.5 * coords[i+1][0] * ( 1.0 - point[1] * point[1] );
      deriv[i+1][1] = - point[1] * ( 1.0 + coords[i+1][0] * point[0] );
    }
  }
  
  // --- Quad 2nd order tensor product ---
  FeH1LagrangeQuad9::FeH1LagrangeQuad9() : FeH1LagrangeQuad(){
    feType_ = Elem::ET_QUAD9;
    completeType_ = TENSOR_TYPE;
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
    completeType_ = SERENDIPITY_TYPE;
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
    
    // From:
    // Zienkiewicz, The Finite Element Method.-  Vol 1, 5th ed., page 185
    UInt i;
    shape.Resize(20);

    //integration points
    const Double & xi   = point[0];
    const Double & eta  = point[1];
    const Double & zeta = point[2];

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
    completeType_ = TENSOR_TYPE;
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
     completeType_ = SERENDIPITY_TYPE;
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
     // From:
     // Zienkiewicz, The Finite Element Method.-  Vol 1, 5th ed., page 182
     
     // symbolically rewritten using sympy (cf. share/python/wedge15.py).
     const Double & xi = point[0];
     const Double & eta = point[1];
     const Double & zeta = point[2];
     
     shape.Resize( 15 );
     shape[0] = (-0.5*xi - 0.5*eta + 0.5)*(zeta*zeta + (zeta - 1)*(2*xi + 2*eta - 1) - 1);
     shape[1] = xi*(-1.0*xi*zeta + 1.0*xi + 0.5*zeta*zeta + 0.5*zeta - 1.0);
     shape[2] = eta*(-1.0*eta*zeta + 1.0*eta + 0.5*zeta*zeta + 0.5*zeta - 1.0);
     shape[3] = (xi + eta - 1)*(-0.5*zeta*zeta + 0.5*(zeta + 1)*(2*xi + 2*eta - 1) + 0.5);
     shape[4] = xi*(1.0*xi*zeta + 1.0*xi + 0.5*zeta*zeta - 0.5*zeta - 1.0);
     shape[5] = eta*(1.0*eta*zeta + 1.0*eta + 0.5*zeta*zeta - 0.5*zeta - 1.0);
     shape[6] = 2*xi*(zeta - 1)*(xi + eta - 1);
     shape[7] = 2*xi*eta*(-zeta + 1);
     shape[8] = 2*eta*(zeta - 1)*(xi + eta - 1);
     shape[9] = -2*xi*(zeta + 1)*(xi + eta - 1);
     shape[10] = 2*xi*eta*(zeta + 1);
     shape[11] = -2*eta*(zeta + 1)*(xi + eta - 1);
     shape[12] = (zeta*zeta - 1)*(xi + eta - 1);
     shape[13] = xi*(-zeta*zeta + 1);
     shape[14] = eta*(-zeta*zeta + 1);
  }

    void FeH1LagrangeWedge2::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                                const Vector<Double>& point,
                                                const Elem* ptElem,
                                                UInt comp )
    {
      const Double & xi = point[0];
      const Double & eta = point[1];
      const Double & zeta = point[2];
      deriv.Resize(15,3);
      deriv.Init();
      // symbolically solved using sympy (cf. share/python/wedge15.py).
      
      deriv[0][0] = -2.0*xi*zeta + 2.0*xi - 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta + 1.5*zeta - 1.0;
      deriv[0][1] = -2.0*xi*zeta + 2.0*xi - 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta + 1.5*zeta - 1.0;
      deriv[0][2] = (xi + eta - 1)*(-1.0*xi - 1.0*eta - 1.0*zeta + 0.5);

      deriv[1][0] = -2.0*xi*zeta + 2.0*xi + 0.5*zeta*zeta + 0.5*zeta - 1.0;
      deriv[1][1] = 0;
      deriv[1][2] = xi*(-1.0*xi + 1.0*zeta + 0.5);

      deriv[2][0] = 0;
      deriv[2][1] = -2.0*eta*zeta + 2.0*eta + 0.5*zeta*zeta + 0.5*zeta - 1.0;
      deriv[2][2] = eta*(-1.0*eta + 1.0*zeta + 0.5);

      deriv[3][0] = 2.0*xi*zeta + 2.0*xi + 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta - 1.5*zeta - 1.0;
      deriv[3][1] = 2.0*xi*zeta + 2.0*xi + 2.0*eta*zeta + 2.0*eta - 0.5*zeta*zeta - 1.5*zeta - 1.0;
      deriv[3][2] = (xi + eta - 1)*(1.0*xi + 1.0*eta - 1.0*zeta - 0.5);

      deriv[4][0] = 2.0*xi*zeta + 2.0*xi + 0.5*zeta*zeta - 0.5*zeta - 1.0;
      deriv[4][1] = 0;
      deriv[4][2] = xi*(1.0*xi + 1.0*zeta - 0.5);

      deriv[5][0] = 0;
      deriv[5][1] = 2.0*eta*zeta + 2.0*eta + 0.5*zeta*zeta - 0.5*zeta - 1.0;
      deriv[5][2] = eta*(1.0*eta + 1.0*zeta - 0.5);

      deriv[6][0] = 4*xi*zeta - 4*xi + 2*eta*zeta - 2*eta - 2*zeta + 2;
      deriv[6][1] = 2*xi*(zeta - 1);
      deriv[6][2] = 2*xi*(xi + eta - 1);

      deriv[7][0] = 2*eta*(-zeta + 1);
      deriv[7][1] = 2*xi*(-zeta + 1);
      deriv[7][2] = -2*xi*eta;

      deriv[8][0] = 2*eta*(zeta - 1);
      deriv[8][1] = 2*xi*zeta - 2*xi + 4*eta*zeta - 4*eta - 2*zeta + 2;
      deriv[8][2] = 2*eta*(xi + eta - 1);

      deriv[9][0] = -4*xi*zeta - 4*xi - 2*eta*zeta - 2*eta + 2*zeta + 2;
      deriv[9][1] = 2*xi*(-zeta - 1);
      deriv[9][2] = 2*xi*(-xi - eta + 1);

      deriv[10][0] = 2*eta*(zeta + 1);
      deriv[10][1] = 2*xi*(zeta + 1);
      deriv[10][2] = 2*xi*eta;

      deriv[11][0] = 2*eta*(-zeta - 1);
      deriv[11][1] = -2*xi*zeta - 2*xi - 4*eta*zeta - 4*eta + 2*zeta + 2;
      deriv[11][2] = 2*eta*(-xi - eta + 1);

      deriv[12][0] = zeta*zeta - 1;
      deriv[12][1] = zeta*zeta - 1;
      deriv[12][2] = 2*zeta*(xi + eta - 1);

      deriv[13][0] = -zeta*zeta + 1;
      deriv[13][1] = 0;
      deriv[13][2] = -2*xi*zeta;

      deriv[14][0] = 0;
      deriv[14][1] = -zeta*zeta + 1;
      deriv[14][2] = -2*eta*zeta;
    }
    
   


   // --- Complete wedge 2nd order ---
    FeH1LagrangeWedge18::FeH1LagrangeWedge18() : FeH1LagrangeWedge() {
      feType_ = Elem::ET_WEDGE18;
      completeType_ = TENSOR_TYPE;
      shape_ = Elem::shapes[feType_];
      actNumFncs_ = 18;
      order_ = 2;
    }

    FeH1LagrangeWedge18::~FeH1LagrangeWedge18() {

    }

    void FeH1LagrangeWedge18::CalcShFnc( Vector<Double>& shape,
                                         const Vector<Double>& point,
                                         const Elem* ptElem,
                                         UInt comp )
    {
      const Double& xi   = point[0];
      const Double& eta  = point[1];
      const Double& zeta = point[2];

      shape.Resize(18);
      
      // These are the standard tensor product basis function obtained
      // from the TRIA6 element  in the xi-eta plane and the LINE3
      // element in zeta-direction. The shape functions have been verified
      // by visualizing them using a calculator filter in ParaView.
      // After that, they have been symbolically rewritten using sympy.
      // (cf. share/python/wedge18.py)

      // Corners
      // shape[0] = -0.5* t * (2*t - 1)* coordsZ * (1 - coordsZ);
      // shape[1] = -0.5* coordsX * (2*coordsX - 1)* coordsZ * (1 - coordsZ);
      // shape[2] = -0.5* coordsY * (2*coordsY - 1)* coordsZ * (1 - coordsZ);      
      // shape[3] = 0.5* t * (2*t - 1)* coordsZ * (1 + coordsZ);
      // shape[4] = 0.5* coordsX * (2*coordsX - 1)* coordsZ * (1 + coordsZ);
      // shape[5] = 0.5* coordsY * (2*coordsY - 1)* coordsZ * (1 + coordsZ);
      shape[0] = 0.5*zeta*(zeta - 1)*(xi + eta - 1)*(2*xi + 2*eta - 1);
      shape[1] = 0.5*xi*zeta*(2*xi - 1)*(zeta - 1);
      shape[2] = 0.5*eta*zeta*(2*eta - 1)*(zeta - 1);
      shape[3] = 0.5*zeta*(zeta + 1)*(xi + eta - 1)*(2*xi + 2*eta - 1);
      shape[4] = 0.5*xi*zeta*(2*xi - 1)*(zeta + 1);
      shape[5] = 0.5*eta*zeta*(2*eta - 1)*(zeta + 1);

      // Mid-sides of quadratic triangles
      // shape[6] = -2*coordsX*t*coordsZ*(1-coordsZ);
      // shape[7] = -2*coordsX*coordsY*coordsZ*(1-coordsZ);
      // shape[8] = -2*coordsY*t*coordsZ*(1-coordsZ);
      // shape[9] = 2*coordsX*t*coordsZ*(1+coordsZ);
      // shape[10] = 2*coordsX*coordsY*coordsZ*(1+coordsZ);
      // shape[11] = 2*coordsY*t*coordsZ*(1+coordsZ);
      shape[6] = -2*xi*zeta*(zeta - 1)*(xi + eta - 1);
      shape[7] = 2*xi*eta*zeta*(zeta - 1);
      shape[8] = -2*eta*zeta*(zeta - 1)*(xi + eta - 1);
      shape[9] = -2*xi*zeta*(zeta + 1)*(xi + eta - 1);
      shape[10] = 2*xi*eta*zeta*(zeta + 1);
      shape[11] = -2*eta*zeta*(zeta + 1)*(xi + eta - 1);

      // Mid-sides of edges between the two triangles
      // shape[12] = t*(2*t-1)*(1+coordsZ)*(1-coordsZ);
      // shape[13] = coordsX*(2*coordsX-1)*(1+coordsZ)*(1-coordsZ);
      // shape[14] = coordsY*(2*coordsY-1)*(1+coordsZ)*(1-coordsZ);
      shape[12] = -(zeta - 1)*(zeta + 1)*(xi + eta - 1)*(2*xi + 2*eta - 1);
      shape[13] = -xi*(2*xi - 1)*(zeta - 1)*(zeta + 1);
      shape[14] = -eta*(2*eta - 1)*(zeta - 1)*(zeta + 1);

      // Centerpoints of the biquadratic quads
      // shape[15] = 4*coordsX*t*(1+coordsZ)*(1-coordsZ)
      // shape[16] = 4*coordsX*coordsY*(1+coordsZ)*(1-coordsZ)
      // shape[17] = 4*coordsY*t*(1+coordsZ)*(1-coordsZ)
      shape[15] = 4*xi*(zeta - 1)*(zeta + 1)*(xi + eta - 1);
      shape[16] = 4*xi*eta*(-zeta*zeta + 1);
      shape[17] = 4*eta*(zeta - 1)*(zeta + 1)*(xi + eta - 1);
    }

    void FeH1LagrangeWedge18::CalcLocDerivShFnc( Matrix<Double> & deriv,
                                                 const Vector<Double>& point,
                                                 const Elem* ptElem,
                                                 UInt comp )
    {
      const Double& xi   = point[0];
      const Double& eta  = point[1];
      const Double& zeta = point[2];

      deriv.Resize(18, 3);

      // symbolically solved using sympy (cf. share/python/wedge18.py).

      deriv[0][0] = zeta*(2.0*xi*zeta - 2.0*xi + 2.0*eta*zeta - 2.0*eta - 1.5*zeta + 1.5);
      deriv[0][1] = zeta*(2.0*xi*zeta - 2.0*xi + 2.0*eta*zeta - 2.0*eta - 1.5*zeta + 1.5);
      deriv[0][2] = 0.5*(2*zeta - 1)*(xi + eta - 1)*(2*xi + 2*eta - 1);
      
      deriv[1][0] = zeta*(2.0*xi*zeta - 2.0*xi - 0.5*zeta + 0.5);
      deriv[1][1] = 0;
      deriv[1][2] = xi*(2.0*xi*zeta - 1.0*xi - 1.0*zeta + 0.5);
      
      deriv[2][0] = 0;
      deriv[2][1] = zeta*(2.0*eta*zeta - 2.0*eta - 0.5*zeta + 0.5);
      deriv[2][2] = eta*(2.0*eta*zeta - 1.0*eta - 1.0*zeta + 0.5);
      
      deriv[3][0] = zeta*(2.0*xi*zeta + 2.0*xi + 2.0*eta*zeta + 2.0*eta - 1.5*zeta - 1.5);
      deriv[3][1] = zeta*(2.0*xi*zeta + 2.0*xi + 2.0*eta*zeta + 2.0*eta - 1.5*zeta - 1.5);
      deriv[3][2] = 0.5*(2*zeta + 1)*(xi + eta - 1)*(2*xi + 2*eta - 1);
      
      deriv[4][0] = zeta*(2.0*xi*zeta + 2.0*xi - 0.5*zeta - 0.5);
      deriv[4][1] = 0;
      deriv[4][2] = xi*(2.0*xi*zeta + 1.0*xi - 1.0*zeta - 0.5);
      
      deriv[5][0] = 0;
      deriv[5][1] = zeta*(2.0*eta*zeta + 2.0*eta - 0.5*zeta - 0.5);
      deriv[5][2] = eta*(2.0*eta*zeta + 1.0*eta - 1.0*zeta - 0.5);
      
      deriv[6][0] = 2*zeta*(-2*xi*zeta + 2*xi - eta*zeta + eta + zeta - 1);
      deriv[6][1] = 2*xi*zeta*(-zeta + 1);
      deriv[6][2] = 2*xi*(-2*xi*zeta + xi - 2*eta*zeta + eta + 2*zeta - 1);
      
      deriv[7][0] = 2*eta*zeta*(zeta - 1);
      deriv[7][1] = 2*xi*zeta*(zeta - 1);
      deriv[7][2] = 2*xi*eta*(2*zeta - 1);
      
      deriv[8][0] = 2*eta*zeta*(-zeta + 1);
      deriv[8][1] = 2*zeta*(-xi*zeta + xi - 2*eta*zeta + 2*eta + zeta - 1);
      deriv[8][2] = 2*eta*(-2*xi*zeta + xi - 2*eta*zeta + eta + 2*zeta - 1);
      
      deriv[9][0] = 2*zeta*(-2*xi*zeta - 2*xi - eta*zeta - eta + zeta + 1);
      deriv[9][1] = 2*xi*zeta*(-zeta - 1);
      deriv[9][2] = 2*xi*(-2*xi*zeta - xi - 2*eta*zeta - eta + 2*zeta + 1);
      
      deriv[10][0] = 2*eta*zeta*(zeta + 1);
      deriv[10][1] = 2*xi*zeta*(zeta + 1);
      deriv[10][2] = 2*xi*eta*(2*zeta + 1);
      
      deriv[11][0] = 2*eta*zeta*(-zeta - 1);
      deriv[11][1] = 2*zeta*(-xi*zeta - xi - 2*eta*zeta - 2*eta + zeta + 1);
      deriv[11][2] = 2*eta*(-2*xi*zeta - xi - 2*eta*zeta - eta + 2*zeta + 1);
      
      deriv[12][0] = -4*xi*zeta*zeta + 4*xi - 4*eta*zeta*zeta + 4*eta + 3*zeta*zeta - 3;
      deriv[12][1] = -4*xi*zeta*zeta + 4*xi - 4*eta*zeta*zeta + 4*eta + 3*zeta*zeta - 3;
      deriv[12][2] = 2*zeta*(-2*xi*xi - 4*xi*eta + 3*xi - 2*eta*eta + 3*eta - 1);
      
      deriv[13][0] = -4*xi*zeta*zeta + 4*xi + zeta*zeta - 1;
      deriv[13][1] = 0;
      deriv[13][2] = 2*xi*zeta*(-2*xi + 1);

      deriv[14][0] = 0;
      deriv[14][1] = -4*eta*zeta*zeta + 4*eta + zeta*zeta - 1;
      deriv[14][2] = 2*eta*zeta*(-2*eta + 1);
      
      deriv[15][0] = 8*xi*zeta*zeta - 8*xi + 4*eta*zeta*zeta - 4*eta - 4*zeta*zeta + 4;
      deriv[15][1] = 4*xi*(zeta*zeta - 1);
      deriv[15][2] = 8*xi*zeta*(xi + eta - 1);
      
      deriv[16][0] = 4*eta*(-zeta*zeta + 1);
      deriv[16][1] = 4*xi*(-zeta*zeta + 1);
      deriv[16][2] = -8*xi*eta*zeta;
      
      deriv[17][0] = 4*eta*(zeta*zeta - 1);
      deriv[17][1] = 4*xi*zeta*zeta - 4*xi + 8*eta*zeta*zeta - 8*eta - 4*zeta*zeta + 4;
      deriv[17][2] = 8*eta*zeta*(xi + eta - 1);
      
    }

   // --- Tetra 1st order ---
    FeH1LagrangeTet1::FeH1LagrangeTet1() : FeH1LagrangeTet() {
      feType_ = Elem::ET_TET4;
      completeType_ = TENSOR_TYPE;
      shape_ = Elem::shapes[feType_];
      actNumFncs_ = 4;
      order_ = 1;
    }

    FeH1LagrangeTet1::~FeH1LagrangeTet1() {

    }

    void FeH1LagrangeTet1::ComputeMonomialCoefficients(Matrix<Integer>& P, Matrix<Double>& C){
      P.Resize(actNumFncs_,3);
      P.Init();
      C.Resize(actNumFncs_,actNumFncs_);
      C.Init();
      P[1][0] = 1;
      P[2][1] = 1;
      P[3][2] = 1;

      C[0][0] = 1;
      C[0][1] = -1;
      C[0][2] = -1;
      C[0][3] = -1;
      C[1][1] = 1;
      C[2][2] = 1;
      C[3][3] = 1;

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
      completeType_ = SERENDIPITY_TYPE;
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
    //  |\\         zeta
    //  | \ \       ^ eta
    //  |  \ +3     |/
    //  |   \ |     0--> xi
    //  |    \ \
    //  |     \|     REFERENCE TETRAHEDRAL ELEMENT
    //  +------+
    //  1      2
    //
    // Node coords: 1:(0,0,0), 2:(1,0,0), 3:(0,1,0), 4:(0,0,1)

    StdVector<UInt> commonIndex(3);
    UInt found = 0;
    UInt indexProduct = 0;
    volIntPoint.coord.Resize(3);
    locNormal.Resize(3);

    // loop over surface connect
    for (UInt iSurf=0; iSurf<3; iSurf++) {
      // loop over volume connect
      for (UInt iVol=0; iVol<4; iVol++) {
        if (surfConnect[iSurf] == volConnect[iVol])
          commonIndex[found++] = iVol+1;
      }
    }

    indexProduct =  commonIndex[0] * commonIndex[1] * commonIndex[2];

    switch(indexProduct) {
      // Surface [1,2,4] is common: front surface, normal in -eta direction
      case 8: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {2, 2}, {4, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = 0.0;
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] =  0.0;
        locNormal[1] = -1.0;
        locNormal[2] =  0.0;
        break;
      }
      // Surface [4,2,3] is common: upper-right surface, normal in +xi+eta+zeta direction
      case 24: {
        // use mapping algorithm of the standard tria element onto a 2D equilateral triangle with coords 1:(-1,0), 2:(1,0), 3:(0,sqrt(3))
        const std::map<UInt, UInt> commonIndexMap{{4, 1}, {2, 2}, {3, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint, true);
        // map back to 3D space
        volIntPoint[0] = (1 + tmpVolIntPoint[0]) / 2.0 - tmpVolIntPoint[1] / sqrt(12.0);
        volIntPoint[1] = tmpVolIntPoint[1] / sqrt(3.0);
        volIntPoint[2] = (1 - tmpVolIntPoint[0]) / 2.0 - tmpVolIntPoint[1] / sqrt(12.0);
        locNormal[0] =  1.0/sqrt(3.0);
        locNormal[1] =  1.0/sqrt(3.0);
        locNormal[2] =  1.0/sqrt(3.0);
        break;
      }
      // Surface[1,3,4] is common: upper-left surface, normal in -xi direction
      case 12: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {3, 2}, {4, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = 0.0;
        volIntPoint[1] = tmpVolIntPoint[0];
        volIntPoint[2] = tmpVolIntPoint[1];
        locNormal[0] = -1.0;
        locNormal[1] =  0.0;
        locNormal[2] =  0.0;
        break;
      }
      // Surface[1,2,3] is common: bottom surface, normal in -zeta direction
      case 6: {
        const std::map<UInt, UInt> commonIndexMap{{1, 1}, {2, 2}, {3, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = tmpVolIntPoint[1];
        volIntPoint[2] = 0.0;
        locNormal[0] =  0.0;
        locNormal[1] =  0.0;
        locNormal[2] = -1.0;
        break;
      }
      default:
        EXCEPTION("FeH1LagrangeTet::GetLocalIntPoints4Surface: surface "
            << "and volume element have not three nodes in common. "
            << "Check your .mesh-file.");
    }
  }

    void FeH1LagrangeTet1::Triangulate(StdVector< StdVector<UInt> > & triConnect){

      triConnect.Resize(1);
      triConnect.Init(StdVector<UInt>(4));

      triConnect[0][0] = 0;
      triConnect[0][1] = 1;
      triConnect[0][2] = 2;
      triConnect[0][3] = 3;
    }

    // --- Pyramid 1st order ---
     FeH1LagrangePyra1::FeH1LagrangePyra1() : FeH1LagrangePyra() {
       feType_ = Elem::ET_PYRA5;
       completeType_ = TENSOR_TYPE;
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
       completeType_ = SERENDIPITY_TYPE;
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

     void FeH1LagrangePyra1::Triangulate(StdVector< StdVector<UInt> > & triConnect){

       triConnect.Resize(2);
       triConnect.Init(StdVector<UInt>(4));

       triConnect[0][0] = 0;
       triConnect[0][1] = 1;
       triConnect[0][2] = 2;
       triConnect[0][3] = 4;

       triConnect[1][0] = 0;
       triConnect[1][1] = 2;
       triConnect[1][2] = 3;
       triConnect[1][3] = 4;

     }

     // --- Pyra 2nd order ---
      FeH1LagrangePyra14::FeH1LagrangePyra14() : FeH1LagrangePyra() {
        feType_ = Elem::ET_PYRA14;
        completeType_ = TENSOR_TYPE;
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
        //WARN("CalcShFnc for ET_PYRA14 implemented but not tested!");

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
        //WARN("CalcLocDerivShFnc for ET_PYRA14 implemented but not tested!");
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
    // Node coords: 1: (1,1,0), 2: (-1,1,0), 3: (-1,-1,0), 4: (1,-1,0), 5: (0,0,1)

    // Check if surface element is triangle or quadrilateral
    UInt numNodes;
    if (surfConnect.GetSize() == 3 || surfConnect.GetSize() == 6)
      numNodes = 3;
    else
      numNodes = 4;

    StdVector<UInt> commonIndex(numNodes);
    UInt found = 0;
    UInt indexProduct = 1;
    volIntPoint.coord.Resize(3);
    locNormal.Resize(3);

    // loop over surface connect
    for (UInt iSurf=0; iSurf<numNodes; iSurf++) {
      // loop over volume connect
      for (UInt iVol=0; iVol<5; iVol++) {
        if (surfConnect[iSurf] == volConnect[iVol]) {
          commonIndex[found] = iVol+1;
          indexProduct *=  commonIndex[found++];
        }
      }
    }

    switch(indexProduct) {
      // Surface[2,1,5] is common: rear triangular surface, normal in +eta+zeta direction
      case 10: {
        const std::map<UInt, UInt> commonIndexMap{{2, 1}, {1, 2}, {5, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint, true);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = 1.0 - tmpVolIntPoint[1] / sqrt(3.0);
        volIntPoint[2] = tmpVolIntPoint[1] / sqrt(3.0);
        locNormal[0] =  0;
        locNormal[1] = sqrt(0.5);
        locNormal[2] = sqrt(0.5);
        break;
      }
      // Surface[2,3,5] is common: left triangular surface , normal in -xi+zeta direction
      case 30: {
        const std::map<UInt, UInt> commonIndexMap{{2, 1}, {3, 2}, {5, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint, true);
        volIntPoint[0] = -1.0 + tmpVolIntPoint[1] / sqrt(3.0);
        volIntPoint[1] = -tmpVolIntPoint[0];
        volIntPoint[2] = tmpVolIntPoint[1] / sqrt(3.0);
        locNormal[0] = -sqrt(0.5);
        locNormal[1] = 0.0;
        locNormal[2] = sqrt(0.5);
        break;
      }
      // Surface[3,4,5] is common: front triangular surface, normal in -eta+zeta direction
      case 60: {
        const std::map<UInt, UInt> commonIndexMap{{3, 1}, {4, 2}, {5, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint, true);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = -1.0 + tmpVolIntPoint[1] / sqrt(3.0);
        volIntPoint[2] = tmpVolIntPoint[1] / sqrt(3.0);
        locNormal[0] = 0.0;
        locNormal[1] = -sqrt(0.5);
        locNormal[2] = sqrt(0.5);
        break;
      }
      // Surface[4,1,5] is common: right triangular surface, normal in +xi+zeta direction
      case 20: {
        const std::map<UInt, UInt> commonIndexMap{{4, 1}, {1, 2}, {5, 3}};
        StdVector<Double> tmpVolIntPoint(2);
        MapTriaSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint, true);
        volIntPoint[0] = 1.0 - tmpVolIntPoint[1] / sqrt(3.0);
        volIntPoint[1] = tmpVolIntPoint[0];
        volIntPoint[2] = tmpVolIntPoint[1] / sqrt(3.0);
        locNormal[0] = sqrt(0.5);
        locNormal[1] =  0;
        locNormal[2] = sqrt(0.5);
        break;
      }
      // Surface[3,4,1,2] is common: bottom quadrilateral surface, normal in zeta direction
      case 24: {
        const std::map<UInt, UInt> commonIndexMap{{3, 1}, {4, 2}, {1, 3}, {2, 4}};
        StdVector<Double> tmpVolIntPoint(2);
        MapQuadSurfOrientation(commonIndexMap, commonIndex, surfIntPoint, tmpVolIntPoint);
        volIntPoint[0] = tmpVolIntPoint[0];
        volIntPoint[1] = tmpVolIntPoint[1];
        volIntPoint[2] = 0.0;
        locNormal[0] =  0;
        locNormal[1] =  0;
        locNormal[2] = -1;
        break;
      }
      default:
        EXCEPTION("FeH1LagrangePyra::GetLocalIntPoints4Surface: surface "
            << "and volume element have not the expected nodes in common. "
            << "Check your mesh file.");
    }
  }
} // namespace CoupledField
