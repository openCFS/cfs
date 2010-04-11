// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "H1ElemsHi.hh"

namespace CoupledField {


  // ========================================================================
  //  FeH1Hi
  // ========================================================================
  
  FeH1Hi::FeH1Hi() {
  }
  
  FeH1Hi::~FeH1Hi() {
  }
  
  void FeH1Hi::GetNumFncs( StdVector<UInt>& numFcns,
                       EntityType entityType,
                       UInt dof ) {
    // this can be implemented in a general way, as the information about
    // the number of function is already stored in entityFncs_
    
    // Note: at the moment we ignore the dof parameter
    
    // numFcns has the length Vertices / Edges / Faces / Inner degrees
    numFcns = entityFncs_[entityType]; 
    
  }

  //void FeH1Hi::EvalPolynom( UInt p, Vector<Double>& vals ) {
    // To be implemented using the recursion formula
  //}
  
  void FeH1Hi::GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                    const Elem* ptElem,
                                    EntityType fctEntityType,
                                    UInt entNumber){
   
    if( fctEntityType == VERTEX ) {
      UInt numFncs = entityFncs_[VERTEX][entNumber];
      fncPermutation.Resize( numFncs );
      for(UInt i = 0; i < numFncs; ++i ) {
        fncPermutation[i] = i;
      }
    }else if( fctEntityType == EDGE ) {
      UInt numFncs = entityFncs_[EDGE][entNumber];
      fncPermutation.Resize( numFncs );
      for(UInt i = 0; i < numFncs; ++i ) {
        fncPermutation[i] = i;
      }
    }else if( fctEntityType == FACE ) {
      UInt numFncs = entityFncs_[FACE][entNumber];
      fncPermutation.Resize( numFncs );
      for(UInt i = 0; i < numFncs; ++i ) {
        fncPermutation[i] = i;
      }
    }else if( fctEntityType == INTERIOR ) {
      UInt numFncs = entityFncs_[INTERIOR][entNumber];
      fncPermutation.Resize( numFncs );
      for(UInt i = 0; i < numFncs; ++i ) {
        fncPermutation[i] = i;
      }
    }
  }

  void FeH1Hi::EvalPolynom( Double& value, Double& deriv,
                            const UInt order, const Double* coeff,
                            const Double xVal ) {

    // Consider the following expression
    // f(xVal) = a0 * (a1*x^order + a2*x^(order-1) + .. + a(order+1))
    // The coefficients a0..a(order+1) are stored in the coeff-array

    value = coeff[1];
    deriv = 0.0;
    for( UInt i = 2; i < order+2; i++ ) {
      deriv = deriv * xVal + value;
      value = value * xVal + coeff[i];
    }
    // Multiply by pre-factor
    deriv *= coeff[0];
    value *= coeff[0];
  }


    // Define coefficients for legendre ansatz functions up to order 8
    Double  FeH1Hi::lCoeff_[9][10] = {
      {0.5                  ,   -1, 1,    0, 0,   0, 0,    0, 0, 0 },
      {0.5                  ,    1, 1,    0, 0,   0, 0,    0, 0, 0 },
      {0.25*sqrt(6.0)       ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
      {0.25*sqrt(10.0)      ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
      {1.0/16.0*sqrt(14.0)  ,    5, 0,   -6, 0,   1, 0,    0, 0, 0 },
      {3.0/16.0*sqrt(2.0)   ,    7, 0,  -10, 0,   3, 0,    0, 0, 0 },
      {1.0/32.0*sqrt(22.0)  ,   21, 0,  -35, 0,  15, 0,   -1, 0, 0 },
      {1.0/32.0*sqrt(26.0)  ,   33, 0,  -63, 0,  35, 0,   -5, 0, 0 },
      {1.0/256.0*sqrt(30.0) ,  429, 0, -924, 0, 630, 0, -140, 0, 5 }
    };
  
  // ========================================================================
  //  FeH1Hi explicit element definition 
  // ========================================================================
  
  FeH1HiLine::FeH1HiLine() {
    feType_ = Elem::ET_LINE2;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 0; // Gets re-setted afterwards
    order_ = 0;
  }
  FeH1HiLine::~FeH1HiLine() {

  }
  void FeH1HiLine::SetIsoOrder(UInt order) {
    
    actNumFncs_ = 0;
    
    // Vertices
    StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
    vertFncs.Resize(shape_.numVertices);
    vertFncs.Init(1); // Vertices have always order 1
    actNumFncs_ += shape_.numVertices * 1;
    
    // Edges
    StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
    edgeFncs.Resize(shape_.numEdges);
    edgeFncs.Init(order-1);
    actNumFncs_ += shape_.numEdges * (order-1);
    
    // Interior
    StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
    innerFncs.Resize(1);
    innerFncs.Init(0);
    
  }
  
  void FeH1HiLine::GetShFnc( Vector<Double>& shape, const LocPoint& lp,
                             const Elem * elem, UInt comp ) {
  }

  void FeH1HiLine::GetDerivShFnc( Matrix<Double> & deriv, const LocPoint& lp,
                                  const Elem * elem,  UInt comp ) {

  }
  
  // --- Quad 1st order ---
   
  FeH1HiQuad::FeH1HiQuad() {
    feType_ = Elem::ET_QUAD4;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 0;
    order_ = 0;
  }
    
  FeH1HiQuad::~FeH1HiQuad() {
    
  }
  
  void FeH1HiQuad::SetIsoOrder(UInt order) {
    
    actNumFncs_ = 0;

    // Vertices
    StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
    vertFncs.Resize(shape_.numVertices);
    vertFncs.Init(1); // Vertices have always order 1
    actNumFncs_ += shape_.numVertices * 1;

    // Edges
    StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
    edgeFncs.Resize(shape_.numEdges);
    edgeFncs.Init(order-1);
    actNumFncs_ += shape_.numEdges * (order-1);
    
    StdVector<UInt>& edgeOrder = entityOrder_[EDGE];
    edgeOrder.Resize(shape_.numEdges);
    edgeOrder.Init(order);
    
    // Faces
    StdVector<UInt>& faceFncs = entityFncs_[FACE];
    faceFncs.Resize(shape_.numFaces);
    faceFncs.Init( (order-1)*(order-1) );
    actNumFncs_ += shape_.numFaces * (order-1)*(order-1) ;
    
    // Interior
    StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
    innerFncs.Resize(1);
    innerFncs.Init(0);
  }
  
  void FeH1HiQuad::GetShFnc( Vector<Double>& shape, const LocPoint& lp,
                             const Elem * elem, UInt comp ) {
    
    // create nodal shape functions
    shape.Resize( actNumFncs_ );
    const Vector<Double> & point = lp.coord;
    // --------------------
    //  a) nodal functions
    // --------------------
    shape[0] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ); 
    shape[1] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] );
    shape[2] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] );
    shape[3] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] );
    
    UInt offset = 4;
    
    // ToDo: The following section is just temporary, as we want to replace
    // it by the recursive calculation of shape functions taken from ngsolve
    // --------------------
    //  b) edge functions
    // --------------------
    // Obtain order of element
    Integer order;
    Double val, factor, deriv;

#define QUAD_EDGE_FCN(edgeNum,  sign_1, dir_1, dir_2 )                  \
    factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
    for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {        \
      EvalPolynom( val, deriv, iDof, lCoeff_[iDof],                   \
                   factor*point[dir_2] );                          \
                   shape[offset] = 0.5 * ( 1 sign_1 point[dir_1] ) * val;       \
    }

    //  EDGE #1
    order = entityOrder_[EDGE][0];
    QUAD_EDGE_FCN( 1, -, 1, 0 );

    //  EDGE #2
    order = entityOrder_[EDGE][1];
    QUAD_EDGE_FCN( 2, +, 0, 1 );

    //  EDGE #3
    order = entityOrder_[EDGE][2];
    QUAD_EDGE_FCN( 3, +, 1, 0 );

    //  EDGE #4
    order = entityOrder_[EDGE][3];
    QUAD_EDGE_FCN( 4, -, 0, 1 );
    
    // --------------------
    //  c) face functions
    // --------------------
    Double val_2, deriv_2, val_3, deriv_3;
    Double sFlag2, sFlag3;
    Integer order_2, order_3;
    UInt d2, d3;

    order_2 = entityOrder_[EDGE][0];
    order_3 = entityOrder_[EDGE][1];
    if( elem->faceFlags[0].test(0) == true) {
      sFlag2 = (elem->faceFlags[0].
          test(2) == true) ? 1.0 : -1.0;
      sFlag3 = (elem->faceFlags[0].
          test(1) == true) ? 1.0 : -1.0;
      d2 = 0;
      d3 = 1;
    } else {
      sFlag3 = (elem->faceFlags[0].test(2) == true) ? 1.0 : -1.0;
      sFlag2 = (elem->faceFlags[0].test(1) == true) ? 1.0 : -1.0;
      d3 = 0;
      d2 = 1;
    }
    for( Integer i = 2; i <= (order_2)-2; i++ ) {
      for( Integer j = 2; j <= (order_3)-2; j++ ) {
        if( (i + j) > std::max(order_2, order_3) ) {
          continue;
        }
        EvalPolynom( val_2, deriv_2, i, lCoeff_[i],
                     sFlag2* point[d2] );
        EvalPolynom( val_3, deriv_3, j, lCoeff_[j],
                     sFlag3*point[d3] );
        shape[offset] = val_2  * val_3;
        offset++;
      }
    }

  }

  void FeH1HiQuad::GetDerivShFnc( Matrix<Double> & deriv, const LocPoint& lp,
                                  const Elem * elem,  UInt comp )  {
    // nodal shape functions
    
   const Vector<Double> & point = lp.coord;
    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
    deriv.Resize( actNumFncs_, 2 );
    for( UInt i = 0; i < 4; i++ ) {
      deriv[i][0] = 0.25 * coords[i][0] * (1 + coords[i][1] * point[1] );
      deriv[i][1] = 0.25 * (1 + coords[i][0] * point[0] ) * coords[i][1];
    }
    
    // -------------------
    //  b) edge functions
    // -------------------
    UInt offset = 4;

    // Obtain order of element
    Integer order_2, order_3;
    Double val, deriv1D, factor;
#define QUAD_EDGE_DERIV( edgeNum, sign_1, dir_1, dir_2 )                \
    factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
    for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {        \
      EvalPolynom( val, deriv1D, iDof, lCoeff_[iDof],                   \
                   factor*point[dir_2] );                          \
                   deriv[offset][dir_1] = sign_1(0.5 * val);                      \
                   deriv[offset][dir_2] =  0.5 * ( 1.0 sign_1 point[dir_1])    \
                   * deriv1D * factor;                  \
    }

    // EDGE #1
    UInt order = entityOrder_[EDGE][0];
    QUAD_EDGE_DERIV( 1, -, 1, 0 );

    // EDGE #2
    order = entityOrder_[EDGE][1];
    QUAD_EDGE_DERIV( 2, +, 0, 1 );

    // EDGE #3
    order = entityOrder_[EDGE][2];
    QUAD_EDGE_DERIV( 3, +, 1, 0 );

    // EDGE #4
    order = entityOrder_[EDGE][3];
    QUAD_EDGE_DERIV( 4, -, 0, 1 );

    // -------------------
    //  c) face functions
    // -------------------

    Double val_2, deriv_2, val_3, deriv_3;
    Double sFlag2, sFlag3;
    UInt d2, d3;

    order_2 = order = entityOrder_[EDGE][1];
    order_3 = order = entityOrder_[EDGE][2];
    if( elem->faceFlags[0].test(0) == true) {
      sFlag2 = (elem->faceFlags[0].
          test(2) == true) ? 1.0 : -1.0;
      sFlag3 = (elem->faceFlags[0].
          test(1) == true) ? 1.0 : -1.0;
      d2 = 0;
      d3 = 1;
    } else {
      sFlag3 = (elem->faceFlags[0].test(2) == true) ? 1.0 : -1.0;
      sFlag2 = (elem->faceFlags[0].test(1) == true) ? 1.0 : -1.0;
      d3 = 0;
      d2 = 1;

    }
//    for( Integer i = 2; i <= (order_2) - 2; i++ ) {
//      for( Integer j = 2; j <= (order_3) - 2; j++ ) {
    for( Integer i = 2; i <= (order_2) ; i++ ) {
      for( Integer j = 2; j <= (order_3) ; j++ ) {
        //if( (i + j) > std::max(order_2,order_3) ) { continue; }
        EvalPolynom( val_2, deriv_2, i, lCoeff_[i],
                     sFlag2* point[d2] );
        EvalPolynom( val_3, deriv_3, j, lCoeff_[j],
                     sFlag3*point[d3] );
        deriv[offset][d2] = deriv_2 * sFlag2 * val_3;
        deriv[offset][d3] = val_2 * deriv_3 * sFlag3;
        offset++;
      }
    }
  }

  
  // --- Hex 1st order ---
  FeH1HiHex::FeH1HiHex() {
    feType_ = Elem::ET_HEXA8;
    shape_ = Elem::shapes[feType_];
    actNumFncs_ = 0;
    order_ = 0;
  }
    
  FeH1HiHex::~FeH1HiHex() {
    
  }

  void FeH1HiHex::SetIsoOrder(UInt order) {
    actNumFncs_ = 0;

    // Vertices
    StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
    vertFncs.Resize(shape_.numVertices);
    vertFncs.Init(1); // Vertices have always order 1
    actNumFncs_ += shape_.numVertices * 1;

    // Edges
    StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
    edgeFncs.Resize(shape_.numEdges);
    edgeFncs.Init(order-1);
    actNumFncs_ += shape_.numEdges * (order-1);

    // Faces
    StdVector<UInt>& faceFncs = entityFncs_[FACE];
    faceFncs.Resize(shape_.numFaces);
    faceFncs.Init( (order-1)*(order-1) );
    actNumFncs_ += shape_.numFaces * (order-1)*(order-1) ;
    
    // Inner
    StdVector<UInt>& intFncs = entityFncs_[INTERIOR];
    intFncs.Resize(1);
    intFncs.Init( (order-1)*(order-1)*(order-1) );
    actNumFncs_ += (order-1)*(order-1)*(order-1) ;

  }
  
  void FeH1HiHex::GetShFnc( Vector<Double>& shape, const LocPoint& lp,
                            const Elem * elem, UInt comp ) {
    
  }
  
  void FeH1HiHex::GetDerivShFnc( Matrix<Double> & deriv, const LocPoint& lp,
                                 const Elem * elem,  UInt comp ) {
  }
  
} // namespace CoupledField
