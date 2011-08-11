// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "H1ElemsHi.hh"
#include "Utils/autodiff.hh"
#include "polynomials.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField {


#define USE_EDGES 1
#define USE_FACES 1
#define USE_INNER 1

// declare class specific logging stream
  DECLARE_LOG(feH1Hi)
  DEFINE_LOG(feH1Hi, "feH1Hi");


  // ========================================================================
  //  FeH1Hi
  // ========================================================================
  
  FeH1Hi::FeH1Hi() {
    updateUnknowns_ = true;
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
    if(updateUnknowns_) CalcNumUnknowns();
    numFcns = entityFncs_[entityType]; 
    
  }

  void FeH1Hi::CalcNumUnknowns() {
    
    LOG_DBG(feH1Hi) << "CalcNumUnknowns for element "
        << Elem::feType.ToString(feType_);
    
    actNumFncs_ = 0;
    
    // Vertices
    StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
    vertFncs.Resize(shape_.numVertices);
    vertFncs.Init(1); // Vertices have always order 1
    actNumFncs_ += shape_.numVertices * 1;
    
    // Edges
    StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
    edgeFncs.Resize(shape_.numEdges);
    edgeFncs.Init(0);
    UInt unknowns = 0;
#ifdef USE_EDGES
    for( UInt i = 0; i < shape_.numEdges; ++i ) {
      unknowns = (orderEdge_[i]-1);
      edgeFncs[i] = unknowns;
      LOG_DBG(feH1Hi) <<   "edge " << i+1 << " has " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
#endif
    
    // Faces
    StdVector<UInt>& faceFncs = entityFncs_[FACE];
    faceFncs.Resize(shape_.numFaces);
    faceFncs.Init(0);
#ifdef USE_FACES
    if( shape_.dim > 1 ) {

      for( UInt i = 0; i < shape_.numFaces; ++i ) {
        unknowns = (orderFace_[i][0]-1) * (orderFace_[i][1]-1);
        faceFncs[i] = unknowns;
        LOG_DBG(feH1Hi) << "face " << i+1 << " has " << unknowns << "unknowns";
        actNumFncs_ += unknowns;
      }
    }
#endif

    // Interior
    StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
    innerFncs.Resize(1);
    innerFncs.Init(0);
#ifdef USE_INNER
    if( shape_.dim > 2 ) {
      unknowns = (orderInner_[0]-1) * (orderInner_[1]-1) 
                                      *(orderInner_[2]-1); 
      innerFncs[0] = unknowns;
      LOG_DBG(feH1Hi) << "interior has " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
#endif

    LOG_DBG(feH1Hi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;
    updateUnknowns_ = false;
  }
  
  void FeH1Hi::SetIsoOrder( UInt order ) {
    
    orderEdge_.Resize(shape_.numEdges);
    orderFace_.Resize(shape_.numFaces);
    
    // set order for edges
    orderEdge_.Init(order);
    
    // set order for faces
    array<UInt,2> faceOrder = {order, order};
    orderFace_.Init(faceOrder);

    // set order for inner
    array<UInt, 3> innerOrder = {order, order, order}; 
    orderInner_ = innerOrder;
    
    updateUnknowns_ = true;
  }
  
  void FeH1Hi::GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                    const Elem* ptElem,
                                    EntityType fctEntityType,
                                    UInt entNumber){
    if (updateUnknowns_) CalcNumUnknowns();
   
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

  // =======================
  //  LINE ELEMENT 
  // =======================

  FeH1HiLine::FeH1HiLine() {
    feType_ = Elem::ET_LINE2;
    shape_ = Elem::shapes[feType_];
  }

  FeH1HiLine::~FeH1HiLine() {

  }

  void FeH1HiLine::CalcShFnc( Vector<Double>& shape, 
                              const Vector<Double>& lp,
                             const Elem * elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();
  }

  void FeH1HiLine::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                      const Vector<Double>& lp,
                                      const Elem * elem,  UInt comp ) {
    if (updateUnknowns_) CalcNumUnknowns();
  }
  
  // =======================
  //  QUADRILATERAL ELEMENT 
  // =======================
   
  FeH1HiQuad::FeH1HiQuad() {
    feType_ = Elem::ET_QUAD4;
    shape_ = Elem::shapes[feType_];
  }
    
  FeH1HiQuad::~FeH1HiQuad() {
    
  }
  
  void FeH1HiQuad::CalcShFnc( Vector<Double>& shape, 
                              const Vector<Double>& lp,
                             const Elem * elem, UInt comp ) {
    if (updateUnknowns_) CalcNumUnknowns();
    _CalcShFnc(lp[0], lp[1], elem, shape);
  }
  
  void FeH1HiQuad::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                      const Vector<Double>& lp,
                                      const Elem * elem,  UInt comp )  {
    if (updateUnknowns_) CalcNumUnknowns();
    AutoDiff<Double,2> x(lp[0],0), y(lp[1],1);
    StdVector<AutoDiff<Double,2> > dShape;
    _CalcShFnc(x,y,elem,dShape);
    UInt size = dShape.GetSize();
    deriv.Resize(size, 2);
    for( UInt i = 0; i < size; ++i ) {
      for(UInt j = 0; j < 2; ++j ) {
        deriv[i][j] = dShape[i].DVal(j);
      }
    }
  }

  template<typename T_SCAL, typename T_VEC>
  void FeH1HiQuad::_CalcShFnc( const T_SCAL x, const T_SCAL y, 
                               const Elem * elem,
                               T_VEC& ret ) {
    
    T_SCAL lambda[4] = {0.25*(1.0-x)*(1.0-y),
                        0.25*(1.0+x)*(1.0-y),
                        0.25*(1.0+x)*(1.0+y),
                        0.25*(1.0-x)*(1.0+y)};
    
    T_SCAL sigma[4]  = {0.5*((1.0-x)+(1.0-y)),
                        0.5*((1.0+x)+(1.0-y)),
                        0.5*((1.0+x)+(1.0+y)),
                        0.5*((1.0-x)+(1.0+y))};
    UInt pos = 0;
    ret.Resize(actNumFncs_);

    // 1) Vertex shape functions
    for( UInt i = 0; i < 4; ++i ) {
      ret[i] = lambda[i];
    }

    pos = 4;
#ifdef USE_EDGES 
    // 2) Edge shape functions
    for( UInt i = 0; i < 4; ++ i ) {
      Double fac = elem->edges[i] < 0 ? -1.0 : 1.0;
      UInt order = orderEdge_[i];
      
      // if order of edge is below two, we leave
      if( order == 1 ) continue;
      UInt index1 = shape_.edgeVertices[i][0]-1;
      UInt index2 = shape_.edgeVertices[i][1]-1;

      // xi: parameterization of edge [-1;+1]
      // eta: parameterization of extension into element [-1;+1]
      T_SCAL xi  =  sigma[index2] -  sigma[index1]; 
      T_SCAL eta = lambda[index1] + lambda[index2];

      T_VEC vals;
      IntLegendrePoly<T_SCAL,T_VEC>( vals, order, fac*xi );

      for( UInt i = 0; i < order-1; ++i ) {
        ret[pos++] = eta * vals[i];
      } 
    }    
#endif
    
#ifdef USE_FACES
    // 3) Inner shape functions
    UInt order1 = orderFace_[0][0];
    UInt order2 = orderFace_[0][1];
    if (order1 >= 2 && order2 >= 2 ) {

      T_SCAL xi = x;
      T_SCAL eta = y;
      // test, if xi and eta get reversed
      if( elem->faceFlags[0][2] == false) {
        std::swap(order1, order2);
        std::swap(xi, eta);
      }
      xi = elem->faceFlags[0][0] ? xi : -xi;
      eta = elem->faceFlags[0][1] ? eta : -eta;
      
      T_VEC xiVals, etaVals;
      IntLegendrePoly<T_SCAL,T_VEC>( xiVals, order1, xi );
      IntLegendrePoly<T_SCAL,T_VEC>( etaVals, order2, eta );
      for( UInt k = 0; k < order1-1; ++k)
        for( UInt j = 0; j < order2-1; ++j)
          ret[pos++] = etaVals[k] * xiVals[j];
    }
#endif
  }
    
  
//  void FeH1HiQuad::CalcShFnc( Vector<Double>& shape, const LocPoint& lp,
//                             const Elem * elem, UInt comp ) {
//    
//    // create nodal shape functions
//    shape.Resize( actNumFncs_ );
//    const Vector<Double> & point = lp.coord;
//    // --------------------
//    //  a) nodal functions
//    // --------------------
//    shape[0] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 - point[1] ); 
//    shape[1] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 - point[1] );
//    shape[2] = 0.25 * ( 1.0 + point[0] ) * ( 1.0 + point[1] );
//    shape[3] = 0.25 * ( 1.0 - point[0] ) * ( 1.0 + point[1] );
//    
//    UInt offset = 4;
//    
//    // ToDo: The following section is just temporary, as we want to replace
//    // it by the recursive calculation of shape functions taken from ngsolve
//    // --------------------
//    //  b) edge functions
//    // --------------------
//    // Obtain order of element
//    Integer order;
//    Double val, factor, deriv;
//
//#define QUAD_EDGE_FCN(edgeNum,  sign_1, dir_1, dir_2 )                  \
//    factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
//    for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {        \
//      EvalPolynom( val, deriv, iDof, lCoeff_[iDof],                   \
//                   factor*point[dir_2] );                          \
//                   shape[offset] = 0.5 * ( 1 sign_1 point[dir_1] ) * val;       \
//    }
//
//    //  EDGE #1
//    order = entityOrder_[EDGE][0];
//    QUAD_EDGE_FCN( 1, -, 1, 0 );
//
//    //  EDGE #2
//    order = entityOrder_[EDGE][1];
//    QUAD_EDGE_FCN( 2, +, 0, 1 );
//
//    //  EDGE #3
//    order = entityOrder_[EDGE][2];
//    QUAD_EDGE_FCN( 3, +, 1, 0 );
//
//    //  EDGE #4
//    order = entityOrder_[EDGE][3];
//    QUAD_EDGE_FCN( 4, -, 0, 1 );
//    
////    // --------------------
////    //  c) face functions
////    // --------------------
////    Double val_2, deriv_2, val_3, deriv_3;
////    Double sFlag2, sFlag3;
////    Integer order_2, order_3;
////    UInt d2, d3;
////
////    order_2 = entityOrder_[EDGE][0];
////    order_3 = entityOrder_[EDGE][1];
////    if( elem->faceFlags[0].test(0) == true) {
////      sFlag2 = (elem->faceFlags[0].
////          test(2) == true) ? 1.0 : -1.0;
////      sFlag3 = (elem->faceFlags[0].
////          test(1) == true) ? 1.0 : -1.0;
////      d2 = 0;
////      d3 = 1;
////    } else {
////      sFlag3 = (elem->faceFlags[0].test(2) == true) ? 1.0 : -1.0;
////      sFlag2 = (elem->faceFlags[0].test(1) == true) ? 1.0 : -1.0;
////      d3 = 0;
////      d2 = 1;
////    }
////    for( Integer i = 2; i <= (order_2)-2; i++ ) {
////      for( Integer j = 2; j <= (order_3)-2; j++ ) {
////        if( (i + j) > std::max(order_2, order_3) ) {
////          continue;
////        }
////        EvalPolynom( val_2, deriv_2, i, lCoeff_[i],
////                     sFlag2* point[d2] );
////        EvalPolynom( val_3, deriv_3, j, lCoeff_[j],
////                     sFlag3*point[d3] );
////        shape[offset] = val_2  * val_3;
////        offset++;
////      }
////    }
//
//  }
//
//  void FeH1HiQuad::GetDerivShFnc( Matrix<Double> & deriv, const LocPoint& lp,
//                                  const Elem * elem,  UInt comp )  {
//    // nodal shape functions
//    
//   const Vector<Double> & point = lp.coord;
//    StdVector<StdVector<Double> >& coords = shape_.nodeCoords;
//    deriv.Resize( actNumFncs_, 2 );
//    for( UInt i = 0; i < 4; i++ ) {
//      deriv[i][0] = 0.25 * coords[i][0] * (1 + coords[i][1] * point[1] );
//      deriv[i][1] = 0.25 * (1 + coords[i][0] * point[0] ) * coords[i][1];
//    }
//    
//    // -------------------
//    //  b) edge functions
//    // -------------------
//    UInt offset = 4;
//
//    // Obtain order of element
//    Integer order_2, order_3;
//    Double val, deriv1D, factor;
//#define QUAD_EDGE_DERIV( edgeNum, sign_1, dir_1, dir_2 )                \
//    factor = elem->edges[edgeNum-1] < 0 ? -1.0 : 1.0;                 \
//    for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {        \
//      EvalPolynom( val, deriv1D, iDof, lCoeff_[iDof],                   \
//                   factor*point[dir_2] );                          \
//                   deriv[offset][dir_1] = sign_1(0.5 * val);                      \
//                   deriv[offset][dir_2] =  0.5 * ( 1.0 sign_1 point[dir_1])    \
//                   * deriv1D * factor;                  \
//    }
//
//    // EDGE #1
//    UInt order = entityOrder_[EDGE][0];
//    QUAD_EDGE_DERIV( 1, -, 1, 0 );
//
//    // EDGE #2
//    order = entityOrder_[EDGE][1];
//    QUAD_EDGE_DERIV( 2, +, 0, 1 );
//
//    // EDGE #3
//    order = entityOrder_[EDGE][2];
//    QUAD_EDGE_DERIV( 3, +, 1, 0 );
//
//    // EDGE #4
//    order = entityOrder_[EDGE][3];
//    QUAD_EDGE_DERIV( 4, -, 0, 1 );
//
////    // -------------------
////    //  c) face functions
////    // -------------------
////
////    Double val_2, deriv_2, val_3, deriv_3;
////    Double sFlag2, sFlag3;
////    UInt d2, d3;
////
////    order_2 = order = entityOrder_[EDGE][1];
////    order_3 = order = entityOrder_[EDGE][2];
////    if( elem->faceFlags[0].test(0) == true) {
////      sFlag2 = (elem->faceFlags[0].
////          test(2) == true) ? 1.0 : -1.0;
////      sFlag3 = (elem->faceFlags[0].
////          test(1) == true) ? 1.0 : -1.0;
////      d2 = 0;
////      d3 = 1;
////    } else {
////      sFlag3 = (elem->faceFlags[0].test(2) == true) ? 1.0 : -1.0;
////      sFlag2 = (elem->faceFlags[0].test(1) == true) ? 1.0 : -1.0;
////      d3 = 0;
////      d2 = 1;
////
////    }
//////    for( Integer i = 2; i <= (order_2) - 2; i++ ) {
//////      for( Integer j = 2; j <= (order_3) - 2; j++ ) {
////    for( Integer i = 2; i <= (order_2) ; i++ ) {
////      for( Integer j = 2; j <= (order_3) ; j++ ) {
////        //if( (i + j) > std::max(order_2,order_3) ) { continue; }
////        EvalPolynom( val_2, deriv_2, i, lCoeff_[i],
////                     sFlag2* point[d2] );
////        EvalPolynom( val_3, deriv_3, j, lCoeff_[j],
////                     sFlag3*point[d3] );
////        deriv[offset][d2] = deriv_2 * sFlag2 * val_3;
////        deriv[offset][d3] = val_2 * deriv_3 * sFlag3;
////        offset++;
////      }
////    }
//  }

  
  // ============================
  //  HEXAHEDRAL ELEMENT ELEMENT 
  // ============================
  
  FeH1HiHex::FeH1HiHex() {
    feType_ = Elem::ET_HEXA8;
    shape_ = Elem::shapes[feType_];
  }
    
  FeH1HiHex::~FeH1HiHex() {
    
  }

  
  void FeH1HiHex::CalcShFnc( Vector<Double>& shape, 
                             const Vector<Double>& lp,
                            const Elem * elem, UInt comp ) {
    if (updateUnknowns_) CalcNumUnknowns();
    _CalcShFnc(lp[0], lp[1], lp[2],  elem, shape);
    
    
  }
  
  void FeH1HiHex::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                     const Vector<Double>& lp,
                                     const Elem * elem,  UInt comp ) {
    if (updateUnknowns_) CalcNumUnknowns();
    AutoDiff<Double,3> x(lp[0],0), y(lp[1],1),z(lp[2],2);
    StdVector<AutoDiff<Double,3> > dShape;
    _CalcShFnc(x,y,z,elem,dShape);
    UInt size = dShape.GetSize();
    deriv.Resize(size, 3);
    for( UInt i = 0; i < size; ++i ) {
      for(UInt j = 0; j < 3; ++j ) {
        deriv[i][j] = dShape[i].DVal(j);
      }
    }
    
    // Debugging stuff
//    StdVector<Double> shape(size);
//    for( UInt i = 0; i < size; ++i ) {
//      shape[i] = dShape[i].Val();
//    }
//    std::cerr << "shape functions are " << shape << std::endl;
  }
  
  template<typename T_SCAL, typename T_VEC>
  void FeH1HiHex::_CalcShFnc( const T_SCAL x,  const T_SCAL y, const T_SCAL z, 
                               const Elem * elem,
                               T_VEC& ret ) {

    T_SCAL lambda[8] = {0.125*(1.0-x)*(1.0-y)*(1.0-z),
                        0.125*(1.0+x)*(1.0-y)*(1.0-z),
                        0.125*(1.0+x)*(1.0+y)*(1.0-z),
                        0.125*(1.0-x)*(1.0+y)*(1.0-z),
                        0.125*(1.0-x)*(1.0-y)*(1.0+z),
                        0.125*(1.0+x)*(1.0-y)*(1.0+z),
                        0.125*(1.0+x)*(1.0+y)*(1.0+z),
                        0.125*(1.0-x)*(1.0+y)*(1.0+z)};

    T_SCAL sigma[8]  = {0.5*((1.0-x)+(1.0-y)+(1.0-z)),
                        0.5*((1.0+x)+(1.0-y)+(1.0-z)),
                        0.5*((1.0+x)+(1.0+y)+(1.0-z)),
                        0.5*((1.0-x)+(1.0+y)+(1.0-z)),
                        0.5*((1.0-x)+(1.0-y)+(1.0+z)),
                        0.5*((1.0+x)+(1.0-y)+(1.0+z)),
                        0.5*((1.0+x)+(1.0+y)+(1.0+z)),
                        0.5*((1.0-x)+(1.0+y)+(1.0+z))};

    UInt pos = 0;
    ret.Resize(actNumFncs_);

    // 1) Vertex shape functions
    for( UInt i = 0; i < 8; ++i ) {
      ret[i] = lambda[i];
    }

    pos = 8;

    // 2) Edge shape functions
    for( UInt i = 0; i < 12; ++ i) {
      
      UInt order = orderEdge_[i];
      // if order of edge is below two, we leave
      if( order == 1 ) continue;

      Double fac = elem->edges[i] < 0 ? -1.0 : 1.0;
      UInt index1 = shape_.edgeVertices[i][0]-1;
      UInt index2 = shape_.edgeVertices[i][1]-1;

      // xi: parameterization of edge [-1;+1]
      // eta: parameterization of extension into element [-1;+1]
      T_SCAL xi  =  sigma[index2] -  sigma[index1]; 
      T_SCAL eta = lambda[index1] + lambda[index2];

      T_VEC vals;
      IntLegendrePoly<T_SCAL,T_VEC>( vals, order, fac*xi );

      for( UInt i = 0; i < order-1; ++i ) {
        ret[pos++] = eta * vals[i];
      } 
    }    
    
    T_VEC xiVals, etaVals;
    // 3) Face shape functions
    for( UInt i = 0; i < 6; ++i ) {
      UInt order1 = orderFace_[i][0];
      UInt order2 = orderFace_[i][1];
      if (order1 >= 2 && order2 >= 2 ) {
        // calculate face extension parameter which is the sum
        // of all lambdas of one face
        T_SCAL sum_lambda = 0.0;
        for( UInt j = 0; j < 4; ++j)
          sum_lambda += lambda[shape_.faceVertices[i][j]-1];
        
        // Parameterization of first edge, connecting the
        // local nodes of the face 1->2
        T_SCAL xi  =  sigma[shape_.faceVertices[i][1]-1]
                    - sigma[shape_.faceVertices[i][0]-1];
        // Parameterization of second edge, connecting the
        // local nodes of the face 1->4
        T_SCAL eta  =  sigma[shape_.faceVertices[i][3]-1]
                     - sigma[shape_.faceVertices[i][0]-1];
        
        // test, if xi and eta get reversed
        if( elem->faceFlags[i][2]) {
          std::swap(order1, order2);
          std::swap(xi, eta);
        }
        xi =  elem->faceFlags[i][0] ? xi : -xi;
        eta = elem->faceFlags[i][2] ? eta : -eta;

        IntLegendrePoly<T_SCAL,T_VEC>( xiVals, order1, xi );
        IntLegendrePoly<T_SCAL,T_VEC>( etaVals, order2, eta );
        for( UInt k = 0; k < order1-1; ++k)
          for( UInt j = 0; j <  order2-1; ++j)
            ret[pos++] = etaVals[k] * xiVals[j] * sum_lambda;
      } // if 
    }
    
    // 4) Inner shape functions
    if( orderInner_[0] >= 2 && 
        orderInner_[1] >= 2 &&
        orderInner_[2] >= 2 ) {
      T_VEC zetaVals;
      IntLegendrePoly<T_SCAL,T_VEC>( xiVals, orderInner_[0], x );
      IntLegendrePoly<T_SCAL,T_VEC>( etaVals, orderInner_[1], y );
      IntLegendrePoly<T_SCAL,T_VEC>( zetaVals, orderInner_[2], z );
      for( UInt i = 0; i < orderInner_[0]-1; ++i ) {
        for( UInt j = 0; j < orderInner_[1]-1; ++j ) {
          T_SCAL temp = xiVals[i] * etaVals[j];
          for( UInt k = 0; k < orderInner_[2]-1; ++k ) {
            ret[pos++] = zetaVals[k] * temp;
          } // k
        } //j
      }//i
    } // if
  }
  
} // namespace CoupledField
