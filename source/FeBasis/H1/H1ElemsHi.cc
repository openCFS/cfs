#include "H1ElemsHi.hh"

#include <algorithm>

#include "Utils/AutoDiff.hh"
#include "Utils/Polynomial.hh"
#include "FeBasis/Polynomials.hh"
#include "Domain/ElemMapping/EdgeFace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {


// the following macros are just used for debugging purposes to
// switch off / on certain degrees of freedoms
#define USE_EDGES 1
#define USE_FACES 1
#define USE_INNER 1

// declare class specific logging stream
DECLARE_LOG(feH1Hi)
DEFINE_LOG(feH1Hi, "feH1Hi")


  // ========================================================================
  //  FeH1Hi
  // ========================================================================
  
  FeH1Hi::FeH1Hi(Elem::FEType feType ) 
  : FeH1(), FeHi(feType )
{
    feType_ = feType;
    shape_ = Elem::shapes[feType];
    updateUnknowns_ = true;
    isIsotropic_ = false;
    isoOrder_ = 0;
    maxOrder_ = 0;
    
    // important: all higher order functions can not
    // pre-compute the shape functions, as the functions depend on the global
    // orientation, i.e. numbering of the nodal connectivity
    preComputShFnc_ = false;
  }
  
  FeH1Hi::~FeH1Hi() {
  }
  
  UInt FeH1Hi::GetNumFncs( ) {
    if(updateUnknowns_) this->CalcNumUnknowns();
    return actNumFncs_;
  }
  
  void FeH1Hi::GetNumFncs( StdVector<UInt>& numFcns,
                           EntityType entityType,
                           UInt dof ) {
    
    // this can be implemented in a general way, as the information about
    // the number of function is already stored in entityFncs_
    
    // Note: at the moment we ignore the dof parameter
    
    // numFcns has the length Vertices / Edges / Faces / Inner degrees
    if(updateUnknowns_) this->CalcNumUnknowns();
    numFcns = entityFncs_[entityType]; 
    
  }

  void FeH1Hi::CalcNumUnknowns() {
    
    // This method calculates the number of unknowns for all
    // elements based on the tensor product, i.e. the 
    // line, quad and hex element.
    
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
      LOG_DBG(feH1Hi) <<   "edge " << i+1 << " has order " <<  orderEdge_[i]-1
            << " and " << unknowns << " unknowns";
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
        LOG_DBG(feH1Hi) << "face " << i+1 << " has " << unknowns << " unknowns";
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
      LOG_DBG(feH1Hi) << "interior has " << unknowns << " unknowns";
      actNumFncs_ += unknowns;
    }
#endif

    LOG_DBG(feH1Hi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;
    updateUnknowns_ = false;
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

  UInt FeH1Hi::GetIsoOrder() const {
    if( isIsotropic_) {
      return isoOrder_;
    } else {
      return 0;  
    }
  }

  UInt FeH1Hi::GetMaxOrder() const {
    return maxOrder_;
  }

  void FeH1Hi::GetAnisoOrder(StdVector<UInt>& order ) const {

    if( isIsotropic_) {
      order.Resize(Elem::shapes[feType_].dim);
      order.Init(isoOrder_);
    } else {
      order = anisoOrder_;  
    }
  }
  
  bool FeH1Hi::operator==( const FeH1Hi& comp) const {
    bool ret = true;
    ret &= this->feType_ == comp.feType_;
    ret &= this->isIsotropic_ == comp.isIsotropic_;
    if( isIsotropic_ ) {
      ret &= this->isoOrder_ == comp.isoOrder_;
    } else {
      ret &= this->anisoOrder_ == comp.anisoOrder_;
      ret &= this->orderEdge_ == comp.orderEdge_;
      ret &= this->orderFace_ == comp.orderFace_;
      ret &= this->orderInner_ == comp.orderInner_;
    }
    return ret;
    
  }
  
  
  
//  void FeH1Hi::EvalPolynom( Double& value, Double& deriv,
//                            const UInt order, const Double* coeff,
//                            const Double xVal ) {
//
//    // Consider the following expression
//    // f(xVal) = a0 * (a1*x^order + a2*x^(order-1) + .. + a(order+1))
//    // The coefficients a0..a(order+1) are stored in the coeff-array
//
//    value = coeff[1];
//    deriv = 0.0;
//    for( UInt i = 2; i < order+2; i++ ) {
//      deriv = deriv * xVal + value;
//      value = value * xVal + coeff[i];
//    }
//    // Multiply by pre-factor
//    deriv *= coeff[0];
//    value *= coeff[0];
//  }
//
//
//  // Define coefficients for legendre ansatz functions up to order 8
//  Double  FeH1Hi::lCoeff_[9][10] = {
//                                    {0.5                  ,   -1, 1,    0, 0,   0, 0,    0, 0, 0 },
//                                    {0.5                  ,    1, 1,    0, 0,   0, 0,    0, 0, 0 },
//                                    {0.25*sqrt(6.0)       ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
//                                    {0.25*sqrt(10.0)      ,    1, 0,   -1, 0,   0, 0,    0, 0, 0 },
//                                    {1.0/16.0*sqrt(14.0)  ,    5, 0,   -6, 0,   1, 0,    0, 0, 0 },
//                                    {3.0/16.0*sqrt(2.0)   ,    7, 0,  -10, 0,   3, 0,    0, 0, 0 },
//                                    {1.0/32.0*sqrt(22.0)  ,   21, 0,  -35, 0,  15, 0,   -1, 0, 0 },
//                                    {1.0/32.0*sqrt(26.0)  ,   33, 0,  -63, 0,  35, 0,   -5, 0, 0 },
//                                    {1.0/256.0*sqrt(30.0) ,  429, 0, -924, 0, 630, 0, -140, 0, 5 }
//  };

  // ========================================================================
  //  FeH1Hi explicit element definition 
  // ========================================================================

  // =======================
  //  LINE ELEMENT 
  // =======================

  FeH1HiLine::FeH1HiLine() 
  : FeH1Hi(Elem::ET_LINE2) 
  {
  }

  FeH1HiLine::~FeH1HiLine() {

  }

  void FeH1HiLine::CalcShFnc( Vector<Double>& shape, 
                              const Vector<Double>& lp,
                             const Elem * elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();
  _CalcShFnc(lp[0], elem, shape);
  }

  void FeH1HiLine::CalcLocDerivShFnc( Matrix<Double> & deriv, 
                                      const Vector<Double>& lp,
                                      const Elem * elem,  UInt comp ) {
    if (updateUnknowns_) CalcNumUnknowns();
    
    AutoDiff<Double,1> x(lp[0],0);
    StdVector<AutoDiff<Double,1> > dShape;
    _CalcShFnc(x,elem,dShape);
    UInt size = dShape.GetSize();
    deriv.Resize(size, 1);
    for( UInt i = 0; i < size; ++i ) {
      deriv[i][0] = dShape[i].DVal(0);
    }
  }
  void FeH1HiLine::GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                        const Elem* ptElem ) {
    if (updateUnknowns_) CalcNumUnknowns();
    Polynomial<Double,1> x(1, 0);
    StdVector<Polynomial<Double,1> > order;
    _CalcShFnc(x,ptElem,order);
    UInt numFncs = order.GetSize();
    polyOrder.Resize( numFncs, 1 );
    for( UInt i = 0; i < numFncs; ++i ) {
      polyOrder[i][0] = order[i].GetMaxOrder(0);
    }
  }

  template<typename T_SCAL, typename T_VEC>
  void FeH1HiLine::_CalcShFnc( const T_SCAL x,
                               const Elem * elem,
                               T_VEC& ret ) {

    Double fac = elem->edges[0] < 0 ? -1.0 : 1.0;
    UInt order = orderEdge_[0];
    ret.Resize(actNumFncs_);
    ret[0] = 0.5 * ( 1.0 - x );
    ret[1] = 0.5 * ( 1.0 + x );
    
    T_VEC vals;
    vals.Init();
    IntLegendreP2<T_SCAL,T_VEC>( vals, order, fac*x );
    UInt pos = 2;
    for( UInt j = 0; j < order-1; ++j ) {
      ret[pos++] = vals[j];
    } 
  }
  
  
  // ====================
  //  TRIANGULAR ELEMENT 
  // ====================

  FeH1HiTria::FeH1HiTria() : FeH1Hi(Elem::ET_TRIA3) {
  }

  FeH1HiTria::~FeH1HiTria() {

  }

  void FeH1HiTria::CalcShFnc( Vector<Double>& shape, 
                              const Vector<Double>& lp,
                              const Elem * elem, UInt comp ) {
    if (updateUnknowns_) CalcNumUnknowns();
    _CalcShFnc(lp[0], lp[1], elem, shape);
  }

  void FeH1HiTria::CalcLocDerivShFnc( Matrix<Double> & deriv, 
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
  
  void FeH1HiTria::GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                          const Elem* ptElem ) {

      
      if (updateUnknowns_) CalcNumUnknowns();
      Polynomial<Double,2> x(1., 0), y(1., 1);
      StdVector<Polynomial<Double,2> > order;
      _CalcShFnc(x,y,ptElem,order);
      UInt numFncs = order.GetSize();
      polyOrder.Resize( numFncs, 2 );
      for( UInt i = 0; i < numFncs; ++i ) {
        for(UInt j = 0; j < 2; ++j ) {
          polyOrder[i][j] = order[i].GetMaxOrder(j);
        }
      }
    }

  template<typename T_SCAL, typename T_VEC>
  void FeH1HiTria::_CalcShFnc( const T_SCAL x, const T_SCAL y, 
                               const Elem * elem,
                               T_VEC& ret ) {

    T_SCAL lambda[3] = { 1.0 - x - y,
                         x,
                         y };

    UInt pos = 0;
    ret.Resize(actNumFncs_);

    // 1) Vertex shape functions
    for( UInt i = 0; i < 3; ++i ) {
      ret[i] = lambda[i];
    }

    pos = 3;
#ifdef USE_EDGES 
    // 2) Edge shape functions
    for( UInt i = 0; i < 3; ++ i ) {
      Double fac = elem->edges[i] < 0 ? -1.0 : 1.0;
      UInt order = orderEdge_[i];

      // if order of edge is below two, we leave
      if( order == 1 ) continue;
      UInt index1 = shape_.edgeVertices[i][0]-1;
      UInt index2 = shape_.edgeVertices[i][1]-1;

      // edge: parameterization of edge [-1;+1]
      // edgeNormal: parameterization of extension into element [-1;+1]
      T_SCAL edge  =  lambda[index2] -  lambda[index1]; 
      T_SCAL edgeNormal = lambda[index2] + lambda[index1];

      T_VEC vals;
      ScaledIntLegendreP2<T_SCAL,T_VEC>( vals, order, edgeNormal, fac*edge );
      for( UInt j = 0; j < order-1; ++j ) {
        ret[pos++] = vals[j];
      } 
    }   
#endif
    //      
#ifdef USE_FACES
    // 3) Inner shape functions
    UInt order = orderFace_[0][0];
    if (order >= 3 ) {
      
      // obtain nodes of face
      const StdVector<UInt>& unsorted = shape_.faceNodes[0];
      StdVector<UInt> ind;
      Face::GetSortedIndices( ind, unsorted, 3, elem->faceFlags[0]);
      // calculate inner shape functions
      UInt nFct =  TriaInnerLegendre( ret, pos, order, 
                                      lambda[ind[0]], lambda[ind[1]],
                                      lambda[ind[2]]);
      pos += nFct;
//        // OLD VERSION WITH EXPLITICLY CODED CALCULATION
//      T_VEC f1, f2;
//      ScaledIntLegendreP2(f1, order-1, lambda[1]+lambda[0],lambda[1]-lambda[0]);
//      Legendre(f2, order-3, lambda[2]*2.0 - T_SCAL(1));
//      for( UInt i = 0; i <= order - 3; ++i ) {
//        for( UInt j = 0; j <= order - 3 - i; ++j ) {
//          ret[pos++] = f1[i] * f2[j] * lambda[2];
//        }
//      }
    }
#endif
  }
  void FeH1HiTria::CalcNumUnknowns() {

    LOG_DBG(feH1Hi) << "CalcNumUnknowns for element "
        << Elem::feType.ToString(feType_);

    actNumFncs_ = 0;

    // Vertices
    StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
    vertFncs.Resize(3);
    vertFncs.Init(1); // Vertices have always order 1
    actNumFncs_ += 3;

    // Edges
    StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
    edgeFncs.Resize(3);
    edgeFncs.Init(0);
    UInt unknowns = 0;
#ifdef USE_EDGES
    for( UInt i = 0; i < 3; ++i ) {
      unknowns = (orderEdge_[i]-1);
      edgeFncs[i] = unknowns;
      LOG_DBG(feH1Hi) <<   "edge " << i+1 << " has order" <<  orderEdge_[i]-1
          << " and " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
#endif

    // Faces
    StdVector<UInt>& faceFncs = entityFncs_[FACE];
    faceFncs.Resize(1);
    faceFncs.Init(0);
#ifdef USE_FACES
    if( orderFace_[0][0] > 2 ) {
      unknowns = (orderFace_[0][0]-2) * (orderFace_[0][0]-1) / 2;
      faceFncs[0] = unknowns;
      LOG_DBG(feH1Hi) << "face 0 has " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
#endif
    // Interior
    StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
    innerFncs.Resize(1);
    innerFncs.Init(0);

    LOG_DBG(feH1Hi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;
    updateUnknowns_ = false;
  }

  // =======================
  //  QUADRILATERAL ELEMENT 
  // =======================
   
  FeH1HiQuad::FeH1HiQuad() : FeH1Hi(Elem::ET_QUAD4) {
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
#ifdef USE_EDGES
    T_SCAL sigma[4]  = {0.5*((1.0-x)+(1.0-y)),
                        0.5*((1.0+x)+(1.0-y)),
                        0.5*((1.0+x)+(1.0+y)),
                        0.5*((1.0-x)+(1.0+y))};
#endif
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
      IntLegendreP2<T_SCAL,T_VEC>( vals, order, fac*xi );

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

      const StdVector<UInt>& unsorted = shape_.faceNodes[0];
      StdVector<UInt> ind;
      Face::GetSortedIndices( ind, unsorted, 4, elem->faceFlags[0]);

      T_SCAL xi  = sigma[ind[1]] - sigma[ind[0]];
      T_SCAL eta = sigma[ind[3]] - sigma[ind[0]];
      // === OLD IMPLEMENTATION ===
//      // test, if xi and eta get reversed
//      if( elem->faceFlags[0][2] == false) {
//        std::swap(order1, order2);
//        std::swap(xi, eta);
//      }
//      xi = elem->faceFlags[0][0] ? xi : -xi;
//      eta = elem->faceFlags[0][1] ? eta : -eta;
      
      T_VEC xiVals, etaVals;
      IntLegendreP2<T_SCAL,T_VEC>( xiVals, order1, xi );
      IntLegendreP2<T_SCAL,T_VEC>( etaVals, order2, eta );
      for( UInt k = 0; k < order1-1; ++k)
        for( UInt j = 0; j < order2-1; ++j)
          ret[pos++] = xiVals[k] * etaVals[j];
    }
#endif
  }
  
  void FeH1HiQuad::GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                        const Elem* ptElem ) {

    
    if (updateUnknowns_) CalcNumUnknowns();
    Polynomial<Double,2> x(1., 0), y(1., 1);
    StdVector<Polynomial<Double,2> > order;
    _CalcShFnc(x,y,ptElem,order);
    UInt numFncs = order.GetSize();
    polyOrder.Resize( numFncs, 2 );
    for( UInt i = 0; i < numFncs; ++i ) {
      for(UInt j = 0; j < 2; ++j ) {
        polyOrder[i][j] = order[i].GetMaxOrder(j);
      }
    }
  }

  
  // ====================
  //  HEXAHEDRAL ELEMENT  
  // ====================
  
  FeH1HiHex::FeH1HiHex() : FeH1Hi(Elem::ET_HEXA8) {
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
  }
  
  void FeH1HiHex::GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                       const Elem* ptElem ) {
    if (updateUnknowns_) CalcNumUnknowns();
    Polynomial<Double,3> x(1., 0), y(1., 1), z(1., 2);
    StdVector<Polynomial<Double,3> > order;
    _CalcShFnc(x,y,z,ptElem,order);
    UInt numFncs = order.GetSize();
    polyOrder.Resize( numFncs, 3 );
    for( UInt i = 0; i < numFncs; ++i ) {
      for(UInt j = 0; j < 3; ++j ) {
        polyOrder[i][j] = order[i].GetMaxOrder(j);
      }
    }
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
      IntLegendreP2<T_SCAL,T_VEC>( vals, order, fac*xi );

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
        
        // get unique sorting of the face
        const StdVector<UInt>& unsorted = shape_.faceNodes[i];
        StdVector<UInt> ind;
        Face::GetSortedIndices( ind, unsorted, 4, elem->faceFlags[i]);
        
//        const StdVector<UInt> & con = elem->connect;

//        if( ! (con[ind[0]] < con[ind[1]] &&
//               con[ind[0]] < con[ind[3]] &&
//               con[ind[1]] < con[ind[3]] ) ) {
//          std::cerr << "\n\n";
//                 std::cerr << "element: " << elem->elemNum << std::endl;
//                 std::cerr << "connectivity: " << elem->connect.Serialize() << std::endl;
//                 std::cerr << "local face# " << i << "\n";
//                 std::cerr << "FACE #"  << elem->faces[i] << std::endl;
//
//                 std::cerr << "local facNodes = " << unsorted.Serialize() << std::endl;
//                 std::cerr << "faceflag is " <<  elem->faceFlags[i] << std::endl;
//                 std::cerr << "Unsorted is:"
//                     << con[unsorted[0]-1] << ", "
//                     << con[unsorted[1]-1] << ", "
//                     << con[unsorted[2]-1] << ", "
//                     << con[unsorted[3]-1] << std::endl;
//                 std::cerr << "ind array is " << ind.Serialize() << std::endl;
//                 std::cerr << "Sorted is :"
//                     << con[ind[0]] << ", "
//                     << con[ind[1]] << ", "
//                     << con[ind[2]] << ", "
//                     << con[ind[3]]  ;
//
//
//            EXCEPTION("Ordering violated in elem #" << elem->elemNum );
//        }

        
        // calculate face extension parameter which is the sum
        // of all lambdas of one face
        T_SCAL sum_lambda = 0.0;
        for( UInt j = 0; j < 4; ++j)
          sum_lambda += lambda[ind[j]];
        
        // Parameterization of first edge, connecting the
        // local nodes of the face 1->2
        T_SCAL xi  =  sigma[ind[1]] - sigma[ind[0]];
        // Parameterization of second edge, connecting the
        // local nodes of the face 1->4
        T_SCAL eta  =  sigma[ind[3]]- sigma[ind[0]];

        IntLegendreP2<T_SCAL,T_VEC>( xiVals, order1, xi );
        IntLegendreP2<T_SCAL,T_VEC>( etaVals, order2, eta );
        for( UInt j = 0; j < order1-1; ++j) {
          for( UInt k = 0; k <  order2-1; ++k)
            ret[pos++] = xiVals[j] * etaVals[k] * sum_lambda;
          }
      } // if 
    }
    
    // 4) Inner shape functions
    if( orderInner_[0] >= 2 && 
        orderInner_[1] >= 2 &&
        orderInner_[2] >= 2 ) {
      T_VEC zetaVals;
      IntLegendreP2<T_SCAL,T_VEC>( xiVals, orderInner_[0], x );
      IntLegendreP2<T_SCAL,T_VEC>( etaVals, orderInner_[1], y );
      IntLegendreP2<T_SCAL,T_VEC>( zetaVals, orderInner_[2], z );
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

  // ===============
  //  WEDGE ELEMENT  
  // ===============
  
  FeH1HiWedge::FeH1HiWedge() : FeH1Hi(Elem::ET_WEDGE6) {
  }
    
  FeH1HiWedge::~FeH1HiWedge() {
    
  }


  void FeH1HiWedge::CalcShFnc( Vector<Double>& shape, 
                               const Vector<Double>& lp,
                               const Elem * elem, UInt comp ) {
    if (updateUnknowns_) CalcNumUnknowns();
    _CalcShFnc(lp[0], lp[1], lp[2],  elem, shape);


  }

  void FeH1HiWedge::CalcLocDerivShFnc( Matrix<Double> & deriv, 
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

  }
  
  void FeH1HiWedge::GetPolyOrderOfNodes( Matrix<UInt>& polyOrder,
                                       const Elem* ptElem ) {
    if (updateUnknowns_) CalcNumUnknowns();
    Polynomial<Double,3> x(1., 0), y(1., 1), z(1., 2);
    StdVector<Polynomial<Double,3> > order;
    _CalcShFnc(x,y,z,ptElem,order);
    UInt numFncs = order.GetSize();
    polyOrder.Resize( numFncs, 3 );
    for( UInt i = 0; i < numFncs; ++i ) {
      for(UInt j = 0; j < 3; ++j ) {
        polyOrder[i][j] = order[i].GetMaxOrder(j);
        
      }
    }
  }
  
  template<typename T_SCAL, typename T_VEC>
  void FeH1HiWedge::_CalcShFnc( const T_SCAL x,  const T_SCAL y, const T_SCAL z, 
                                const Elem * elem,
                                T_VEC& ret ) {

//    std::cerr << "\n-------------------\n"
//              << " ELEM " << elem->elemNum << "\n"
//              << "-------------------\n";
    
    T_SCAL lambda[6] = { 1.0 - x - y, x,  y, 
                         1.0 - x - y, x,  y };
    T_SCAL mu[6]     = { 0.5 * (1.0-z), 0.5 * (1.0-z), 0.5 * (1.0-z),
                         0.5 * (1.0+z), 0.5 * (1.0+z), 0.5 * (1.0+z) };  
    
    UInt pos = 0;
    ret.Resize(actNumFncs_);

    // 1) Vertex shape functions
    for( UInt i = 0; i < 6; ++i ) {
      ret[i] = lambda[i] * mu[i];
    }

    pos = 6;
    
#ifdef USE_EDGES
    // 2) Edge shape functions

    // a) horizontal edges (triangular edges)
    for( UInt i = 0; i < 6; ++i ) {
      Double fac = elem->edges[i] < 0 ? -1.0 : 1.0;
      UInt order = orderEdge_[i];

      // if order of edge is below two, we leave
      if( order == 1 ) continue;
      UInt index1 = shape_.edgeVertices[i][0]-1;
      UInt index2 = shape_.edgeVertices[i][1]-1;

      // edge: parameterization of edge [-1;+1]
      // edgeNormal: parameterization of extension into element [-1;+1]
      T_SCAL edge  =  lambda[index2] -  lambda[index1]; 
      T_SCAL edgeNormal = lambda[index2] + lambda[index1];

      T_VEC vals;
      ScaledIntLegendreP2<T_SCAL,T_VEC>( vals, order, edgeNormal, fac*edge );
      for( UInt j = 0; j < order-1; ++j ) {
        ret[pos++] = vals[j] * mu[index1];
      } 
    }
    
    // b) vertical edges (quadrilateral related edges)
    for( UInt i = 6; i < 9; ++i ) {
      Double fac = elem->edges[i] < 0 ? -1.0 : 1.0;
      UInt order = orderEdge_[i];

      // if order of edge is below two, we leave
      if( order == 1 ) continue;
      UInt index1 = shape_.edgeVertices[i][0]-1;
      UInt index2 = shape_.edgeVertices[i][1]-1;

      // zeta: parameterization of edge [-1;+1]
      T_SCAL zeta  =  mu[index2] -  mu[index1]; 

      T_VEC vals;
      IntLegendreP2<T_SCAL,T_VEC>( vals, order, fac*zeta );

      for( UInt i = 0; i < order-1; ++i ) {
        ret[pos++] =  lambda[index1]* vals[i];
      } 
    }
//    std::cerr << "pos = " << pos << std::endl;
//    std::cerr << "size of functions: " << actNumFncs_;
#endif 
    // 3) Face shape functions
#ifdef USE_FACES
    // a) horizontal faces (triangular faces top/bottom)
    for( UInt i = 0; i < 2; ++i ) {
      UInt order = orderFace_[i][0]; 
      if( order < 3) continue;
      
      // obtain nodal permutation, i.e. sort the indices in ascending order
      const StdVector<UInt>& unsorted = shape_.faceNodes[i];
      StdVector<UInt> ind;
      Face::GetSortedIndices( ind, unsorted, 3, elem->faceFlags[i]);
      // take inner shape functions of triangle and
      // extend them via mu-lifting function
      UInt nFct =  TriaInnerLegendre( ret, pos, order, 
                                      lambda[ind[0]], lambda[ind[1]],
                                      lambda[ind[2]]);
      for( UInt j = 0; j < nFct; ++j ) {
        ret[pos] = ret[pos] * mu[ind[0]];
        pos++;
      }
    }
    
    // b) vertical faces (quadrilateral faces)
    for( UInt i = 2; i < 5; ++i ) {
     
      if( orderFace_[i][0] < 2 || orderFace_[i][1] < 2) continue;
      
      const StdVector<UInt>& unsorted = shape_.faceNodes[i];
      
      StdVector<UInt> ind;
      
      // obtain nodal permutation, i.e. sort the indices in ascending order
      Face::GetSortedIndices( ind, unsorted, 4, elem->faceFlags[i]);

      T_VEC horiz, vert;
      // we have to determine 2 things:
      // 1) Determine, if first edge in sorted array is horizontal or vertical one
      if( shape_.nodeCoords[ind[0]][2] == shape_.nodeCoords[ind[1]][2] ) {
        // edge [0] -> [1] is the horizontal one with order p[0] 
        // edge [0] -> [3] is the vertical one  with order p[1]

        // triangular shape function on horizontal edge 
        ScaledIntLegendreP2( horiz, orderFace_[i][0], 
                             lambda[ind[1]]+lambda[ind[0]],
                             lambda[ind[1]]-lambda[ind[0]] );

        // exension on vertical edge with p[1]
        IntLegendreP2( vert, orderFace_[i][1], mu[ind[0]]*2.0-1.0 );
        
        for( UInt j = 0;  j < orderFace_[i][0]-1; ++j ) {
          for( UInt k = 0;  k < orderFace_[i][1]-1; ++k ) {
            ret[pos++] = horiz[j] * vert[k];
          }
        }

      } else {
        // edge [0] -> [1] is the vertical one with order p[0]
        // edge [0] -> [3] is the horizontal one with order p[1]
        
        // triangular shape function on horizontal edge
        ScaledIntLegendreP2( horiz, orderFace_[i][1], 
                             lambda[ind[3]]+lambda[ind[0]],
                             lambda[ind[3]]-lambda[ind[0]] );

        // extension on vertical edge with p[0]
        IntLegendreP2( vert, orderFace_[i][0], mu[ind[0]]*2.0-1.0 );
        
        for( UInt k = 0;  k < orderFace_[i][0]-1; ++k ) {
          for( UInt j = 0;  j < orderFace_[i][1]-1; ++j ) {
            ret[pos++] = vert[k] * horiz[j];
          }
        }
      }
    } // faced (triangular)
#endif

    // 4) Inner shape functions
#ifdef USE_INNER
    if( orderInner_[0] > 2 && 
        orderInner_[2] > 1 ) {
      
      T_VEC triaVals(orderInner_[0]*orderInner_[0]);
      UInt nFct =  TriaInnerLegendre( triaVals, 0, orderInner_[0], 
                                      lambda[0], lambda[1],
                                      lambda[2]);
      T_VEC zetaVals;
      IntLegendreP2( zetaVals, orderInner_[2], mu[0]*2.0-1.0 );
      for( UInt i = 0; i < nFct; ++i ) {
        for( UInt k = 0; k < orderInner_[2]-1; ++k ) {
          ret[pos++] = triaVals[i] * zetaVals[k];
        } //i 
      } //j 
    } // if
#endif
    
    // additional check
    if( ret.GetSize() != pos ) {
      EXCEPTION("Should have " << ret.GetSize()
                << " unknowns but just got " << pos << " functions");
    }
  }

  void FeH1HiWedge::CalcNumUnknowns( ) {
    LOG_DBG(feH1Hi) << "CalcNumUnknowns for element "
        << Elem::feType.ToString(feType_);

    actNumFncs_ = 0;

    // Vertices
    StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
    vertFncs.Resize(6);
    vertFncs.Init(1); // Vertices have always order 1
    actNumFncs_ += 6;

    // Edges
    StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
    edgeFncs.Resize(9);
    edgeFncs.Init(0);
    UInt unknowns = 0;
#ifdef USE_EDGES
    for( UInt i = 0; i < 9; ++i ) {
      unknowns = (orderEdge_[i]-1);
      edgeFncs[i] = unknowns;
      LOG_DBG(feH1Hi) <<   "edge " << i+1 << " has order" <<  orderEdge_[i]-1
          << " and " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
#endif

    // Faces
    StdVector<UInt>& faceFncs = entityFncs_[FACE];
    faceFncs.Resize(5);
    faceFncs.Init(0);
#ifdef USE_FACES
    if( orderFace_[0][0] > 1 ) {
      for( UInt i = 0; i < 5; ++i ) {
        // check for triangular face
        if( i < 2 ) {
          unknowns = (orderFace_[i][0]-2) * (orderFace_[0][0]-1) / 2;
        } else {
          unknowns = (orderInner_[0]-1) * (orderInner_[1]-1);
//          unknowns = 0;
         
        }
        faceFncs[i] = unknowns;
        LOG_DBG(feH1Hi) << "face 0 has " << unknowns << "unknowns";
        actNumFncs_ += unknowns;
      }
    }
#endif
    // Interior
    StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
    innerFncs.Resize(1);
    innerFncs.Init(0);
#ifdef USE_INNER
    if( orderInner_[0] > 2 && orderInner_[2] > 1 ) {
      UInt unknowns = UInt ((orderInner_[0]-1) * (orderInner_[0]-2) * 0.5 
                    * (orderInner_[2]-1) ); 
    innerFncs.Init( unknowns );
    LOG_DBG(feH1Hi) << "interior has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
    }
#endif

    LOG_DBG(feH1Hi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;
    updateUnknowns_ = false;
  }
  
} // namespace CoupledField
