// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode c
#include "HCurlElemsHi.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Utils/autodiff.hh"
#include "polynomials.hh"

namespace CoupledField {

#define USE_FACES 1
#define USE_INNER 1

// declare class specific logging stream
DECLARE_LOG(feHCurlHi)
DEFINE_LOG(feHCurlHi, "feHCurlHi")


//! This method makes uses of the following vector identity
//!
//!     curl(b*F) = b x curl(F) + grad(b) x F
//!
//! If we further assume, that F=grad(b), we arrive at
//!
//!     curl(b*grad(a)) = b x curl(grad(a)) + grad(b) x grad(a)
//!
//! As the first term is 0, we get
//!
//!     curl(b*(grad(a)) = grad(b) x grad(a) = Cross(b,a)
//!
//! This second term is calculated within this method.

inline AutoDiff<Double,3> Cross (const AutoDiff<Double,3> & u,
                                 const AutoDiff<Double,3> & v)
  {
    AutoDiff<Double,3> hv;
    hv.Val() = 0.0;
    hv.DVal(0) = u.DVal(1)*v.DVal(2)-u.DVal(2)*v.DVal(1);
    hv.DVal(1) = u.DVal(2)*v.DVal(0)-u.DVal(0)*v.DVal(2);
    hv.DVal(2) = u.DVal(0)*v.DVal(1)-u.DVal(1)*v.DVal(0);
    return hv;
  }


FeHCurlHi::FeHCurlHi() {
  updateUnknowns_ = true;
  onlyLowestOrder_ = false;
  isoOrder_ = 0;
  
  // disable by default use of gradient functions
  useGrad_[EDGE]     = false;
  useGrad_[FACE]     = false;
  useGrad_[INTERIOR] = false;
}
FeHCurlHi::~FeHCurlHi() {

}

void FeHCurlHi::GetNumFncs( StdVector<UInt>& numFcns,
                              EntityType entityType,
                              UInt dof) {
  if(updateUnknowns_) CalcNumUnknowns();
  numFcns = entityFncs_[entityType]; 
}
  
void FeHCurlHi::UseGradient(EntityType entity, bool usage) {
  useGrad_[entity] = usage;
  updateUnknowns_ = true;
}

void FeHCurlHi::SetIsoOrder(UInt order) {
  orderEdge_.Resize(shape_.numEdges);
  orderFace_.Resize(shape_.numFaces);

  // set order for edges
  orderEdge_.Init(order);

  // set order for faces
  array<UInt,2> faceOrder = {{order, order}};
  orderFace_.Init(faceOrder);

  // set order for inner
  array<UInt, 3> innerOrder = {{order, order, order}}; 
  orderInner_ = innerOrder;

  updateUnknowns_ = true;
  isIsotropic_ = true;
  isoOrder_ = order;
}

void FeHCurlHi::GetNodalPermutation( StdVector<UInt>& fncPermutation,
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

UInt FeHCurlHi::GetIsoOrder() const {
    if( isIsotropic_) {
      return isoOrder_;
    } else {
      EXCEPTION("Implement me");
      return 0;
    }
}

UInt FeHCurlHi::GetMaxOrder() const {
  if( isIsotropic_) {
    return isoOrder_;
  } else {
    EXCEPTION("Implement me");
    return 0;
  }
}

void FeHCurlHi::GetMaxOrderLocDir(StdVector<UInt>& order ) {
  if( isIsotropic_ ) {
    order.Resize( Elem::shapes[feType_].dim);
    order.Init(isoOrder_);
  } else {
    EXCEPTION("Implement me");
  }
}

// ========================================================================
//  FeHCurlHi explicit element definition 
// ========================================================================

// =======================
//  QUADRILATERAL ELEMENT 
// =======================

FeHCurlHiQuad::FeHCurlHiQuad() {
  feType_ = Elem::ET_QUAD4;
  shape_ = Elem::shapes[feType_];
}


FeHCurlHiQuad::~FeHCurlHiQuad() {
  
}

void FeHCurlHiQuad::CalcNumUnknowns() {
  LOG_DBG(feHCurlHi) << "CalcNumUnknowns for element "
      << Elem::feType.ToString(feType_);

  actNumFncs_ = 0;

  // Vertices 
  StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
  vertFncs.Resize(shape_.numVertices);
  vertFncs.Init(0); // -> no unknowns on vertices
  
  // Edges
  StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
  edgeFncs.Resize(shape_.numEdges);
  UInt unknowns = 0;
  for( UInt i = 0; i < shape_.numEdges; ++i ) {
    unknowns = 1; // Lowest order Nedelc functions
    if( useGrad_[EDGE])
      unknowns += orderEdge_[i];
    edgeFncs[i] = unknowns;
    LOG_DBG(feHCurlHi) <<   "edge " << i+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
  }

  // Faces
  StdVector<UInt>& faceFncs = entityFncs_[FACE];
  faceFncs.Resize(shape_.numFaces);
  faceFncs.Init(0);
  for( UInt i = 0; i < shape_.numFaces; ++i ) {
    if( orderFace_[i][0] > 0 &&
        orderFace_[i][1] > 0 ) {
      unknowns = orderFace_[i][0] * orderFace_[i][1] // face functions of 1st kind
                + orderFace_[i][0] + orderFace_[i][1];
      if( useGrad_[FACE])
        unknowns +=  orderFace_[i][0] * orderFace_[i][1];
      faceFncs[i] = unknowns;
      LOG_DBG(feHCurlHi) << "face " << i+1 << " has " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
  }

  // Interior
  StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
  innerFncs.Resize(1);
  innerFncs.Init(0);
  LOG_DBG(feHCurlHi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;
  updateUnknowns_ = false;

}
  
void FeHCurlHiQuad::CalcLocShFnc( Matrix<Double>& shape, LocPointMapped& lp,
                              const Elem* elem, UInt comp  ) {
  if (updateUnknowns_) CalcNumUnknowns();
  EXCEPTION("Implement me");
  
}

void FeHCurlHiQuad::CalcLocCurlShFnc( Matrix<Double>& curl, LocPointMapped& lp,
                                     const Elem* elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();
  EXCEPTION("Implement me");
}


// =======================
//  HEXAHEDRAL ELEMENT 
// =======================
FeHCurlHiHex::FeHCurlHiHex() {
  feType_ = Elem::ET_HEXA8;
  shape_ = Elem::shapes[feType_];
}


//! Destructor
FeHCurlHiHex::~FeHCurlHiHex() {

}

void FeHCurlHiHex::CalcNumUnknowns() {
  LOG_DBG(feHCurlHi) << "CalcNumUnknowns for element "
      << Elem::feType.ToString(feType_);

  actNumFncs_ = 0;

  // Vertices 
  StdVector<UInt>& vertFncs = entityFncs_[VERTEX];
  vertFncs.Resize(shape_.numVertices);
  vertFncs.Init(0); // -> no unknowns on vertices

  // Edges
  StdVector<UInt>& edgeFncs = entityFncs_[EDGE];
  edgeFncs.Resize(shape_.numEdges);
  UInt unknowns = 0;
  for( UInt i = 0; i < shape_.numEdges; ++i ) {
    unknowns = 1; // Lowest order Nedelc functions
    if( useGrad_[EDGE])
      unknowns += orderEdge_[i];
    edgeFncs[i] = unknowns;
    LOG_DBG(feHCurlHi) <<   "edge " << i+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
  }

  // Faces
  StdVector<UInt>& faceFncs = entityFncs_[FACE];
  faceFncs.Resize(shape_.numFaces);
  faceFncs.Init(0);
#ifdef USE_FACES
  for( UInt i = 0; i < shape_.numFaces; ++i ) {
    if( orderFace_[i][0] > 0 &&
        orderFace_[i][1] > 0 ) {
      unknowns = orderFace_[i][0] * orderFace_[i][1] // face functions of 1st kind
                + orderFace_[i][0] + orderFace_[i][1];
      if( useGrad_[FACE])
        unknowns +=  orderFace_[i][0] * orderFace_[i][1];
      faceFncs[i] = unknowns;
      LOG_DBG(feHCurlHi) << "face " << i+1 << " has " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
  }
#endif

  // Interior
  StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
  innerFncs.Resize(1);
  innerFncs.Init(0);

  #ifdef USE_INNER
  if( orderInner_[0] > 0 && 
      orderInner_[1] > 0 && 
      orderInner_[2] > 0 ) {

    unknowns = 2 * (orderInner_[0] * orderInner_[1] * orderInner_[2]) 
                   + orderInner_[1] * orderInner_[2] 
                   + orderInner_[0] * (orderInner_[2] + orderInner_[1]);
    if( useGrad_[INTERIOR] ) { 
      unknowns += orderInner_[0] * orderInner_[1] * orderInner_[2];
    }
    actNumFncs_ += unknowns;
    innerFncs[0] = unknowns;
    LOG_DBG(feHCurlHi) << "interior has " << unknowns << "unknowns";
  }
#endif

  LOG_DBG(feHCurlHi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;
  updateUnknowns_ = false;
}


void FeHCurlHiHex::CalcLocShFnc( Matrix<Double>& shape, LocPointMapped& lpm,
                             const Elem* elem, UInt comp  ) {
  if (updateUnknowns_) CalcNumUnknowns();
  
  AutoDiff<Double, 3> x (lpm.lp.coord[0],0);
  AutoDiff<Double, 3> y (lpm.lp.coord[1],1);
  AutoDiff<Double, 3> z (lpm.lp.coord[2],2);
  
  AutoDiff<Double, 3> lambda[8] = {0.125*(1.0-x)*(1.0-y)*(1.0-z),
                                   0.125*(1.0+x)*(1.0-y)*(1.0-z),
                                   0.125*(1.0+x)*(1.0+y)*(1.0-z),
                                   0.125*(1.0-x)*(1.0+y)*(1.0-z),
                                   0.125*(1.0-x)*(1.0-y)*(1.0+z),
                                   0.125*(1.0+x)*(1.0-y)*(1.0+z),
                                   0.125*(1.0+x)*(1.0+y)*(1.0+z),
                                   0.125*(1.0-x)*(1.0+y)*(1.0+z)};

  AutoDiff<Double, 3> sigma[8]  = {0.5*((1.0-x)+(1.0-y)+(1.0-z)),
                                   0.5*((1.0+x)+(1.0-y)+(1.0-z)),
                                   0.5*((1.0+x)+(1.0+y)+(1.0-z)),
                                   0.5*((1.0-x)+(1.0+y)+(1.0-z)),
                                   0.5*((1.0-x)+(1.0-y)+(1.0+z)),
                                   0.5*((1.0+x)+(1.0-y)+(1.0+z)),
                                   0.5*((1.0+x)+(1.0+y)+(1.0+z)),
                                   0.5*((1.0-x)+(1.0+y)+(1.0+z))};
  UInt pos = 12;
  shape.Resize(3,actNumFncs_);
  shape.Init();
  
  StdVector<AutoDiff<Double, 3> > xiVals, etaVals, zetaVals;
  
  // ------------------------
  // 1) Edge shape functions
  // ------------------------
  for( UInt i = 0; i < 12; ++i ) {

    UInt order = orderEdge_[i];
    Double fac = elem->edges[i] < 0 ? -1.0 : 1.0;
    UInt index1 = shape_.edgeVertices[i][0]-1;
    UInt index2 = shape_.edgeVertices[i][1]-1;

    // xi: parameterization of edge [-1;+1]
    // eta: parameterization of extension into element [-1;+1]
    AutoDiff<Double, 3> xi  =  sigma[index2] -  sigma[index1]; 
    AutoDiff<Double, 3> eta = lambda[index1] + lambda[index2];

    // a) standard Nedelec shape functions
    for( UInt k = 0; k < 3; ++k ) {
      shape[k][i] = xi.DVal(k)*eta.Val()*fac;
    }
    // b) gradient functions
    if( useGrad_[EDGE] && !onlyLowestOrder_ ) {
      IntLegendreP2(xiVals, order+1, fac*xi );
      
      for( UInt j = 0; j < order; ++j ) {
        for( UInt k = 0; k < 3; ++k ) {
          shape[k][pos] = eta.DVal(k)* xiVals[j].Val()
                          + eta.Val() * xiVals[j].DVal(k);
        }
        pos++;
      }
    }
  }
  if(onlyLowestOrder_) return;
  
  // ------------------------
  // 2) Face shape functions
  // ------------------------
#ifdef USE_FACES
  for( UInt i = 0; i < 6; ++i ) {
    UInt order1 = orderFace_[i][0];
    UInt order2 = orderFace_[i][1];
    if (order1 > 0 && order2 > 0 ) {
      // calculate face extension parameter which is the sum
             // of all lambdas of one face
      AutoDiff<Double,3> sum_lambda = 0.0;
      for( UInt j = 0; j < 4; ++j)
        sum_lambda += lambda[shape_.faceVertices[i][j]-1];
      
#define NETGEN_VERSION
#ifdef NETGEN_VERSION
      int qmin = 0;
      for (UInt p = 1; p < 4; p++)
        if ( elem->connect[shape_.faceVertices[i][p]-1] <
             elem->connect[shape_.faceVertices[i][qmin]-1] ) {
          qmin = p;
        }

      int q1 = (qmin+3)%4; 
      int q2 = (qmin+1)%4; 

      if(elem->connect[shape_.faceVertices[i][q2]-1] <
          elem->connect[shape_.faceVertices[i][q1]-1] )
        std::swap(q1,q2);  // fmax > f1 > f2

      // horiz = 1 if sigma_{fmax}-sigma_{f1} is horizontal coordinate 
      // for anisotropic meshes (vertical(esp. z) is the longer one) 
      double horiz=1.; 
      if( (qmin & 2) == (q2 & 2)) 
        horiz = -1.; 

      int fmax = shape_.faceVertices[i][qmin]-1; 
      int f1 = shape_.faceVertices[i][q1]-1; 
      int f2 = shape_.faceVertices[i][q2]-1; 

      AutoDiff<Double,3> xi = sigma[fmax]-sigma[f1]; 
      AutoDiff<Double,3> eta = sigma[fmax]-sigma[f2];
#endif
//#define OWN_VERSION
#ifdef OWN_VERSION
      AutoDiff<Double,3> xi  =   sigma[shape_.faceVertices[i][1]-1]
                               - sigma[shape_.faceVertices[i][0]-1];
      AutoDiff<Double,3> eta  =  sigma[shape_.faceVertices[i][3]-1]
                               - sigma[shape_.faceVertices[i][0]-1];
      // Test for orientation:
      // 1) Check, if local surface directions are rotated (90degree or 
      //    270 degree) by testing the MSB bit of the faceFlags.
      // 2) Sine we always mutltiply two functions for a surface function, we
      //    just have to consider one sign, which is the product of the signs of
      //    both surface directions. This sign is incorporated into the 
      //    xi-variable.
      if( elem->faceFlags[i].test(2) == false) {
        std::swap(order1, order2);
        std::swap(xi, eta);
      }
      // Model a XOR B = ((a || b) && !(a && b) 
//      if ( (elem->faceFlags[i].test(0) || (elem->faceFlags[i].test(1))) && 
//          !(elem->faceFlags[i].test(0) && (elem->faceFlags[i].test(1))) ) {
//        xi *= -1.0;
//      }
      if ( (elem->faceFlags[i].test(0)) && 
          !(elem->faceFlags[i].test(0) && (elem->faceFlags[i].test(1))) ) {
        xi *= -1.0;
      }
      if ( (elem->faceFlags[i].test(1)) && 
          !(elem->faceFlags[i].test(0) && (elem->faceFlags[i].test(1))) ) {
        eta *= -1.0;
      }

#endif
      
      IntLegendreP2(xiVals, order1+1, xi );
      IntLegendreP2(etaVals, order2+1, eta );
      
      // a) gradient fields
      if( useGrad_[FACE])
        WARN("Calculation of face gradient fields not yet implemented");
      
      // b) curl of gradient fields
      for( UInt k = 0; k < order1; ++k ) {
        for( UInt j = 0; j < order2; ++j ) {
          for( UInt l = 0; l < 3; ++l) {
//            shape[l][pos] =
//                sum_lambda * (xiVals[k].DVal(l) * etaVals[j].Val()* -
//                              xiVals[k].Val() * etVals[j]*DVal(l));
            // alternative variant:
            shape[l][pos] =
                sum_lambda.Val() * (xiVals[k].DVal(l) * etaVals[j].Val()-
                                    xiVals[k].Val() * etaVals[j].DVal(l))
                                    + horiz * sum_lambda.DVal(l)*xiVals[k].Val()*etaVals[j].Val();
          }
          pos++;
        }
      }
        
      // c) remaining functions
      for( UInt k = 0; k < order1; ++k ) {
        for( UInt l = 0; l < 3; ++l) {
          shape[l][pos] = eta.DVal(l) * xiVals[k].Val() * sum_lambda.Val();
        }
        pos++;
      }
      
      for( UInt k = 0; k < order1; ++k ) {
        for( UInt l = 0; l < 3; ++l) {
          shape[l][pos] =  xi.DVal(l) * etaVals[k].Val() * sum_lambda.Val();
        }
        pos++;
      }
    }
  }
#endif
  // ----------------------------
  // 3) Interior shape functions
  // ----------------------------
#ifdef USE_INNER
  IntLegendreP2(xiVals,   orderInner_[0]+1, x );
  IntLegendreP2(etaVals,  orderInner_[1]+1, y );
  IntLegendreP2(zetaVals, orderInner_[2]+1, z );

  // a) gradient fiels
  if( useGrad_[INTERIOR]) {
    for( UInt i = 0; i < orderInner_[0]; ++i ) {
      for( UInt j = 0; j < orderInner_[1]; ++j ) {
        for( UInt k = 0; k < orderInner_[2]; ++k ) {
          shape[0][pos] = xiVals[i].DVal(0) * etaVals[j].Val()   * zetaVals[k].Val();
          shape[1][pos] = xiVals[i].Val()   * etaVals[j].DVal(1) * zetaVals[k].Val();
          shape[2][pos] = xiVals[i].Val()   * etaVals[j].Val()   * zetaVals[k].DVal(2);
        }
        pos++;
      }
    }
  } // useGrad_

  // b) rotation of gradient fields
  for( UInt i = 0; i < orderInner_[0]; ++i ) {
    for( UInt j = 0; j < orderInner_[1]; ++j ) {
      for( UInt k = 0; k < orderInner_[2]; ++k ) {
        shape[0][pos] =  xiVals[i].DVal(0) * etaVals[j].Val()   * zetaVals[k].Val();
        shape[1][pos] = -xiVals[i].Val()   * etaVals[j].DVal(1) * zetaVals[k].Val();
        shape[2][pos] =  xiVals[i].Val()   * etaVals[j].Val()   * zetaVals[k].DVal(2);
        pos++;
        shape[0][pos] =  xiVals[i].DVal(0) * etaVals[j].Val()   * zetaVals[k].Val();
        shape[1][pos] =  xiVals[i].Val()   * etaVals[j].DVal(1) * zetaVals[k].Val();
        shape[2][pos] = -xiVals[i].Val()   * etaVals[j].Val()   * zetaVals[k].DVal(2);
        pos++;
      }
    }
  }

  // c) remaining functions
  for( UInt i = 0; i < orderInner_[0]; ++i ) {
    for( UInt j = 0; j < orderInner_[1]; ++j ) {
      shape[2][pos++] = xiVals[i].Val() * etaVals[j].Val();
    }
  }
  for( UInt i = 0; i < orderInner_[0]; ++i ) {
    for( UInt k = 0; k < orderInner_[2]; ++k ) {
      shape[1][pos++] = xiVals[i].Val() * zetaVals[k].Val();
    }
  }
  for( UInt j = 0; j < orderInner_[1]; ++j ) {
    for( UInt k = 0; k < orderInner_[2]; ++k ) {
      shape[0][pos++] = etaVals[j].Val() * zetaVals[k].Val();
    }
  }

#endif
}

void FeHCurlHiHex::CalcLocCurlShFnc( Matrix<Double>& curl, LocPointMapped& lpm,
                                     const Elem* elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();

  AutoDiff<Double, 3> x (lpm.lp.coord[0],0);
  AutoDiff<Double, 3> y (lpm.lp.coord[1],1);
  AutoDiff<Double, 3> z (lpm.lp.coord[2],2);
  
  AutoDiff<Double, 3> lambda[8] = {0.125*(1.0-x)*(1.0-y)*(1.0-z),
                                   0.125*(1.0+x)*(1.0-y)*(1.0-z),
                                   0.125*(1.0+x)*(1.0+y)*(1.0-z),
                                   0.125*(1.0-x)*(1.0+y)*(1.0-z),
                                   0.125*(1.0-x)*(1.0-y)*(1.0+z),
                                   0.125*(1.0+x)*(1.0-y)*(1.0+z),
                                   0.125*(1.0+x)*(1.0+y)*(1.0+z),
                                   0.125*(1.0-x)*(1.0+y)*(1.0+z)};

  AutoDiff<Double, 3> sigma[8]  = {0.5*((1.0-x)+(1.0-y)+(1.0-z)),
                                   0.5*((1.0+x)+(1.0-y)+(1.0-z)),
                                   0.5*((1.0+x)+(1.0+y)+(1.0-z)),
                                   0.5*((1.0-x)+(1.0+y)+(1.0-z)),
                                   0.5*((1.0-x)+(1.0-y)+(1.0+z)),
                                   0.5*((1.0+x)+(1.0-y)+(1.0+z)),
                                   0.5*((1.0+x)+(1.0+y)+(1.0+z)),
                                   0.5*((1.0-x)+(1.0+y)+(1.0+z))};

  StdVector<AutoDiff<Double, 3> > xiVals, etaVals, zetaVals;
  
  UInt pos = 12;
  curl.Resize(3,actNumFncs_);
  curl.Init();
  
  // -------------------------
  // 1) Edge shape functions
  // -------------------------
  for( UInt i = 0; i < 12; ++i) {
    UInt order = orderEdge_[i];
    Double fac = elem->edges[i] < 0 ? -1.0 : 1.0;
    UInt index1 = shape_.edgeVertices[i][0]-1;
    UInt index2 = shape_.edgeVertices[i][1]-1;

    // xi: parameterization of edge [-1;+1]
    // eta: parameterization of extension into element [-1;+1]
    AutoDiff<Double, 3> xi  =  sigma[index2] -  sigma[index1]; 
    AutoDiff<Double, 3> eta = lambda[index1] + lambda[index2];

    // a) standard Nedelec shape functions
    AutoDiff<Double, 3> temp =  Cross(eta, fac * xi);
    for( UInt j = 0; j < 3; ++j ) {
      curl[j][i] = temp.DVal(j);
    }
    // b) gradient functions -> get skipped
    if( useGrad_[EDGE] ) {
      pos += order;
    }
  }
  
  if(onlyLowestOrder_) return;
  // -------------------------
  // 2) Face shape functions
  // -------------------------

  bool print = false;
 //if ( elem->elemNum == 20 ) print = true;

  #ifdef USE_FACES
  for( UInt f = 0; f < 6; ++f ) {
    
    if(print) {
      std::cerr << "\n===========================\n"
                 <<" Face #" << f << " of elem# " << elem->elemNum << std::endl
                 << "===========================\n";
    }
    UInt order1 = orderFace_[f][0];
    UInt order2 = orderFace_[f][1];
    if (order1 > 0 && order2 > 0 ) {
      // calculate face extension parameter which is the sum
      // of all lambdas of one face
      AutoDiff<Double,3> sum_lambda = 0.0;
      for( UInt j = 0; j < 4; ++j)
        sum_lambda += lambda[shape_.faceVertices[f][j]-1];
      

#ifdef OWN_VERSION
      AutoDiff<Double,3> xi  =  sigma[shape_.faceVertices[f][1]-1]
                              - sigma[shape_.faceVertices[f][0]-1];
      AutoDiff<Double,3> eta  =  sigma[shape_.faceVertices[f][3]-1]
                               - sigma[shape_.faceVertices[f][0]-1];

//      // test, if xi and eta get reversed
      if( print ) {
        std::cerr << "nodes = "  
            << elem->connect[shape_.faceVertices[f][0]-1] << ", "
            << elem->connect[shape_.faceVertices[f][1]-1] << ", "
            << elem->connect[shape_.faceVertices[f][2]-1] << ", "
            << elem->connect[shape_.faceVertices[f][3]-1] << "\n";
        std::cerr << "local coords: " << lpm.lp.coord.ToString() << std::endl;
        std::cerr << "faceFlags are " << elem->faceFlags[f] << std::endl;; 
        std::cerr << "I = " << xi << std::endl;
        std::cerr << "II = " << eta << std::endl << std::endl;
      }
      
      // Test for orientation:
      // 1) Check, if local surface directions are rotated (90degree or 
      //    270 degree) by testing the MSB bit of the faceFlags.
      // 2) Sine we always mutltiply two functions for a surface function, we
      //    just have to consider one sign, which is the product of the signs of
      //    both surface directions. This sign is incorporated into the 
      //    xi-variable.
      if( elem->faceFlags[f].test(2) == false) {
        std::swap(order1, order2);
        std::swap(xi, eta);
      }
      // Model a XOR B = ((a || b) && !(a && b) 
//      if ( (elem->faceFlags[f].test(0) || (elem->faceFlags[f].test(1))) && 
//          !(elem->faceFlags[f].test(0) && (elem->faceFlags[f].test(1))) ) {
//        xi *= -1.0;
//      }
      if ( (elem->faceFlags[f].test(0) ) && 
          !(elem->faceFlags[f].test(0) && (elem->faceFlags[f].test(1))) ) {
        xi *= -1.0;
      }
      if ( (elem->faceFlags[f].test(1)) && 
          !(elem->faceFlags[f].test(0) && (elem->faceFlags[f].test(1))) ) {
        eta *= -1.0;
      }
      if(print )  {
        std::cerr << "xi    is " << xi << std::endl;
        std::cerr << "eta    is " << eta << std::endl;
      }
#endif

#define NETGEN_VERSION
#ifdef NETGEN_VERSION
      int qmin = 0;
      for (UInt p = 0; p < 4; p++)
        if ( elem->connect[shape_.faceVertices[f][p]-1] <
             elem->connect[shape_.faceVertices[f][qmin]-1] ) {
          qmin = p;
        }

      int q1 = (qmin+3)%4; 
      int q2 = (qmin+1)%4; 

      if(elem->connect[shape_.faceVertices[f][q2]-1] <
          elem->connect[shape_.faceVertices[f][q1]-1] )
        std::swap(q1,q2);  // fmin < f1 < f2

      // horiz = 1 if sigma_{fmax}-sigma_{f1} is horizontal coordinate 
      // for anisotropic meshes (vertical(esp. z) is the longer one) 
      double horiz=1.; 
      if( (qmin & 2) == (q2 & 2)) 
        horiz = -1.; 

      int fmax = shape_.faceVertices[f][qmin]-1; 
      int f1 = shape_.faceVertices[f][q1]-1; 
      int f2 = shape_.faceVertices[f][q2]-1; 

      //AutoDiff<Double,3> xi = sigma[fmax]-sigma[f1]; 
      //AutoDiff<Double,3> eta = sigma[fmax]-sigma[f2];
      AutoDiff<Double,3> xi = sigma[fmax]-sigma[f1];
      AutoDiff<Double,3> eta = sigma[fmax]-sigma[f2];
      if(print )  {
        std::cerr << "xi_ng is " << xi << std::endl;
        std::cerr << "eta_ng is " << eta << std::endl;
      }
#endif

      
      IntLegendreP2(xiVals, order1+1, xi );
      IntLegendreP2(etaVals, order2+1,eta );

      // a) gradient fields
      if( useGrad_[FACE])
        WARN("Calculation of face gradient fields not yet implemented");

      // b) curl of gradient fields
      StdVector<AutoDiff<Double, 3> > xiLambdaVals(xiVals.GetSize()), 
          etaLambdaVals(etaVals.GetSize());
      
      for( UInt j = 0; j < order1; ++j ) {
        xiLambdaVals[j]= sum_lambda * xiVals[j];
      }
      for( UInt j = 0; j < order2; ++j ) {
        etaLambdaVals[j] = sum_lambda * etaVals[j];
      }

//      if( print ) {
//        std::cerr << "xiLambdaVals = " << xiLambdaVals << std::endl;
//        std::cerr << "etaLambdaVals = " << etaLambdaVals << std::endl;
//      }
      
    
//        for( UInt i = 0; i < order1; ++i ) {
//          for( UInt j = 0; j < order2; ++j ) {
//            AutoDiff<Double,3> temp1 = Cross(etaLambdaVals[j], xiVals[i]);
//            AutoDiff<Double,3> temp2 = Cross(xiLambdaVals[i], etaVals[j]);
//            for( UInt l = 0; l < 3; ++l) {  
//              curl[l][pos] = temp1.DVal(l) - temp2.DVal(l);
//            }
//            pos++;
//          }
//        } 
      if(horiz == 1 ) {
        for( UInt i = 0; i < order1; ++i ) {
          for( UInt j = 0; j < order2; ++j ) {
            AutoDiff<Double,3> temp = Cross(etaVals[j], xiLambdaVals[i]);
            for( UInt l = 0; l < 3; ++l) {  
              curl[l][pos] = temp.DVal(l);
            }
            pos++;
          }
        }
      } else {
        for( UInt i = 0; i < order1; ++i ) {
          for( UInt j = 0; j < order2; ++j ) {
            AutoDiff<Double,3> temp = Cross(etaLambdaVals[j], xiVals[i]);
            for( UInt l = 0; l < 3; ++l) {  
              curl[l][pos] = temp.DVal(l);
            }
            pos++;
          }
        }
      }
        

      // c) remaining functions
      for( UInt i = 0; i < order1; ++i ) {
        AutoDiff<Double,3> hd =  Cross(xiLambdaVals[i],eta);
        for( UInt l = 0; l < 3; ++l) {
          curl[l][pos] = hd.DVal(l);
        }
        pos++;
      }
      for( UInt j = 0; j < order2; ++j ) {
        AutoDiff<Double,3> hd =  Cross(etaLambdaVals[j],xi);
        for( UInt l = 0; l < 3; ++l) {
          curl[l][pos] = hd.DVal(l);
        }
        pos++;
      }
    } // if order
  } // loop over faces
#endif

  // ----------------------------
  // 3) Interior shape functions
  // ----------------------------

  
#ifdef USE_INNER
  IntLegendreP2(xiVals,   orderInner_[0]+1, x );
  IntLegendreP2(etaVals,  orderInner_[1]+1, y );
  IntLegendreP2(zetaVals, orderInner_[2]+1, z );

  // a) gradient fiels
  if( useGrad_[INTERIOR]) {
    pos += orderInner_[0] * orderInner_[1] * orderInner_[2];
  }

  // b) curl of gradient fields
  for( UInt i = 0; i < orderInner_[0]; ++i ) {
    for( UInt j = 0; j < orderInner_[1]; ++j ) {
      for( UInt k = 0; k < orderInner_[2]; ++k ) {

        AutoDiff<Double, 3> xiZetaVal = xiVals[i] * zetaVals[k];
        AutoDiff<Double, 3> tempCurl  =  Cross( etaVals[j], xiZetaVal);
        tempCurl *= 2;  
        for( UInt l = 0; l < 3; ++l ) {
          curl[l][pos] = tempCurl.DVal(l);
        }

        pos++;
        AutoDiff<Double, 3> xiEtaVal = xiVals[i] * etaVals[j];
        tempCurl = Cross(zetaVals[k], xiEtaVal);
        tempCurl *= 2;
        for( UInt l = 0; l < 3; ++l ) {
          curl[l][pos] = tempCurl.DVal(l);
        }
        pos++;
      }
    }
  }

  // c) remaining functions
  for( UInt i = 0; i < orderInner_[0]; ++i ) {
    for( UInt j = 0; j < orderInner_[1]; ++j ) {
      AutoDiff<Double, 3> xiEtaVal = xiVals[i] * etaVals[j];
      AutoDiff<Double, 3> tempCurl = Cross(xiEtaVal, z);
      for( UInt l = 0; l < 3; ++l ) {
        curl[l][pos] = tempCurl.DVal(l);
      }
      pos++;
    }
  }
  for( UInt i = 0; i < orderInner_[0]; ++i ) {
    for( UInt k = 0; k < orderInner_[2]; ++k ) {
      AutoDiff<Double, 3> xiZetaVal = xiVals[i] * zetaVals[k];
      AutoDiff<Double, 3> tempCurl = Cross(xiZetaVal, y);
      for( UInt l = 0; l < 3; ++l ) {
        curl[l][pos] = tempCurl.DVal(l);
      }
      pos++;
    }
  }
  for( UInt j = 0; j < orderInner_[1]; ++j ) {
    for( UInt k = 0; k < orderInner_[2]; ++k ) {
      AutoDiff<Double, 3> etaZetaVal = etaVals[j] * zetaVals[k];
      AutoDiff<Double, 3> tempCurl = Cross(etaZetaVal, x);
      for( UInt l = 0; l < 3; ++l ) {
        curl[l][pos] = tempCurl.DVal(l);
      }
      pos++;
    }
  }
#endif
  
//std::cerr << "local curl of matrix is\n" << curl << std::endl;
//std::cerr << "size of matrix is 3 x " << curl.GetNumCols() << std::endl;
//std::cerr << "\t-> last entry was " << pos-1 << std::endl;
}

}// end of namespace
