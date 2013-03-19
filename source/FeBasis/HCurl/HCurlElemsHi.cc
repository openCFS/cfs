#include "HCurlElemsHi.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/AutoDiff.hh"
#include "FeBasis/Polynomials.hh"
#include "Domain/ElemMapping/EdgeFace.hh"

namespace CoupledField {

#define USE_FACES 1
#define USE_INNER 1

// declare class specific logging stream
DECLARE_LOG(feHCurlHi)
DEFINE_LOG(feHCurlHi, "feHCurlHi")


// ===============================================================
//  VECTOR IDENTITIES NEEDED FOR THE CALCULATION OF CURL MATRICES
// ===============================================================
//! This method makes use of the following vector identity
//! (assuming b is a scalar and F a vector function):
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


// -------------------
//  Expr = grad(U)
// -------------------
template<int D, FeHCurlHi::DiffType DIFF>
class Xpr_GradU {
};

template<int D>
class Xpr_GradU<D, FeHCurlHi::ID> {
  public:
    Xpr_GradU(const AutoDiff<Double,D>& u) {
      for( UInt i = 0; i < D; ++i ) {
        val[i] = u.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return val[i];}
  private:
    Double val[D];
};

template<int D>
class Xpr_GradU<D, FeHCurlHi::CURL> {
  public:
    Xpr_GradU(const AutoDiff<Double,D>& u) {
      for( UInt i = 0; i < D; ++i ) {
        curl[i] = 0.0;
      }
    }
    Double operator[] (UInt i ) {return curl[i];}
  private:
    Double curl[D];
};

// -------------------
//  Expr = s*grad(U)
// -------------------

//! Auxiliary class for evaluating expressions and curls of the type
//!    s * grad(U)
template<int D, FeHCurlHi::DiffType DIFF>
class Xpr_SGradU {
};

template<int D>
class Xpr_SGradU<D, FeHCurlHi::ID> {
  public:
    Xpr_SGradU(const AutoDiff<Double,D>& s, 
                 const AutoDiff<Double,D>& u ) {
      for( UInt i = 0; i < D; ++i ) {
        val[i] = s.Val() * u.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return val[i];}
  private:
    Double val[D];
};

template<int D>
class Xpr_SGradU<D, FeHCurlHi::CURL> {
  public:
  Xpr_SGradU(const AutoDiff<Double,D>& s, 
               const AutoDiff<Double,D>& u ) {
    AutoDiff<Double,D> c = Cross(s,u);
      for( UInt i = 0; i < D; ++i ) {
        curl[i] = c.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return curl[i];}
  private:
    Double curl[D];
};


// -------------------------------------
//  Expr = grad(U)*V - U*grad(V)
// -------------------------------------

template<int D, FeHCurlHi::DiffType DIFF>
class Xpr_Diff_VGradU {
};

template<int D>
class Xpr_Diff_VGradU<D, FeHCurlHi::ID> {
  public:
  Xpr_Diff_VGradU(const AutoDiff<Double,D>& u, 
                  const AutoDiff<Double,D>& v ) {
      for( UInt i = 0; i < D; ++i ) {
        val[i] = u.Val() * v.DVal(i) - v.Val() * u.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return val[i];}
  private:
    Double val[D];
};

template<int D>
class Xpr_Diff_VGradU<D, FeHCurlHi::CURL> {
  public:
  Xpr_Diff_VGradU(const AutoDiff<Double,D>& u, 
                  const AutoDiff<Double,D>& v ) {
    AutoDiff<Double,D> c = Cross(u,v);
      for( UInt i = 0; i < D; ++i ) {
        curl[i] = 2 * c.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return curl[i];}
  private:
    Double curl[D];
};


// -------------------------------------
//  Expr = s( grad(U)*V - U*grad(V))
// -------------------------------------

template<int D, FeHCurlHi::DiffType DIFF>
class Xpr_Diff_SVGradU {
};

template<int D>
class Xpr_Diff_SVGradU<D, FeHCurlHi::ID> {
  public:
  Xpr_Diff_SVGradU(const AutoDiff<Double,D>& u, 
                   const AutoDiff<Double,D>& v,
                   const AutoDiff<Double,D>& s) {
      for( UInt i = 0; i < D; ++i ) {
        val[i] = s.Val() * (u.Val() * v.DVal(i) - v.Val() * u.DVal(i) );
      }
    }
    Double operator[] (UInt i ) {return val[i];}
  private:
    Double val[D];
};

template<int D>
class Xpr_Diff_SVGradU<D, FeHCurlHi::CURL> {
  public:
  Xpr_Diff_SVGradU(const AutoDiff<Double,D>& u, 
                   const AutoDiff<Double,D>& v,
                   const AutoDiff<Double,D>& s) {
    AutoDiff<Double,D> c = Cross(u*s, v) + Cross(u, s*v);
      for( UInt i = 0; i < D; ++i ) {
        curl[i] = c.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return curl[i];}
  private:
    Double curl[D];
};

// ===============================================================
// ELEMENTS SECTION
// ===============================================================



FeHCurlHi::FeHCurlHi(Elem::FEType feType ) 
 :  FeHCurl(), FeHi(feType ){
  feType_ = feType;
  shape_ = Elem::shapes[feType];
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

UInt FeHCurlHi::GetNumFncs( ) {
  if(updateUnknowns_) this->CalcNumUnknowns();
  return actNumFncs_;
}

void FeHCurlHi::GetNumFncs( StdVector<UInt>& numFcns,
                              EntityType entityType,
                              UInt dof) {
  if(updateUnknowns_) this->CalcNumUnknowns();
  numFcns = entityFncs_[entityType]; 
}
  
void FeHCurlHi::UseGradient(EntityType entity, bool usage) {
  useGrad_[entity] = usage;
  updateUnknowns_ = true;
}

void FeHCurlHi::GetNodalPermutation( StdVector<UInt>& fncPermutation,
                                  const Elem* ptElem,
                                  EntityType fctEntityType,
                                  UInt entNumber){
  if(updateUnknowns_) this->CalcNumUnknowns();

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
  return maxOrder_;
}

void FeHCurlHi::GetAnisoOrder(StdVector<UInt>& order ) const {

  if( isIsotropic_) {
    order.Resize(Elem::shapes[feType_].dim);
    order.Init(isoOrder_);
  } else {
    order = anisoOrder_;  
  }
}



// ========================================================================
//  FeHCurlHi explicit element definition 
// ========================================================================

// =======================
//  QUADRILATERAL ELEMENT 
// =======================

FeHCurlHiQuad::FeHCurlHiQuad() : FeHCurlHi( Elem::ET_QUAD4 ) {
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

}
  
void FeHCurlHiQuad::CalcLocShFnc( Matrix<Double>& shape, const LocPointMapped& lp,
                              const Elem* elem, UInt comp  ) {
  EXCEPTION("Implement me");
  
}

void FeHCurlHiQuad::CalcLocCurlShFnc( Matrix<Double>& curl, const LocPointMapped& lp,
                                     const Elem* elem, UInt comp ) {
  EXCEPTION("Implement me");
}


// =======================
//  HEXAHEDRAL ELEMENT 
// =======================
FeHCurlHiHex::FeHCurlHiHex() : FeHCurlHi( Elem::ET_HEXA8) {
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
}


void FeHCurlHiHex::CalcLocShFnc( Matrix<Double>& shape, const LocPointMapped& lpm,
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
    //fac *= 0.5;
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


// define small macro to ease assignment of xpression to shape functon
#define COPYSHFNC                \
for( UInt n = 0; n < 3; ++n ) {  \
  shape[n][pos] = xpr[n];        \
  }                              \
pos++;

template<FeHCurlHiHex::DiffType DIFF_TYPE>
void FeHCurlHiHex::CalcLocShFnc2( Matrix<Double>& shape, const LocPointMapped& lpm,
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

      // === a) standard Nedelec shape functions ===
      Xpr_SGradU<3,DIFF_TYPE> xpr(eta,xi*fac);
      for( UInt k = 0; k < 3; ++k ) {
        shape[k][i] =  xpr[k];
      }
      
      // ===  b) gradient functions ===
      if( useGrad_[EDGE] && !onlyLowestOrder_ ) {
        IntLegendreP2(xiVals, order+1, fac*xi );

        for( UInt j = 0; j < order; ++j ) {
          Xpr_GradU<3,DIFF_TYPE> xpr(xiVals[j]*eta);
          COPYSHFNC
        }
      }
    }
    if(onlyLowestOrder_) return;
    
    // -------------------------
    // 2) Face shape functions
    // -------------------------
#ifdef USE_FACES
    for( UInt iFace = 0; iFace < 6; ++iFace ) {
      UInt order1 = orderFace_[iFace][0];
      UInt order2 = orderFace_[iFace][1];
      if (order1 >0 && order2 >0 ) {

        // get unique sorting of the face
        const StdVector<UInt>& unsorted = shape_.faceNodes[iFace];
        StdVector<UInt> ind;
        Face::GetSortedIndices( ind, unsorted, 4, elem->faceFlags[iFace]);

        // calculate face extension parameter which is the sum
        // of all lambdas of one face
        AutoDiff<Double,3> sum_lambda = 0.0;
        for( UInt i = 0; i < 4; ++i)
          sum_lambda += lambda[ind[i]];

        // Parameterization of first edge, connecting the
        // local nodes of the face 1->2
        AutoDiff<Double,3> xi  =  sigma[ind[1]] - sigma[ind[0]];
        // Parameterization of second edge, connecting the
        // local nodes of the face 1->4
        AutoDiff<Double,3> eta  =  sigma[ind[3]]- sigma[ind[0]];

        IntLegendreP2(xiVals,  order1+1, xi );
        IntLegendreP2(etaVals, order2+1,eta );
        
        // === a) type 1: gradient fields ===
        if( useGrad_[FACE]) {
          for( UInt i = 0; i < order1; ++i ) {
            for( UInt j = 0; j < order2; ++j ) {
              Xpr_GradU<3,DIFF_TYPE> xpr( xiVals[i] * etaVals[j] * sum_lambda);
              COPYSHFNC
            }
          }
        }
        
        // === b) type 2: face functions ===
        for( UInt i = 0; i < order1; ++i ) {
          for( UInt j = 0; j < order2; ++j ) {
            Xpr_Diff_VGradU<3,DIFF_TYPE> xpr(etaVals[j], xiVals[i] * sum_lambda);
            COPYSHFNC
          }
        }
        
        // === c) type 3 face functions ===
        for( UInt i = 0; i < order1; ++i ) {
          Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(1.0, eta, xiVals[i] * sum_lambda );
          COPYSHFNC
        }
        
        for( UInt j = 0; j < order2; ++j ) {
          Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(1.0, xi, etaVals[j] * sum_lambda );
          COPYSHFNC
        }
      } // if order
    }  // loop faces
#endif
    
    // -------------------------
    // 3) Interior shape functions
    // -------------------------
#ifdef USE_INNER
    if( orderInner_[0] > 0 && 
        orderInner_[1] > 0 &&
        orderInner_[2] > 0 ) {

      IntLegendreP2(xiVals,   orderInner_[0]+1, x );
      IntLegendreP2(etaVals,  orderInner_[1]+1, y );
      IntLegendreP2(zetaVals, orderInner_[2]+1, z );

      // === a) type 1: gradient fields ===
      if( useGrad_[INTERIOR] ) {
        for( UInt i = 0; i < orderInner_[0]; ++i ) {
          for( UInt j = 0; j < orderInner_[1]; ++j ) {
            for( UInt k = 0; k < orderInner_[2]; ++k ) {
              Xpr_GradU<3, DIFF_TYPE> xpr( xiVals[i] * etaVals[j] * zetaVals[k] );
              COPYSHFNC
            }
          }
        }
      }

      // === b) type 2 volume functions ===
      for( UInt i = 0; i < orderInner_[0]; ++i ) {
        for( UInt j = 0; j < orderInner_[1]; ++j ) {
          for( UInt k = 0; k < orderInner_[2]; ++k ) {
            Xpr_Diff_VGradU<3, DIFF_TYPE> xpr1( xiVals[i] * etaVals[j], zetaVals[k] );
            Xpr_Diff_VGradU<3, DIFF_TYPE> xpr2( xiVals[i],  etaVals[j]* zetaVals[k] );
            for( UInt n = 0; n < 3; ++n ) {
              shape[n][pos]   = xpr1[n];
              shape[n][pos+1] = xpr2[n];
            }
            pos+=2;
          }
        }
      }

      // === c) type 3 volume functions ===
      for( UInt i = 0; i < orderInner_[0]; ++i ) {
        for( UInt j = 0; j < orderInner_[1]; ++j ) {
          Xpr_Diff_SVGradU<3, DIFF_TYPE> xpr(z, AutoDiff<Double,3>(1)-z, 
                                             xiVals[i] * etaVals[j] );
          COPYSHFNC
        }
      }
      for( UInt i = 0; i < orderInner_[0]; ++i ) {
        for( UInt k = 0; k < orderInner_[2]; ++k ) {
          Xpr_Diff_SVGradU<3, DIFF_TYPE> xpr(y, AutoDiff<Double,3>(1)-y, 
                                             xiVals[i] * zetaVals[k] );
          COPYSHFNC
        }
      }
      for( UInt j = 0; j < orderInner_[1]; ++j ) {
        for( UInt k = 0; k < orderInner_[2]; ++k ) {
          Xpr_Diff_SVGradU<3, DIFF_TYPE> xpr(x, AutoDiff<Double,3>(1)-x, 
                                             etaVals[j] * zetaVals[k] );
          COPYSHFNC
        }
      }
    } // if order > 1 
#endif


}

void FeHCurlHiHex::CalcLocCurlShFnc( Matrix<Double>& curl, const LocPointMapped& lpm,
                                     const Elem* elem, UInt comp ) {

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
    //fac *= 0.5;
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

  #ifdef USE_FACES
  for( UInt f = 0; f < 6; ++f ) {
    
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


//! Return HCurl shape functions 
void FeHCurlHiHex::GetShFnc( Matrix<Double>& shape, const LocPointMapped& lpm,
                             const Elem* elem, UInt comp ) {
  // Perform local->global gradient transformation
  Matrix<Double> locShape;
  //this->CalcLocShFnc( locShape, lpm, elem, comp );
  //std::cerr << "Old local shape\n" << locShape << std::endl;
  CalcLocShFnc2<ID>(locShape, lpm, elem, comp);
  //std::cerr << "New local shape\n" << locShape << std::endl;
  shape =  Transpose(lpm.jacInv) * locShape;
}

//! Return global curl of shape functions
void FeHCurlHiHex::GetCurlShFnc( Matrix<Double>& curl, const LocPointMapped& lpm,
                                 const Elem* elem, UInt comp ) {
  // Perform local->global curl transformation
  Matrix<Double> locCurl;    
  
  //this->CalcLocCurlShFnc( locCurl, lpm, elem, comp );
  //std::cerr << "Old local curl\n" << locCurl << std::endl;
  CalcLocShFnc2<CURL>( locCurl, lpm, elem, comp );
  //std::cerr << "New local curl\n" << locCurl << std::endl;
  
  curl = lpm.jac * locCurl;
  curl *= ( 1.0 / std::abs(lpm.jacDet) );
  
}

}// end of namespace
