#include "HCurlElemsHi.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/AutoDiff.hh"
#include "FeBasis/Polynomials.hh"
#include "Domain/ElemMapping/EdgeFace.hh"

namespace CoupledField {

#define USE_FACES 1
#define USE_INNER 1

// declare class specific logging stream
DEFINE_LOG(feHCurlHi, "feHCurlHi")


#ifndef PYRA_DOF_PRINT
#define PYRA_DOF_PRINT(elemName, onlyLO, what, idx, val) \
    //std::cout << "[DBG] " << elemName << " onlyLO=" << onlyLO \
    //          << "  " << what << "[" << idx << "] = " << val << '\n';
#endif


// ===============================================================
//  VECTOR IDENTITIES NEEDED FOR THE CALCULATION OF CURL MATRICES
// ===============================================================
//! This method makes use of the following vector identity
//! (assuming b is a scalar and F a vector function):
//!
//!     curl(b*F) = b x curl(F) + grad(b) x F
//!
//! If we further assume, that F = grad(a), we arrive at
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


// -------------------------------------
//  Expr = U*grad(V) - W*grad(X)
// -------------------------------------

template<int D, FeHCurlHi::DiffType DIFF>
class Xpr_Diff_UGradV_min_WGradX {
};

template<int D>
class Xpr_Diff_UGradV_min_WGradX<D, FeHCurlHi::ID> {
  public:
  Xpr_Diff_UGradV_min_WGradX(const AutoDiff<Double,D>& u, 
                             const AutoDiff<Double,D>& v,
                             const AutoDiff<Double,D>& w,
                             const AutoDiff<Double,D>& x) {
    for( UInt i = 0; i < D; ++i ) {
        val[i] = u.Val() * v.DVal(i) - w.Val() * x.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return val[i];}
  private:
    Double val[D];
};

template<int D>
class Xpr_Diff_UGradV_min_WGradX<D, FeHCurlHi::CURL> {
  public:
  Xpr_Diff_UGradV_min_WGradX(const AutoDiff<Double,D>& u, 
                             const AutoDiff<Double,D>& v,
                             const AutoDiff<Double,D>& x,
                             const AutoDiff<Double,D>& w) {
    AutoDiff<Double,D> c = Cross(u,v) - Cross(x,w);
      for( UInt i = 0; i < D; ++i ) {
        curl[i] = c.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return curl[i];}
  private:
    Double curl[D];
};

// -------------------------------------
//  Expr = V*W*Grad(U)
// -------------------------------------
template<int D, FeHCurlHi::DiffType DIFF>
class Xpr_Diff_VWGradU {
};

template<int D>
class Xpr_Diff_VWGradU<D, FeHCurlHi::ID> {
  public:
  Xpr_Diff_VWGradU(const AutoDiff<Double,D>& u,
                             const AutoDiff<Double,D>& v,
                             const AutoDiff<Double,D>& w) {
    for( UInt i = 0; i < D; ++i ) {
        val[i] = u.DVal(i) * v.Val() * w.Val();
      }
    }
    Double operator[] (UInt i ) {return val[i];}
  private:
    Double val[D];
};

template<int D>
class Xpr_Diff_VWGradU<D, FeHCurlHi::CURL> {
  public:
  Xpr_Diff_VWGradU(const AutoDiff<Double,D>& u,
                             const AutoDiff<Double,D>& v,
                             const AutoDiff<Double,D>& w) {


  AutoDiff<Double,D> vGradw = Cross(v,w); //from Xpr_SGradU(s,u)
  AutoDiff<Double,D> wGradv = Cross(w,v); //from Xpr_SGradU(s,u)
  AutoDiff<Double,D> vGradwPwGradv= vGradw + wGradv;


    AutoDiff<Double,D> c = Cross(vGradwPwGradv, u);
      for( UInt i = 0; i < D; ++i ) {
        curl[i] = c.DVal(i);
      }
    }
    Double operator[] (UInt i ) {return curl[i];}
  private:
    Double curl[D];
};

// -------------------------------------
//  Expr = U*V
// -------------------------------------
template<int D, FeHCurlHi::DiffType DIFF>
class Xpr_Diff_UV {
};

template<int D>
class Xpr_Diff_UV<D, FeHCurlHi::ID> {
  public:
  Xpr_Diff_UV(const AutoDiff<Double,D>& u,
                             const AutoDiff<Double,D>& v) {
    for( UInt i = 0; i < D; ++i ) {
        val[i] = u.Val() * v.Val();
      }
    }
    Double operator[] (UInt i ) {return val[i];}
  private:
    Double val[D];
};

template<int D>
class Xpr_Diff_UV<D, FeHCurlHi::CURL> {
  public:
  Xpr_Diff_UV(const AutoDiff<Double,D>& u,
                             const AutoDiff<Double,D>& v) {


  AutoDiff<Double,D> c = Cross(u,v);
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
  
  // initialize useage of gradients
  useEdgeGrad_.Resize( shape_.numEdges );
  useFaceGrad_.Resize( shape_.numFaces );
  useEdgeGrad_.Init( false );
  useFaceGrad_.Init( false );
  useInteriorGrad_ = false;
  
}

FeHCurlHi::FeHCurlHi(const FeHCurlHi & other)
          : FeHCurl(other), FeHi(other){
  this->onlyLowestOrder_ = other.onlyLowestOrder_;
  this->useEdgeGrad_ = other.useEdgeGrad_;
  this->useFaceGrad_ = other.useFaceGrad_;
  this->useInteriorGrad_ = other.useInteriorGrad_;
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
  
void FeHCurlHi::SetOnlyLowestOrder( bool flag ) {
  onlyLowestOrder_ = flag;
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

bool FeHCurlHi::operator==( const FeHCurlHi& comp) const {
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

void FeHCurlHi::SetUseGradients(bool useGrad) {
  useInteriorGrad_ = useGrad;
  useEdgeGrad_.Init( useGrad );
  useFaceGrad_.Init( useGrad );

  updateUnknowns_ = true;
}

void FeHCurlHi::SetEdgeGradient(UInt edgeNum, bool useGrad) {
  assert( edgeNum <= elemShape_.numEdges);
  useEdgeGrad_[edgeNum] = useGrad;
  updateUnknowns_ = true;
}

void FeHCurlHi::SetFaceGradient(UInt faceNum, bool useGrad) {
  assert( faceNum <= elemShape_.numFaces);
  useFaceGrad_[faceNum] = useGrad;
  updateUnknowns_ = true;
}


// ========================================================================
//  FeHCurlHi explicit element definition 
// ========================================================================

// =======================
//  TRIANGULAR ELEMENT 
// =======================

FeHCurlHiTria::FeHCurlHiTria() : FeHCurlHi( Elem::ET_TRIA3 ) {
}


FeHCurlHiTria::~FeHCurlHiTria() {
  
}

void FeHCurlHiTria::CalcNumUnknowns() {
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
    if( useEdgeGrad_[i])
      unknowns += orderEdge_[i];
    edgeFncs[i] = unknowns;
    LOG_DBG(feHCurlHi) <<   "edge " << i+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
  }

  // Faces
  StdVector<UInt>& faceFncs = entityFncs_[FACE];
  faceFncs.Resize(shape_.numFaces);
  faceFncs.Init(0);
  
  
  // Note: Curtently we do not support higher order shape functions for
  // the triangular element
  
//  for( UInt i = 0; i < shape_.numFaces; ++i ) {
//    if( orderFace_[i][0] > 0 &&
//        orderFace_[i][1] > 0 ) {
//      unknowns = orderFace_[i][0] * orderFace_[i][1] // face functions of 1st kind
//                + orderFace_[i][0] + orderFace_[i][1];
//      if( useFaceGrad_[i])
//        unknowns +=  orderFace_[i][0] * orderFace_[i][1];
//      faceFncs[i] = unknowns;
//      LOG_DBG(feHCurlHi) << "face " << i+1 << " has " << unknowns << "unknowns";
//      actNumFncs_ += unknowns;
//    }
//  }

  // Interior
  StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
  innerFncs.Resize(1);
  innerFncs.Init(0);
  LOG_DBG(feHCurlHi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;

}

void FeHCurlHiTria::GetShFnc( Matrix<Double>& shape,
                             const LocPointMapped& lpm,
                             const Elem* elem, UInt comp ) {

  Matrix<Double> locShape;
  CalcLocShFnc2<ID>(locShape, lpm, elem, comp);

  // Perform local->global gradient transformation
    if (lpm.isSurface) {
      // result taking from PhD of S. Zaglmayr, p.44, based on (4.20)
      shape = locShape *(1.0/ lpm.lpmVol->jacDet);
      //shape = locShape;
    } else {
  shape =  Transpose(lpm.jacInv) * locShape;
    }
}

void FeHCurlHiTria::GetCurlShFnc( Matrix<Double>& curl,
                                 const LocPointMapped& lpm,
                                 const Elem* elem, UInt comp ) {

  Matrix<Double> locCurl;
  CalcLocShFnc2<CURL>( locCurl, lpm, elem, comp );

  // Perform local->global curl transformation
  curl = lpm.jac * locCurl;
  curl *= ( 1.0 / std::fabs(lpm.jacDet) );
}

template<FeHCurlHi::DiffType DIFF_TYPE>
void FeHCurlHiTria::CalcLocShFnc2( Matrix<Double>& shape,
                                  const LocPointMapped& lpm,
                                  const Elem* elem, UInt comp ) {

  if (updateUnknowns_) CalcNumUnknowns();
  AutoDiff<Double, 2> x (lpm.lp.coord[0],0);
  AutoDiff<Double, 2> y (lpm.lp.coord[1],1);
  AutoDiff<Double, 2> lambda[3] = {1.0 - x - y, x, y};
  shape.Resize(2,actNumFncs_);
  shape.Init();

  StdVector<AutoDiff<Double, 2> > Vals;
  // ------------------------
  // 1) Edge shape functions
  // ------------------------
  for( UInt i = 0; i < 3; ++i ) {
    UInt order = orderEdge_[i];
    if(order > 1) EXCEPTION("HCurl TRIA shape functions only defined for order 0 and 1!");
    UInt index1 = shape_.edgeVertices[i][0]-1;
    UInt index2 = shape_.edgeVertices[i][1]-1;
    if ( elem->extended->edges[i] < 0 ) std::swap(index1, index2);  // fmax > f1 > f2

    // === a) standard Nedelec shape functions ===
    //Xpr_Diff_SVGradU<2,DIFF_TYPE> xpr( lambda[index1], lambda[index2], 1.0);
    for( UInt k = 0; k < 2; ++k ){
      shape[k][i] = lambda[index1].DVal(k) * lambda[index2].Val() + lambda[index2].DVal(k) * lambda[index1].Val();
    }

    // === b) gradient functions
    if( useEdgeGrad_[i] ) {
      WARN("Gradient fields for TRIA elements not thoroughly tested!!!");
      if (!onlyLowestOrder_){
        ScaledIntLegendreP2(Vals, order + 1, lambda[index2]+lambda[index1], lambda[index2]-lambda[index1]);
        for (UInt j = 0; j < order; ++j) {
          Xpr_GradU<2,DIFF_TYPE> xpr(Vals[j]);
          for( UInt k = 0; k < 2; ++k ) {
            shape[k][j] = xpr[k];
          }
        }
      } //if: edgeGrad
    }
    if(onlyLowestOrder_) return;
  }
  return;
}

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
    if( useEdgeGrad_[i])
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
      if( useFaceGrad_[i])
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
  
void FeHCurlHiQuad::GetShFnc( Matrix<Double>& shape, 
                             const LocPointMapped& lpm,
                             const Elem* elem, UInt comp ) {
  
  Matrix<Double> locShape;
  CalcLocShFnc2<ID>(locShape, lpm, elem, comp);

  // Perform local->global gradient transformation
    if (lpm.isSurface) {
      // result taking from PhD of S. Zaglmayr, p.44, based on (4.20)
      shape = locShape *(1.0/ lpm.lpmVol->jacDet); //This simple transformation may not work
      //shape = locShape;
    } else {
      shape =  Transpose(lpm.jacInv) * locShape;
    }
}

void FeHCurlHiQuad::GetCurlShFnc( Matrix<Double>& curl, 
                                 const LocPointMapped& lpm,
                                 const Elem* elem, UInt comp ) {

  Matrix<Double> locCurl;    
  CalcLocShFnc2<CURL>( locCurl, lpm, elem, comp );
  
  // Perform local->global curl transformation
  curl = lpm.jac * locCurl;
  curl *= ( 1.0 / std::fabs(lpm.jacDet) );
}

template<FeHCurlHi::DiffType DIFF_TYPE>
void FeHCurlHiQuad::CalcLocShFnc2( Matrix<Double>& shape, 
                                  const LocPointMapped& lpm,
                                  const Elem* elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();

  AutoDiff<Double, 2> x (lpm.lp.coord[0],0);
  AutoDiff<Double, 2> y (lpm.lp.coord[1],1);

  AutoDiff<Double, 2> lambda[4] = {0.25*(1.0-x)*(1.0-y),
                                   0.25*(1.0+x)*(1.0-y),
                                   0.25*(1.0+x)*(1.0+y),
                                   0.25*(1.0-x)*(1.0+y)};

  AutoDiff<Double, 2> sigma[4]  = {0.5*((1.0-x)+(1.0-y)),
                                   0.5*((1.0+x)+(1.0-y)),
                                   0.5*((1.0+x)+(1.0+y)),
                                   0.5*((1.0-x)+(1.0+y))};
  shape.Resize(2,actNumFncs_);
  shape.Init();

  StdVector<AutoDiff<Double, 2> > xiVals, etaVals, zetaVals;

  // ------------------------
  // 1) Edge shape functions
  // ------------------------
  for( UInt i = 0; i < 4; ++i ) {

    Double fac = elem->extended->edges[i] < 0 ? -1.0 : 1.0;
    fac *= 0.5;
    UInt index1 = shape_.edgeVertices[i][0]-1;
    UInt index2 = shape_.edgeVertices[i][1]-1;

    // xi: parameterization of edge [-1;+1]
    // eta: parameterization of extension into element [-1;+1]
    AutoDiff<Double, 2> xi  =  sigma[index2] -  sigma[index1]; 
    AutoDiff<Double, 2> eta = lambda[index1] + lambda[index2];

    // a) standard Nedelec shape functions
    for( UInt k = 0; k < 2; ++k ) {
      shape[k][i] = xi.DVal(k)*eta.Val()*fac;
    }
    // b) gradient functions
    if( useEdgeGrad_[i] ) {
      EXCEPTION("Not implemented");
    }
  }
  return;

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
    if( useEdgeGrad_[i])
      unknowns += orderEdge_[i];
    edgeFncs[i] = unknowns;
    LOG_DBG(feHCurlHi) <<   "edge " << i+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
    // std::cout << "[HEX] edge " << i
    //   << "  order=" << orderEdge_[i]
    //   << "  grad="  << (useEdgeGrad_[i] ? "yes" : "no")
    //   << "  -> "    << i << " DOFs\n";

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
      if( useFaceGrad_[i])
        unknowns +=  orderFace_[i][0] * orderFace_[i][1];
      faceFncs[i] = unknowns;
      LOG_DBG(feHCurlHi) << "face " << i+1 << " has " << unknowns << "unknowns";
      actNumFncs_ += unknowns;
    }
    // std::cout << "[HEX] face " << i
    // << "  orders=(" << orderFace_[i][0] << "," << orderFace_[i][1] << ")"
    // << "  grad="    << (useFaceGrad_[i] ? "yes" : "no")
    // << "  -> "      << unknowns << " DOFs\n";

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
    if( useInteriorGrad_ ) { 
      unknowns += orderInner_[0] * orderInner_[1] * orderInner_[2];
    }
    actNumFncs_ += unknowns;
    innerFncs[0] = unknowns;
    LOG_DBG(feHCurlHi) << "interior has " << unknowns << "unknowns";
  }
#endif

  LOG_DBG(feHCurlHi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;

  for (UInt f = 0; f < shape_.numFaces; ++f)
    PYRA_DOF_PRINT("HEX ", onlyLowestOrder_, "face", f, entityFncs_[FACE][f]);
for (UInt e = 0; e < shape_.numEdges; ++e)
    PYRA_DOF_PRINT("HEX ", onlyLowestOrder_, "edge", e, entityFncs_[EDGE][e]);
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
    Double fac = elem->extended->edges[i] < 0 ? -1.0 : 1.0;
    fac *= 0.5;
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
    if( useEdgeGrad_[i] ) {
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
//      if( useGrad_[FACE])
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
  if( useInteriorGrad_ ) {
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
    UInt pos = 0;
    shape.Resize(3,actNumFncs_);
    shape.Init();
    
    StdVector<AutoDiff<Double, 3> > xiVals, etaVals, zetaVals;
    // ------------------------
    // 1) Edge shape functions
    // ------------------------
    for( UInt i = 0; i < 12; ++i ) {

      UInt order = orderEdge_[i];
      UInt index1 = shape_.edgeVertices[i][0]-1;
      UInt index2 = shape_.edgeVertices[i][1]-1;
      if ( elem->extended->edges[i] < 0 ) {
        std::swap(index1, index2);  // fmax > f1 > f2
      }

      // xi: parameterization of edge [-1;+1]
      // eta: parameterization of extension into element [-1;+1]
      AutoDiff<Double, 3> xi  =  sigma[index2] -  sigma[index1]; 
      AutoDiff<Double, 3> eta = lambda[index1] + lambda[index2];

      // === a) standard Nedelec shape functions ===
      Xpr_SGradU<3,DIFF_TYPE> xpr(eta*0.5,xi);
      COPYSHFNC
      
      // ===  b) gradient functions ===
      if( useEdgeGrad_[i]) {
        if (onlyLowestOrder_) {
          for( UInt j = 0; j < order; ++j ) {
            pos++;
          }
        } else {
          IntLegendreP2(xiVals, order+1, xi );

          for( UInt j = 0; j < order; ++j ) {
            Xpr_GradU<3,DIFF_TYPE> xpr(xiVals[j]*eta);
            COPYSHFNC
          }
        }// if: lowestOrder
      } //if: edgeGrad
    } //loop: edges
    
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
        Face::GetSortedIndices( ind, unsorted, 4, elem->extended->faceFlags[iFace]);

        // calculate face extension parameter which is the sum
        // of all lambdas of one face
        AutoDiff<Double,3> sum_lambda = 0.0;
        for( UInt i = 0; i < 4; ++i)
        {
          sum_lambda += lambda[ind[i]];
        }

        // Parameterization of first edge, connecting the
        // local nodes of the face 1->2
        AutoDiff<Double,3> xi  =  sigma[ind[1]] - sigma[ind[0]];
        // Parameterization of second edge, connecting the
        // local nodes of the face 1->4
        AutoDiff<Double,3> eta  =  sigma[ind[3]]- sigma[ind[0]];

        IntLegendreP2(xiVals,  order1+1, xi );
        IntLegendreP2(etaVals, order2+1,eta );
        
        // === a) type 1: gradient fields ===
        if( useFaceGrad_[iFace]) {
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
        
        // === c) type 3: face functions ===
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
      if( useInteriorGrad_ ) {
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
    Double fac = elem->extended->edges[i] < 0 ? -1.0 : 1.0;
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
    if( useEdgeGrad_[i] ) {
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
      if( elem->extended->faceFlags[f].test(2) == false) {
        std::swap(order1, order2);
        std::swap(xi, eta);
      }
      // Model a XOR B = ((a || b) && !(a && b) 
//      if ( (elem->faceFlags[f].test(0) || (elem->faceFlags[f].test(1))) && 
//          !(elem->faceFlags[f].test(0) && (elem->faceFlags[f].test(1))) ) {
//        xi *= -1.0;
//      }
      if ( (elem->extended->faceFlags[f].test(0) ) &&
          !(elem->extended->faceFlags[f].test(0) && (elem->extended->faceFlags[f].test(1))) ) {
        xi *= -1.0;
      }
      if ( (elem->extended->faceFlags[f].test(1)) &&
          !(elem->extended->faceFlags[f].test(0) && (elem->extended->faceFlags[f].test(1))) ) {
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
      if( useFaceGrad_[f])
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
  if( useInteriorGrad_) {
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
  curl = lpm.jac * locCurl;
  curl *= ( 1.0 / std::fabs(lpm.jacDet) );

}


// ===========================
//  Wedge / Prismatic element
// ===========================

FeHCurlHiWedge::FeHCurlHiWedge() : FeHCurlHi( Elem::ET_WEDGE6) {

}

FeHCurlHiWedge::~FeHCurlHiWedge() {

}

void FeHCurlHiWedge::GetShFnc( Matrix<Double>& shape, 
                             const LocPointMapped& lpm,
                             const Elem* elem, UInt comp ) {
  
  Matrix<Double> locShape;
  CalcLocShFnc2<ID>(locShape, lpm, elem, comp);

  // Perform local->global gradient transformation
  shape =  Transpose(lpm.jacInv) * locShape;
}

void FeHCurlHiWedge::GetCurlShFnc( Matrix<Double>& curl, 
                                 const LocPointMapped& lpm,
                                 const Elem* elem, UInt comp ) {

  Matrix<Double> locCurl;    
  CalcLocShFnc2<CURL>( locCurl, lpm, elem, comp );
  
  // Perform local->global curl transformation
  curl = lpm.jac * locCurl;
  curl *= ( 1.0 / std::fabs(lpm.jacDet) );
}

template<FeHCurlHi::DiffType DIFF_TYPE>
void FeHCurlHiWedge::CalcLocShFnc2( Matrix<Double>& shape, 
                                  const LocPointMapped& lpm,
                                  const Elem* elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();

  AutoDiff<Double, 3> x (lpm.lp.coord[0],0);
  AutoDiff<Double, 3> y (lpm.lp.coord[1],1);
  AutoDiff<Double, 3> z (lpm.lp.coord[2],2);

  AutoDiff<Double, 3> lambda[6] = { 1.0 - x - y, x,  y, 
                                    1.0 - x - y, x,  y };
  AutoDiff<Double, 3> mu[6]     = { 0.5 * (1.0-z), 0.5 * (1.0-z), 0.5 * (1.0-z),
                                    0.5 * (1.0+z), 0.5 * (1.0+z), 0.5 * (1.0+z) };  
  UInt pos = 0;
  shape.Resize(3, actNumFncs_);
  shape.Init();
  
  // ------------------------
  // 1) Edge shape functions
  // ------------------------
  StdVector<AutoDiff<Double, 3> > Vals;

  // a) horizontal edges (= triangular)
  for( UInt i = 0; i < 6; ++i ) {

    UInt order = orderEdge_[i];
    UInt index1 = shape_.edgeVertices[i][0]-1;
    UInt index2 = shape_.edgeVertices[i][1]-1;
    if ( elem->extended->edges[i] < 0 ) {
      std::swap(index1, index2);  // fmax > f1 > f2
    }
    
    // === a) standard Nedelec shape functions ===
    Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr( lambda[index1], lambda[index2], mu[index1] );
    COPYSHFNC

    // ===  b) gradient functions ===
    if( useEdgeGrad_[i]) {
      if (onlyLowestOrder_) {
        for( UInt j = 0; j < order; ++j ) {
          pos++;
        }
      } else {
        ScaledIntLegendreP2(Vals, order + 1, lambda[index2]+lambda[index1], lambda[index2]-lambda[index1]);
        for (UInt j = 0; j < order; ++j) {
            Xpr_GradU<3,DIFF_TYPE> xpr(Vals[j] * mu[index1]);
            COPYSHFNC
          }
        }// if: lowestOrder
      } //if: edgeGrad
    } //loop: horizontal edges
     
  // b) vertical edges (= quad)
  for( UInt i = 6; i < 9; ++i ) {

    UInt order = orderEdge_[i];
    UInt index1 = shape_.edgeVertices[i][0]-1;
    UInt index2 = shape_.edgeVertices[i][1]-1;
    if ( elem->extended->edges[i] < 0 ) {
      std::swap(index1, index2);  // fmax > f1 > f2
    }

    // === a) standard Nedelec shape functions ===
    Xpr_SGradU<3,DIFF_TYPE> xpr(  lambda[index1], mu[index2]  );
    COPYSHFNC

    // ===  b) gradient functions ===
    if( useEdgeGrad_[i]) {
      if (onlyLowestOrder_) {
        for( UInt j = 0; j < order; ++j ) {
          pos++;
        }
      } else {
    IntLegendreP2(Vals, order + 1, 2.0 * mu[index1] - 1.0);
    for (UInt j = 0; j < order; ++j) {
       Xpr_GradU<3,DIFF_TYPE> xpr(Vals[j] * (lambda[index1] + lambda[index2]));
       COPYSHFNC
        }
       }// if: lowestOrder
      } //if: edgeGrad
  } //loop: vertical edges

  if(onlyLowestOrder_) return;


  // -------------------------
  // 2) Face shape functions
  // -------------------------

#ifdef USE_FACES


  // === triangular faces
  StdVector<AutoDiff<Double, 3> > ui;
  StdVector<AutoDiff<Double, 3> > vj;
  StdVector<AutoDiff<Double, 3> > temp_vj;

  for( UInt iFace = 0; iFace < 2; ++iFace ) {
   //only valid for isotropic polynomial order!!

   UInt order = orderFace_[iFace][0];
   if (order > 1) {
        // get unique sorting of the face
        const StdVector<UInt>& unsorted = shape_.faceNodes[iFace];
        StdVector<UInt> ind;
        Face::GetSortedIndices( ind, unsorted, 3, elem->extended->faceFlags[iFace]);

        //definition of ui and vj according to PHD thesis of Sabine Zaglmayr p.92
        ScaledIntLegendreP2(ui, order+1, lambda[ind[0]]+lambda[ind[1]], lambda[ind[1]]-lambda[ind[0]]);
        Legendre(temp_vj, order+1, 2.0 * lambda[ind[2]] - 1.0);
        vj.Init();
        vj.Resize(order);
         for (UInt j=0; j<order; ++j){
          vj[j] = lambda[ind[2]] * temp_vj[j];
         }

         // === a) type 1: gradient fields ===
         if( useFaceGrad_[iFace]) {
             for( UInt i = 0; i <= order - 2; ++i ) {
               for( UInt j = 0; j <= order - 2 - i; ++j ) {
               Xpr_GradU<3,DIFF_TYPE> xpr(ui[i] * vj[j] * mu[ind[0]]);
               COPYSHFNC
             }
           }
         }

         // === b) type 2: face functions ===
         for( UInt i = 0; i <= order - 2; ++i ) {
           for( UInt j = 0; j <= order - 2 - i; ++j ) {
             Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(ui[i], vj[j], mu[ind[0]]);
             COPYSHFNC
           }
         }

         // === c) type 3: face functions ===
         for( UInt j = 0; j <= order - 2; ++j ) {
           Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(lambda[ind[0]], lambda[ind[1]], vj[j] * mu[ind[0]]);
           COPYSHFNC
         }
   } // if: order > 1
  } // loop over all triangular faces



  // === quadrilateral faces
  StdVector<AutoDiff<Double, 3> > wj;
  for( UInt iFace = 2; iFace < 5; ++iFace ) {
    //only valid for isotropic polynomial order!!
  // this means orderFace_[iFace][0] = orderFace_[iFace][1] = orderFace_[iFace][2]

    UInt order = orderFace_[iFace][0];
    if (order > 0) {
        // get unique sorting of the face
        const StdVector<UInt>& unsorted = shape_.faceNodes[iFace];
        StdVector<UInt> ind;
        Face::GetSortedIndices( ind, unsorted, 4, elem->extended->faceFlags[iFace]);

        // extension on vertical edge
        IntLegendreP2( wj, orderFace_[iFace][1] + 1, 2.0 * mu[ind[0]] - 1.0 );

        // we have to determine 2 things:
        // 1) Determine, if first edge in sorted array is horizontal or vertical one
        if( shape_.nodeCoords[ind[0]][2] == shape_.nodeCoords[ind[1]][2] ) {
          // edge [0] -> [1] is the horizontal one
          // edge [0] -> [3] is the vertical one
          // this means, this case is the horizontal case

            //definition of ui and wj according to PHD thesis of Sabine Zaglmayr p.92
            ScaledIntLegendreP2( ui, orderFace_[iFace][0]+1,
                                 lambda[ind[1]]+lambda[ind[0]],
                                 lambda[ind[1]]-lambda[ind[0]] );

            // === a) type 1: gradient fields ===
            if( useFaceGrad_[iFace]) {
                for( UInt i = 0; i <= order - 1; ++i ) {
                  for( UInt j = 0; j <= order - 1; ++j ) {
                  Xpr_GradU<3,DIFF_TYPE> xpr(ui[i] * wj[j]);
                  COPYSHFNC
                }
              }
            }
            // === b) type 2: face functions ===
            for( UInt i = 0; i <= order - 1; ++i ) {
              for( UInt j = 0; j <= order - 1; ++j ) {
                Xpr_Diff_VGradU<3,DIFF_TYPE> xpr(ui[i], wj[j]);
                COPYSHFNC
              }
            }
            // === c) type 3: face functions ===
            for( UInt j = 0; j <= order - 1; ++j ) {
              Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(lambda[ind[0]], lambda[ind[1]], wj[j]);
              COPYSHFNC
            }
            for( UInt j = 0; j <= order - 1; ++j ) {
              Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(lambda[ind[0]], lambda[ind[3]], wj[j]);
              COPYSHFNC
            }

        } else{
          // this now is the vertical case: edge[0] -> edge[3]

            ScaledIntLegendreP2( ui, orderFace_[iFace][1]+1,
                                 lambda[ind[3]]+lambda[ind[0]],
                                 lambda[ind[3]]-lambda[ind[0]] );

            // === a) type 1: gradient fields ===
            if( useFaceGrad_[iFace]) {
                for( UInt i = 0; i <= order - 1; ++i ) {
                  for( UInt j = 0; j <= order - 1; ++j ) {
                  Xpr_GradU<3,DIFF_TYPE> xpr(ui[i] * wj[j]);
                  COPYSHFNC
                }
              }
            }
            // === b) type 2: face functions ===
            for( UInt i = 0; i <= order - 1; ++i ) {
              for( UInt j = 0; j <= order - 1; ++j ) {
                Xpr_Diff_VGradU<3,DIFF_TYPE> WGradU(ui[i], wj[j]);
                for( UInt n = 0; n < 3; ++n ) {
                  shape[n][pos]   = WGradU[n] * (-1.0);
                }
                pos++;
              }
            }

            // === c) type 3: face functions ===
            for( UInt j = 0; j <= order - 1; ++j ) {
              Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(lambda[ind[0]], lambda[ind[3]], wj[j]);
              COPYSHFNC
            }
            for( UInt j = 0; j <= order - 1; ++j ) {
              Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(lambda[ind[0]], lambda[ind[1]], wj[j]);
              COPYSHFNC
            }
        } // horizontal or vertical edge selection
    } // if: order > 0
  } // loop over all quadrilateral faces
#endif



  // -------------------------
  // 3) Interior shape functions
  // -------------------------
#ifdef USE_INNER
  StdVector<AutoDiff<Double, 3> > wk;
  // only isotropic polynomial order
  if( orderInner_[0] > 1) {

    //definition of ui, vj and wk according to PHD thesis of Sabine Zaglmayr p.93
    ScaledIntLegendreP2(ui, orderInner_[0]+1, lambda[0]+lambda[1], lambda[1]-lambda[0]);
    Legendre(temp_vj, orderInner_[0]+1, 2.0 * lambda[2] - 1.0);
    IntLegendreP2( wk, orderInner_[0]+1, 2.0 * mu[0] - 1.0);
    vj.Resize(orderInner_[0]);
    for (UInt j=0; j<orderInner_[0]; ++j){
     vj[j] = lambda[2] * temp_vj[j];
    }


    // === a) type 1: gradient fields ===
    if( useInteriorGrad_ ) {
    for( UInt i = 0; i <= orderInner_[0] - 2; ++i ) {
      for( UInt j = 0; j <= orderInner_[0] - 2 - i; ++j ) {
        for ( UInt k = 0; k <= orderInner_[0] - 1; ++k){
          Xpr_GradU<3, DIFF_TYPE> xpr( ui[i] * vj[j] * wk[k] );
          COPYSHFNC
        }
      }
    }
    }
    // === b) type 2 volume functions ===
    for( UInt i = 0; i <= orderInner_[0] - 2; ++i ) {
      for( UInt j = 0; j <= orderInner_[0] - 2 - i; ++j ) {
        for ( UInt k = 0; k <= orderInner_[0] - 1; ++k){
          Xpr_Diff_VGradU<3, DIFF_TYPE> xpr1( ui[i], vj[j] * wk[k] );
          Xpr_Diff_VGradU<3, DIFF_TYPE> xpr2( vj[j], ui[i] * wk[k] );
          for( UInt n = 0; n < 3; ++n ) {
            shape[n][pos]   = xpr1[n];
            shape[n][pos+1] = xpr2[n];
          }
          pos+=2;
        }
      }
    }
    // === b) type 3 volume functions ===
    for( UInt i = 0; i <= orderInner_[0] - 2; ++i ) {
      for( UInt j = 0; j <= orderInner_[0] - 2 - i; ++j ) {
          Xpr_Diff_SVGradU<3, DIFF_TYPE> xpr (z, AutoDiff<Double,3>(1)-z, ui[i] * vj[j]);
          COPYSHFNC
      }
    }
      for( UInt j = 0; j <= orderInner_[0] - 2; ++j ) {
        for ( UInt k = 0; k <= orderInner_[0] - 1; ++k){
          Xpr_Diff_SVGradU<3, DIFF_TYPE> xpr( lambda[0], lambda[1], vj[j] * wk[k] );
          COPYSHFNC
        }
      }
  } // if order > 1
#endif

return;
}

void FeHCurlHiWedge::CalcNumUnknowns() {
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

    if( useEdgeGrad_[i]) {
      if (onlyLowestOrder_) {return;}
      else{
        unknowns += orderEdge_[i];
      }
    }
    edgeFncs[i] = unknowns;
    LOG_DBG(feHCurlHi) <<   "edge " << i+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
  }


  // Faces
  StdVector<UInt>& faceFncs = entityFncs_[FACE];
  faceFncs.Resize(shape_.numFaces);
  faceFncs.Init(0);

#ifdef USE_FACES

  // === triangular faces
  for( UInt iFace = 0; iFace < 2; ++iFace ) {
    unknowns = 0;
    UInt order = orderFace_[iFace][0];
    if (order >1) {

         // === a) type 1: gradient fields ===
         if( useFaceGrad_[iFace]) {
             for( UInt i = 0; i <= order - 2; ++i ) {
               for( UInt j = 0; j <= order - 2 - i; ++j ) {
                 unknowns+=1;
             }
           }
         }

         // === b) type 2: face functions ===
         for( UInt i = 0; i <= order - 2; ++i ) {
           for( UInt j = 0; j <= order - 2 - i; ++j ) {
             unknowns+=1;
           }
         }

         // === c) type 3: face functions ===
         for( UInt j = 0; j <= order - 2; ++j ) {
           unknowns+=1;
         }
    } //if order > 0
    faceFncs[iFace] = unknowns;
    LOG_DBG(feHCurlHi) << "face " << iFace+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
  } //loop over all triangular faces



    // === quadrilateral faces
  for( UInt iFace = 2; iFace < 5; ++iFace ) {
    unknowns = 0;
    UInt order = orderFace_[iFace][0];
    if (order > 0) {
      if( useFaceGrad_[iFace]) {
          for( UInt i = 0; i <= order - 1; ++i ) {
            for( UInt j = 0; j <= order - 1 ; ++j ) {
            unknowns+=1;
          }
        }
      }


    // === b) type 2: face functions ===
    for( UInt i = 0; i <= order - 1; ++i ) {
      for( UInt j = 0; j <= order - 1; ++j ) {
        unknowns+=1;
      }
    }

    // === c) type 3: face functions ===
    for( UInt j = 0; j <= order - 1; ++j ) {
      unknowns+=2;
    }
    } //if order > 0
    faceFncs[iFace] = unknowns;
    LOG_DBG(feHCurlHi) << "face " << iFace+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
  } //loop over all quadrilateral faces
#endif

  // Interior
  StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
  innerFncs.Resize(1);
  innerFncs.Init(0);

  #ifdef USE_INNER
  unknowns=0;
  // only isotropic polynomial order
  if( orderInner_[0] > 1) {
    // === a) type 1: gradient fields ===
    if( useInteriorGrad_ ) {
    for( UInt i = 0; i <= orderInner_[0] - 2; ++i ) {
      for( UInt j = 0; j <= orderInner_[0] - 2 - i; ++j ) {
        for ( UInt k = 0; k <= orderInner_[0] - 1; ++k){
          unknowns+=1;
        }
      }
    }
    }

    // === b) type 2 volume functions ===
    for( UInt i = 0; i <= orderInner_[0] - 2; ++i ) {
      for( UInt j = 0; j <= orderInner_[0] - 2 - i; ++j ) {
        for ( UInt k = 0; k <= orderInner_[0] - 1; ++k){
          unknowns+=2;
        }
      }
    }

    // === c) type 3 volume functions ===
    for( UInt i = 0; i <= orderInner_[0] - 2; ++i ) {
      for( UInt j = 0; j <= orderInner_[0] - 2 - i; ++j ) {
          unknowns+=1;
      }
    }

      for( UInt j = 0; j <= orderInner_[0] - 2; ++j ) {
        for ( UInt k = 0; k <= orderInner_[0] - 1; ++k){
          unknowns+=1;
        }
      }

  } //if order > 1
  actNumFncs_ += unknowns;
  innerFncs[0] = unknowns;
  LOG_DBG(feHCurlHi) << "interior has " << unknowns << "unknowns";
  #endif

  LOG_DBG(feHCurlHi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;
}







// =======================
//  Tetrahedral element
// =======================

FeHCurlHiTet::FeHCurlHiTet() : FeHCurlHi( Elem::ET_TET4) {

}

FeHCurlHiTet::~FeHCurlHiTet() {

}

void FeHCurlHiTet::GetShFnc( Matrix<Double>& shape, 
                             const LocPointMapped& lpm,
                             const Elem* elem, UInt comp ) {
  
  Matrix<Double> locShape;
  CalcLocShFnc2<ID>(locShape, lpm, elem, comp);

  // Perform local->global gradient transformation
  shape =  Transpose(lpm.jacInv) * locShape;
}

void FeHCurlHiTet::GetCurlShFnc( Matrix<Double>& curl, 
                                 const LocPointMapped& lpm,
                                 const Elem* elem, UInt comp ) {

  Matrix<Double> locCurl;    
  CalcLocShFnc2<CURL>( locCurl, lpm, elem, comp );
  
  // Perform local->global curl transformation
  curl = lpm.jac * locCurl;
  curl *= ( 1.0 / std::fabs(lpm.jacDet) );
}

template<FeHCurlHi::DiffType DIFF_TYPE>
void FeHCurlHiTet::CalcLocShFnc2( Matrix<Double>& shape, 
                                  const LocPointMapped& lpm,
                                  const Elem* elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();

  AutoDiff<Double, 3> x (lpm.lp.coord[0],0);
  AutoDiff<Double, 3> y (lpm.lp.coord[1],1);
  AutoDiff<Double, 3> z (lpm.lp.coord[2],2);

  AutoDiff<Double, 3> lambda[4] = { (1.0 - x- y- z),
                                 x,
                                 y,
                                 z};
  UInt pos = 0;
  shape.Resize(3,actNumFncs_);
  shape.Init();
  
  StdVector<AutoDiff<Double, 3> > Vals;
  // ------------------------
  // 1) Edge shape functions
  // ------------------------
  for( UInt i = 0; i < 6; ++i ) {

    //UInt order = orderEdge_[i];
  UInt order = orderEdge_[i];
    UInt ind1 = shape_.edgeVertices[i][0]-1;
    UInt ind2 = shape_.edgeVertices[i][1]-1;
    if ( elem->extended->edges[i] < 0 ) {
      std::swap(ind1, ind2);  // fmax > f1 > f2
    }
    
    // === a) standard Nedelec shape functions ===
    Xpr_Diff_VGradU<3,DIFF_TYPE> xpr( lambda[ind1], lambda[ind2] );
    COPYSHFNC
    

  // === b) gradient functions ===
    if( useEdgeGrad_[i]) {
      if (onlyLowestOrder_) {
        for( UInt j = 0; j < order; ++j ) {
          pos++;
        }
      } else {
        ScaledIntLegendreP2(Vals, order+1, lambda[ind1]+lambda[ind2], lambda[ind2]-lambda[ind1]);

          for( UInt j = 0; j < order; ++j ) {
            Xpr_GradU<3,DIFF_TYPE> xpr(Vals[j]);
            COPYSHFNC
          }
      } // if: lowestOrder
    } //if: edgeGrad
  } //loop: edges
     
  if(onlyLowestOrder_) return;
  
  StdVector<AutoDiff<Double, 3> > ui;
  StdVector<AutoDiff<Double, 3> > vj;
  StdVector<AutoDiff<Double, 3> > temp_vj;
  // -------------------------
  // 2) Face shape functions
  // -------------------------

#ifdef USE_FACES
  for( UInt iFace = 0; iFace < 4; ++iFace ) {
    //only valid for isotropic polynomial order!! Is there a way
    //to use anisotropic order? Maybe via a transformation to hexahedral
    //there we could apply anisotropy and then transform it back to tet.
    //The problem is, that the polynomial-"isosurfaces" are then warped
    //in the tetrahedron?!
    UInt order = orderFace_[iFace][0];
    if (order >1) {
        // get unique sorting of the face
        const StdVector<UInt>& unsorted = shape_.faceNodes[iFace];
        StdVector<UInt> ind;
        Face::GetSortedIndices( ind, unsorted, 3, elem->extended->faceFlags[iFace]);

        // calculate face extension parameter which is the sum
        // of all lambdas of one face
        AutoDiff<Double,3> sum_lambda = 0.0;
        for( UInt i = 0; i < 3; ++i){
          sum_lambda += lambda[ind[i]];
        }

        //definition of ui and vj according to PHD thesis of Sabine Zaglmayr p.103
         ScaledIntLegendreP2(ui, order+1, lambda[ind[0]]+lambda[ind[1]], lambda[ind[1]]-lambda[ind[0]]);
         ScaledLegendre(temp_vj, order+1, sum_lambda, 2.0 * lambda[ind[2]] - sum_lambda);
         vj.Init();
         vj.Resize(order);
          for (UInt j=0; j<order; ++j){
           vj[j] = lambda[ind[2]] * temp_vj[j];
          }


         // === a) type 1: gradient fields ===
         if( useFaceGrad_[iFace]) {
             for( UInt i = 0; i <= order - 2; ++i ) {
               for( UInt j = 0; j <= order - 2 - i; ++j ) {
               Xpr_GradU<3,DIFF_TYPE> xpr(ui[i]*vj[j]);
               COPYSHFNC
             }
           }
         }

         // === b) type 2: face functions ===
         for( UInt i = 0; i <= order - 2; ++i ) {
           for( UInt j = 0; j <= order - 2 - i; ++j ) {
             Xpr_Diff_VGradU<3,DIFF_TYPE> xpr(ui[i], vj[j]);
             COPYSHFNC
           }
         }

         // === c) type 3: face functions ===
         for( UInt j = 0; j <= order - 2; ++j ) {
           Xpr_Diff_SVGradU<3,DIFF_TYPE> xpr(lambda[ind[0]], lambda[ind[1]], vj[j] );
           COPYSHFNC
         }

    } //if order >0

  } //loop over all faces
#endif



  StdVector<AutoDiff<Double, 3> > wk;
  StdVector<AutoDiff<Double, 3> > temp_wk;
  // -------------------------
  // 3) Interior shape functions
  // -------------------------
#ifdef USE_INNER
  // only isotropic polynomial order
    if( orderInner_[0] > 2) {

        //definition of ui, vj and wk according to PHD thesis of Sabine Zaglmayr p.103
      ScaledIntLegendreP2(ui, orderInner_[0]+1, lambda[0]+lambda[1], lambda[1]-lambda[0]);
        ScaledLegendre(temp_vj, orderInner_[0]+1, 1.0-lambda[3], 2.0*lambda[2]-(1.0-lambda[3]));
        Legendre( temp_wk, orderInner_[0]+1, 2.0*lambda[3]-1.0);
        vj.Resize(orderInner_[0]);
        wk.Resize(orderInner_[0]);
        for (UInt j=0; j<orderInner_[0]; ++j){
         vj[j] = lambda[2] * temp_vj[j];
        }
        for (UInt k=0; k<orderInner_[0]; ++k){
         wk[k] = lambda[3] * temp_wk[k];
        }



        // === a) type 1: gradient fields ===
        if( useInteriorGrad_ ) {
        for( UInt i = 0; i <= orderInner_[0] - 3; ++i ) {
          for( UInt j = 0; j <= orderInner_[0] - 3 - i; ++j ) {
            for ( UInt k = 0; k <= orderInner_[0] - 3 - i - j; ++k){
              Xpr_GradU<3, DIFF_TYPE> xpr( ui[i] * vj[j] * wk[k] );
              COPYSHFNC
            }
          }
        }
        }


        // === b) type 2 volume functions ===
        for( UInt i = 0; i <= orderInner_[0] - 3; ++i ) {
          for( UInt j = 0; j <= orderInner_[0] - 3 - i; ++j ) {
            for( UInt k = 0; k < orderInner_[0] - 3 - i - j; ++k ) {

              Xpr_Diff_VWGradU<3, DIFF_TYPE> VWGradU( ui[i], vj[j], wk[k] );
              Xpr_Diff_VWGradU<3, DIFF_TYPE> UWGradV( vj[j], ui[i], wk[k] );
              Xpr_Diff_VWGradU<3, DIFF_TYPE> UVGradW( wk[k], ui[i], vj[j] );
              for( UInt n = 0; n < 3; ++n ) {
                shape[n][pos]   = VWGradU[n] - UWGradV[n] + UVGradW[n];
                shape[n][pos+1] = VWGradU[n] + UWGradV[n] - UVGradW[n];
              }
              pos+=2;
            }
          }
        }

        // === c) type 3 volume functions ===
        for( UInt j = 0; j < orderInner_[0] - 3; ++j ) {
          for( UInt k = 0; k < orderInner_[0] - 3 - j; ++k ) {
            Xpr_Diff_SVGradU<3, DIFF_TYPE> xpr(lambda[0], lambda[1], vj[j]*wk[k] );
            COPYSHFNC
          }
        }
    }
#endif
}

void FeHCurlHiTet::CalcNumUnknowns() {
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

     if( useEdgeGrad_[i]) {
       if (onlyLowestOrder_) {return;}
       else{
         unknowns += orderEdge_[i];
       }
     }
     edgeFncs[i] = unknowns;
     LOG_DBG(feHCurlHi) <<   "edge " << i+1 << " has " << unknowns << "unknowns";
     actNumFncs_ += unknowns;
    //  std::cout << "[TET] edge " << i
    //     << "  order=" << orderEdge_[i]
    //     << "  grad="  << (useEdgeGrad_[i] ? "yes" : "no")
    //     << "  -> "    << i << " DOFs\n";

   }

   // Faces
   StdVector<UInt>& faceFncs = entityFncs_[FACE];
   faceFncs.Resize(shape_.numFaces);
   faceFncs.Init(0);


#ifdef USE_FACES
  for( UInt iFace = 0; iFace < 4; ++iFace ) {
    unknowns = 0;
    UInt order = orderFace_[iFace][0];
    if (order >1) {
      unknowns=0;

         // === a) type 1: gradient fields ===
         if( useFaceGrad_[iFace]) {
             for( UInt i = 0; i <= order - 2; ++i ) {
               for( UInt j = 0; j <= order - 2 - i; ++j ) {
                 unknowns+=1;
             }
           }
         }

         // === b) type 2: face functions ===
         for( UInt i = 0; i <= order - 2; ++i ) {
           for( UInt j = 0; j <= order - 2 - i; ++j ) {
             unknowns+=1;
           }
         }

         // === c) type 3: face functions ===
         for( UInt j = 0; j <= order - 2; ++j ) {
           unknowns+=1;
         }
    } //if order >0

    faceFncs[iFace] = unknowns;
    LOG_DBG(feHCurlHi) << "face " << iFace+1 << " has " << unknowns << "unknowns";
    actNumFncs_ += unknowns;
    // std::cout << "[TET] face " << iFace
    // << "  orders=(" << orderFace_[iFace][0] << "," << orderFace_[iFace][1] << ")"
    // << "  grad="    << (useFaceGrad_[iFace] ? "yes" : "no")
    // << "  -> "      << unknowns << " DOFs\n";

  } //loop over all faces
#endif


   // Interior
   StdVector<UInt>& innerFncs = entityFncs_[INTERIOR];
   innerFncs.Resize(1);
   innerFncs.Init(0);

#ifdef USE_INNER
    if( orderInner_[0] > 2) {
      unknowns=0;

        // === a) type 1: gradient fields ===
        if( useInteriorGrad_ ) {
        for( UInt i = 0; i <= orderInner_[0] - 3; ++i ) {
          for( UInt j = 0; j <= orderInner_[0] - 3 - i; ++j ) {
            for ( UInt k = 0; k <= orderInner_[0] - 3 - i - j; ++k){
              unknowns+=1;
            }
          }
        }
        }

        // === b) type 2 volume functions ===
        for( UInt i = 0; i <= orderInner_[0] - 3; ++i ) {
          for( UInt j = 0; j <= orderInner_[0] - 3 - i; ++j ) {
            for( UInt k = 0; k < orderInner_[0] - 3 - i - j; ++k ) {
              unknowns+=1;
            }
          }
        }

        // === c) type 3 volume functions ===
        for( UInt j = 0; j < orderInner_[0] - 3; ++j ) {
          for( UInt k = 0; k < orderInner_[0] - 3 - j; ++k ) {
            unknowns+=1;
          }
        }

       actNumFncs_ += unknowns;
       innerFncs[0] = unknowns;
       LOG_DBG(feHCurlHi) << "interior has " << unknowns << "unknowns";
    }

#endif

   LOG_DBG(feHCurlHi) <<  "totalUnknowns: " << actNumFncs_  << std::endl;

   for (UInt f = 0; f < shape_.numFaces; ++f)
    PYRA_DOF_PRINT("TET ", onlyLowestOrder_, "face", f, entityFncs_[FACE][f]);
for (UInt e = 0; e < shape_.numEdges; ++e)
    PYRA_DOF_PRINT("TET ", onlyLowestOrder_, "edge", e, entityFncs_[EDGE][e]);
}



// =======================
//  Pyramidal element
// =======================

FeHCurlHiPyra::FeHCurlHiPyra() : FeHCurlHi( Elem::ET_PYRA5) {

}

FeHCurlHiPyra::~FeHCurlHiPyra() {

}

void FeHCurlHiPyra::GetShFnc( Matrix<Double>& shape, 
                             const LocPointMapped& lpm,
                             const Elem* elem, UInt comp ) {
  
  Matrix<Double> locShape;
  CalcLocShFnc2<ID>(locShape, lpm, elem, comp);

  // Perform local->global gradient transformation
  shape =  Transpose(lpm.jacInv) * locShape;
}

void FeHCurlHiPyra::GetCurlShFnc( Matrix<Double>& curl, 
                                 const LocPointMapped& lpm,
                                 const Elem* elem, UInt comp ) {

  Matrix<Double> locCurl;    
  CalcLocShFnc2<CURL>( locCurl, lpm, elem, comp );
  
  // Perform local->global curl transformation
  curl = lpm.jac * locCurl;
  curl *= ( 1.0 / std::fabs(lpm.jacDet) );
}
zu checken:
- jacobi matrix wird ja mit H1 elementen berechnet. passt da vl was ned?
- vgl der massenmatrix einträge von hex und tet und pyra
template<FeHCurlHi::DiffType DIFF_TYPE>
void FeHCurlHiPyra::CalcLocShFnc2( Matrix<Double>& shape, 
                                  const LocPointMapped& lpm,
                                  const Elem* elem, UInt comp ) {
  if (updateUnknowns_) CalcNumUnknowns();

    AutoDiff<Double,3>  x (lpm.lp.coord[0],0);
  AutoDiff<Double,3>  y (lpm.lp.coord[1],1);
  AutoDiff<Double,3>  z (lpm.lp.coord[2],2);

  // guarded against (1-z)->0
  const double EPS = 1e-12;
  AutoDiff<Double,3> one(1.0);
  AutoDiff<Double,3> rho = one - z;                 // rho = 1 - z
  double rho_val = rho.Val();
  AutoDiff<Double,3> invrho = (rho_val < EPS) ? AutoDiff<Double,3>(0.0) : (1.0 / rho);
  AutoDiff<Double,3> frac   = (rho_val < EPS) ? AutoDiff<Double,3>(0.0) : (x*y*z) * invrho;

  AutoDiff<Double,3>  lam[5] =
  {
    0.25*((one+x)*(one+y) - z + frac),   // lambda_1
    0.25*((one-x)*(one+y) - z - frac),   // lambda_2
    0.25*((one-x)*(one-y) - z + frac),   // lambda_3
    0.25*((one+x)*(one-y) - z - frac),   // lambda_4
    z                                    // lambda_5 (apex)
  };

  //---------------------------------------------------------------------------
  // 1) allocate output & keep running index
  //---------------------------------------------------------------------------
  if (updateUnknowns_) CalcNumUnknowns();
  UInt pos = 0;
  shape.Resize(3, actNumFncs_);
  shape.Init();

  //---------------------------------------------------------------------------
  // 2) EDGE functions  (horizontal edges first, then vertical)
  //---------------------------------------------------------------------------
  StdVector<AutoDiff<Double,3> > Vals;
  // a) base-quad edges 0...3
  for (UInt e = 0; e < 4; ++e)
  {
    UInt i1 = shape_.edgeVertices[e][0]-1;         // starting vertex
    UInt i2 = shape_.edgeVertices[e][1]-1;         // end vertex
    UInt i3 = (i2+1) % 4;                          // next CCW base vertex
    UInt i4 = (i1+3) % 4;                          // previous base vertex
    if (elem->extended->edges[e] < 0) { std::swap(i1,i2); std::swap(i3,i4); }

    /* ---- 2a-1  lowest-order Nédélec (edge-tangential) */
    {
      Xpr_Diff_UGradV_min_WGradX<3,DIFF_TYPE> xpr(lam[i1], lam[i2]+lam[i3],
                                                  lam[i2], lam[i1]+lam[i4]);
      COPYSHFNC
    }

    /* ---- 2a-2  higher-order edge-gradients   (Zaglmayr §4.1.2) */
    UInt p  = orderEdge_[e];
    if (p > 1 && useEdgeGrad_[e] && !onlyLowestOrder_)
    {
      ScaledIntLegendreP2(Vals, p, lam[i2]+lam[i1], lam[i2]-lam[i1]);
      for (UInt k = 0; k < p; ++k)
      {
        Xpr_GradU<3,DIFF_TYPE>  xpr( Vals[k] );
        COPYSHFNC
      }
    }
  }

  // b) vertical pyramid edges 4...7
  for (UInt e = 4; e < 8; ++e)
  {
    UInt i1 = shape_.edgeVertices[e][0]-1;
    UInt i2 = shape_.edgeVertices[e][1]-1;
    if (elem->extended->edges[e] < 0) std::swap(i1,i2);

    /* ---- 2b-1  lowest order */
    {
      Xpr_Diff_VGradU<3,DIFF_TYPE> xpr(lam[i1], lam[i2]);
      COPYSHFNC
    }

    /* ---- 2b-2  HO gradient bubble  */
    UInt p = orderEdge_[e];
    if (p > 1 && useEdgeGrad_[e] && !onlyLowestOrder_)
    {
      ScaledIntLegendreP2(Vals, p, lam[i2]+lam[i1], lam[i2]-lam[i1]);
      for (UInt k = 0; k < p; ++k)
      {
        Xpr_GradU<3,DIFF_TYPE> xpr( Vals[k]*lam[i2] ); // bubble factor
        COPYSHFNC
      }
    }
  }

  //---------------------------------------------------------------------------
  // 3)  FACE functions   (4 triangular + 1 quadrilateral)
  //---------------------------------------------------------------------------
  //
  //  3.1  Triangular faces (apex + two base vertices)  –  iFace=0...3
  //
  for (UInt f = 1; f < 5; ++f)   // faces 1...4 are the triangular ones
  {
    UInt ind[3] = { shape_.faceVertices[f][0]-1,
                    shape_.faceVertices[f][1]-1,
                    shape_.faceVertices[f][2]-1 };

    UInt p  = orderFace_[f][0];          // Zaglmayr uses (p,q) but we only
    UInt q  = orderFace_[f][1];          //   support p==q => isotropic

    if (!useFaceGrad_[f] || p < 2 || onlyLowestOrder_)
    {   /* skip – nothing to add */   }
    else
    {
      StdVector<AutoDiff<Double,3> > ui, vj;

      // u_i  : integrated Legendre in edge-direction lambda_ind0–lambda_ind1
      ScaledIntLegendreP2(ui, p+1,
                          lam[ind[1]]+lam[ind[0]],
                          lam[ind[1]]-lam[ind[0]]);

      // v_j  : scalar Legendre in “vertical” lambda_ind2
      StdVector<AutoDiff<Double,3> > temp;
      Legendre(temp, q+1, 2.0*lam[ind[2]] - 1.0);
      vj.Resize(q);
      for (UInt j=0;j<q;++j) vj[j] = lam[ind[2]] * temp[j];

      /* ---- Type-I surface gradient fields  */
      for (UInt i=0;i<=p-2;++i)
        for (UInt j=0;j<=q-2-i;++j)
        {
          Xpr_GradU<3,DIFF_TYPE> xpr(ui[i]*vj[j]);
          COPYSHFNC
        }

      /* ---- Type-II  tangential curls (u_i × grad(v_j))  */
      for (UInt i=0;i<=p-2;++i)
        for (UInt j=0;j<=q-1-i;++j)
        {
          Xpr_Diff_UGradV_min_WGradX<3,DIFF_TYPE>
             xpr( ui[i+1], vj[j], ui[i], vj[j] );
          COPYSHFNC
        }
    }
  }

  //
  //  3.2  Base quad face (iFace == 0)
  //
  {
    UInt f = 0;   // face‑0 is the quad base in OpenCFS
    UInt p = orderFace_[f][0];
    UInt q = orderFace_[f][1];

    if (p > 0 && q > 0 && !onlyLowestOrder_)
    {
      //------------------------------------------------------------------
      //  Same coordinate pull-back NGSolve uses:
      //        zeta = x/(1.0-z)   ,   eta = y/(1.0-z)
      //------------------------------------------------------------------
      // Reuse rho, rho_val, invrho from above to avoid division by ~0
      AutoDiff<Double,3> xi, eta;
      if (rho_val < EPS) {
        // limit: (1-z)^2 multiplies all these anyway -> set to 0 safely
        xi  = AutoDiff<Double,3>(0.0);
        eta = AutoDiff<Double,3>(0.0);
      } else {
        xi  = x * invrho;
        eta = y * invrho;
      }

      // Legendre bases
      StdVector<AutoDiff<Double,3>> lp, lq, ip, iq;
      Legendre(lp,  p+1, xi);
      Legendre(lq,  q+1, eta);
      IntLegendre(ip, p+1, xi);
      IntLegendre(iq, q+1, eta);

      // ---------- FAMILY A: surface gradients (add ONLY if useFaceGrad_[0]) ----------
      if (useFaceGrad_[0]) {
          for (UInt i=1; i<=p; ++i)
            for (UInt j=0; j<=q; ++j) {
              AutoDiff<Double,3> phi = ip[i-1]*lq[j]*rho*rho;
              Xpr_GradU<3,DIFF_TYPE> xpr(phi);
              COPYSHFNC
            }

          for (UInt i=0; i<=p; ++i)
            for (UInt j=1; j<=q; ++j) {
              AutoDiff<Double,3> phi = lp[i]*iq[j-1]*rho*rho;
              Xpr_GradU<3,DIFF_TYPE> xpr(phi);
              COPYSHFNC
            }
      }

      // ---------- FAMILY B: “type-2” (p*q)  -- ALWAYS ----------
      for (UInt i=0; i<p; ++i)
        for (UInt j=0; j<q; ++j) {
          AutoDiff<Double,3> u = lp[i+1]*lq[j]*rho*rho;   // NOTE: i+1 exists since i<p
          AutoDiff<Double,3> v = lp[i]*lq[j+1]*rho*rho;
          Xpr_Diff_UGradV_min_WGradX<3,DIFF_TYPE> xpr(u,v,u,v);
          COPYSHFNC
        }

      // ---------- FAMILY C: “type-3” (p + q)  -- ALWAYS ----------
      for (UInt i=0; i<p; ++i) {
        AutoDiff<Double,3> u = ip[i]*rho*rho;          // or lp[i+1]*rho*rho ... keep consistent
        Xpr_GradU<3,DIFF_TYPE> xpr(u);                 // choose the same expr you used on HEX
        COPYSHFNC
      }
      for (UInt j=0; j<q; ++j) {
        AutoDiff<Double,3> v = iq[j]*rho*rho;
        Xpr_GradU<3,DIFF_TYPE> xpr(v);
        COPYSHFNC
      }
    }
  }

  //---------------------------------------------------------------------------
  // 4)  INTERIOR bubble functions  
  //---------------------------------------------------------------------------
  if (useInteriorGrad_ &&
      orderInner_[0]>0 && orderInner_[1]>0 && orderInner_[2]>0 &&
      !onlyLowestOrder_)
  {
    UInt p = orderInner_[0];
    UInt q = orderInner_[1];
    UInt r = orderInner_[2];

    StdVector<AutoDiff<Double,3> >  up,vq,wr;
    IntLegendre(up, p+1, lam[0]-lam[2]);          // just one possible set
    IntLegendre(vq, q+1, lam[1]-lam[3]);
    IntLegendre(wr, r+1, 2.0*lam[4]-1.0);

    for (UInt i=0;i<p;++i)
      for (UInt j=0;j<q;++j)
        for (UInt k=0;k<r;++k)
        {
          AutoDiff<Double,3> phi = up[i]*vq[j]*wr[k]*lam[4];   // bubble
          Xpr_GradU<3,DIFF_TYPE> xpr(phi);                     // grad(phi)
          COPYSHFNC
        }
  }
}

void FeHCurlHiPyra::CalcNumUnknowns() {
  
  actNumFncs_ = 0;

  /* ---- EDGES ---------------------------------------------------------- */
  entityFncs_[EDGE].Resize(shape_.numEdges);
  for (UInt e = 0; e < shape_.numEdges; ++e)
  {
      UInt n = 1;                            // always the lowest-order tangential mode
      if (useEdgeGrad_[e])                   // add all gradient modes if requested
          n += orderEdge_[e];                //  -> total = 1 + p  (p = orderEdge_[e])

      entityFncs_[EDGE][e] = n;
      actNumFncs_         += n;
      // std::cout << "[Pyra] edge " << e
      //     << "  order=" << orderEdge_[e]
      //     << "  grad="  << (useEdgeGrad_[e] ? "yes" : "no")
      //     << "  -> "    << n << " DOFs\n";

  }

  /* ---- TRIANGULAR FACES (1...4) ---------------------------------------- */
  entityFncs_[FACE].Resize(shape_.numFaces);
  entityFncs_[FACE].Init(0);
  for (UInt f = 1; f < 5; ++f) {
    UInt p = orderFace_[f][0];  // isotropic
    UInt q = orderFace_[f][1];

    UInt n = 0;
    if (!onlyLowestOrder_ && p > 1 && q > 1) {
      // Type-I (surface gradients) only if useFaceGrad_
      if (useFaceGrad_[f])
        n += (p-1)*p/2;               // sum_{i=0}^{p-2} (p-1-i)

      // Type-II tangential curls (always)
      n += (p-1)*p/2;

      // Type-III edge-blends (always)
      n += (p-1);
    }

    entityFncs_[FACE][f] = n;
    actNumFncs_         += n;

    // std::cout <<"[Pyra] tria face " << f
    //       << "  orders=(" << orderFace_[f][0] << "," << orderFace_[f][1] << ")"
    //       << "  grad="    << (useFaceGrad_[f] ? "yes" : "no")
    //       << "  -> "      << n << " DOFs\n";

  }

  /* ---- QUAD BASE FACE (index 0) ------------------------------------ */
  {
    const UInt f = 0;
    const UInt p = orderFace_[f][0];
    const UInt q = orderFace_[f][1];

    UInt n = 0;
    if (!onlyLowestOrder_ && p > 0 && q > 0) {
        // mandatory families (must match HEX):
        // type-2 (p*q)  + type-3 (p + q)
        n = p*q + p + q;

        // optional surface-gradient family (HEX adds when useFaceGrad=true)
        if (useFaceGrad_[f])
            n += p*q;     // add the extra p*q block

        // your additional curl-type block ONLY if you really generate it
        // and ONLY for p,q >= 2 (HEX doesn't have that for p=1)
        // If you keep it, count it here as well:
        // if (p > 1 && q > 1) n += p*q_extra;   // decide if you keep it
    }

    entityFncs_[FACE][f] = n;
    actNumFncs_         += n;

    // std::cout <<"[Pyra] quad face " << f
    //       << "  orders=(" << orderFace_[f][0] << "," << orderFace_[f][1] << ")"
    //       << "  grad="    << (useFaceGrad_[f] ? "yes" : "no")
    //       << "  -> "      << n << " DOFs\n";
  }

  /* ---- INTERIOR ------------------------------------------------------ */
  entityFncs_[INTERIOR].Resize(1);
  {
    UInt n = 0;
    if (useInteriorGrad_ && !onlyLowestOrder_)
        n = orderInner_[0] * orderInner_[1] * orderInner_[2];
    entityFncs_[INTERIOR][0] = n;
    actNumFncs_             += n;
    // std::cout <<"[Pyra] interior " 
    //       << "  orders=(" << orderInner_[0]
    //       << "  grad="    << (useInteriorGrad_ ? "yes" : "no")
    //       << "  -> "      << n << " DOFs\n";
  }

  
}

}// end of namespace
