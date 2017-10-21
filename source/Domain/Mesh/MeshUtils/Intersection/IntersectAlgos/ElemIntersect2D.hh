// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     PlaneIntersect.hh
 *       \brief    <Description>
 *
 *       \date     Apr 13, 2016
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_INTERSECTALGOS_ELEMINTERSECT2D_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_INTERSECTALGOS_ELEMINTERSECT2D_HH_

#include "ElemIntersect.hh"

namespace CoupledField{

class ElemIntersect2D : public ElemIntersect{

public:


  //! \copydoc ElemIntersect::ElemIntersect
  ElemIntersect2D(Grid* targetGrid, Grid* sourceGrid, bool isRegular = false)
      : ElemIntersect(targetGrid,sourceGrid){
    scaleFac_ = 1e4;
    tElemNum_ = 0;
    sElemNum_ = 0;
  }

  //! copy constructor
  ElemIntersect2D(const ElemIntersect2D& inter)
   : ElemIntersect(inter){

    this->sGrid_ = inter.sGrid_;
    this->tGrid_ = inter.tGrid_;
    tElemNum_ = inter.tElemNum_;
    sElemNum_ = inter.sElemNum_;
    lastIntersection_ = inter.lastIntersection_;
  }

  //! assignment
  ElemIntersect2D & operator=(const ElemIntersect2D& inter){
    this->sGrid_ = inter.sGrid_;
    this->tGrid_ = inter.tGrid_;
    tElemNum_ = inter.tElemNum_;
    sElemNum_ = inter.sElemNum_;
    lastIntersection_ = inter.lastIntersection_;
    return *this;
  }

  //! deep pointer copy
  virtual ElemIntersect2D* Clone(){
     return new ElemIntersect2D(*this);
  }


  virtual ~ElemIntersect2D(){

  }

  //! \copydoc ElemIntersect::SetTElem
  virtual void SetTElem(UInt tNum);

  //! \copydoc ElemIntersect::Intersect
  virtual bool Intersect(UInt sNum);

  //! \copydoc ElemIntersect::GetIntersectionElems
  virtual void GetIntersectionElems(StdVector<IntersectionElem*>& interElems);

  //! \copydoc ElemIntersect::GetVolumeAndCenters
  virtual void GetVolumeAndCenters(StdVector<VolCenterInfo>& infos);

protected:

private:

  //! intersection of two arbitrary 2D elements
  Double IntersectPolys( const Elem * trgElem , const Elem *srcElem);

  //! CGAL version of polygon intersection
  Double CutPolysCGAL(StdVector<Vector<Double> > &p1,
      StdVector<Vector<Double> > &p2, const bool coplanar,
      StdVector<Vector<Double> > &r);

  UInt TriangulatePoly( const StdVector< Vector<Double> > &p,
                        StdVector<MortarNcSurfElem*> &tri,
                        StdVector<UInt> &newNodes );

  void PolyCentroid(const StdVector< Vector<Double> > &p,
                    Vector<Double> &c);

  Double PolyVolume(const StdVector< Vector<Double> > &p);

  //! storage of target element number
  UInt tElemNum_;

  //! storage of source element number
  UInt sElemNum_;

  //! scale factor for element intersections and volume conputations
  Double scaleFac_; // = 1e4;

  // dummy intersection grid structures
  StdVector< Vector<Double> > lastIntersection_;

  //! flag if 2D element is surface
  bool is2dTrg_;

  //cachings
  Double lastVolume_;
  StdVector< Vector<Double> > tPoly_,sPoly_;
  StdVector< Vector<Double> > p1Poly_,p2Poly_,rPoly_,oldPoly_;
  Vector<Double> temp1_,temp2_;
  Vector<Double> e_,c1_,c2_;
  Matrix<Double> rMat_,rMatTrans_;
  StdVector<Vector<Double> > p1Rot_,p2Rot_;


};

}


#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_INTERSECTALGOS_ELEMINTERSECT2D_HH_ */
