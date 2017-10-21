// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CutTetra.hh
 *       \brief    <Description>
 *
 *       \date     Dez 11, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_INTERSECTALGOS_COORDTETRA_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_INTERSECTALGOS_COORDTETRA_HH_

#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"

namespace CoupledField{

struct CoordTetra{
  //cutplane defined by vector and magnitude
  typedef std::pair<Vector<Double>, Double> CutPlane;

  CoordTetra() :
    points(4),
    clipPlanes_(4),
    tetSize_(0.0),
    planesComputed_(false),
    e10(3),e20(3),e30(3),e21(3),e31(3)
  {

    points.Init(Vector<Double>(3));
    AllocateClipPlaneMem();

  }

  CoordTetra(const CoordTetra & origTet)  :
      points(4),
      clipPlanes_(4),
      tetSize_(0.0),
      planesComputed_(false),
      e10(3),e20(3),e30(3),e21(3),e31(3)
  {
    points.Init(Vector<Double>(3));
    for(UInt i=0;i<4;++i){
      for(UInt d=0;d<3;++d){
        points[i][d] = origTet.points[i][d];
      }
    }
    planesComputed_ = false;
    AllocateClipPlaneMem();
    if(origTet.planesComputed_){
      clipPlanes_ = origTet.clipPlanes_;
      planesComputed_ = true;
    }
    tetSize_ = origTet.tetSize_;
  }


  inline Vector<Double>& operator[](UInt i){
    assert(i<4);
    return points[i];
  }

  inline const Vector<Double>& operator[](UInt i) const {
    assert(i<4);
    return points[i];
  }

  inline CutPlane& GetPlane(UInt i){
    assert(i<4);
    if(!planesComputed_){
      ComputeClipPlanes();
    }
    return clipPlanes_[i];
  }

  inline void ScaleTet(Double factor){
    for(UInt i=0;i<4;++i){
      for(UInt d=0;d<3;++d){
        points[i][d] *= factor;
      }
    }
    InitClipPlanes();
    planesComputed_ = false;
  }

  void ReversePoints(){
    //WARNING THIS RENDERS THE CUT PLANES INVALID
    InitClipPlanes();
    planesComputed_ = false;
    Vector<Double> dummy = points[3];
    points[3] = points[2];
    points[2] = points[1];
    points[1] = points[0];
    points[0] = dummy;
  }

private:


  inline void AllocateClipPlaneMem(){
    for(UInt i =0;i<4;++i){
      clipPlanes_[i].first.Resize(3);
    }
  }
  inline void InitClipPlanes(){
    for(UInt i =0;i<4;++i){
      clipPlanes_[i].first.Init();
      clipPlanes_[i].second = 0;
    }
  }

  void ComputeClipPlanes();

  StdVector< Vector<Double> > points;

  StdVector<CutPlane> clipPlanes_;

  Double tetSize_;

  bool planesComputed_;

  //=============================================================
  // temporary storage
  //=============================================================
  Vector<Double> e10,e20,e30,e21,e31;


};

}
#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_INTERSECTALGOS_COORDTETRA_HH_ */
