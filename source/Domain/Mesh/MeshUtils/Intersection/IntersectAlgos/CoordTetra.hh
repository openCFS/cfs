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

  CoordTetra(){
    points.Resize(4);
    for(UInt i=0;i<4;++i){
      points[i].Resize(3);
      points[i].Init();
    }
    planesComputed_ = false;
    tetSize_ = 0;
  }

  CoordTetra(const CoordTetra & origTet)
  {
    for(UInt i=0;i<4;++i){
      points[i] = origTet.points[i];
    }
    planesComputed_ = false;
    if(origTet.planesComputed_){
      clipPlanes_ = origTet.clipPlanes_;
      planesComputed_ = true;
    }
    tetSize_ = 0;
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

  void ReversePoints(){
    //WARNING THIS RENDERS THE CUT PLANES INVALID
    clipPlanes_.Clear(true);
    planesComputed_ = false;
    Vector<Double> dummy = points[3];
    points[3] = points[2];
    points[2] = points[1];
    points[1] = points[0];
    points[0] = dummy;
  }

private:

  void ComputeClipPlanes();

  StdVector< Vector<Double> > points;

  StdVector<CutPlane> clipPlanes_;

  Double tetSize_;

  bool planesComputed_;


};

}


#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_INTERSECTALGOS_COORDTETRA_HH_ */
