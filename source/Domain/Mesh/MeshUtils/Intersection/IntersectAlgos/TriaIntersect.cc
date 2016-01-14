// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TriaIntersect.cc
 *       \brief    interset tetrahedrons based on
 *                 3D Game Engine Design. A Practical Approach to Real-Time Computer Graphics
 *                 (Morgan Kaufmann Series in Interactive 3D Technology)
 *                 and code fragments found in OpenFOAM which is based on the same citation
 *
 *       \date     Dec 10, 2015
 *       \author   ahueppe
 */
//================================================================================================

#include "TriaIntersect.hh"
#include "Domain/ElemMapping/Elem.hh"


namespace CoupledField{

void TriaIntersect::SetTElem( UInt tNum ){
  const Elem* newTElem = g1_->GetElem(tNum);
  if(tElem_->elemNum != newTElem->elemNum){
    tElem_ = newTElem;
    tTets_ = GetTetsFromElem(tElem_, g1_);
  }
}

bool TriaIntersect::Intersect(UInt sNum){
  const Elem* elGrid2 = g2_->GetElem(sNum);
  //1. Get tetrahedral information
  StdVector<CoordTetra> srcTets = GetTetsFromElem(elGrid2,g2_);
  //this is a n^2 approach. Could be less overhead
  UInt numtTets = tTets_.GetSize();
  UInt numSTets = srcTets.GetSize();
  intersectingTets_.Clear(true);
  for(UInt trgT = 0;trgT < numtTets; ++trgT){
    for(UInt aPlane = 0; aPlane < 4; ++aPlane){
      //now cut each src to cut plane
      for(UInt srcT = 0; srcT < numSTets ; ++srcT){
        SplitAndDecompose(trgT, aPlane, srcTets[srcT], intersectingTets_);
      }
    }
  }
  return (intersectingTets_.GetSize()!=0);
}

void TriaIntersect::GetVolumeAndCenters(StdVector<VolCenterInfo>& infos){
  UInt numTets = intersectingTets_.GetSize();
  infos.Resize(numTets);
  Vector<Double> tC(3);

  for(UInt aTet = 0;aTet < numTets; ++aTet){
    Double& volume = infos[aTet].volume;
    Vector<Double>& center = infos[aTet].center;
    center.Resize(3);
    for(unsigned int tetI=0; tetI < numTets; tetI++) {
      const CoordTetra& t = intersectingTets_[tetI];
      // Calculate volume (no check for orientation)
      //TODO needs better support from vector class!
      Double a1 = ( (t[1][1] - t[0][1]) * (t[2][2] - t[0][2]) ) - ( (t[1][2] - t[0][2]) * (t[2][1] - t[0][1]) );
      a1 *= (t[3][0] - t[0][0]);

      Double a2 = ( (t[1][2] - t[0][2]) * (t[2][0] - t[0][0]) ) - ( (t[1][0] - t[0][0]) * (t[2][2] - t[0][2]) );
      a2 *= (t[3][1] - t[0][1]);

      Double a3 = ( (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) ) - ( (t[1][1] - t[0][1]) * (t[2][0] - t[0][0]) );
      a3 *= (t[3][2] - t[0][2]);

      Double tV = (1.0/6.0) * std::sqrt(a1*a1+a2*a2+a3*a3);
      //      Vector<Double> sub1 = (t[1] - t[0]);
      //      Vector<Double> sub2 = (t[2] - t[0]);
      //      Vector<Double> sub3 = (t[3] - t[0]);
      //      Vector<Double> cross;
      //      sub1.CrossProduct(sub2,cross);
      //      sub1 = cross * sub3;
      //      Double tV = (1.0/6.0)*sub1.NormL2();
      // Calculate centroid
      //tC = (0.25 * (t[0] + t[1] + t[2] + t[3]));

      for(UInt d=0;d<3;++d){
        tC[d] = 0.25 * (t[0][d] + t[1][d] + t[2][d] + t[3][d]);
        center[d] += (tV * tC[d]);
      }
      volume += tV;

    }
    center /= (volume + CoordTetra::EPS);
  }
}

inline StdVector<CoordTetra> TriaIntersect::GetTetsFromElem(const Elem* newTElem, Grid* aGrid){
    StdVector<CoordTetra> genTets;

    //now get triangular information
    shared_ptr<ElemShapeMap> esm = aGrid->GetElemShapeMap(newTElem,true);
    BaseFE* feElement = esm->GetBaseFE();
    feElement->Triangulate(lastTetIdx_);
    //loop over each Tet
    UInt numTets = lastTetIdx_.GetSize();
    genTets.Resize(numTets);
    for(UInt aT = 0;aT<numTets;++aT){
      const StdVector<UInt>& aTetIdx = lastTetIdx_[aT];
      CoordTetra & aTet = genTets[aT];
      for(UInt aNode =0;aNode<4;++aNode){
        UInt nodeNum = tElem_->connect[aTetIdx[aNode]];
        aGrid->GetNodeCoordinate3D( aTet[aNode], nodeNum,true);
      }
    }
    return genTets;
  }

void TriaIntersect::GetIntersectionElems(StdVector<IntersectionElem*>& interElems){
  EXCEPTION("Not yet implemented");
}

inline void TriaIntersect::SplitAndDecompose(UInt tetIdx, UInt planeIdx, CoordTetra& tetra,
                              StdVector<CoordTetra>& genTets){

  StdVector<Double> C(4);
  CoordTetra tmpTetra;
  //flags for indicating the halfspace location
  StdVector<UInt> pos(4), neg(4), zero(4);
  UInt i = 0, nPos = 0, nNeg = 0, nZero = 0;
  // Fetch reference to plane
  CoordTetra::CutPlane& tetPlane = tTets_[tetIdx].GetPlane(planeIdx);
  for (i = 0; i < 4; ++i) {
    // Compute distance to plane
    C[i] = (tetra[i] * tetPlane.first) - tetPlane.second;
    if (C[i] > 0.0) {
      pos[nPos++] = i;
    }
    else{
      if (C[i] < 0.0){
        neg[nNeg++] = i;
      } else {
        zero[nZero++] = i;
      }
    }
  }
  if (nNeg == 0){
    //tetrahedron is not intersected
    return;
  }
  if (nPos == 0){
    //tetrahedron is completely embedded
    genTets.Push_back(tetra);
    return;
  }

  //perform the intersection
  Double w0, w1, invCDiff;
  StdVector< Vector<Double> > intp(4);
  if (nPos == 3) {
    // +++-
    for (i = 0; i < nPos; ++i) {
      invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));
      w0 = -C[neg[0]] * invCDiff;
      w1 = +C[pos[i]] * invCDiff;
      tetra[pos[i]] = (tetra[pos[i]] * w0) + (tetra[neg[0]] * w1);
    }
    genTets.Push_back(tetra);
  } else {
    if (nPos == 2) {
      if (nNeg == 2) {
        // ++--
        for (i = 0; i < nPos; ++i) {
          invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));
          w0 = -C[neg[0]] * invCDiff;
          w1 = +C[pos[i]] * invCDiff;
          intp[i] = (tetra[pos[i]] * w0) + (tetra[neg[0]] * w1);
        }
        for (i = 0; i < nNeg; ++i) {
          invCDiff = (1.0 / (C[pos[i]] - C[neg[1]]));
          w0 = -C[neg[1]] * invCDiff;
          w1 = +C[pos[i]] * invCDiff;
          intp[i+2] = (tetra[pos[i]] * w0) + (tetra[neg[1]] * w1);
        }
        tetra[pos[0]] = intp[2];
        tetra[pos[1]] = intp[1];
        genTets.Push_back(tetra);
        tmpTetra[0] = tetra[neg[1]];
        tmpTetra[1] = intp[3];
        tmpTetra[2] = intp[2];
        tmpTetra[3] = intp[1];
        genTets.Push_back(tmpTetra);
        tmpTetra[0] = tetra[neg[0]];
        tmpTetra[1] = intp[0];
        tmpTetra[2] = intp[1];
        tmpTetra[3] = intp[2];
        genTets.Push_back(tmpTetra);
      } else {
        // ++-0
        for (i = 0; i < nPos; ++i) {
          invCDiff = (1.0 / (C[pos[i]] - C[neg[0]]));
          w0 = -C[neg[0]] * invCDiff;
          w1 = +C[pos[i]] * invCDiff;
          tetra[pos[i]] = (tetra[pos[i]] * w0) + (tetra[neg[0]] * w1);
        }
        genTets.Push_back(tetra);
      }
    } else {
      if (nPos == 1) {
        if (nNeg == 3) {
          // +---
          for (i = 0; i < nNeg; ++i) {
            invCDiff = (1.0 / (C[pos[0]] - C[neg[i]]));
            w0 = -C[neg[i]] * invCDiff;
            w1 = +C[pos[0]] * invCDiff;
            intp[i] = (tetra[pos[0]] * w0) + (tetra[neg[i]] * w1);
          }
          tetra[pos[0]] = intp[0];
          genTets.Push_back(tetra);
          tmpTetra[0] = intp[0];
          tmpTetra[1] = tetra[neg[1]];
          tmpTetra[2] = tetra[neg[2]];
          tmpTetra[3] = intp[1];
          genTets.Push_back(tmpTetra);
          tmpTetra[0] = tetra[neg[2]];
          tmpTetra[1] = intp[1];
          tmpTetra[2] = intp[2];
          tmpTetra[3] = intp[0];
          genTets.Push_back(tmpTetra);
        } else {
          if (nNeg == 2){
            // +--0
            for (i = 0; i < nNeg; ++i) {
              invCDiff = (1.0 / (C[pos[0]] - C[neg[i]]));
              w0 = -C[neg[i]] * invCDiff;
              w1 = +C[pos[0]] * invCDiff;
              intp[i] = (tetra[pos[0]] * w0) + (tetra[neg[i]] * w1);
            }
            tetra[pos[0]] = intp[0];
            genTets.Push_back(tetra);
            tmpTetra[0] = intp[1];
            tmpTetra[1] = tetra[zero[0]];
            tmpTetra[2] = tetra[neg[1]];
            tmpTetra[3] = intp[0];
            genTets.Push_back(tmpTetra);
          }else {
            // +-00
            invCDiff = (1.0 / (C[pos[0]] - C[neg[0]]));
            w0 = -C[neg[0]] * invCDiff;
            w1 = +C[pos[0]] * invCDiff;
            tetra[pos[0]] = (tetra[pos[0]] * w0) + (tetra[neg[0]] * w1);
            genTets.Push_back(tetra);
          }
        }
      }
    }
  }
}

}
