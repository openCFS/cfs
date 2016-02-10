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
#include "FeBasis/H1/H1ElemsLagExpl.hh"
#include <fstream>

namespace CoupledField{

void TriaIntersect::InitElemMap(){
  refFeMap[Elem::ET_LINE2]   = new FeH1LagrangeLine1();
  refFeMap[Elem::ET_LINE3]   = new FeH1LagrangeLine2();
  refFeMap[Elem::ET_TRIA3]   = new FeH1LagrangeTria1();
  refFeMap[Elem::ET_TRIA6]   = new FeH1LagrangeTria2();
  refFeMap[Elem::ET_QUAD4]   = new FeH1LagrangeQuad1();
  refFeMap[Elem::ET_QUAD8]   = new FeH1LagrangeQuad2();
  refFeMap[Elem::ET_QUAD9]   = new FeH1LagrangeQuad9();
  refFeMap[Elem::ET_TET4]    = new FeH1LagrangeTet1();
  refFeMap[Elem::ET_TET10]   = new FeH1LagrangeTet2();
  refFeMap[Elem::ET_HEXA8]   = new FeH1LagrangeHex1();
  refFeMap[Elem::ET_HEXA20]  = new FeH1LagrangeHex2();
  refFeMap[Elem::ET_HEXA27]  = new FeH1LagrangeHex27();
  refFeMap[Elem::ET_WEDGE6]  = new FeH1LagrangeWedge1();
  refFeMap[Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
  refFeMap[Elem::ET_WEDGE18] = new FeH1LagrangeWedge18();
  refFeMap[Elem::ET_PYRA5]   = new FeH1LagrangePyra1();
  refFeMap[Elem::ET_PYRA13]  = new FeH1LagrangePyra2();
  refFeMap[Elem::ET_PYRA14]  = new FeH1LagrangePyra14();

}

void TriaIntersect::SetTElem( UInt tNum ){
  const Elem* newTElem = tGrid_->GetElem(tNum);
  if(tElemNum_ != newTElem->elemNum){
    GetTetsFromElem(newTElem, tGrid_,tTets_);
  }
}

bool TriaIntersect::Intersect(UInt sNum){
  /*
   * We intersect here all tetras from triangulation
   * the key in this method is the cut against the four planes
   * of each tetrahedron.
   * we initially cover the current pair of tetras and begin with
   * cutting the source tetra against the first plane
   * the rsult of this is used as a new input for the cut against the
   * next tetra plane. and so on. only the remaining tetras are stored
   */

  //clear last run but retain storage to avoid realloc
  lastSrcTets_.Clear(true);
  intersectingTets_.Clear(true);
  currentTets.Clear(true);
  tmpTets.Clear(true);

  //1. Get tetrahedral information
  const Elem* elGrid2 = sGrid_->GetElem(sNum);
  sElemNum_ = sNum;
  GetTetsFromElem(elGrid2,sGrid_,lastSrcTets_);

  //Scale every tatrahedron
  ScaleTetras(lastSrcTets_,scaleFac_);
  ScaleTetras(tTets_,scaleFac_);

  //export source and target tetras
  //ExportTetras(lastSrcTets_,"sources");
  //ExportTetras(tTets_,"targets");

  //this is a n^2 approach. Could be less overhead
  UInt numtTets = tTets_.GetSize();
  UInt numSTets = lastSrcTets_.GetSize();

  for(UInt trgT = 0;trgT < numtTets; ++trgT){
    for(UInt srcT = 0; srcT < numSTets ; ++srcT){
      currentTets.Resize(1);
      currentTets[0] = lastSrcTets_[srcT];

      tmpTets.Clear(true);
      for(UInt aPlane = 0; aPlane < 4; ++aPlane){
        //now cut each src to cut plane
        for(UInt cutT = 0; cutT < currentTets.GetSize() ; ++cutT){
          SplitAndDecompose(trgT, aPlane, currentTets[cutT] , tmpTets);
        }
        currentTets = tmpTets;
        tmpTets.Clear(true);      
      }
      //ExportTetras(currentTets,"Intersect");
      //store the remaining tets
      for(UInt cutT = 0; cutT < currentTets.GetSize() ; ++cutT){
        intersectingTets_.Push_back(currentTets[cutT]);
      }
    }
  }
  //scale back
  Double invScale = 1.0/scaleFac_;
  ScaleTetras(lastSrcTets_,invScale);
  ScaleTetras(tTets_,invScale);
  ScaleTetras(intersectingTets_,invScale);
  //ExportTetras(intersectingTets_,"Intersect");
  return (intersectingTets_.GetSize()!=0);
}

void TriaIntersect::GetVolumeAndCenters(StdVector<VolCenterInfo>& infos){
  UInt numTets = intersectingTets_.GetSize();
  if(numTets == 0){
    infos.Resize(0);
    return;
  }
  infos.Resize(1);
  Vector<Double> tC(3);
  
  //scale it
  ScaleTetras(intersectingTets_,scaleFac_);

  Double& volume = infos[0].volume;
  Vector<Double>& center = infos[0].center;
  infos[0].targetElemNum = tElemNum_;
  infos[0].sourceElemNum = sElemNum_;
  center.Resize(3);
  center.Init();
  volume = 0.0;
  UInt reversedCounter = 0;
  for(unsigned int tetI=0; tetI < numTets; tetI++) {
    CoordTetra& t = intersectingTets_[tetI];
    // Calculate volume
    //basic formula : volume = ((p1 - p0)\times(p2-p0))\cdot(p3-p0)
    //in case of degenerated elements, we only try to correct the connectivity
    //if the volume is smaller than -1e-10 otherwise, the element will be ignored
    // Question: could we also take the absolute value for now?
    //TODO needs better support from vector class to improve performance!
    Double a1 = ( (t[1][1] - t[0][1]) * (t[2][2] - t[0][2]) ) - ( (t[1][2] - t[0][2]) * (t[2][1] - t[0][1]) );
    a1 *= (t[3][0] - t[0][0]);

    Double a2 = ( (t[1][2] - t[0][2]) * (t[2][0] - t[0][0]) ) - ( (t[1][0] - t[0][0]) * (t[2][2] - t[0][2]) );
    a2 *= (t[3][1] - t[0][1]);

    Double a3 = ( (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) ) - ( (t[1][1] - t[0][1]) * (t[2][0] - t[0][0]) );
    a3 *= (t[3][2] - t[0][2]);

    Double tV = (1.0/6.0) * (a1+a2+a3);
    //in case the volume is smaller than zero, maybe the tetra was wrongly
    //created during intersection. for now we ignore it
    if(tV < 1e-32){
      if(reversedCounter < 4){
        t.ReversePoints();
        tetI--;
        reversedCounter++;
      }else{
        reversedCounter = 0;
      }    
      continue;
    }else if(tV == 0){
      continue;
    }
    reversedCounter = 0;
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
    //
  }
  if(volume>0){
    center /= (volume);
  }else{
    volume = 0;
  }
  Double invScale = 1.0/scaleFac_;
  ScaleTetras(intersectingTets_,invScale);
  center *= invScale;
  volume *= invScale*invScale*invScale;

}

inline void TriaIntersect::GetTetsFromElem(const Elem* newTElem, Grid* aGrid, StdVector<CoordTetra> & genTets){
    //now get triangular information
    refFeMap[newTElem->type]->Triangulate(lastTetIdx_);
    //loop over each Tet
    UInt numTets = lastTetIdx_.GetSize();
    genTets.Resize(numTets);
    for(UInt aT = 0;aT<numTets;++aT){
      const StdVector<UInt>& aTetIdx = lastTetIdx_[aT];
      //CoordTetra & aTet = genTets[aT];
      for(UInt aNode =0;aNode<4;++aNode){
        UInt nodeNum = newTElem->connect[aTetIdx[aNode]];
        aGrid->GetNodeCoordinate( genTets[aT][aNode], nodeNum,true);
      }
    }
    return;
  }

void TriaIntersect::GetIntersectionElems(StdVector<IntersectionElem*>& interElems){
  EXCEPTION("Not yet implemented");
}

inline void TriaIntersect::SplitAndDecompose(UInt tetIdx, UInt planeIdx, CoordTetra& tetra,
                              StdVector<CoordTetra>& genTets){

  StdVector<Double> C(4);

  //flags for indicating the halfspace location
  pos.Init(0);
  neg.Init(0);
  zero.Init(0);
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
        tmpTetra_[0] = tetra[neg[1]];
        tmpTetra_[1] = intp[3];
        tmpTetra_[2] = intp[2];
        tmpTetra_[3] = intp[1];
        genTets.Push_back(tmpTetra_);
        tmpTetra_[0] = tetra[neg[0]];
        tmpTetra_[1] = intp[0];
        tmpTetra_[2] = intp[1];
        tmpTetra_[3] = intp[2];
        genTets.Push_back(tmpTetra_);
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
          tmpTetra_[0] = intp[0];
          tmpTetra_[1] = tetra[neg[1]];
          tmpTetra_[2] = tetra[neg[2]];
          tmpTetra_[3] = intp[1];
          genTets.Push_back(tmpTetra_);
          tmpTetra_[0] = tetra[neg[2]];
          tmpTetra_[1] = intp[1];
          tmpTetra_[2] = intp[2];
          tmpTetra_[3] = intp[0];
          genTets.Push_back(tmpTetra_);
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
            tmpTetra_[0] = intp[1];
            tmpTetra_[1] = tetra[zero[0]];
            tmpTetra_[2] = tetra[neg[1]];
            tmpTetra_[3] = intp[0];
            genTets.Push_back(tmpTetra_);
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

void TriaIntersect::ExportTetras(StdVector<CoordTetra> tetList,std::string baseFName){
  for(UInt i=0;i<tetList.GetSize();++i ){
    std::stringstream aFil;
    aFil << baseFName << i << ".off";
    std::ofstream aStream(aFil.str().c_str(), std::ios::trunc | std::ios::out);
    aStream << "OFF" << std::endl;
    aStream << "4 4 6" << std::endl;
    aStream << tetList[i][0][0] << "\t" << tetList[i][0][1] << "\t" << tetList[i][0][2] << std::endl;
    aStream << tetList[i][1][0] << "\t" << tetList[i][1][1] << "\t" << tetList[i][1][2] << std::endl;
    aStream << tetList[i][2][0] << "\t" << tetList[i][2][1] << "\t" << tetList[i][2][2] << std::endl;
    aStream << tetList[i][3][0] << "\t" << tetList[i][3][1] << "\t" << tetList[i][3][2] << std::endl;

    aStream << "3 0 1 2" << std::endl;
    aStream << "3 0 1 3" << std::endl;
    aStream << "3 0 2 3" << std::endl;
    aStream << "3 1 2 3" << std::endl;
    aStream.close();
  }

}

}
