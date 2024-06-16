// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SimInputVTKBased.cc
 *       \brief    <Description>
 *
 *       \date     28.02.2014
 *       \author   ahueppe
 */
//================================================================================================

#include "SimInputVTKBased.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/ElemMapping/Elem.hh"

//boost includes
#include "boost/regex.hpp"

//VTK Includes

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Winconsistent-missing-override"
#pragma clang diagnostic ignored "-Warray-parameter"

#include "vtkMultiBlockDataSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkMergePoints.h"
#include "vtkSMPMergePoints.h"
#include "vtkInformation.h"

#pragma clang diagnostic pop

#include <unordered_map>

// This is due to the fucking OLAS New operator!!!
#undef New

namespace CoupledField {

DEFINE_LOG(siminputVTK, "simInput.vtk")

//! This is a special class for the identification of unqie finite volume faces.
// It is intended to find the two cells adjacent to one face. Therefore it provides 
// a hash value und a == operator. 
// As long as slave and master element are equal, it is a boundary face
VTKFVFace::VTKFVFace(UInt masterElement, StdVector<UInt>& points) {
  masterElement_ = masterElement;
  slaveElement_ = masterElement;
  points_ = points;
  
  StdVector<UInt> sortedPoints = points;
  std::sort(sortedPoints.begin(), sortedPoints.end());

  UInt Isize = sortedPoints.GetSize();
  hash_ = 11;
  for (UInt i = 0; i < Isize; i++) {
    hash_ += 5 * sortedPoints[i];
  }
}

VTKFVFace::VTKFVFace() {
  masterElement_ = 0;
  slaveElement_ = 0;
  hash_ = 0;
}

UInt VTKFVFace::GetMasterElement() {
  return masterElement_;
}
  
UInt VTKFVFace::GetSlaveElement() {
  return slaveElement_;
}
  
void VTKFVFace::SetSlaveElement(UInt slaveElement) {
  slaveElement_ = slaveElement;
}

bool VTKFVFace::HasSlaveElement() {
 return masterElement_ != slaveElement_;
}
  
UInt VTKFVFace::GetPointCount() {
  return points_.GetSize();
}
  
void VTKFVFace::GetPoints(StdVector<UInt>& points) {
  points = points_;
}

bool VTKFVFace::HasPoint(UInt point) {
  for (UInt i = 0; i < points_.GetSize(); i++) {
    if (points_[i] == point) {
      return true;
    }
  }
  return false;
}

size_t VTKFVFace::GetHash() const {
  return hash_;
}

bool VTKFVFace::operator==(const VTKFVFace &other) const {
  return EqualTest(other) != 0;
}
  
//! This is for testing equality to other faces.
//! returns 0 if the faces are not equal
//! returns 1 if the points are equal and in same orientation
//! returns -1 if the points are equal and in reverse orientation
Integer VTKFVFace::EqualTest(const VTKFVFace &other) const {
  const UInt size = points_.GetSize();
  if (size != other.points_.GetSize()) {
    return 0;
  }
  if (size == 0) {
    return 0;
  }
  UInt firstPoint = points_[0];
  UInt firstOtherIndex = 0;
  bool notFound = true;
  for (UInt i = 0; i < size && notFound; i++) {
    if (other.points_[i] == firstPoint) {
      firstOtherIndex = i;
      notFound = false;
    }
  }
  if (notFound) {
    return 0;
  }
  if (size == 1) {
    return -1;
  }
  
  // forward check
  bool forwardEqual = true;
  UInt oI = firstOtherIndex;
  for (UInt i = 0; i < size && forwardEqual; i++) {
    forwardEqual = points_[i] == other.points_[oI];
    oI++;
    if (oI == size) {
      oI = 0;
    }
  }
  if (forwardEqual) {
    return 1;
  }
  // backward check
  bool backwardEqual = true;
  oI = firstOtherIndex;
  for (UInt i = 0; i < size && backwardEqual; i++) {
    backwardEqual = points_[i] == other.points_[oI];
    if (oI == 0) {
      oI = size - 1;
    } else {
      oI--;
    }
  }
  if (backwardEqual) {
    return -1;
  }
  return 0;
}

VTKFVFace::~VTKFVFace() {
}

void SimInputVTKBased::InitModule(){

  CreateReader();
  dim_ = GetDim();

  //make this variable...
  if(false)
    reader_->DebugOn();

  SetTimeValue(-1.0);
  reader_->Modified();
  reader_->Update();

  GetTimeValues();

  EnableRegions();

  SetTimeValue(-1);
  reader_->Modified();
  reader_->Update();



}

void SimInputVTKBased::ReadMesh( Grid *mi ){
  mi_ = mi;

  LOG_TRACE(siminputVTK) << "Generating VTK Grid";

  std::map<std::string, StdVector<UInt> > pRegMap;
  GetPointIDsByRegion(pRegMap);

  ReadElemData(pRegMap);

}

UInt SimInputVTKBased::ComputeNumElems(vtkCompositeDataIterator* iter){
  vtkCell* cell;
  UInt numElems = 0;
  vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());

  if(ugrid){
     vtkCellArray* cells = ugrid->GetCells();
     vtkIdType npts;
     vtkIdType *pts;
     vtkIdType cellId = 0;

     cells->InitTraversal();

     while(cells->GetNextCell(npts, pts)) {
       cell = ugrid->GetCell(cellId);
       Elem::FEType curT = this->VTKCellTypeToFEType(cell->GetCellType());
       if(curT == Elem::ET_POLYGON || curT == Elem::ET_POLYHEDRON){
         vtkIdList *   polyPtsIds = vtkIdList::New();
         vtkPoints *   polyPts = vtkPoints::New();
         UInt result = cell->Triangulate(0,polyPtsIds,polyPts);
         UInt eDim = 0;
         if(curT == Elem::ET_POLYGON){
           eDim = 2; 
         }else{
           eDim = 3; 
         }
         UInt numTets = polyPtsIds->GetNumberOfIds()/(eDim+1);
         if(!result && eDim == 3){
           EXCEPTION("Triangulation of polyhedral cell of dim " << eDim << " returned degenerated simplices. Return value " << result);
         }else if(!result){
           WARN("Triangulation of polygon cell returned degenerated simplices. As this may not be fatal, we continue..." );
         
         }

         numElems += numTets;
         cellId++;
       }else{
         numElems++;
         cellId++;
       }
     }
  }else{
    //polyhedra are not very common in structured grids
    vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    UInt numCells = ds->GetNumberOfCells();
    numElems += numCells;
  }
  return numElems;
}

void SimInputVTKBased::GetPointIDsByRegion(std::map<std::string, StdVector<UInt> >& pointMap){
  LOG_TRACE(siminputVTK) << "Enter SimInputVTKBased::GetPointIDsByRegion";
  //Determine number of Regions
  vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
  iter->GoToFirstItem();
  /* find out the number of partitions */
  //store region names and their dimension
  UInt numRegions  = 0;
  UInt numNodes = 0;
  UInt numElems = 0;
  std::string regionName;
  vtkDataSet* ds;
  Double pt[3] = {0.0,0.0,0.0};
  std::unordered_map<Point, UInt> coordNodeMap;
  
  //determine dimension of the Grid to distinguish between surface and volume regions
  Integer gDim = mi_->GetDim();
  while (! iter->IsDoneWithTraversal()){
    numRegions++;
    GetRegionName(regionName,iter,numRegions);

    ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    UInt curNumPoints = ds->GetNumberOfPoints();
    UInt curNumElems = ComputeNumElems(iter);
    numElems += curNumElems;
    if(curNumPoints == 0 && curNumElems == 0){
      iter->GoToNextItem();
      continue;
    }

    //determine surface or volume region
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());
    int regDim  = ugrid->GetCell(0)->GetCellDimension();
    if(regDim == gDim-1){
      RegionIdType curReg = mi_->AddRegion(regionName);
      regionAssoc_[curReg] = numRegions-1;
    }else if (regDim == gDim){
      RegionIdType curReg = mi_->AddRegion(regionName);
      regionAssoc_[curReg] = numRegions-1;
    }else{
      EXCEPTION("Got region with unexpected dimension " << regDim << " for a " << gDim << "d Grid.");
    }

    regionNamesOfDim_[regDim].Push_back(regionName);

    //mi_->AddNodes(curNumPoints);
    mi_->AddElems(curNumElems);

    pointMap[regionName].Reserve(curNumPoints);

    //make points unique if we need this
    //right now, every region consists of a disconnected set of points
    for (UInt j = 0; j < curNumPoints; ++j){
      ds->GetPoint(j, pt);
      Point ip(pt[0], pt[1], pt[2]);
      if (coordNodeMap.find(ip) == coordNodeMap.end()) {
        numNodes++;
        coordNodeMap[ip] = numNodes;
      }
      pointMap[regionName].Push_back(coordNodeMap[ip]);
      /**
      numNodes++;
      ds->GetPoint(j, pt);
      Vector<Double> aNodeC;
      aNodeC.Fill(pt,gDim);
      mi_->SetNodeCoordinate(numNodes,aNodeC);
      pointMap[regionName].Push_back(numNodes);
      **/

    }
    iter->GoToNextItem();
  }
  iter->Delete();
  
  // inserting node coordinates_
  mi_->AddNodes(numNodes);
  for (std::unordered_map<Point, UInt>::iterator it = coordNodeMap.begin(); it != coordNodeMap.end(); it++) {
    mi_->SetNodeCoordinate(it->second,it->first.data);
  }
  
  numNodes_ = numNodes;
  numRegions_ = numRegions;
  numElems_ = numElems;
  globElemLocElem_.reserve(numElems);

  LOG_DBG(siminputVTK) << "Number of regions in VTK grid: " << numRegions;
  LOG_DBG(siminputVTK) << "Number of nodes in VTK grid: " << numNodes;
  LOG_DBG(siminputVTK) << "Number of volume elements in VTK grid: " << mi_->GetNumVolElems();
  LOG_DBG(siminputVTK) << "Number of surface elements in VTK grid: " << mi_->GetNumSurfElems();
  LOG_TRACE(siminputVTK) << "Exit SimInputVTKBased::GetPointIDsByRegion";
}

void SimInputVTKBased::ReadElemData(std::map<std::string, StdVector<UInt> >& pointMap){
  vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
  iter->GoToFirstItem();
  UInt idx = 1;
  UInt elemIdx = 1;
  
  // FV stuff
  UInt numFVPoly = 0;
  UInt numFVHexa = 0;
  UInt numFVWedge = 0;
  UInt numFVPyra = 0;
  UInt numFVTetra = 0;
  Grid::FiniteVolumeRepresentation& fvrep = mi_->GetFiniteVolumeRepresentation();
  std::vector<bool>& isVolumeElement = fvrep.isVolumeElement;
  StdVector<UInt>& masterElement = fvrep.masterElement;
  if (readFVMesh_) {
    fvrep.isSet = true;
    fvrep.grid = mi_;
    isVolumeElement.clear();
    isVolumeElement.resize(numElems_ + 1, false);
    masterElement.Resize(numElems_ + 1);
    masterElement.Init(0);
  }
  // FV stuff end
  
  while (! iter->IsDoneWithTraversal()){
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());
    vtkCell* cell;

    UInt curNumPoints = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject())->GetNumberOfPoints();
    UInt curNumElems = ComputeNumElems(iter);
    if(curNumPoints == 0 && curNumElems == 0){
      iter->GoToNextItem();
      continue;
    }

    std::string curRegName;
    GetRegionName(curRegName,iter,idx);
    RegionIdType curReg = mi_->GetRegion().Parse(curRegName);
    StdVector<UInt>& regPointMap = pointMap[curRegName];

    if(ugrid){
      vtkCellArray* cells = ugrid->GetCells();
      vtkIdType npts;
      vtkIdType *pts;
      vtkIdType cellId = 0;

      cells->InitTraversal();

      while(cells->GetNextCell(npts, pts)) {
        cell = ugrid->GetCell(cellId);
        Elem::FEType curT = this->VTKCellTypeToFEType(cell->GetCellType());
        if(curT == Elem::ET_POLYGON || curT == Elem::ET_POLYHEDRON){
          vtkIdList *   polyPtsIds = vtkIdList::New();
          vtkPoints *   polyPts = vtkPoints::New();
          UInt result = cell->Triangulate(0,polyPtsIds,polyPts);
          UInt masterFVElement = elemIdx;  // FV stuff
          bool isVolume = curT == Elem::ET_POLYHEDRON;
          if(result || curT == Elem::ET_POLYGON){
            Elem::FEType newType;
            UInt eDim = 0;
            if(curT == Elem::ET_POLYGON){
              newType = Elem::ET_TRIA3;
              eDim = 2; 
            }else{
              newType = Elem::ET_TET4;
              eDim = 3; 
            }

            UInt numTets = polyPtsIds->GetNumberOfIds()/(eDim+1);
            UInt * con = new UInt[eDim+1];
            for(UInt aTet = 0; aTet < numTets; aTet++){
              for(UInt node = 0; node < eDim+1; node++){
                UInt locId = polyPtsIds->GetId(aTet*(eDim+1)+node);
                con[node] = regPointMap[locId];
              }
              mi_->SetElemData(elemIdx,newType,curReg,con);
              globElemLocElem_[elemIdx] = cellId;
              if (readFVMesh_) { // FV stuff
                masterElement[elemIdx] = masterFVElement;
                isVolumeElement[elemIdx] = isVolume;
              } // FV stuff end
              elemIdx++;
            }
            delete[] con;
          }else{
            EXCEPTION("Triangulation of polygon/polyhedral cell return degenerate simplices.")
          }
          polyPtsIds->Delete();
          polyPts->Delete();
          if (readFVMesh_ && curT == Elem::ET_POLYHEDRON) { // FV stuff
            numFVPoly++;
            vtkIdList* faceStream = vtkIdList::New();
            ugrid->GetFaceStream(cellId, faceStream);
            vtkIdType idx = 0;
            UInt numFaces = faceStream->GetId(idx);
            idx++;
            StdVector<UInt> facePoints;
            for (UInt iFace = 0; iFace < numFaces; iFace++) {
              UInt numFacePoints = faceStream->GetId(idx);
              facePoints.Resize(numFacePoints);
              idx++;
              for (UInt iPoint = 0; iPoint < numFacePoints; iPoint++) {
                facePoints[iPoint] = regPointMap[faceStream->GetId(idx)];
                idx++;
              }
              AddFVFace(masterFVElement, facePoints);
            }
            faceStream->Delete();
          } // FV stuff end
        }else{
          UInt * con = new UInt[npts];
          for (UInt k = 0; k < npts; ++k){
            UInt locId = pts[k];
            con[k] = regPointMap[locId];
          }
          mi_->SetElemData(elemIdx,curT,curReg,con);
          globElemLocElem_[elemIdx] = cellId;
          if (readFVMesh_) { // FV stuff
            masterElement[elemIdx] = elemIdx;
            isVolumeElement[elemIdx] = true;
            if (curT == Elem::ET_HEXA8) {
              AddFVQuadFace(elemIdx, con[0], con[4], con[7], con[3]);
              AddFVQuadFace(elemIdx, con[1], con[2], con[6], con[5]);
              AddFVQuadFace(elemIdx, con[0], con[1], con[5], con[4]);
              AddFVQuadFace(elemIdx, con[3], con[7], con[6], con[2]);
              AddFVQuadFace(elemIdx, con[0], con[3], con[2], con[1]);
              AddFVQuadFace(elemIdx, con[4], con[5], con[6], con[7]);
              numFVHexa++;
            } else if (curT == Elem::ET_WEDGE6) {
              AddFVTriFace(elemIdx, con[0], con[1], con[2]);
              AddFVTriFace(elemIdx, con[3], con[5], con[4]);
              AddFVQuadFace(elemIdx, con[0], con[3], con[4], con[1]);
              AddFVQuadFace(elemIdx, con[1], con[4], con[5], con[2]);
              AddFVQuadFace(elemIdx, con[2], con[5], con[3], con[0]);
              numFVWedge++;
            } else if (curT == Elem::ET_PYRA5) {
              AddFVTriFace(elemIdx, con[1], con[4], con[0]);
              AddFVTriFace(elemIdx, con[2], con[4], con[1]);
              AddFVTriFace(elemIdx, con[3], con[4], con[2]);
              AddFVTriFace(elemIdx, con[0], con[4], con[3]);
              AddFVQuadFace(elemIdx, con[3], con[2], con[1], con[0]);
              numFVPyra++;
            } else if (curT == Elem::ET_TET4) {
              if (elemIdx == 13517665 || elemIdx == 13519266 || elemIdx == 13519271) {
                std::cout << " real tetra " << elemIdx << std::endl;
              }
              AddFVTriFace(elemIdx, con[0], con[1], con[3]);
              AddFVTriFace(elemIdx, con[1], con[2], con[3]);
              AddFVTriFace(elemIdx, con[2], con[0], con[3]);
              AddFVTriFace(elemIdx, con[0], con[2], con[1]);
              numFVTetra++;
            }
          } // FV stuff end
          elemIdx++;
          delete[] con;
        }
        cellId++;
      }

    }else{
      // We have structured or even uniform grids
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      UInt numCells = ds->GetNumberOfCells();
      for (UInt j = 0; j < numCells; ++j) {
        cell = ds->GetCell(j);
        UInt numPoints = cell->GetNumberOfPoints();
        Elem::FEType curT = VTKCellTypeToFEType(cell->GetCellType());

        switch(ds->GetDataObjectType()) {
          case VTK_UNIFORM_GRID: // 10
          case VTK_RECTILINEAR_GRID: // 3
          case VTK_IMAGE_DATA: // 6
            if (numPoints == 4) {
              curT = Elem::ET_QUAD4;
            } else {
              curT = Elem::ET_HEXA8;
            }
            break;
          case VTK_STRUCTURED_GRID: // 2
            break;
        }

        UInt * con = new UInt[numPoints];
        for (UInt k = 0; k < numPoints; ++k){
          UInt locId = cell->GetPointId(k);
          con[k] = regPointMap[locId];
        }
        mi_->SetElemData(elemIdx,curT,curReg,con);
        globElemLocElem_[elemIdx] = j;
        if (readFVMesh_ && curT == Elem::ET_HEXA8) { // FV stuff
          masterElement[elemIdx] = elemIdx;
          isVolumeElement[elemIdx] = true;
          AddFVQuadFace(elemIdx, con[0], con[4], con[7], con[3]);
          AddFVQuadFace(elemIdx, con[1], con[2], con[6], con[5]);
          AddFVQuadFace(elemIdx, con[0], con[1], con[5], con[4]);
          AddFVQuadFace(elemIdx, con[3], con[7], con[6], con[2]);
          AddFVQuadFace(elemIdx, con[0], con[3], con[2], con[1]);
          AddFVQuadFace(elemIdx, con[4], con[5], con[6], con[7]);
          numFVHexa++;
        } // FV stuff end
        elemIdx++;
      }
    }
    idx++;
    iter->GoToNextItem();
  }
  iter->Delete();
  numElems_ = elemIdx-1;
  
  if (readFVMesh_) { // FV stuff
    std::cout << " Finite volume tested faces: " << numFVTestedFaces_ << std::endl;
    std::cout << " Finite volume internal faces " << numFVInternalFaces_ << std::endl;
    std::cout << " Finite volume boundary faces " << numFVBoundaryFaces_ << std::endl;
    std::cout << " Finite volume polyhedrons " << numFVPoly << std::endl;
    std::cout << " Finite volume hexahedrons " << numFVHexa << std::endl;
    std::cout << " Finite volume wedges " << numFVWedge << std::endl;
    std::cout << " Finite volume pyramids " << numFVPyra << std::endl;
    std::cout << " Finite volume tetrahedrons " << numFVTetra << std::endl;
    UInt sumFV = numFVPoly + numFVHexa + numFVWedge + numFVPyra + numFVTetra;
    std::cout << " Finite volume cells " <<  sumFV << std::endl;

    UInt numBoundaryFaces = numFVBoundaryFaces_;
    UInt numFaces = numFVBoundaryFaces_ + numFVInternalFaces_;
    
    if (fixFVPyramids_ && numFVPyra > 0) {
      // creating a vector of boundary faces
      StdVector<VTKFVFace> allFaces(numBoundaryFaces);
      UInt bIdx = 0;
      for (std::unordered_set<VTKFVFace>::iterator itr = FVFaces_.begin(); itr != FVFaces_.end(); ++itr) {
        VTKFVFace face = *itr;
        if (!face.HasSlaveElement()) {
          allFaces[bIdx] = face;
          bIdx++;
        }
      }
  
      // clearing and saving FVFaces
      std::unordered_set<VTKFVFace> FVFaces;
      FVFaces.swap(FVFaces_);
      FVFaces_.clear();
      std::unordered_set<VTKFVFace>().swap(FVFaces_);
      
      // adding edges as faces
      numFVTestedFaces_ = 0;
      numFVInternalFaces_ = 0;
      numFVBoundaryFaces_ = 0;
      StdVector<UInt> fPoints;
      for (UInt iFace = 0; iFace < numBoundaryFaces; iFace++) {
        allFaces[iFace].GetPoints(fPoints);
        UInt nPoints = fPoints.GetSize();
        UInt pointA = fPoints[nPoints - 1];
        UInt pointB;
        for (UInt ip = 0; ip < nPoints; ip++) {
          pointB = fPoints[ip];
          AddFVEdge(iFace, pointA, pointB);
          pointA = pointB;
        }
      }
      std::cout << " Finite volume tested boundary edges " << numFVTestedFaces_ << std::endl;
      std::cout << " Finite volume internal boundary edges " << numFVInternalFaces_ << std::endl;
      std::cout << " Finite volume boundary boundary edges " << numFVBoundaryFaces_ << std::endl;

      // copying  edges into a vector
      StdVector<VTKFVFace> allEdges;
      UInt edgeCount = 0;
      for (std::unordered_set<VTKFVFace>::iterator itr = FVFaces_.begin(); itr != FVFaces_.end(); ++itr) {
        VTKFVFace face = *itr;
        allEdges.Push_back(face);
        edgeCount++;
      }
      FVFaces_.clear();
      std::unordered_set<VTKFVFace>().swap(FVFaces_);

      // creating boundary face connectivity
      StdVector<UInt> bFaceFaceIndex(numBoundaryFaces + 1);
      bFaceFaceIndex.Init(0);
      for (UInt iE = 0; iE < edgeCount; iE++) {
        bFaceFaceIndex[allEdges[iE].GetMasterElement()]++;
        bFaceFaceIndex[allEdges[iE].GetSlaveElement()]++;
      }
      UInt iIdx = 0;
      for (UInt ibf = 0; ibf < numBoundaryFaces; ibf++) {
        UInt idxP = bFaceFaceIndex[ibf];
        bFaceFaceIndex[ibf] = iIdx;
        iIdx += idxP;
      }
      bFaceFaceIndex[numBoundaryFaces] = iIdx;
      const UInt dummyIdx = numBoundaryFaces + 2;
      StdVector<UInt> bFaceFace(iIdx);
      bFaceFace.Init(dummyIdx);
      for (UInt iE = 0; iE < edgeCount; iE++) {
        UInt masterCell = allEdges[iE].GetMasterElement();
        UInt slaveCell = allEdges[iE].GetSlaveElement();
        UInt* bFaces = &bFaceFace[bFaceFaceIndex[masterCell]];
        UInt iSF;
        for (iSF = 0; bFaces[iSF] != dummyIdx; iSF++) {}
        //bFaces[iSF] = iE;
        bFaces[iSF] = slaveCell;
        if (allEdges[iE].HasSlaveElement()) {
          
          bFaces = &bFaceFace[bFaceFaceIndex[slaveCell]];
          for (iSF = 0; bFaces[iSF] != dummyIdx; iSF++) {}
          //bFaces[iSF] = iE;
          bFaces[iSF] = masterCell;
        }
      }
       
       
      // sorting boundary faces into adjacent groups
      UInt markedBFaces = 0;
      StdVector<UInt> markings(numBoundaryFaces);
      markings.Init(0);
      std::vector<bool> checkings(numBoundaryFaces,false);
      UInt marker = 1;
      std::vector<UInt> facesPerMarker(2,0);
      while (markedBFaces < numBoundaryFaces) {
        UInt startIdx = 0;
        while (markings[startIdx] != 0) {
          startIdx++;
        }
        
        markings[startIdx] = marker;
        markedBFaces++;
        facesPerMarker[marker]++;
        std::unordered_set<UInt> toCheckSet;
        toCheckSet.insert(startIdx);
        while (toCheckSet.size() > 0) {
          UInt checkFaceIdx = *toCheckSet.begin();
          toCheckSet.erase(checkFaceIdx);
          UInt maxIndex = bFaceFaceIndex[checkFaceIdx + 1];
          for (UInt nIndex = bFaceFaceIndex[checkFaceIdx]; nIndex < maxIndex; nIndex++) {
            UInt nFaceIndex = bFaceFace[nIndex];
            if (markings[nFaceIndex] == 0) {
              markings[nFaceIndex] = marker;
              markedBFaces++;
              facesPerMarker[marker]++;
              toCheckSet.insert(nFaceIndex);
            }
          }
          
        }
        marker++;
        facesPerMarker.push_back(0);
      }
      marker--;
      UInt fixedPyramids = 0;
      //std::cout << " Number of markers " << marker << std::endl;
      for (UInt iM = 1; iM <= marker; iM++) {
        //std::cout << "    marker " << iM << ": " << facesPerMarker[iM] << " bFaces" << std::endl;
        if (facesPerMarker[iM] == 4) { // case for patcheble pyramid
          bool isOk = true;
          StdVector<VTKFVFace> bFaces(4);
          StdVector<UInt> bFaceMasterElements(4);
          std::vector<bool> bFaceIsPyra(4,false);
          UInt bwIdx = 0;
          for (UInt ibf = 0; ibf < numBoundaryFaces && isOk; ibf++) {
            if (markings[ibf] == iM) {
              if (allFaces[ibf].GetPointCount() == 3) {
                bFaces[bwIdx] = allFaces[ibf];
                bFaceMasterElements[bwIdx] = bFaces[bwIdx].GetMasterElement();
                const Elem* feElem = mi_->GetElem(bFaceMasterElements[bwIdx]);
                bFaceIsPyra[bwIdx] = feElem->type == Elem::ET_PYRA5;
                bwIdx++;
              } else {
                isOk = false;
              }
            }
          }
          if (!isOk) {
            continue;
          }
          // searching for master pyramid
          StdVector<UInt> appears(4);
          appears.Init(0);
          for (UInt i = 0; i < 4; i++) {
            if (bFaceIsPyra[i]) {
              for (UInt j = 0; j < 4; j++) {
                if (bFaceMasterElements[j] == bFaceMasterElements[i]) {
                  appears[j]++;
                }
              }            
            }
          }
          UInt mIdxA = 6;
          UInt mIdxB = 6;
          for (UInt i = 0; i < 4; i++) {
            if (appears[i] == 2) {
              if (mIdxA == 6) {
                mIdxA = i;
              } else {
                mIdxB = i;
              }
            }
          }
          if (mIdxA == 6 || mIdxB == 6) {
            continue;
          }
          const Elem* mPyra = mi_->GetElem(bFaceMasterElements[mIdxA]);
          
          // searching neighbour face connected to top node
          UInt topNode = mPyra->connect[4];
          UInt tIdx = 6;
          for (UInt i = 0; i < 4 && tIdx == 6; i++) {
            if (i != mIdxA && i != mIdxB) {
              if (bFaces[i].HasPoint(topNode)) {
                tIdx = i;
              }
            }
          }
          if (tIdx == 6) {
            continue;
          }
          // getting left of connect neighbour
          UInt bIdx = 6;
          for (UInt i = 0; i < 4 && bIdx == 6; i++) {
            if (i != mIdxA && i != mIdxB && i != tIdx) {
              bIdx = i;
            }
          }
          // getting bottom face of pyramid
          StdVector<UInt> points(4);
          points[0] = mPyra->connect[3];          
          points[1] = mPyra->connect[2];          
          points[2] = mPyra->connect[1];          
          points[3] = mPyra->connect[0];          
          VTKFVFace botFace(bFaceMasterElements[mIdxA], points);
          bool hasBotElement = false;
          UInt botElement = bFaceMasterElements[mIdxA];
          std::unordered_set<VTKFVFace>::const_iterator got = FVFaces.find(botFace);
          if (got == FVFaces.end()) {
            continue;
          }
          VTKFVFace initBotFace = *got;
          if (initBotFace.GetMasterElement() != bFaceMasterElements[mIdxA]) {
            hasBotElement = true;
            botElement = initBotFace.GetMasterElement();
          }
          if (initBotFace.GetSlaveElement() != bFaceMasterElements[mIdxA]) {
            hasBotElement = true;
            botElement = initBotFace.GetSlaveElement();
          }
          // evaluate replacement for bottom face
          bool found = false;
          UInt fNiPIdx = 0;
          for (UInt i = 0; i < 4 && !found; i++) {
            if (bFaces[bIdx].HasPoint(points[i]) && !bFaces[tIdx].HasPoint(points[i])) {
              found = true;
              fNiPIdx = i;
            }
          }
          if (!found) {
            continue;
          }
          StdVector<UInt> npoints(3);
          UInt addForNewBot = 0;
          for (UInt i = 0; i < 3; i++) {
            if (i == fNiPIdx) {
              addForNewBot++;
            }
            npoints[i] = points[i + addForNewBot];
          }
          VTKFVFace newBotFace(bFaceMasterElements[mIdxA], npoints);
          if (hasBotElement) {
            newBotFace.SetSlaveElement(botElement);
          }
          
          // erasing old faces
          for (UInt i = 0; i < 4; i++) {
            FVFaces.erase(bFaces[i]);
          }
          FVFaces.erase(botFace);
          
          // inserting face replacements;
          FVFaces.insert(newBotFace);
          bFaces[tIdx].SetSlaveElement(bFaceMasterElements[mIdxA]);
          FVFaces.insert(bFaces[tIdx]);
          if (hasBotElement) {
            bFaces[bIdx].SetSlaveElement(botElement);
          }
          FVFaces.insert(bFaces[bIdx]);
          fixedPyramids++;
          numFaces -= 2;
          numBoundaryFaces -= 4;
          
          // inserting element replacement
          RegionIdType regionId = mPyra->regionId;
          StdVector<UInt> tetCon = npoints;
          tetCon.Push_back(topNode);
          mi_->SetElemData(bFaceMasterElements[mIdxA], Elem::ET_TET4, regionId, &tetCon[0]);
          //std::cout << " patched one pyra " << std::endl;
        }
        
      }
      std::cout << " Finite volume fixed " << fixedPyramids << " pyramids" << std::endl;
      
      // restoring FVFaces
      FVFaces_.clear();
      std::unordered_set<VTKFVFace>().swap(FVFaces_);
      FVFaces.swap(FVFaces_);
      
      std::cout << " Final finite volume internal faces " << (numFaces - numBoundaryFaces) << std::endl;
      std::cout << " Final finite volume boundary faces " << numBoundaryFaces << std::endl;
      std::cout << " Final finite volume faces " << numFaces << std::endl;
    }
      
    // copy faces to finite volume representation
    fvrep.isFaceBoundary.clear();
    fvrep.isFaceBoundary.resize(numFaces,false);
    fvrep.faceMasterElement.Resize(numFaces);
    fvrep.faceSlaveElement.Resize(numFaces);
    fvrep.facePointIndex.Resize(numFaces + 1);
    fvrep.faceCount = numFaces;
    
    // copy main face interformation
    UInt fIdx = 0;
    UInt pointIdx = 0;
    for (std::unordered_set<VTKFVFace>::iterator itr = FVFaces_.begin(); itr != FVFaces_.end(); ++itr) {
      VTKFVFace face = *itr;
      fvrep.isFaceBoundary[fIdx] = !face.HasSlaveElement();
      fvrep.faceMasterElement[fIdx] = face.GetMasterElement();
      fvrep.faceSlaveElement[fIdx] = face.GetSlaveElement();
      fvrep.facePointIndex[fIdx] = pointIdx;
      pointIdx += face.GetPointCount();
      fIdx++;
    }
    fvrep.facePointIndex[fIdx] = pointIdx;
    
    // copy face points
    fvrep.facePoint.Resize(pointIdx);
    pointIdx = 0;
    StdVector<UInt> facePoints;
    for (std::unordered_set<VTKFVFace>::iterator itr = FVFaces_.begin(); itr != FVFaces_.end(); ++itr) {
      VTKFVFace face = *itr;
      face.GetPoints(facePoints);
      for(UInt ip = 0; ip < facePoints.GetSize(); ip++) {
        fvrep.facePoint[pointIdx] = facePoints[ip];
        pointIdx++;
      }
    }
  } // FV stuff end
}

void SimInputVTKBased::AddFVEdge(UInt masterElement, UInt pA, UInt pB) { // FV stuff
  StdVector<UInt> points(2);
  points[0] = pA;
  points[1] = pB;
  AddFVFace(masterElement, points);
} // FV stuff end

void SimInputVTKBased::AddFVTriFace(UInt masterElement, UInt pA, UInt pB, UInt pC) { // FV stuff
  StdVector<UInt> points(3);
  points[0] = pA;
  points[1] = pB;
  points[2] = pC;
  AddFVFace(masterElement, points);
} // FV stuff end

void SimInputVTKBased::AddFVQuadFace(UInt masterElement, UInt pA, UInt pB, UInt pC, UInt pD) { // FV stuff
  StdVector<UInt> points(4);
  points[0] = pA;
  points[1] = pB;
  points[2] = pC;
  points[3] = pD;
  AddFVFace(masterElement, points);
} // FV stuff end

void SimInputVTKBased::AddFVFace(UInt masterElement, StdVector<UInt>& points) { // FV stuff
  numFVTestedFaces_++;
  VTKFVFace fvFace(masterElement, points);
  std::unordered_set<VTKFVFace>::const_iterator got = FVFaces_.find(fvFace);
  if (got == FVFaces_.end()) {
    FVFaces_.insert(fvFace);
    numFVBoundaryFaces_++;
  } else {
    VTKFVFace fvFaceFound = *got;
    if (fvFaceFound.EqualTest(fvFace) == 1 && points.GetSize() > 2) {
      WARN("Equal finite volume faces are forward equal");
    }
    FVFaces_.erase(fvFaceFound);
    if (fvFaceFound.HasSlaveElement()) {
      WARN("Internal finite volume faces appeared more than two times");
    }
    fvFaceFound.SetSlaveElement(masterElement);
    FVFaces_.insert(fvFaceFound);
    numFVInternalFaces_++;
    numFVBoundaryFaces_--;
  }
} // FV stuff end


UInt SimInputVTKBased::GetDim(){
  Integer dim = 0;
  vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
  iter->GoToFirstItem();
  while (! iter->IsDoneWithTraversal()){
    //determine the dimension of the region by
    //checking the first element works only if the regions
    //are defined consistently
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());
    if(ugrid->GetNumberOfCells() != 0){
      vtkCell* cell = ugrid->GetCell(1);
      dim = std::max(cell->GetCellDimension(),dim);
    }
    iter->GoToNextItem();
  }
  if(dim<0){
    EXCEPTION("Detected grid dimension of " << dim << ". Check your input Grid");
  }
  iter->Delete();
  return (UInt)dim;
}

void SimInputVTKBased::GetRegionName(std::string& name, vtkCompositeDataIterator* iter, UInt idx){
  std::stringstream tmp;
  //first determine the region name
  bool regNameRead = false;
  if(iter->HasCurrentMetaData() ){
    vtkInformation* info = iter->GetCurrentMetaData();
    if(info->Has(vtkCompositeDataSet::NAME())){
      name = iter->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME());
      regNameRead = true;
    }
  }
  if(!regNameRead){
    tmp << "Region_" << idx;
    name = tmp.str();
  }else{
    const boost::regex datExp("[,;/ ]");
    const std::string name_format("_");

    tmp << boost::regex_replace(name, datExp, name_format, boost::match_default | boost::format_sed);
    name = tmp.str();
  }
  bool removed = true;
  while (removed) {
    removed = false;
    if (name.size() > 0) {
      if (name.back() == '_') {
        name.pop_back();
        removed = true;
      }
    }
  }
}

Elem::FEType SimInputVTKBased::VTKCellTypeToFEType(UInt cellType){
    switch(cellType)
    {
    case VTK_VERTEX:
      return Elem::ET_POINT;
    case VTK_LINE:
      return Elem::ET_LINE2;
    case VTK_TRIANGLE:
      return Elem::ET_TRIA3;
    case VTK_QUAD:
      return Elem::ET_QUAD4;
    case VTK_TETRA:
      return Elem::ET_TET4;
    case VTK_HEXAHEDRON:
      return Elem::ET_HEXA8;
    case VTK_WEDGE:
      return Elem::ET_WEDGE6;
    case VTK_PYRAMID:
      return Elem::ET_PYRA5;
    case VTK_POLYGON:
      return Elem::ET_POLYGON;
    case VTK_POLYHEDRON:
      return Elem::ET_POLYHEDRON;

    // Quadratic, isoparametric cells
    case VTK_QUADRATIC_EDGE:
      return Elem::ET_LINE3;
    case VTK_QUADRATIC_TRIANGLE:
      return Elem::ET_TRIA6;
    case VTK_QUADRATIC_QUAD:
      return Elem::ET_QUAD8;
    case VTK_QUADRATIC_TETRA:
      return Elem::ET_TET10;
    case VTK_QUADRATIC_HEXAHEDRON:
      return Elem::ET_HEXA20;
    case VTK_QUADRATIC_WEDGE:
      return Elem::ET_WEDGE15;
    case VTK_QUADRATIC_PYRAMID:
      return Elem::ET_PYRA13;
    case VTK_BIQUADRATIC_QUAD:
      return Elem::ET_QUAD9;
    case VTK_TRIQUADRATIC_HEXAHEDRON:
      return Elem::ET_HEXA27;
    case VTK_QUADRATIC_LINEAR_QUAD:
      return Elem::ET_QUAD4;
    case VTK_QUADRATIC_LINEAR_WEDGE:
      return Elem::ET_WEDGE6;
    case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
      return Elem::ET_WEDGE15;
    case VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON:
      return Elem::ET_HEXA20;

    default:
      return Elem::ET_UNDEF;
    }
  }

}

