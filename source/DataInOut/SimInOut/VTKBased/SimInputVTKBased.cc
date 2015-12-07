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

#include "vtkMultiBlockDataSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkMergePoints.h"
#include "vtkSMPMergePoints.h"
#include "vtkInformation.h"

// This is due to the fucking OLAS New operator!!!
#undef New

namespace CoupledField {

DECLARE_LOG(siminputVTK)
DEFINE_LOG(siminputVTK, "simInput.vtk")


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
  Double pt[3];

  //determine dimension of the Grid to distinguish between surface and volume regions
  Integer gDim = mi_->GetDim();
  while (! iter->IsDoneWithTraversal()){
    numRegions++;
    GetRegionName(regionName,iter,numRegions);

    ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    UInt curNumPoints = ds->GetNumberOfPoints();
    UInt curNumElems = ComputeNumElems(iter);
    numElems += curNumElems;
    if(curNumPoints == 0 && curNumElems == 0)
      continue;

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

    mi_->AddNodes(curNumPoints);
    mi_->AddElems(curNumElems);

    pointMap[regionName].Reserve(curNumPoints);

    //make points unique if we need this
    //right now, every region consists of a disconnected set of points
    for (UInt j = 0; j < curNumPoints; ++j){
      numNodes++;
      ds->GetPoint(j, pt);
      Vector<Double> aNodeC;
      aNodeC.Fill(pt,gDim);
      mi_->SetNodeCoordinate(numNodes,aNodeC);
      pointMap[regionName].Push_back(numNodes);
    }
    iter->GoToNextItem();
  }
  iter->Delete();
  numNodes_ = numNodes;
  numRegions_ = numRegions;
  numElems_ = numElems;
  globElemLocElem_.Reserve(numElems);

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
  while (! iter->IsDoneWithTraversal()){
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());
    vtkCell* cell;

    std::string curRegName;
    GetRegionName(curRegName,iter,idx);
    RegionIdType curReg = mi_->GetRegion().Parse(curRegName);

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
                con[node] = pointMap[curRegName][locId];
              }
              mi_->SetElemData(elemIdx,newType,curReg,con);
              globElemLocElem_.Push_back(cellId);
              elemIdx++;
            }
            delete[] con;
          }else{
            EXCEPTION("Triangulation of polygon/polyhedral cell return degenerate simplices.")
          }
          polyPtsIds->Delete();
          polyPts->Delete();
        }else{
          UInt * con = new UInt[npts];
          for (UInt k = 0; k < npts; ++k){
            UInt locId = pts[k];
            con[k] = pointMap[curRegName][locId];
          }
          mi_->SetElemData(elemIdx,curT,curReg,con);
          globElemLocElem_.Push_back(cellId);
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
          con[k] = pointMap[curRegName][locId];
        }
        mi_->SetElemData(elemIdx,curT,curReg,con);
        globElemLocElem_.Push_back(j);
        elemIdx++;
      }
    }
    idx++;
    iter->GoToNextItem();
  }
  iter->Delete();
  numElems_ = elemIdx-1;
}

UInt SimInputVTKBased::GetDim(){
  Integer dim = 0;
  vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
  iter->GoToFirstItem();
  while (! iter->IsDoneWithTraversal()){
    //determine the dimension of the region by
    //checking the first element works only if the regions
    //are defined consistently
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());
    vtkCell* cell = ugrid->GetCell(1);
    dim = std::max(cell->GetCellDimension(),dim);
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

