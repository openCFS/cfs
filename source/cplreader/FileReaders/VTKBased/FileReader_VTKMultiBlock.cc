// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <iomanip>
#include <sstream>
#include <set>
#include "sys/stat.h"

#include "boost/regex.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
namespace fs=boost::filesystem;

#include "Domain/resultInfo.hh"

#include "cplreader/Settings.hh"
#include "FileReader_VTKMultiBlock.hh"


// This is due to the fucking OLAS New operator!!!
#undef New

#include "vtkMultiBlockDataSetAlgorithm.h"

#include "vtkDataArrayCollection.h"
#include "vtkDataArrayCollectionIterator.h"
#include "vtkInformation.h"

#include "vtkCellArray.h"
#include "vtkCellType.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMergePoints.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"

#include "vtkMultiBlockDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkCellDataToPointData.h"
#include "vtkPointSet.h"
#include "vtkCell.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkCellDataToPointData.h"
#include "vtkArrayIteratorTemplate.h"

namespace CoupledField
{

  FileReader_VTKMultiBlock::FileReader_VTKMultiBlock(const std::string& name,
                                                     const UInt dim,
                                                     const UInt numFiles) :
    FileReader(name, dim, numFiles, 1),
    reader_(NULL),
    pts_(NULL),
    numElems_(0)
  {
  }

  FileReader_VTKMultiBlock::~FileReader_VTKMultiBlock()
  {
    pts_->Delete();
    reader_->Delete();
  }

  void FileReader_VTKMultiBlock::Init()
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_VTKMultiBlock::Init" << std::endl;
    std::stringstream sstr;
    UInt numPoints;
    UInt numElems;

    // Let subclasses generate our reader
    CreateReader();

    if(settings.GetInt("verbose"))
      reader_->DebugOn();

    SetTimeValue(-1.0);
    reader_->Modified();
    reader_->Update();

    if (settings.GetInt("verbose"))
    {
        reader_->PrintSelf(std::cout, (vtkIndent)0);
    }

    GetTimeValues();
    
    UInt numTSFromFile = 0;
    if(timeValues_.size()) 
    {
      numTSFromFile = timeValues_.size();
    } else 
    {
      numTSFromFile = 1;
      timeValues_.push_back( 0 );
    }
    
    if (settings.GetInt("numsteps"))
    {
      numSteps_ = (UInt) settings.GetInt("numsteps");
      /* don not exceed the maximal number of timesteps possible */
      if (numSteps_ >= numTSFromFile)
      {
          numSteps_ = numTSFromFile;
      }
    } else {
      numSteps_ = numTSFromFile;
    }

    EnableRegions();

    SetTimeValue(-1);
    reader_->Modified();
    reader_->Update();

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    /* find out the number of partitions */
    numRegions_ = 0;
    while (! iter->IsDoneWithTraversal())
    {
      ++numRegions_;
      iter->GoToNextItem();
      //      std::cout << "region: " << numRegions_ << std::endl;      
    }

    if(settings.GetInt("verbose"))
    {
      std::cout << " Name: " << name_ << std::endl
                << " Dim: " << dim_ << std::endl
                << " numfiles: " << numSteps_ << std::endl
                << " numPartitions: " << numRegions_ << std::endl
                << " timeStep: " << (timeValues_[1]-timeValues_[0]) << std::endl;
    }

    numNodesPerRegion_.resize(numRegions_);
    numElemsPerRegion_.resize(numRegions_);

    pts_ = vtkUnstructuredGrid::New();
    pts_->Initialize();

    UInt p_cnt = 0;
    iter->GoToFirstItem();
    vtkDataSet* ds;

    // Compute bounding box  of whole multi block dataset in  order to be able
    // to partition this 3D space for merging points.
    Double bounds[6];
    std::fill(bounds, bounds+6, 0.0);
    while (! iter->IsDoneWithTraversal())
    {
      ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      Double* bnds = ds->GetBounds();
      
      for(UInt i=0; i<3; i++) 
      {
        UInt idx = i*2;
        bounds[idx+0] = bnds[idx+0] < bounds[idx+0] ? bnds[idx+0] : bounds[idx+0];
        bounds[idx+1] = bnds[idx+1] > bounds[idx+1] ? bnds[idx+1] : bounds[idx+1];
      }
      iter->GoToNextItem();      
    }
    
    // Fill pts_ dataset with points of first block.
    iter->GoToFirstItem();
    vtkDataSet* ds1;
    ds1 = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());

    UInt n=ds1->GetNumberOfPoints();
    vtkPoints* ps = vtkPoints::New();
    pts_->SetPoints(ps);
    ps->SetNumberOfPoints(n);    

    for(UInt i=0; i<n; i++) 
    {
      ps->InsertPoint(i, ds1->GetPoint(i));
    }
    ps->Delete();

    //    std::cout << "Building vtkMergePoints locator..." << std::endl;
    vtkMergePoints* merger = vtkMergePoints::New();
    merger->SetDataSet(pts_);
    //    merger->BuildLocator();
    merger->InitPointInsertion(pts_->GetPoints(), bounds);
    Double pt[3];
    //    std::cout << "Finished building vtkMergePoints locator..." << std::endl;

    nodeMap_.resize(numRegions_);
    
    while (! iter->IsDoneWithTraversal())
    {
      ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      numPoints = ds->GetNumberOfPoints();
      numNodesPerRegion_[p_cnt] = numPoints;
      /* cell values will be interpolated to element values */
      numElems = ds->GetNumberOfCells();
      numElemsPerRegion_[p_cnt] = numElems;
      numElems_ += numElems;

      /* Find out the maximum number of element nodes in region */
      numPoints = ds->GetMaxCellSize();
      maxNumElemNodes_ = numPoints > maxNumElemNodes_ ? numPoints : maxNumElemNodes_;
      
      if(settings.GetInt("verbose"))
      {
        std::cout << "Partition " << (p_cnt+1)
                  << " nodes: " << numNodesPerRegion_[p_cnt]
                  << " elems: " << numElemsPerRegion_[p_cnt]
                  << std::endl;
      }

      // Map all points to nodes in the inner mesh (partitionIdx=0 -> ds1)
      for (UInt j = 0, n= numNodesPerRegion_[p_cnt]; j < n; ++j)
      {
        ds->GetPoint(j, pt);
          
        vtkIdType ptId;
        UInt unique;
        unique = merger->InsertUniquePoint(pt, ptId);

        nodeMap_[p_cnt][j] = ptId+1;

        if(settings.GetInt("verbose"))
        {
          std::cout << "unique: " << unique << ", ptId: " << std::setw(10) << ptId
                    << ", nodeMap_[" << p_cnt << "][" << j << "]: " << nodeMap_[p_cnt][j] << std::endl;
        }
        
      }
      

      ++p_cnt;
      iter->GoToNextItem();
    }
    merger->Delete();
    iter->Delete();

    InitElemNodeMapping();
    
    std::cout << "Exiting FileReader_VTKMultiBlock::Init" << std::endl;
  }


  void FileReader_VTKMultiBlock::ReadNodalCoords(std::vector<Double>& NODECOORD)
  {
    std::cout << "Entering FileReader_VTKMultiBlock::ReadNodalCoords" << std::endl;
    /* number of coordinates. One Tupel of 3D coordinate (x,y,z) */

    const UInt numCoords = 3;

    /* just read coords for internal mesh partitionIdx = 0 */
    const UInt& numPoints = pts_->GetNumberOfPoints();
    NODECOORD.resize(numCoords * numPoints);

    /* get the tuples and write them into NODECOORD */
    UInt tmp_ij = 0;
    Double* point_coords;
    for (UInt i = 0; i < numPoints; ++i)
    {
      point_coords = pts_->GetPoint(i);
      /* tmp_ij = numCoords * i; */
      for (UInt j = 0; j < numCoords; ++j)
      {
        NODECOORD[tmp_ij] = point_coords[j];
        ++tmp_ij;
      }
    }
  }

  void FileReader_VTKMultiBlock::ReadTopology(std::vector<UInt>& TOPOLOGYDATA,
                                         std::vector<UInt>& elemTypes)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_VTKMultiBlock::ReadTopology" << std::endl;

    UInt numPoints;
    UInt numCells;
    UInt elemIdx = 0;
    UInt topoIdx = 0;
    std::vector<UInt> strucElemNodeMap;
    strucElemNodeMap.resize(maxNumElemNodes_);
    for(UInt i=0; i<maxNumElemNodes_; i++) 
    {
      strucElemNodeMap[i] = i;      
    }
    
    

    elemTypes.resize(numElems_);
    TOPOLOGYDATA.resize(numElems_ * maxNumElemNodes_);

    //    std::cout << "maxNumElemNodes_: " << maxNumElemNodes_ << std::endl;
    

    /* disregard pointZones and faceZones */
    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    //    vtkDataSet* ds1 = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    //    std::cout << "Building Kd tree point locator..." << std::endl;
    //    vtkKdTreePointLocator* locator = vtkKdTreePointLocator::New();
    //    locator->SetDataSet(ds1);
    //    locator->BuildLocator();
    //    std::cout << "Finished building Kd tree point locator..." << std::endl;

    actualRegionNodes_.resize(numRegions_);
    actualRegionNodesToNewIdx_.resize(numRegions_);
    
    /* goto to the desired partition */
    for (UInt i = 0; i < numRegions_; ++i)
    {
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());
      vtkCell* cell;
      std::set<UInt> regionNodes;
      std::map< UInt, UInt > regionNodeMap;
      std::vector<UInt> &elemNodeMapping = unstrucElemNodeMapping_[ Elem::POINT ];
  
      if(ugrid)
      {
        vtkCellArray* cells = ugrid->GetCells();
        vtkIdType cellId = 0;
	vtkIdType npts;
        vtkIdType *pts;
  
        //        std::cout << "Unstructured grid!" << std::endl;
        
        cells->InitTraversal();
        while(cells->GetNextCell(npts, pts)) 
        {
          regionNodes.insert(pts, pts+npts);
          cell = ugrid->GetCell(cellId);
          elemTypes[elemIdx] = VTKCellTypeToFEType(cell->GetCellType());

          elemNodeMapping = unstrucElemNodeMapping_[ elemTypes[elemIdx] ];
            //strucElemNodeMap;
            //unstrucElemNodeMapping_[ elemTypes[elemIdx] ];

          for (UInt k = 0; k < npts; ++k) 
          {
            UInt id = pts[ elemNodeMapping[k] ];
            TOPOLOGYDATA[topoIdx + k] = nodeMap_[i][id];
            regionNodeMap[ nodeMap_[i][id] ] = id;
          }
          
          elemIdx++;
          topoIdx += maxNumElemNodes_;
          cellId++;
        }
      } else 
      {
        // We have structured or even uniform grids
        vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
        numCells = ds->GetNumberOfCells();
        numPoints = ds->GetNumberOfPoints();

        for (UInt j = 0; j < numCells; ++j)
        {
          cell = ds->GetCell(j);
          numPoints = cell->GetNumberOfPoints();
          elemTypes[elemIdx] = VTKCellTypeToFEType(cell->GetCellType());

          elemNodeMapping = unstrucElemNodeMapping_[ elemTypes[elemIdx] ];

          //std::cout << "data object type " << ds->GetDataObjectType() << std::endl;
          switch(ds->GetDataObjectType()) {
          case VTK_UNIFORM_GRID: // 10
          case VTK_RECTILINEAR_GRID: // 3
          case VTK_IMAGE_DATA: // 6
            if (numPoints == 4) {
              elemTypes[elemIdx] = Elem::QUAD4;
              elemNodeMapping = uniformElemNodeMapping_[ elemTypes[elemIdx] ];
            } else {
              elemTypes[elemIdx] = Elem::HEXA8;
              elemNodeMapping = uniformElemNodeMapping_[ elemTypes[elemIdx] ];
            }            
            break;
          case VTK_STRUCTURED_GRID: // 2
            elemNodeMapping = strucElemNodeMap;
            break;
          }

          UInt localId, id;
          
          for (UInt k = 0; k < numPoints; ++k) 
          {
            localId = cell->GetPointId(elemNodeMapping[k]);
            regionNodes.insert(localId);

            id = nodeMap_[i][localId];
            TOPOLOGYDATA[topoIdx + k] = id;

            regionNodeMap[ id ] = localId;
          }

          elemIdx++;
          topoIdx += maxNumElemNodes_;
        }
      }

      //      std::cout << "#### Num nodes in region " << i << ": " << regionNodes.size() << std::endl;
      std::copy(regionNodes.begin(), regionNodes.end(),
                std::back_inserter(actualRegionNodes_[i]));
      
      if(settings.GetInt("verbose")) 
      {
        std::cout << "actualRegionNodes[" << i << "]:" << std::endl;       
        std::copy(regionNodes.begin(), regionNodes.end(),
                  std::ostream_iterator<UInt>(std::cout, " "));
        std::cout << std::endl;
      }
      

      std::map< UInt, UInt >::iterator it, end;
      it = regionNodeMap.begin();
      end = regionNodeMap.end();
      
      UInt idx = 0;
      for( ; it != end; it++ ) {
        UInt regionNodeIdx = it->second;

        //        std::cout << "global node: " << it->first << ", region node: " << regionNodeIdx << " -> " << idx << std::endl;

        actualRegionNodesToNewIdx_[i][regionNodeIdx] = idx++;
      }

      iter->GoToNextItem();
    }

    iter->Delete();
  }

  void FileReader_VTKMultiBlock::GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx)
  {
    UInt elemOffset = 0;

    for(UInt i=0; i < regionIdx; i++)
      elemOffset += numElemsPerRegion_[i];

    regionElements.resize(numElemsPerRegion_[regionIdx]);

    for(UInt i=0; i < numElemsPerRegion_[regionIdx]; i++)
      regionElements[i] = elemOffset + i + 1;
  }

#if 0
  /* get nodal values from the corresponding fluid datafile the new way */
  void FileReader_VTKMultiBlock::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const std::vector<bool>& activeParts,
                                            const UInt timeStepIdx)
  {

    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_VTKMultiBlock::ReadNodalValues..." << std::endl;

    SetTimeValue(timeValues_[timeStepIdx]);
    reader_->Update();

    if(settings.GetInt("verbose"))
      reader_->PrintSelf(std::cerr, vtkIndent());

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();

    iter->GoToFirstItem();

    for (UInt actRegion=0; actRegion < numRegions_; ++actRegion)
    {
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      vtkInformation* info = iter->GetCurrentMetaData();
      std::cout << "##### Region: " << info->Get(vtkCompositeDataSet::NAME())<< std::endl; 

      c2p->SetInput(ds);
      c2p->Update();

      vtkDataSet* ds_point = c2p->GetOutput();

      int nvx = ds->GetNumberOfPoints();
      UInt numTuples = actualRegionNodes_[actRegion].size();

      FlowDataType& fd = nodalFlowData[actRegion];
      UInt numDOFs = 0;

      FlowDataPartStruct* fdps = NULL;
      vtkArrayIteratorTemplate<float>* floatIt = NULL;

      vtkPointData* pointData = ds_point->GetPointData();
      UInt numArrays = pointData->GetNumberOfArrays();
      /* go through the solutions and store them in nodalFlowData */
      for (UInt array=0; array < numArrays; array++)
      {

        //        std::cout << "Array name: " << pointData->GetArrayName(array) << std::endl;
        vtkDataArray* data = pointData->GetArray(array);
        vtkArrayIterator* dataIt = data->NewIterator();
        /*vtkArrayIteratorTemplate<float>* doubleIt = NULL;*/
        std::string dsName = pointData->GetArrayName(array);
        UInt numComps = data->GetNumberOfComponents();


        //        std::cout << dsName << " " << numComps << " " << numTuples << std::endl;

        /* fluidVel_array = pointData->GetScalars(&u_char); <-- does not work
         * because the data is stored inside the array-variable not inside the
         * scalar- or vector-variable of the vtk class. */
        /* Get access to the fluid velocity data */
        /* check if the first array is fluid velocity data */
        if (dsName == "Velocity" &&
            (requiredResults_[FLUIDMECH_VELOCITY] ||
             requiredResults_[NO_SOLUTION_TYPE]) )
        {

          /* copy the fluid velocity values */
          fdps = &fd[FLUIDMECH_VELOCITY];
          fdps->isActive = 1; // all partitions have results

          if (fdps->dofNames.empty())
          {
            fdps->definedOn = ResultInfo::NODE; // nodes
            fdps->entryType = ResultInfo::VECTOR;
            fdps->dofNames.push_back("x");
            fdps->dofNames.push_back("y");
            if(dim_ == 3)
            {
              fdps->dofNames.push_back("z");
            }

            fdps->unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
            fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
          }
        }

        /* check if the array is fluid pressure data */
        if (dsName == "DENS" &&
            (requiredResults_[FLUIDMECH_PRESSURE] ||
             requiredResults_[NO_SOLUTION_TYPE]))
        {
          fdps = &fd[FLUIDMECH_PRESSURE];
          fdps->isActive = 1; // all partitions have results
          if (fdps->dofNames.empty())
          {
            fdps->definedOn = ResultInfo::NODE; // nodes
            fdps->entryType = ResultInfo::SCALAR;
            fdps->dofNames.push_back("-");
            fdps->unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
            fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE);
            fdps->data.resize(numDOFs * numTuples);
          }
        }

        if ((dsName == "Velocity" &&
             (requiredResults_[FLUIDMECH_VELOCITY] ||
              requiredResults_[NO_SOLUTION_TYPE])) ||
            (dsName == "DENS" &&
             (requiredResults_[FLUIDMECH_PRESSURE] ||
              requiredResults_[NO_SOLUTION_TYPE])))
        {
          numDOFs = fdps->dofNames.size();
          fdps->data.resize(numDOFs * numTuples);

          std::cout << dsName << " " << numComps << " " << numTuples << " nvx " << nvx << std::endl;


          /* copy the fluid velocity values */
          for(UInt i=0; i<numTuples; i++)
          {
            //            UInt node = actualRegionNodes_[actRegion][i];
            
            switch(data->GetDataType())
            {
            case VTK_FLOAT:
              floatIt = static_cast< vtkArrayIteratorTemplate<float>* >(dataIt);
              for(UInt j=0; j<numComps; j++) {
                float val = floatIt->GetValue(i*numComps+j);
                fdps->data[i*numDOFs+j] = val;
              }
              break;
              
            case VTK_DOUBLE:
              vtkArrayIteratorTemplate<Double>* doubleIt;
              doubleIt = static_cast< vtkArrayIteratorTemplate<Double>* >(dataIt);
              for(UInt j=0; j<numComps; j++) 
              {
                double val = doubleIt->GetValue(i*numComps+j);
                fdps->data[i*numDOFs+j] = val;
              }
              break;
            }
          }
        }
        
        dataIt->Delete();
      }

      iter->GoToNextItem();
    }

    iter->Delete();
    c2p->Delete();
  }
#endif

  Double FileReader_VTKMultiBlock::GetTimeStep(UInt stepNumber)
  {
    Settings& settings = Settings::Instance();
    Double timestep = settings.GetDouble("timestep");
    // check if timestep has been set as an argument
    if (timestep != -1)
    {
      return timestep * stepNumber;
    }

    return timeValues_[stepNumber];
  }

  std::string FileReader_VTKMultiBlock::GetRegionName(const UInt partitionIdx)
  {
    std::ostringstream sstr;
    std::string blockName;
    
    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    UInt i = 0;
    while (i <= partitionIdx && !iter->IsDoneWithTraversal())
    {
      sstr << "partition" << (i+1);
      blockName = sstr.str();
      sstr.clear(); sstr.str("");

      if(iter->HasCurrentMetaData()) 
      {
        vtkInformation* info = iter->GetCurrentMetaData();
        if(info->Has(vtkCompositeDataSet::NAME()))
          blockName = iter->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME());
      }
      
      //      std::cout << "region: " << blockName << std::endl;
      iter->GoToNextItem(); i++;
    }
    iter->Delete();
    
    const boost::regex datExp("[,;/ ]");
    const std::string name_format("_");
    
    sstr << regex_replace(blockName, datExp, name_format, boost::match_default | boost::format_sed);
    
    return sstr.str();
  }

  Elem::FEType FileReader_VTKMultiBlock::VTKCellTypeToFEType(UInt cellType)
  {
    switch(cellType) 
    {
    case VTK_VERTEX:
      return Elem::POINT;
    case VTK_LINE:
      return Elem::LINE2;
    case VTK_TRIANGLE:
      return Elem::TRIA3;
    case VTK_QUAD:
      return Elem::QUAD4;
    case VTK_TETRA:
      return Elem::TET4;
    case VTK_HEXAHEDRON:
      return Elem::HEXA8;
    case VTK_WEDGE:
      return Elem::WEDGE6;
    case VTK_PYRAMID:
      return Elem::PYRA5;

    // Quadratic, isoparametric cells
    case VTK_QUADRATIC_EDGE:
      return Elem::LINE3;
    case VTK_QUADRATIC_TRIANGLE:
      return Elem::TRIA6;
    case VTK_QUADRATIC_QUAD:
      return Elem::QUAD8;
    case VTK_QUADRATIC_TETRA:
      return Elem::TET10;
    case VTK_QUADRATIC_HEXAHEDRON:
      return Elem::HEXA20;
    case VTK_QUADRATIC_WEDGE:
      return Elem::WEDGE15;
    case VTK_QUADRATIC_PYRAMID:
      return Elem::PYRA13;
    case VTK_BIQUADRATIC_QUAD:
      return Elem::QUAD9;
    case VTK_TRIQUADRATIC_HEXAHEDRON:
      return Elem::HEXA27;
    case VTK_QUADRATIC_LINEAR_QUAD:
      return Elem::QUAD4;
    case VTK_QUADRATIC_LINEAR_WEDGE:
      return Elem::WEDGE6;
    case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
      return Elem::WEDGE15;
    case VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON:
      return Elem::HEXA20;

    default:
      return Elem::UNDEF;
    }
  }

  //BE warned this is now overloaded in Ensight file reader!!!!
  void FileReader_VTKMultiBlock::InitElemNodeMapping()
  {
    EnumMap::iterator it, end;

    it = Elem::feType.map.begin();
    end = Elem::feType.map.end();

    for( ; it != end; it++ ) 
    {
      UInt et = it->first;
      UInt numElemNodes = Elem::GetNumElemNodes((Elem::FEType)et);
      
      unstrucElemNodeMapping_[et].resize(numElemNodes);
      uniformElemNodeMapping_[et].resize(numElemNodes);
      
      for(UInt i=0; i<numElemNodes; i++) 
      {
        unstrucElemNodeMapping_[et][i] = i;
        uniformElemNodeMapping_[et][i] = i;
      }

      switch((Elem::FEType)et)
      {
      case Elem::QUAD4:
        uniformElemNodeMapping_[et][0] = 0;
        uniformElemNodeMapping_[et][1] = 1;
        uniformElemNodeMapping_[et][2] = 3;
        uniformElemNodeMapping_[et][3] = 2;
        break;
        
      case Elem::WEDGE6:
        unstrucElemNodeMapping_[et][0] = 3;
        unstrucElemNodeMapping_[et][1] = 4;
        unstrucElemNodeMapping_[et][2] = 5;
        unstrucElemNodeMapping_[et][3] = 0;
        unstrucElemNodeMapping_[et][4] = 1;
        unstrucElemNodeMapping_[et][5] = 2;
        break;
      case Elem::HEXA8:
        unstrucElemNodeMapping_[et][0] = 4;
        unstrucElemNodeMapping_[et][1] = 5;
        unstrucElemNodeMapping_[et][2] = 6;
        unstrucElemNodeMapping_[et][3] = 7;
        unstrucElemNodeMapping_[et][4] = 0;
        unstrucElemNodeMapping_[et][5] = 1;
        unstrucElemNodeMapping_[et][6] = 2;
        unstrucElemNodeMapping_[et][7] = 3;

        uniformElemNodeMapping_[et][0] = 0;
        uniformElemNodeMapping_[et][1] = 1;
        uniformElemNodeMapping_[et][2] = 3;
        uniformElemNodeMapping_[et][3] = 2;
        uniformElemNodeMapping_[et][4] = 4;
        uniformElemNodeMapping_[et][5] = 5;
        uniformElemNodeMapping_[et][6] = 7;
        uniformElemNodeMapping_[et][7] = 6;
        break;

      default:
        break;
      }
    }
    
  }
  
} // end of namespace
