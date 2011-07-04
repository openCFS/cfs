// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
#include <set>
#include <sys/stat.h>

#include <boost/regex.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;

#include "Domain/resultInfo.hh"

#include "cplreader/Settings.hh"
#include "FileReader_EnSight.hh"


// This is due to the fucking OLAS New operator!!!
#undef New

#include <vtkEnSight6Reader.h>
#include <vtkEnSight6BinaryReader.h>
#include <vtkEnSightGoldReader.h>
#include <vtkEnSightGoldBinaryReader.h>

#include <vtkDataArrayCollection.h>
#include <vtkDataArrayCollectionIterator.h>
#include <vtkInformation.h>

#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMergePoints.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointSet.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellDataToPointData.h>
#include <vtkArrayIteratorTemplate.h>

namespace CoupledField
{

  FileReader_EnSight::FileReader_EnSight(const std::string& name,
                                           const UInt dim,
                                           const UInt numFiles) :
    FileReader(name, dim, numFiles, 1),
    reader_(NULL),
    pts_(NULL),
    numElems_(0)
  {
  }

  FileReader_EnSight::~FileReader_EnSight()
  {
  }

  void FileReader_EnSight::Init()
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_EnSight::Init" << std::endl;
    std::stringstream sstr;
    UInt numPoints;
    UInt numElems;
    UInt numResults;
    std::stringstream caseFileName;

    /* TODO: check if file even exists, otherwise vtkEnSightGoldReader will hang */
    caseFileName << settings.GetString("basedir") << "/"
                    << name_.c_str() << "/" << name_.c_str() << ".case";

    // First try to determine the version of the EnSight file
    reader_ = vtkEnSightReader::New();
    if(settings.GetInt("verbose"))
      reader_->DebugOn();

    reader_->SetCaseFileName(caseFileName.str().c_str());
#if 0
    int fileVersion = reader_->DetermineEnSightVersion(0);
    reader_->Delete();
    reader_ = NULL;
    
    switch(fileVersion) 
    {
    case vtkGenericEnSightReader::ENSIGHT_6:
      reader_ = vtkEnSight6Reader::New();
      break;
    case vtkGenericEnSightReader::ENSIGHT_6_BINARY:
      reader_ = vtkEnSight6BinaryReader::New();
      break;
    case vtkGenericEnSightReader::ENSIGHT_GOLD:
      reader_ = vtkEnSightGoldReader::New();
      break;
    case vtkGenericEnSightReader::ENSIGHT_GOLD_BINARY:
      reader_ = vtkEnSightGoldBinaryReader::New();
      break;
    case vtkGenericEnSightReader::ENSIGHT_MASTER_SERVER:
      EXCEPTION("Cannot read EnSight master server files at the moment.");
      break;
    default:
      EXCEPTION("A problem occurred while determining EnSight file type of '" << caseFileName << "'");
      break;
    }
#endif    
    if(settings.GetInt("verbose"))
      reader_->DebugOn();

    reader_->SetCaseFileName(caseFileName.str().c_str());
    reader_->Update();

    if (settings.GetInt("verbose"))
    {
        reader_->PrintSelf(std::cout, (vtkIndent)0);
    }

    // Get the time sets from the reader.
    std::set<Double> timeSteps;
    std::set<Double>::const_iterator tsIt, tsEnd;
    vtkDataArrayCollection* timeSets = reader_->GetTimeSets();

    // Iterate through the time sets.
    vtkDataArrayCollectionIterator* daciter = vtkDataArrayCollectionIterator::New();
    daciter->SetCollection(timeSets);
    for(daciter->GoToFirstItem(); !daciter->IsDoneWithTraversal();
        daciter->GoToNextItem())
    {
      // Each time set is stored in one message.
      vtkDataArray* da = daciter->GetDataArray();
      for(int i=0; i < da->GetNumberOfTuples(); ++i)
      {
        timeSteps.insert(da->GetTuple1(i));
      }
    }
    daciter->Delete();

    for(tsIt = timeSteps.begin(), tsEnd = timeSteps.end(); tsIt != tsEnd; tsIt++) 
    {
      timeValues_.push_back( *tsIt );      
      std::cout << "ts: " << (*tsIt) << std::endl;
    }
    
    UInt numTSFromFile = 0;
    if(timeSteps.size()) 
    {
      numTSFromFile = timeSteps.size();
    } else 
    {
      numTSFromFile = 1;
      timeValues_.push_back( 0 );
    }
    
    reader_->SetTimeValue(timeValues_[0]);
    reader_->Update();

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

    // TODO: numResults richtig setzen
    numResults = reader_->GetNumberOfCellArrays();

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    /* find out the number of partitions */
    numRegions_ = 0;
    while (! iter->IsDoneWithTraversal())
    {
      ++numRegions_;
      iter->GoToNextItem();
      std::cout << "region: " << numRegions_ << std::endl;      
    }

    if(settings.GetInt("verbose"))
    {
      std::cout << " Name: " << name_ << std::endl
                << " Dim: " << dim_ << std::endl
                << " numfiles: " << numSteps_ << std::endl
                << " numPartitions: " << numRegions_ << std::endl
                << " numResults: " << numResults << std::endl
                << " timeStep: " << (timeValues_[1]-timeValues_[0]) << std::endl

                << "NumberOfScalarsPerNode: " << reader_->GetNumberOfScalarsPerNode() << std::endl
                << "NumberOfVectorsPerNode: " << reader_->GetNumberOfVectorsPerNode() << std::endl
                << "NumberOfTensorsSymmPerNode: " << reader_->GetNumberOfTensorsSymmPerNode() << std::endl  
                << "NumberOfScalarsPerElement: " << reader_->GetNumberOfScalarsPerElement() << std::endl
                << "NumberOfVectorsPerElement: " << reader_->GetNumberOfVectorsPerElement() << std::endl
                << "NumberOfTensorsSymmPerElement: " << reader_->GetNumberOfTensorsSymmPerElement() << std::endl
                << "NumberOfScalarsPerMeasuredNode: " << reader_->GetNumberOfScalarsPerMeasuredNode() << std::endl
                << "NumberOfVectorsPerMeasuredNode: " << reader_->GetNumberOfVectorsPerMeasuredNode() << std::endl
                << "NumberOfComplexScalarsPerNode: " << reader_->GetNumberOfComplexScalarsPerNode() << std::endl
                << "NumberOfComplexVectorsPerNode: " << reader_->GetNumberOfComplexVectorsPerNode() << std::endl
                << "NumberOfComplexScalarsPerElement: " << reader_->GetNumberOfComplexScalarsPerElement() << std::endl
                << "NumberOfComplexVectorsPerElement: " << reader_->GetNumberOfComplexVectorsPerElement() << std::endl;
    }

    numNodesPerRegion_.resize(numRegions_);
    numElemsPerRegion_.resize(numRegions_);

    pts_ = vtkUnstructuredGrid::New();
    pts_->Initialize();

    UInt p_cnt = 0;
    iter->GoToFirstItem();
    vtkDataSet* ds;
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


    std::cout << "Building vtkMergePoints locator..." << std::endl;
    vtkMergePoints* merger = vtkMergePoints::New();
    merger->SetDataSet(pts_);
    //    merger->BuildLocator();
    merger->InitPointInsertion(pts_->GetPoints(), pts_->GetBounds());
    Double pt[3];
    std::cout << "Finished building vtkMergePoints locator..." << std::endl;

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

        //        std::cout << "unique: " << unique << ", ptId: " << ptId << std::endl;
          
        nodeMap_[p_cnt][j] = ptId+1;
      }
      

      ++p_cnt;
      iter->GoToNextItem();
    }
    iter->Delete();

    //    numNodesPerRegion_[0] = ds1->GetNumberOfPoints();


#if 0
    if(settings.GetInt("verbose"))
    {
      std::cout << "Number of boundaries: "
                << reader_->GetNumBoundaries() << std::endl;
      std::cout << "Number of point zones: "
                << reader_->GetNumPointZones() << std::endl;
      std::cout << "Number of face zones: "
                << reader_->GetNumFaceZones() << std::endl;
      std::cout << "Number of cell zones: "
                << reader_->GetNumCellZones() << std::endl;
    }
#endif

    std::cout << "Exiting FileReader_EnSight::Init" << std::endl;
    /* nodalCoords_ should store the first mesh, which may be needed if we have
     * a moving mesh*/
    this->ReadNodalCoords(nodalCoords_);
  }


  void FileReader_EnSight::ReadNodalCoords(std::vector<Double>& NODECOORD)
  {
    std::cout << "Entering FileReader_EnSight::ReadNodalCoords" << std::endl;
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

  void FileReader_EnSight::ReadTopology(std::vector<UInt>& TOPOLOGYDATA,
                                         std::vector<UInt>& elemTypes)
  {
    std::cout << "Entering FileReader_EnSight::ReadTopology" << std::endl;

    UInt numPoints;
    UInt numCells;
    UInt elemIdx = 0;
    UInt topoIdx = 0;
    //    Double pt[3];

    elemTypes.resize(numElems_);
    TOPOLOGYDATA.resize(numElems_ * maxNumElemNodes_);

    std::cout << "maxNumElemNodes_: " << maxNumElemNodes_ << std::endl;
    

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

    /* goto to the desired partition */
    for (UInt i = 0; i < numRegions_; ++i)
    {
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(iter->GetCurrentDataObject());
      vtkCell* cell;
      std::set<UInt> regionNodes;

      if(ugrid)
      {
        vtkCellArray* cells = ugrid->GetCells();
        vtkIdType cellId = 0;
	vtkIdType npts;
        vtkIdType *pts;
        
        cells->InitTraversal();
        while(cells->GetNextCell(npts, pts)) 
        {
          regionNodes.insert(pts, pts+npts);
          cell = ugrid->GetCell(cellId);
          elemTypes[elemIdx] = VTKCellTypeToFEType(cell->GetCellType());

          for (UInt k = 0; k < npts; ++k)
            TOPOLOGYDATA[topoIdx + k] = nodeMap_[i][ pts[k] ];

          elemIdx++;
          topoIdx += maxNumElemNodes_;
          cellId++;
        }

        std::cout << "#### Num nodes in region " << i << ": " << regionNodes.size() << std::endl;
        std::copy(regionNodes.begin(), regionNodes.end(),
                  std::back_inserter(actualRegionNodes_[i]));
        
      } else 
      {
        // We have structured or even uniform grids
        vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
        numCells = ds->GetNumberOfCells();
        numPoints = ds->GetNumberOfPoints();
        std::vector<UInt> elemNodeMapping(maxNumElemNodes_);
        
        for (UInt j = 0; j < numCells; ++j)
        {
          cell = ds->GetCell(j);
          numPoints = cell->GetNumberOfPoints();
          elemTypes[elemIdx] = VTKCellTypeToFEType(cell->GetCellType());

          switch(elemTypes[elemIdx])
          {
          case Elem::WEDGE6:
            elemNodeMapping[0] = 3;
            elemNodeMapping[1] = 4;
            elemNodeMapping[2] = 5;
            elemNodeMapping[3] = 0;
            elemNodeMapping[4] = 1;
            elemNodeMapping[5] = 2;
            break;
          case Elem::HEXA8:
            elemNodeMapping[0] = 4;
            elemNodeMapping[1] = 5;
            elemNodeMapping[2] = 6;
            elemNodeMapping[3] = 7;
            elemNodeMapping[4] = 0;
            elemNodeMapping[5] = 1;
            elemNodeMapping[6] = 2;
            elemNodeMapping[7] = 3;
            break;
          default:
            for (UInt k = 0; k < numPoints; ++k) 
            {
              elemNodeMapping[k] = k;
            }
            break;
          }
          
          std::cout << "data object type " << ds->GetDataObjectType() << std::endl;
          switch(ds->GetDataObjectType()) {
          case VTK_UNIFORM_GRID:
            //          case VTK_STRUCTURED_GRID:
          case VTK_IMAGE_DATA:
            if (numPoints == 4) {
              elemTypes[elemIdx] = Elem::QUAD4;
              elemNodeMapping[0] = 0;
              elemNodeMapping[1] = 1;
              elemNodeMapping[2] = 3;
              elemNodeMapping[3] = 2;
            } else {
              elemTypes[elemIdx] = Elem::HEXA8;
              elemNodeMapping[0] = 0;
              elemNodeMapping[1] = 1;
              elemNodeMapping[2] = 3;
              elemNodeMapping[3] = 2;
              elemNodeMapping[4] = 4;
              elemNodeMapping[5] = 5;
              elemNodeMapping[6] = 7;
              elemNodeMapping[7] = 6;
            }            
            break;
          }

          UInt id;
          
          for (UInt k = 0; k < numPoints; ++k) 
          {
            id = nodeMap_[i][cell->GetPointId(elemNodeMapping[k])];
            TOPOLOGYDATA[topoIdx + k] = id;
            regionNodes.insert(id);
          }

          elemIdx++;
          topoIdx += maxNumElemNodes_;
        }

        std::cout << "#### Num nodes in region " << i << ": " << regionNodes.size() << std::endl;
        std::copy(regionNodes.begin(), regionNodes.end(),
                  std::back_inserter(actualRegionNodes_[i]));

      }
      
      iter->GoToNextItem();
    }

    iter->Delete();
  }

  void FileReader_EnSight::GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx)
  {
    UInt elemOffset = 0;

    for(UInt i=0; i < regionIdx; i++)
      elemOffset += numElemsPerRegion_[i];

    regionElements.resize(numElemsPerRegion_[regionIdx]);

    for(UInt i=0; i < numElemsPerRegion_[regionIdx]; i++)
      regionElements[i] = elemOffset + i + 1;
  }

  /* get nodal values from the corresponding fluid datafile the new way */
  void FileReader_EnSight::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const std::vector<bool>& activeParts,
                                            const UInt timeStepIdx)
  {

    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_EnSight::ReadNodalValues..." << std::endl;

    reader_->SetTimeValue(timeValues_[timeStepIdx]);
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

  Double FileReader_EnSight::GetTimeStep(UInt stepNumber)
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

  std::string FileReader_EnSight::GetRegionName(const UInt partitionIdx)
  {
    std::ostringstream sstr;
    std::string partId;
    
    partId = reader_->GetOutput()->GetMetaData(partitionIdx)->Get(vtkCompositeDataSet::NAME());
    
    const boost::regex datExp("[,;/ ]");
    const std::string name_format("_");

    sstr << regex_replace(partId, datExp, name_format, boost::match_default | boost::format_sed);

    return sstr.str();
  }

  Elem::FEType FileReader_EnSight::VTKCellTypeToFEType(UInt cellType)
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

  //! get user data from file reader
  void FileReader_EnSight::GetUserData(std::map<std::string, std::string>& userData)
  {
    std::ifstream fin;
    std::stringstream sstr;

    sstr << reader_->GetFilePath() << reader_->GetCaseFileName();
      fin.open( sstr.str().c_str(), std::ios::binary );

      if(fin.fail())
        EXCEPTION("Cannot open file '" << sstr.str()
            <<"' to dump into HDF5!");

      // seek to the end of the file
      fin.seekg (0, std::ios::end);
      UInt numBytes = fin.tellg();
      fin.seekg (0, std::ios::beg);

      std::string str;
      str.resize(numBytes);
      fin.read(&str[0], numBytes);
      sstr.clear(); sstr.str("");
      sstr << reader_->GetCaseFileName();
      userData[sstr.str()] = str;
      fin.close();
  }

} // end of namespace
