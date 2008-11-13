// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;

#include "Domain/resultInfo.hh"

#include "cplreader/Settings.hh"
#include "FileReader_OPENFOAM.hh"


// This is due to the fucking OLAS New operator!!!
#undef New
#include <vtkOpenFOAMReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataSet.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellDataToPointData.h>
#include <vtkArrayIteratorTemplate.h>

namespace CoupledField
{

  FileReader_OPENFOAM::FileReader_OPENFOAM(const std::string& name,
                                           const UInt dim,
                                           const UInt numFiles) :
    FileReader(name, dim, numFiles),
    reader_(NULL),
    numElems_(0)
  {
  }

  FileReader_OPENFOAM::~FileReader_OPENFOAM()
  {
  }

  void FileReader_OPENFOAM::Init()
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::Init" << std::endl;
    std::stringstream sstr;
    UInt numPoints;
    UInt numElems;
    UInt numResults;
    std::stringstream controlDictName;

    /* TODO: check if file even exists, otherwise vtkOpenFOAMReader will hang */
    controlDictName << settings.GetString("basedir") << "/"
                    << name_.c_str() << "/system/controlDict";
    reader_ = vtkOpenFOAMReader::New();
    if(settings.GetInt("verbose"))
      reader_->DebugOn();

    reader_->SetFileName(controlDictName.str().c_str());
    reader_->SetTimeStep(0);
    reader_->Update();

    if (settings.GetInt("verbose"))
    {
        reader_->PrintSelf(std::cout, (vtkIndent)0);
    }

    if (settings.GetInt("numsteps"))
    {
      numSteps_ = (UInt) settings.GetInt("numsteps");
      /* don not exceed the maximal number of timesteps possible */
      if (numSteps_ >= (UInt) reader_->GetNumberOfTimeSteps())
      {
          numSteps_ = (UInt) reader_->GetNumberOfTimeSteps()-1;
      }
    } else {
      numSteps_ = (UInt) reader_->GetNumberOfTimeSteps()-1;
    }
    numResults = reader_->GetNumberOfCellArrays();

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    /* find out the number of partitions */
    numRegions_ = 0;
    while (! iter->IsDoneWithTraversal())
    {
      ++numRegions_;
      iter->GoToNextItem();
    }

    if(settings.GetInt("verbose"))
    {
      std::cout << " Name: " << name_ << std::endl
                << " Dim: " << dim_ << std::endl
                << " numfiles: " << numSteps_ << std::endl
                << " numPartitions: " << numRegions_ << std::endl
                << " numResults: " << numResults << std::endl
                << " timeStep: " << settings.GetDouble("timestep") << std::endl;
    }

    numNodesPerRegion_.resize(numRegions_);
    numElemsPerRegion_.resize(numRegions_);

    UInt p_cnt = 0;
    iter->GoToFirstItem();
    vtkDataSet* ds;
    while (! iter->IsDoneWithTraversal())
    {
      ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      numPoints = ds->GetNumberOfPoints();
      numNodesPerRegion_[p_cnt] = numPoints;
      /* cell values will be interpolated to element values */
      numElems = ds->GetNumberOfCells();
      numElemsPerRegion_[p_cnt] = numElems;
      numElems_ += numElems;

      /* just get the first cell, the number of points per element has to stay
       * the same inside a partition */
      vtkCell* cell = ds->GetCell(0);
      numPoints = cell->GetNumberOfPoints();
      maxNumElemNodes_ = numPoints > maxNumElemNodes_ ? numPoints : maxNumElemNodes_;

      if(settings.GetInt("verbose"))
      {
        std::cout << "Partition " << (p_cnt+1)
                  << " nodes: " << numNodesPerRegion_[p_cnt]
                  << " elems: " << numElemsPerRegion_[p_cnt]
                  <<std::endl;
      }

      ++p_cnt;
      iter->GoToNextItem();
    }
    iter->Delete();

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

    std::cout << "Exiting FileReader_OPENFOAM::Init" << std::endl;
    /* nodalCoords_ should store the first mesh, which may be needed if we have
     * a moving mesh*/
    this->ReadNodalCoords(nodalCoords_);
  }


  void FileReader_OPENFOAM::ReadNodalCoords(std::vector<Double>& NODECOORD)
  {
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalCoords" << std::endl;
    /* number of coordinates. One Tupel of 3D coordinate (x,y,z) */
    const UInt numCoords = 3;

    /* just read coords for internal mesh partitionIdx = 0 */
    const UInt& numPoints = numNodesPerRegion_[0];
    NODECOORD.resize(numCoords * numPoints);

    /* Goto first dataset (internal mesh) */
    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());

    /* get the tuples and write them into NODECOORD */
    UInt tmp_ij = 0;
    Double* point_coords;
    for (UInt i = 0; i < numPoints; ++i)
    {
      point_coords = ds->GetPoint(i);
      /* tmp_ij = numCoords * i; */
      for (UInt j = 0; j < numCoords; ++j)
      {
        NODECOORD[tmp_ij] = point_coords[j];
        ++tmp_ij;
      }
    }
  }

  void FileReader_OPENFOAM::ReadTopology(std::vector<UInt>& TOPOLOGYDATA,
                                         std::vector<UInt>& elemTypes)
  {
    std::cout << "Entering FileReader_OPENFOAM::ReadTopology" << std::endl;
    UInt numPoints;
    UInt numCells;
    UInt elemIdx = 0;
    UInt topoIdx = 0;
    Double pt[3];
    std::map<UInt, UInt> nodeMap;

    elemTypes.resize(numElems_);
    TOPOLOGYDATA.resize(numElems_ * maxNumElemNodes_);

    /* disregard pointZones and faceZones */
    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    vtkDataSet* ds1 = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());

    /* goto to the desired partition */
    for (UInt i = 0; i < numRegions_; ++i)
    {
      vtkCell* cell;
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      numCells = ds->GetNumberOfCells();
      numPoints = ds->GetNumberOfPoints();

      // Map all points to nodes in the inner mesh (partitionIdx=0 -> ds1)
      for (UInt j = 0; j < numPoints; ++j)
      {
        ds->GetPoint(j, pt);

        vtkIdType ptId = ds1->FindPoint(pt);

        nodeMap[j] = ptId+1;
      }

      for (UInt j = 0; j < numCells; ++j)
      {
        cell = ds->GetCell(j);
        numPoints = cell->GetNumberOfPoints();
        elemTypes[elemIdx] = VTKCellTypeToFEType(cell->GetCellType());

        switch(elemTypes[elemIdx])
        {
        case ET_WEDGE6:
          TOPOLOGYDATA[topoIdx+0] = nodeMap[cell->GetPointId(3)];
          TOPOLOGYDATA[topoIdx+1] = nodeMap[cell->GetPointId(4)];
          TOPOLOGYDATA[topoIdx+2] = nodeMap[cell->GetPointId(5)];
          TOPOLOGYDATA[topoIdx+3] = nodeMap[cell->GetPointId(0)];
          TOPOLOGYDATA[topoIdx+4] = nodeMap[cell->GetPointId(1)];
          TOPOLOGYDATA[topoIdx+5] = nodeMap[cell->GetPointId(2)];
          break;
        case ET_HEXA8:
          TOPOLOGYDATA[topoIdx+0] = nodeMap[cell->GetPointId(4)];
          TOPOLOGYDATA[topoIdx+1] = nodeMap[cell->GetPointId(5)];
          TOPOLOGYDATA[topoIdx+2] = nodeMap[cell->GetPointId(6)];
          TOPOLOGYDATA[topoIdx+3] = nodeMap[cell->GetPointId(7)];
          TOPOLOGYDATA[topoIdx+4] = nodeMap[cell->GetPointId(0)];
          TOPOLOGYDATA[topoIdx+5] = nodeMap[cell->GetPointId(1)];
          TOPOLOGYDATA[topoIdx+6] = nodeMap[cell->GetPointId(2)];
          TOPOLOGYDATA[topoIdx+7] = nodeMap[cell->GetPointId(3)];
          break;
        default:
          for (UInt k = 0; k < numPoints; ++k)
            TOPOLOGYDATA[topoIdx + k] = nodeMap[cell->GetPointId(k)];
          break;
        }

        elemIdx++;
        topoIdx += maxNumElemNodes_;
      }

      iter->GoToNextItem();
    }

    iter->Delete();
  }

  void FileReader_OPENFOAM::GetRegionElements(std::vector<UInt> & regionElements,
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
  void FileReader_OPENFOAM::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const std::vector<bool>& activeParts,
                                            const UInt timeStepIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalValues..." << std::endl;

    reader_->SetTimeStep(timeStepIdx + 1);
    reader_->Update();

    if(settings.GetInt("verbose"))
      reader_->PrintSelf(std::cerr, vtkIndent());

    /* store gap to jump over pointZones and faceZones */
    /* const UInt idx_gap = reader_->GetNumPointZones() + reader_->GetNumFaceZones();*/

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();

    iter->GoToFirstItem();

    for (UInt actRegion=0; actRegion < 1; ++actRegion)
    {
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      c2p->SetInput(ds);
      c2p->Update();

      vtkDataSet* ds_point = c2p->GetOutput();

      int nvx = ds_point->GetNumberOfPoints();
      FlowDataType& fd = nodalFlowData[actRegion];
      UInt numDOFs;

      FlowDataPartStruct* fdps;
      vtkArrayIteratorTemplate<float>* floatIt = NULL;

      /*  Helper string which should help check if a seperate polyMesh in each timestep exists.
       *  This would indicate moving mesh */
      std::stringstream path_mechDispl;
      path_mechDispl << settings.GetString("basedir") << "/"
                     << settings.GetString("name") << "/"
                     << reader_->GetTimeStepValue(timeStepIdx + 1) <<"/polyMesh";
      /* if a new mesh for this timestep exists, get it an store the
       * mechDisplacement */
      if (fs::exists(path_mechDispl.str()) &&
          fs::is_directory(path_mechDispl.str()))
      {
        fdps = &fd[MECH_DISPLACEMENT];
        fdps->isActive = !actRegion; // all partitions have results
        fdps->definedOn = ResultInfo::NODE; // nodes
        fdps->entryType = ResultInfo::VECTOR;

        if (fdps->dofNames.empty())
        {
          fdps->dofNames.push_back("x");
          fdps->dofNames.push_back("y");
          if(dim_ == 3)
          {
            fdps->dofNames.push_back("z");
          }
        }

        fdps->unit = "m";
        Enum2String(MECH_DISPLACEMENT, fdps->resultName);
        numDOFs = fdps->dofNames.size();
        fdps->data.resize(numDOFs * nvx);
        if (!actRegion)
        {
          this->ReadNodalCoords(fdps->data);
          this->calcMechDisplacement(nodalCoords_, fdps->data);
        }
      }

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
        UInt numTuples = data->GetNumberOfTuples();

        /* fluidVel_array = pointData->GetScalars(&u_char); <-- does not work
         * because the data is stored inside the array-variable not inside the
         * scalar- or vector-variable of the vtk class. */
        /* Get access to the fluid velocity data */
        /* check if the first array is fluid velocity data */
        if (dsName == "U" &&
            (requiredResults_[FLUIDMECH_VELOCITY] ||
             requiredResults_[NO_SOLUTION_TYPE]) )
        {

          /* copy the fluid velocity values */
          fdps = &fd[FLUIDMECH_VELOCITY];
          fdps->isActive = !actRegion; // all partitions have results

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
            Enum2String(FLUIDMECH_VELOCITY, fdps->resultName);
          }
        }

        /* check if the array is fluid pressure data */
        if (dsName == "p" &&
            (requiredResults_[FLUIDMECH_PRESSURE] ||
             requiredResults_[NO_SOLUTION_TYPE]))
        {
          fdps = &fd[FLUIDMECH_PRESSURE];
          fdps->isActive = !actRegion; // all partitions have results
          if (fdps->dofNames.empty())
          {
            fdps->definedOn = ResultInfo::NODE; // nodes
            fdps->entryType = ResultInfo::SCALAR;
            fdps->dofNames.push_back("-");
            fdps->unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
            Enum2String(FLUIDMECH_PRESSURE, fdps->resultName);
            fdps->data.resize(numDOFs * nvx);
          }
        }

#if 0
        if ((dsName == "U" || dsName == "p") &&
            (requiredResults_[FLUIDMECH_PRESSURE] ||
             requiredResults_[FLUIDMECH_VELOCITY] ||
             requiredResults_[NO_SOLUTION_TYPE]))
#endif
          if ((dsName == "U" &&
                (requiredResults_[FLUIDMECH_VELOCITY] ||
                 requiredResults_[NO_SOLUTION_TYPE])) ||
              (dsName == "p" &&
               (requiredResults_[FLUIDMECH_PRESSURE] ||
                requiredResults_[NO_SOLUTION_TYPE])))
          {
            numDOFs = fdps->dofNames.size();
            fdps->data.resize(numDOFs * nvx);

            /* copy the fluid velocity values */
            for(UInt i=0; i<numTuples; i++)
            {
              switch(data->GetDataType())
              {
                case VTK_FLOAT:
                  floatIt = static_cast< vtkArrayIteratorTemplate<float>* >(dataIt);
                  for(UInt j=0; j<numComps; j++)
                    fdps->data[i*numDOFs+j] = floatIt->GetValue(i*numComps+j);
                  break;

                case VTK_DOUBLE:
                  vtkArrayIteratorTemplate<Double>* doubleIt;
                  doubleIt = static_cast< vtkArrayIteratorTemplate<Double>* >(dataIt);
                  for(UInt j=0; j<numComps; j++)
                    fdps->data[i*numDOFs+j] = doubleIt->GetValue(i*numComps+j);
                  break;
              }
            }
          }

        dataIt->Delete();
      }
    }

    iter->Delete();
    c2p->Delete();
  }

  Double FileReader_OPENFOAM::GetTimeStep(UInt stepNumber)
  {
    Settings& settings = Settings::Instance();
    Double timestep = settings.GetDouble("timestep");
    // check if timestep has been set as an argument
    if (timestep != -1)
    {
      return timestep * stepNumber;
    }

    return reader_->GetTimeStepValue(stepNumber);
  }

  std::string FileReader_OPENFOAM::GetRegionName(const UInt partitionIdx)
  {
    std::ostringstream sstr;
    const UInt idx_tmp = reader_->GetNumBoundaries();

    if (partitionIdx == 0)
    {
      sstr << "InternalMesh";
    } else {
      if (partitionIdx <= idx_tmp)
      {
        sstr << reader_->GetBoundaryName(partitionIdx - 1);
      } else {
        sstr << reader_->GetCellZoneName(partitionIdx - idx_tmp - 1);
      }
    }

    return sstr.str();
  }

  FEType FileReader_OPENFOAM::VTKCellTypeToFEType(UInt cellType)
  {
    static std::map<UInt, FEType> elemTypeMap;

    if(elemTypeMap.empty())
    {
      elemTypeMap[VTK_LINE] = ET_LINE2;
      elemTypeMap[VTK_TRIANGLE] = ET_TRIA3;
      elemTypeMap[VTK_QUAD] = ET_QUAD4;
      elemTypeMap[VTK_TETRA] = ET_TET4;
      elemTypeMap[VTK_HEXAHEDRON] = ET_HEXA8;
      elemTypeMap[VTK_WEDGE] = ET_WEDGE6;
      elemTypeMap[VTK_PYRAMID] = ET_PYRA5;
    }

    if(elemTypeMap.find(cellType) == elemTypeMap.end())
      return ET_UNDEF;
    else
      return elemTypeMap[cellType];
  }

  //! get user data from file reader
  void FileReader_OPENFOAM::GetUserData(std::map<std::string, std::string>& userData)
  {
    Settings& settings = Settings::Instance();
    std::vector<std::string> fileNames;
    std::vector<std::string> dataSetNames;
    std::ifstream fin;
    std::stringstream sstr;



    sstr << settings.GetString("basedir") << "/" << name_.c_str();
    fs::path foamDir = sstr.str();
    fs::directory_iterator end_iter;
    for ( fs::directory_iterator dir_itr( foamDir );
        dir_itr != end_iter;
        ++dir_itr )
    {
      if ( !fs::is_directory( *dir_itr ) )
      {
        std::string fn = dir_itr->leaf();
        if(fn.substr(0, 4) == "log.")
        {
          sstr.clear(); sstr.str("");
          sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/" << fn;
          fileNames.push_back(sstr.str());

          sstr.clear(); sstr.str("");
          sstr << name_.c_str() << "_" << fn;
          dataSetNames.push_back(sstr.str());
        }
      }
    }

    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/0";
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
        dir_itr != end_iter;
        ++dir_itr )
    {
      if ( !fs::is_directory( *dir_itr ) )
      {
        std::string fn = dir_itr->leaf();
        sstr.clear(); sstr.str("");
        sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/0/" << fn;
        fileNames.push_back(sstr.str());

        sstr.clear(); sstr.str("");
        sstr << name_.c_str() << "_0_" << fn;
        dataSetNames.push_back(sstr.str());
      }
    }

    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/constant";
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
        dir_itr != end_iter;
        ++dir_itr )
    {
      if ( !fs::is_directory( *dir_itr ) )
      {
        std::string fn = dir_itr->leaf();
        sstr.clear(); sstr.str("");
        sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/constant/" << fn;
        fileNames.push_back(sstr.str());

        sstr.clear(); sstr.str("");
        sstr << name_.c_str() << "_constant_" << fn;
        dataSetNames.push_back(sstr.str());
      }
    }

    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/constant/polyMesh";
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
        dir_itr != end_iter;
        ++dir_itr )
    {
      if ( !fs::is_directory( *dir_itr ) )
      {
        std::string fn = dir_itr->leaf();
        sstr.clear(); sstr.str("");
        sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/constant/polyMesh/" << fn;
        fileNames.push_back(sstr.str());

        sstr.clear(); sstr.str("");
        sstr << name_.c_str() << "_constant_polyMesh_" << fn;
        dataSetNames.push_back(sstr.str());
      }
    }

    sstr.clear(); sstr.str("");
    sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/system";
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
        dir_itr != end_iter;
        ++dir_itr )
    {
      if ( !fs::is_directory( *dir_itr ) )
      {
        std::string fn = dir_itr->leaf();
        sstr.clear(); sstr.str("");
        sstr << settings.GetString("basedir") << "/" << name_.c_str() << "/system/" << fn;
        fileNames.push_back(sstr.str());

        sstr.clear(); sstr.str("");
        sstr << name_.c_str() << "_system_" << fn;
        dataSetNames.push_back(sstr.str());
      }
    }

    for(UInt i=0; i<fileNames.size(); i++)
    {
      fin.open( fileNames[i].c_str(), std::ios::binary );

      if(fin.fail())
        EXCEPTION("Cannot open file '" << fileNames[i]
            <<"' to dump into HDF5!");

      // seek to the end of the file
      fin.seekg (0, std::ios::end);
      UInt numBytes = fin.tellg();
      fin.seekg (0, std::ios::beg);

      std::string str;
      str.resize(numBytes);
      fin.read(&str[0], numBytes);
      userData[dataSetNames[i]] = str;
      fin.close();
    }
  }

  void FileReader_OPENFOAM::calcMechDisplacement(const std::vector<Double>& origin, \
      std::vector<Double>& newCoords) const
  {
    const UInt sizeVec = newCoords.size();
    for (UInt i = 0; i < sizeVec; ++i)
    {
      newCoords[i] -= origin[i];
    }
  }

} // end of namespace
