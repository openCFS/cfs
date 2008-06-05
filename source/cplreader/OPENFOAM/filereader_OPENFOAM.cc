#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream> 

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;

#include "Domain/resultInfo.hh"

#include "../params.hh"
#include "../settings.hh"
#include "filereader_OPENFOAM.hh"


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
    reader_(NULL)
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
    std::stringstream controlDictName;

    /* TODO: check if file even exists, otherwise vtkOpenFOAMReader will hang */
    controlDictName << name_.c_str() << "/system/controlDict";
    reader_ = vtkOpenFOAMReader::New();
    if(settings.GetInt("verbose"))
      reader_->DebugOn();
    
    reader_->SetFileName(controlDictName.str().c_str());
    reader_->Update();

    if(settings.GetInt("verbose"))
      reader_->PrintSelf(std::cout, (vtkIndent)0);

    numFiles_ = reader_->GetNumberOfTimeSteps()-1;
    numResults_ = reader_->GetNumberOfCellArrays();

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    /* find out the number of partitions */
    numPartitions_ = 0;
    while (! iter->IsDoneWithTraversal())
    {
      ++numPartitions_;
      iter->GoToNextItem();
    }

    if(settings.GetInt("verbose")) 
    {
      std::cout << " Name: " << name_ << std::endl
                << " Dim: " << dim_ << std::endl
                << " numfiles: " << numFiles_ << std::endl
                << " numPartitions: " << numPartitions_ << std::endl
                << " numResults: " << numResults_ << std::endl
                << " timeStep: " << settings.GetDouble("timeStep") << std::endl;
    }
    
    elsize_.resize(numPartitions_);
    MpCCInodes_.resize(numPartitions_);
    MpCCIelems_.resize(numPartitions_);

    UInt p_cnt = 0;
    iter->GoToFirstItem();
    vtkDataSet* ds;
    while (! iter->IsDoneWithTraversal())
    {
      ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      numPoints = ds->GetNumberOfPoints();
      MpCCInodes_[p_cnt] = numPoints;
      /* cell values will be interpolated to element values */
      numElems = ds->GetNumberOfCells();
      MpCCIelems_[p_cnt] = numElems;

      /* just get the first cell, the number of points per element has to stay
       * the same inside a partition */
      vtkCell* cell = ds->GetCell(0);
      numPoints = cell->GetNumberOfPoints();
      elsize_[p_cnt] = numPoints;

      if(settings.GetInt("verbose"))
      {
        std::cout << "Partition " << (p_cnt+1)
                  << " nodes: " << MpCCInodes_[p_cnt]
                  << " elems: " << MpCCIelems_[p_cnt]
                  << " elsize: " << elsize_[p_cnt]
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
  }


  void FileReader_OPENFOAM::ReadNodalCoords(std::vector<Double>& NODECOORD,
                                            const UInt partitionIdx)
  {
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalCoords" << std::endl;
    /* number of coordinates. One Tupel of 3D coordinate (x,y,z) */
    const UInt numCoords = dim_;

    const UInt& numPoints = MpCCInodes_[partitionIdx];
    NODECOORD.resize(numCoords * numPoints);

    /* just read coords for internal mesh partitionIdx = 0 */
    if(partitionIdx)
      return;

    /* Goto first dataset (internal mesh) */
    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());

    /* get the tuples and write them into NODECOORD */
    UInt tmp_ij = 0;
    double* point_coords;
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
                                         std::vector<UInt>& numNodesPerElem,
                                         std::vector<UInt>& elemTypes,
                                         const UInt partitionIdx)
  {
    std::cout << "Entering FileReader_OPENFOAM::ReadTopology" << std::endl;
    UInt numPoints;
    UInt numCells;
    UInt topo_size = 0;
    double pt[3];
    std::map<UInt, UInt> nodeMap;

    /* disregard pointZones and faceZones */
    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    vtkDataSet* ds1 = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());

    /* goto to the desired partition */
    for (UInt i = 0; i < partitionIdx; ++i)
    {
        iter->GoToNextItem();
    }

    vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    iter->Delete();
    vtkCell* cell;
    numCells = ds->GetNumberOfCells();
    numPoints = ds->GetNumberOfPoints();

    // Map all points to nodes in the inner mesh (partitionIdx=0 -> ds1)
    for (UInt i = 0; i < numPoints; ++i)
    {
      ds->GetPoint(i, pt);

      vtkIdType ptId = ds1->FindPoint(pt);

      nodeMap[i] = ptId+1;
    }
    
    numNodesPerElem.resize(numCells);
    elemTypes.resize(numCells);
    for (UInt i = 0; i < numCells; ++i)
    {
      cell = ds->GetCell(i);
      numPoints = cell->GetNumberOfPoints();
      topo_size += numPoints;
      numNodesPerElem[i] = numPoints;
      elemTypes[i] = VTKCellTypeToFEType(cell->GetCellType());
    }

    TOPOLOGYDATA.resize(topo_size);
    UInt cnt = 0;
    for (UInt i = 0; i < numCells; ++i)
    {
      cell = ds->GetCell(i);
      numPoints = cell->GetNumberOfPoints();
      
      for (UInt j = 0; j < numPoints; ++j)
      {
        TOPOLOGYDATA[cnt] = nodeMap[cell->GetPointId(j)];
        ++cnt;
      }
    }


    // Hier werden die Elemente für eine Partition/Subdomain eingelesen.
    // partitionIdx kann hierbei von 0 ... n-1 gehen.
    // Für eine Übersicht über die Elementtypen siehe elementtypes.pdf

    // elemTypes[0] = ET_QUAD4
    // numNodesPerElem[0] = 4
    // TOPOLOGYDATA[0] = 1
    // TOPOLOGYDATA[1] = 2
    // TOPOLOGYDATA[2] = 3
    // TOPOLOGYDATA[3] = 4

    // elemTypes[1] = ET_TRIA3
    // numNodesPerElem[0] = 3
    // TOPOLOGYDATA[4] = 3
    // TOPOLOGYDATA[5] = 4
    // TOPOLOGYDATA[6] = 5
    //...
  }

  /* get nodal values from the corresponding fluid datafile the new way */
  void FileReader_OPENFOAM::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const UInt timeStepIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalValues..." << std::endl;

    reader_->SetTimeStep(timeStepIdx+1);
    reader_->Update();

    if(settings.GetInt("verbose"))
      reader_->PrintSelf(std::cerr, vtkIndent());

    /* store gap to jump over pointZones and faceZones */
    const UInt idx_gap = reader_->GetNumPointZones() + reader_->GetNumFaceZones();

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();
    
    iter->GoToFirstItem();

    for (UInt actPart=0; actPart < 1; actPart++)
    {
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      c2p->SetInput(ds);
      c2p->Update();

      vtkDataSet* ds_point = c2p->GetOutput();

      int nvx = ds_point->GetNumberOfPoints();
      FlowDataType& fd = nodalFlowData[actPart];
      UInt numDOFs;
      vtkPointData* pointData = ds_point->GetPointData();
      UInt numArrays = pointData->GetNumberOfArrays();


      for (UInt array=0; array < numArrays; array++)
      {

        //        std::cout << "Array name: " << pointData->GetArrayName(array) << std::endl;
        vtkDataArray* data = pointData->GetArray(array);
        vtkArrayIterator* dataIt = data->NewIterator();
        vtkArrayIteratorTemplate<float>* floatIt = NULL;
        vtkArrayIteratorTemplate<float>* doubleIt = NULL;
        std::string dsName = pointData->GetArrayName(array);
        UInt numComps = data->GetNumberOfComponents();
        UInt numTuples = data->GetNumberOfTuples();
        FlowDataPartStruct* fdps;
        
        /* fluidVel_array = pointData->GetScalars(&u_char); <-- does not work 
         * because the data is stored inside the array-variable not inside the
         * scalar- or vector-variable of the vtk class. */
        /* Get access to the fluid velocity data */
        /* check if the first array is fluid velocity data */
        if (dsName == "U")
        {
          /* copy the fluid velocity values */
          fdps = &fd[FLUIDMECH_VELOCITY];
          fdps->isActive = !actPart; // all partitions have results
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

          fdps->unit = "m s^-1";
          Enum2String(FLUIDMECH_VELOCITY, fdps->resultName);
          numDOFs = fdps->dofNames.size();
          fdps->data.resize(numDOFs * nvx);
        }

        /* check if the array is fluid pressure data */
        if (dsName == "p")
        {
          fdps = &fd[FLUIDMECH_PRESSURE];
          fdps->isActive = !actPart; // all partitions have results
          fdps->definedOn = ResultInfo::NODE; // nodes
          fdps->entryType = ResultInfo::SCALAR;
          if (!fdps->dofNames.size())
            fdps->dofNames.push_back("-");
          fdps->unit = "Pa";
          Enum2String(FLUIDMECH_PRESSURE, fdps->resultName);
          numDOFs = fdps->dofNames.size();
          fdps->data.resize(numDOFs * nvx);
        }

        if (dsName == "U" || dsName == "p") 
        {
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
              vtkArrayIteratorTemplate<double>* doubleIt;
              doubleIt = static_cast< vtkArrayIteratorTemplate<double>* >(dataIt);
              for(UInt j=0; j<numComps; j++)
                fdps->data[i*numDOFs+j] = doubleIt->GetValue(i*numComps+j);
              break;
            }
          }

          dataIt->Delete();
        }
      }

      iter->GoToNextItem();
    }

    iter->Delete();
    c2p->Delete();
  }

  double FileReader_OPENFOAM::GetTimeStep(UInt t)
  {
    return reader_->GetTimeStepValue(t);
  }

  std::string FileReader_OPENFOAM::GetPartitionName(const UInt partitionIdx)
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
    std::vector<std::string> fileNames;
    std::vector<std::string> dataSetNames;
    std::ifstream fin;
    std::stringstream sstr;

    fs::path foamDir( name_.c_str() );
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
           sstr << name_.c_str() << "/" << fn;
           fileNames.push_back(sstr.str());

           sstr.clear(); sstr.str("");
           sstr << name_.c_str() << "_" << fn;
           dataSetNames.push_back(sstr.str());
         }
       }
    }

    sstr.clear(); sstr.str("");
    sstr << name_.c_str() << "/" << 0;
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
          dir_itr != end_iter;
          ++dir_itr ) 
    {
       if ( !fs::is_directory( *dir_itr ) )
       {
         std::string fn = dir_itr->leaf();
         sstr.clear(); sstr.str("");
         sstr << name_.c_str() << "/0/" << fn;
         fileNames.push_back(sstr.str());

         sstr.clear(); sstr.str("");
         sstr << name_.c_str() << "_0_" << fn;
         dataSetNames.push_back(sstr.str());
       }
    }

    sstr.clear(); sstr.str("");
    sstr << name_.c_str() << "/constant";
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
          dir_itr != end_iter;
          ++dir_itr ) 
    {
       if ( !fs::is_directory( *dir_itr ) )
       {
         std::string fn = dir_itr->leaf();
         sstr.clear(); sstr.str("");
         sstr << name_.c_str() << "/constant/" << fn;
         fileNames.push_back(sstr.str());

         sstr.clear(); sstr.str("");
         sstr << name_.c_str() << "_constant_" << fn;
         dataSetNames.push_back(sstr.str());
       }
    }

    sstr.clear(); sstr.str("");
    sstr << name_.c_str() << "/constant/polyMesh";
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
          dir_itr != end_iter;
          ++dir_itr ) 
    {
       if ( !fs::is_directory( *dir_itr ) )
       {
         std::string fn = dir_itr->leaf();
         sstr.clear(); sstr.str("");
         sstr << name_.c_str() << "/constant/polyMesh/" << fn;
         fileNames.push_back(sstr.str());

         sstr.clear(); sstr.str("");
         sstr << name_.c_str() << "_constant_polyMesh_" << fn;
         dataSetNames.push_back(sstr.str());
       }
    }

    sstr.clear(); sstr.str("");
    sstr << name_.c_str() << "/system";
    foamDir = sstr.str();

    for ( fs::directory_iterator dir_itr( foamDir );
          dir_itr != end_iter;
          ++dir_itr ) 
    {
       if ( !fs::is_directory( *dir_itr ) )
       {
         std::string fn = dir_itr->leaf();
         sstr.clear(); sstr.str("");
         sstr << name_.c_str() << "/system/" << fn;
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

} // end of namespace
