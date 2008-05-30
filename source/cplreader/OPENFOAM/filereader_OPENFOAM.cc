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
    reader_->SetFileName(controlDictName.str().c_str());
    reader_->Update();
    numFiles_ = reader_->GetNumberOfTimeSteps();
    numResults_ = reader_->GetNumberOfCellArrays();

    vtkMultiBlockDataSet* readOuput = reader_->GetOutput();

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    /* find out the number of partitions */
    numPartitions_ = 0;
    while (! iter->IsDoneWithTraversal())
    {
      ++numPartitions_;
      iter->GoToNextItem();
    }

    std::cout << " Name: " << name_ << std::endl
              << " Dim: " << dim_ << std::endl
              << " numfiles: " << numFiles_ << std::endl
              << " numPartitions: " << numPartitions_ << std::endl
              << " numResults: " << numResults_ << std::endl
              << " timeStep: " << settings.GetDouble("timeStep") << std::endl;
    

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

      /* jsut get the first cell, the number of points per element has to stay
       * the same*/
      /* FIXME: is this correct? */
      vtkCell* cell = ds->GetCell(0);
      numPoints = cell->GetNumberOfPoints();
      elsize_[p_cnt] = numPoints;

      std::cout << "Partition " << (p_cnt+1)
        << " nodes: " << MpCCInodes_[p_cnt]
        << " elems: " << MpCCIelems_[p_cnt]
        << " elsize: " << elsize_[p_cnt]
        <<std::endl;
      ++p_cnt;
      iter->GoToNextItem();
    }
    iter->Delete();
    std::cout << "Exiting FileReader_OPENFOAM::Init" << std::endl;
  }


  void FileReader_OPENFOAM::ReadNodalCoords(std::vector<Double>& NODECOORD,
                                            const UInt partitionIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalCoords" << std::endl;
    /* number of coordinates. One Tupel of 3D coordinate (x,y,z) */
    const UInt numCoords = dim_;

    const UInt& numPoints = MpCCInodes_[partitionIdx];
    NODECOORD.resize(numCoords * numPoints);

    /* goto to the desired partition */
    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    for (UInt i = 0; i < partitionIdx; ++i)
    {
        iter->GoToNextItem();
    }
    vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());

    /* get the tupels and write them into NODECOORD */
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

    // Hier werden die Koordinaten für eine Partition/Subdomain eingelesen.
    // partitionIdx kann hierbei von 0 ... n-1 gehen.
    // NODECOORD[0] = x1
    // NODECOORD[1] = y1
    // NODECOORD[2] = z1

    // NODECOORD[3] = x2
    // NODECOORD[4] = y2
    // NODECOORD[5] = z2

    // In ReadNodalValues müssen die Werte in der gleichen Reihenfolge
    // eingelesen werden.
  }

  void FileReader_OPENFOAM::ReadTopology(std::vector<UInt>& TOPOLOGYDATA,
                                         std::vector<UInt>& numNodesPerElem,
                                         std::vector<UInt>& elemTypes,
                                         const UInt partitionIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadTopology" << std::endl;
    UInt numPoints;
    UInt numCells;
    UInt topo_size = 0;

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    /* goto to the desired partition */
    for (UInt i = 0; i < partitionIdx; ++i)
    {
        iter->GoToNextItem();
    }

    vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    iter->Delete();
    vtkCell* cell;
    numCells = ds->GetNumberOfCells();
    numNodesPerElem.resize(numCells);
    elemTypes.resize(numCells);
    for (UInt i = 0; i < numCells; ++i)
    {
      cell = ds->GetCell(i);
      numPoints = cell->GetNumberOfPoints();
      topo_size += numPoints;
      numNodesPerElem[i] = numPoints;
      /* TODO: a map should be made which reads the type from vtk and maps it to
       * the CFS type */
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
        TOPOLOGYDATA[cnt] = cell->GetPointId(j)+1;
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

  void FileReader_OPENFOAM::ReadNodalValues(std::vector<double>& flowdata,
                                            const UInt partitionIdx,
                                            const UInt timeStepIdx)
  {
    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalValues..." << std::endl;
    reader_->SetTimeStep(timeStepIdx);
    reader_->Update();

    flowdata.resize(7 * MpCCInodes_[partitionIdx]);

    char u_char = 'U';
    char p_char = 'p';

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    iter->GoToFirstItem();
    /* goto to the desired partition */
    for (UInt i = 0; i < partitionIdx; ++i)
    {
        iter->GoToNextItem();
    }
    
    vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();
    vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
    c2p->SetInput(ds);
    c2p->Update();
    vtkDataSet* ds_point = c2p->GetOutput();
    vtkPointData* pointData = ds_point->GetPointData();

    vtkDataArray* tupel_vel;
    vtkDataArray* scalar_pres;
    double* velo_vals;
    std::vector<Double> tempVec;
    tempVec.resize(numResults_);
    UInt tmp_ij;
    for (UInt i = 0; i < MpCCInodes_[partitionIdx]; ++i)
    {
      tmp_ij = 7 * i;
      std::fill(&flowdata[tmp_ij + 0], &flowdata[tmp_ij + 7], 0.0);

      tupel_vel = pointData->GetScalars(&u_char);
      velo_vals = tupel_vel->GetTuple3(i);
      tempVec[dataColumns_[1]] = velo_vals[1];
      tempVec[dataColumns_[2]] = velo_vals[2];
      tempVec[dataColumns_[3]] = velo_vals[3];
      scalar_pres = pointData->GetScalars(&p_char);
      tempVec[dataColumns_[4]] = scalar_pres->GetTuple1(i);

      for (UInt j = 0; j < dataColumns_.size(); ++j)
      {
        if (dataColumns_[j] > -1)
        {
            flowdata[tmp_ij] = tempVec[dataColumns_[j]];
        }
        ++tmp_ij;
      }
    }

    // flowdata layout:
    // flowdata[0] = Lighthill Quellterm 1
    // flowdata[1] = u1
    // flowdata[2] = v1 
    // flowdata[3] = w1
    // flowdata[4] = Druck1
    // flowdata[5] = N/A
    // flowdata[6] = N/A

    // flowdata[7] = Lighthill Quellterm 2
    // flowdata[8] = u2
    // flowdata[9] = v2
    // flowdata[10] = w2
    // flowdata[11] = Druck2
    // ...

    iter->Delete();
  }

  double FileReader_OPENFOAM::GetTimeStep(UInt t)
  {
    double ts;

    ts = reader_->GetTimeStepValue(t);
    return ts;
  }

  std::string FileReader_OPENFOAM::GetPartitionName(const UInt partitionIdx)
  {
    if(partitionIdx == 0)
      return "fluid";
    else
      return reader_->GetBoundaryName(partitionIdx-1);
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
