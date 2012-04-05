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

#include "boost/iostreams/filtering_streambuf.hpp"
#include "boost/iostreams/copy.hpp"
#include "boost/iostreams/filter/gzip.hpp"

#include "boost/algorithm/string/predicate.hpp"
namespace algo=boost::algorithm;

namespace fs=boost::filesystem;
namespace bios=boost::iostreams;

#include "Domain/resultInfo.hh"

#include "cplreader/Settings.hh"
#include "FileReader_FLUENT.hh"


// This is due to the fucking OLAS New operator!!!
#undef New

#include "vtkFLUENTReader.h"

#include "vtkDataArrayCollection.h"
#include "vtkDataArrayCollectionIterator.h"
#include "vtkInformation.h"

#include "vtkCellArray.h"
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

  FileReader_FLUENT::FileReader_FLUENT(const std::string& name,
                                           const UInt dim,
                                           const UInt numFiles) :
    FileReader_VTKMultiBlock(name, dim, numFiles),
    oneCasToManyDat_(false)
  {
  }

  FileReader_FLUENT::~FileReader_FLUENT()
  {
    std::string fn;
    
    if(!fs::exists(workDir_))
       return;

    fn = workDir_ + "/cplreader_fluent_current.cas";
    if(fs::exists(fn))
      fs::remove(fn);

    fn = workDir_ + "/cplreader_fluent_current.dat";
    if(fs::exists(fn))
      fs::remove(fn);
  }

  void FileReader_FLUENT::CreateReader()
  {
    Settings& settings = Settings::Instance();
    std::stringstream sstr;

    sstr << settings.GetString("basedir") << "/"
         << name_.c_str();

    workDir_ = sstr.str();

    fs::path casDir( workDir_ );
    fs::directory_iterator end_iter;
    std::set<std::string> casFileNames;
    std::set<std::string> datFileNames;

    std::set<std::string>::const_iterator it, end, casIt;
    std::string fn;

    // Get all .cas[.gz]  and .dat[.gz] file names. We  either accept one .cas
    // and n .dat files or n .cas and n .dat files.
    for ( fs::directory_iterator dir_itr( casDir );
          dir_itr != end_iter;
          ++dir_itr )
    {
      if ( fs::is_directory( *dir_itr ) )
        continue;
      
      fn = dir_itr->path().filename().string();

      if(!algo::starts_with(fn, name_))
        continue;

      if(algo::ends_with(fn, ".cas") || algo::ends_with(fn, ".cas.gz"))
      {
        casFileNames.insert(fn);
        // std::cout << fn << std::endl;
      }

      if(algo::ends_with(fn, ".dat") || algo::ends_with(fn, ".dat.gz"))
      {
        datFileNames.insert(fn);
        // std::cout << fn << std::endl;
      }
    }

    UInt numCasFiles = casFileNames.size();
    UInt numDatFiles = datFileNames.size();

    if(!numDatFiles) 
    {
      EXCEPTION("No .dat[.gz] file could be found");
    }

    if(numCasFiles != 1) 
    {
      if(!numCasFiles) 
      {
        EXCEPTION("No .cas[.gz] file could be found");
      }
      
      if(numCasFiles != numDatFiles) 
      {
        EXCEPTION("The association of " << numCasFiles << 
                  " .cas[.gz] files with " << numDatFiles <<
                  " .dat[.gz] files is not one-to-one!");
      }

      it = datFileNames.begin();
      end = datFileNames.end();
      casIt = casFileNames.begin();
      UInt timeidx = 0;
      for(; it != end; it++, timeidx++, casIt++) 
      {
        tsCasFiles_[timeidx] = *casIt;
        tsDatFiles_[timeidx]  = *it;
      }      

      oneCasToManyDat_ = false;
    } else 
    {
      it = datFileNames.begin();
      end = datFileNames.end();
      UInt timeidx = 0;
      for(; it != end; it++, timeidx++) 
      {
        tsCasFiles_[timeidx] = *casFileNames.begin();
        tsDatFiles_[timeidx]  = *it;

        //        std::cout << "ts[" << timeidx << "]: " << tsCasFiles_[timeidx] <<
        //          " -> " << tsDatFiles_[timeidx] << std::endl;        
      }

      oneCasToManyDat_ = true;
    }

    GetTimeValues();
    SetTimeValue(timeValues_[0]);
  }
  
  void FileReader_FLUENT::EnableRegions() 
  {
  }

  void FileReader_FLUENT::GetTimeValues() 
  {
    Settings& settings = Settings::Instance();
    Double ts = settings.GetDouble("timestep");
    UInt n=tsDatFiles_.size();
    timeValues_.resize(n);
    
    for(UInt i=0; i<n; i++) 
    {
      timeValues_[i] = i*ts;
    }
  }

  void FileReader_FLUENT::SetTimeValue(Double val) 
  {
    std::stringstream sstr;
    UInt timeidx = 0;
    std::vector<Double>::iterator it;
    vtkFLUENTReader* reader = NULL;
    std::string currCasFile;
    std::string currDatFile;

    it = std::find(timeValues_.begin(), timeValues_.end(), val);
    timeidx = std::distance(timeValues_.begin(), it);

    if(val < 0.0) 
    {
      timeidx = 0;
    }
    
    currCasFile = tsCasFiles_[timeidx];
    currDatFile = tsDatFiles_[timeidx];

    if(oneCasToManyDat_) 
    {
      if(timeidx == 0)
        GUnzipFile(currCasFile, "cplreader_fluent_current.cas");
    } else 
    {
      GUnzipFile(currCasFile, "cplreader_fluent_current.cas");
    }
    GUnzipFile(currDatFile, "cplreader_fluent_current.dat");

    if(reader_)
      reader_->Delete();

    reader = vtkFLUENTReader::New();
    reader_ = reader;

    sstr << workDir_ << "/" << "cplreader_fluent_current.cas";

    reader->SetFileName(sstr.str().c_str());

    // Activate only arrays which are needed.
    reader->DisableAllCellArrays();
    if(requiredResults_[FLUIDMECH_VELOCITY]) 
    {
      reader->SetCellArrayStatus("X_VELOCITY", 1);
      reader->SetCellArrayStatus("Y_VELOCITY", 1);
      reader->SetCellArrayStatus("Z_VELOCITY", 1);
    }
    if(requiredResults_[FLUIDMECH_PRESSURE]) 
    {
      reader->SetCellArrayStatus("PRESSURE", 1);
    }    

    reader->Modified();
    reader->Update();
  }
  
  void FileReader_FLUENT::GUnzipFile(const std::string& fn,
                                     const std::string& outFN)
  {
    std::stringstream sstr;
    
    sstr << workDir_ << "/" << fn;

    std::ifstream file(sstr.str().c_str(),
                       std::ios_base::in | std::ios_base::binary);
    bios::filtering_streambuf<bios::input> in;
    if( algo::ends_with(fn, ".gz") ) 
    {
      in.push(bios::gzip_decompressor());
    }    
    in.push(file);

    sstr.clear(); sstr.str("");
    sstr << workDir_ << "/" << outFN;

    std::ofstream out(sstr.str().c_str(),
                      std::ios_base::out |
                      std::ios_base::binary |
                      std::ios_base::trunc );

    bios::copy(in, out);

    out.close();
    file.close();
  }

  /* get nodal values from the corresponding fluid datafile the new way */
  void FileReader_FLUENT::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const std::vector<bool>& activeParts,
                                            const UInt timeStepIdx)
  {

    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_FLUENT::ReadNodalValues..." << std::endl;

    SetTimeValue(timeValues_[timeStepIdx]);
    reader_->Modified();
    reader_->Update();

    if(settings.GetInt("verbose"))
      reader_->PrintSelf(std::cerr, vtkIndent());

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();

    iter->GoToFirstItem();

    UInt nodeOffset = 0;

    for (UInt actRegion=0; actRegion < numRegions_; ++actRegion)
    {
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      //      std::cout << "##### Region: " << actRegion << std::endl; 

      c2p->SetInput(ds);
      c2p->Update();

      vtkDataSet* ds_point = c2p->GetOutput();

      UInt numNodes = ds->GetNumberOfPoints();
      UInt numElems = ds->GetNumberOfCells();

      if( (numNodesPerRegion_[actRegion] != numNodes) ||
          (numElemsPerRegion_[actRegion] != numElems) ) 
      {
        EXCEPTION("Number of nodes or elements changed with time:\n"
                  << "numNodes at t = " << timeValues_[0] << "s: " 
                  << numNodesPerRegion_[actRegion] << "\n"
                  << "numNodes at t = " << timeValues_[timeStepIdx] << "s: " 
                  << numNodes << "\n"
                  << "numElems at t = " << timeValues_[0] << "s: " 
                  << numElemsPerRegion_[actRegion] << "\n"
                  << "numElems at t = " << timeValues_[timeStepIdx] << "s: " 
                  << numElems << "\n" <<
                  "Time dependent meshes are not supported by cplreader at the moment!");
      }
      

      UInt numTuples = actualRegionNodes_[actRegion].size();

      FlowDataType& fd = nodalFlowData[actRegion];
      UInt numDOFs = 0;

      FlowDataPartStruct* fdps = NULL;

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
        // UInt numTuples = data->GetNumberOfTuples();
        /* fluidVel_array = pointData->GetScalars(&u_char); <-- does not work
         * because the data is stored inside the array-variable not inside the
         * scalar- or vector-variable of the vtk class. */
        /* Get access to the fluid velocity data */
        /* check if the first array is fluid velocity data */
        if ((((dsName == "X_VELOCITY") || (dsName == "Y_VELOCITY") || (dsName == "Z_VELOCITY")) &&
             (requiredResults_[FLUIDMECH_VELOCITY] || requiredResults_[NO_SOLUTION_TYPE])))
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
        if (dsName == "PRESSURE" &&
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

        if ((((dsName == "X_VELOCITY") || (dsName == "Y_VELOCITY") || (dsName == "Z_VELOCITY")) &&
            (requiredResults_[FLUIDMECH_VELOCITY] || requiredResults_[NO_SOLUTION_TYPE])) ||
            (dsName == "PRESSURE" &&
             (requiredResults_[FLUIDMECH_PRESSURE] ||
              requiredResults_[NO_SOLUTION_TYPE])))
        {
          numDOFs = fdps->dofNames.size();
          fdps->data.resize(numDOFs * numTuples);

          //          std::cout << dsName << " " << numComps << " " << numTuples << " numNodes " << numNodes << std::endl;

          UInt offset = 0;
          if (dsName == "X_VELOCITY")
            offset = 0;
          if (dsName == "Y_VELOCITY")
            offset = 1;
          if (dsName == "Z_VELOCITY")
            offset = 2;

          /* copy the fluid velocity values */
          for(UInt i=0; i<numTuples; i++)
          {
            UInt in = actualRegionNodes_[actRegion][i];
            UInt out = actualRegionNodesToNewIdx_[actRegion][in];

            //            std::cout << "in: " << in << " -> " << out << std::endl;

#define EXTRACT_DATA_ARRAY(TYPE) \
            { \
              vtkArrayIteratorTemplate<TYPE>* TYPEIt = NULL; \
              TYPEIt = static_cast< vtkArrayIteratorTemplate<TYPE>* >(dataIt); \
              for(UInt j=0; j<numComps; j++) { \
                TYPE val = TYPEIt->GetValue(in*numComps+j); \
                fdps->data[out*numDOFs+offset] = val; \
              } \
            }   
            
            switch(data->GetDataType())
            {
            case VTK_FLOAT:
              EXTRACT_DATA_ARRAY(float);
              break;
              
            case VTK_DOUBLE:
              EXTRACT_DATA_ARRAY(double);
              break;
            }
          }

        }
        
        dataIt->Delete();
      }
      nodeOffset += numTuples;

      iter->GoToNextItem();
    }

    iter->Delete();
    c2p->Delete();
  }

  //! get user data from file reader
  void FileReader_FLUENT::GetUserData(std::map<std::string, std::string>& userData)
  {
    std::stringstream sstr;

    for(UInt timeidx=0; timeidx<tsDatFiles_.size(); timeidx++) 
    {
      sstr << "ts[" << timeidx << "] (t=" << timeValues_[timeidx] << 
        "s): " << tsCasFiles_[timeidx] <<
        " -> " << tsDatFiles_[timeidx] << std::endl;        
    }

    userData["CasToDatMapping"] = sstr.str();

    #if 0
    std::ifstream fin;
    std::stringstream sstr;
    vtkFLUENTReader* reader = dynamic_cast<vtkFLUENTReader*>(reader_);

    sstr << reader->GetFilePath() << reader->GetCaseFileName();
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
      sstr << reader->GetCaseFileName();
      userData[sstr.str()] = str;
      fin.close();

      #endif
  }

} // end of namespace
