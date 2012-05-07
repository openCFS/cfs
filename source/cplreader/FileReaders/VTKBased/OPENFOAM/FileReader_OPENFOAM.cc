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
#include "FileReader_OPENFOAM.hh"


// This is due to the fucking OLAS New operator!!!
#undef New
#include "vtkOpenFOAMReader.h"
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
#include "vtkDataSet.h"
#include "vtkCell.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkCellDataToPointData.h"
#include "vtkArrayIteratorTemplate.h"
#include "vtkDoubleArray.h"

namespace CoupledField
{

  FileReader_OPENFOAM::FileReader_OPENFOAM(const std::string& name,
                                           const UInt dim,
                                           const UInt numFiles) :
    FileReader_VTKMultiBlock(name, dim, numFiles)
  {
  }

  FileReader_OPENFOAM::~FileReader_OPENFOAM()
  {
  }

  void FileReader_OPENFOAM::CreateReader()
  {
    Settings& settings = Settings::Instance();
    std::stringstream controlDictName;

    /* TODO: check if file even exists, otherwise vtkOpenFOAMReader will hang */
    controlDictName << settings.GetString("basedir") << "/"
                    << name_.c_str() << "/system/controlDict";

    vtkOpenFOAMReader* reader = vtkOpenFOAMReader::New();
    reader->SetFileName(controlDictName.str().c_str());
    reader_ = reader;
  }
  
  void FileReader_OPENFOAM::EnableRegions() 
  {
    vtkOpenFOAMReader* reader = dynamic_cast<vtkOpenFOAMReader*>(reader_);
    reader->EnableAllPatchArrays();
  }

  void FileReader_OPENFOAM::GetTimeValues() 
  {
    // Get the time sets from the reader.
    std::set<Double> timeSteps;
    std::set<Double>::const_iterator tsIt, tsEnd;
    vtkOpenFOAMReader* reader = dynamic_cast<vtkOpenFOAMReader*>(reader_);

    for (int bla = 0; bla<reader->GetTimeValues()->GetSize(); bla++)
        timeSteps.insert(reader->GetTimeValues()->GetValue(bla));

    for(tsIt = timeSteps.begin(), tsEnd = timeSteps.end(); tsIt != tsEnd; tsIt++) 
    {
      // Do not add time value 0.0 to array of available time values since 
      // it contains OpenFOAM's initial value fields.
      if( (*tsIt) > 0.0 ) 
      {
        timeValues_.push_back( *tsIt );      
        //        std::cout << "ts: " << (*tsIt) << std::endl;
      }
      
    }
  }

  void FileReader_OPENFOAM::SetTimeValue(Double val) 
  {
    vtkOpenFOAMReader* reader = dynamic_cast<vtkOpenFOAMReader*>(reader_);
    
    if(val < 0.0) 
    {
      reader->SetTimeValue(0.0);
    } else
    {
      reader->SetTimeValue(val);
    }    
  }
  
  /* get nodal values from the corresponding fluid datafile the new way */
  void FileReader_OPENFOAM::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const std::vector<bool>& activeParts,
                                            const UInt timeStepIdx)
  {

    Settings& settings = Settings::Instance();
    std::cout << "Entering FileReader_OPENFOAM::ReadNodalValues..." << std::endl;
    
    
    
    SetTimeValue(timeValues_[timeStepIdx]);
    reader_->Modified();
    reader_->Update();

    if(settings.GetInt("verbose"))
      reader_->PrintSelf(std::cerr, vtkIndent());

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();

    iter->GoToFirstItem();

    for (UInt actRegion=0; actRegion < 1; ++actRegion)
    {
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      //      vtkInformation* info = iter->GetCurrentMetaData();
      //      std::cout << "##### Region: " << info->Get(vtkCompositeDataSet::NAME())<< std::endl; 

      c2p->SetInput(ds);
      c2p->Update();

      vtkDataSet* ds_point = c2p->GetOutput();

      int nvx = ds_point->GetNumberOfPoints();
      FlowDataType& fd = nodalFlowData[actRegion];
      UInt numDOFs = 0;

      FlowDataPartStruct* fdps = NULL;
      vtkArrayIteratorTemplate<float>* floatIt = NULL;

      /*  Helper string which should help check if a separate polyMesh in each timestep exists.
       *  This would indicate moving mesh */
      std::stringstream path_mechDispl;
      path_mechDispl << settings.GetString("basedir") << "/"
                     << settings.GetString("name") << "/"
                     << timeValues_[timeStepIdx]
                     << "/polyMesh";
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
        fdps->resultName = SolutionTypeEnum.ToString(MECH_DISPLACEMENT);
        numDOFs = fdps->dofNames.size();
        fdps->data.resize(numDOFs * nvx);
        if (!actRegion)
        {
          /* Goto first dataset (internal mesh) */
          vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
          iter->GoToFirstItem();
          vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
          Double *pt1;
          Double *pt2;

          /* Calculate the displacements in respect to the first time step */
          for (UInt i = 0, n=ds->GetNumberOfPoints(); i < n; ++i)
          {
            pt1 = pts_->GetPoint(i);
            pt2 = ds->GetPoint(i);
            for (UInt j = 0; j < 3; ++j)
            {
              fdps->data[i*3+j] = pt2[j] - pt1[j];
            }
          }
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
            fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
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
            fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE);
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
                    {
                    fdps->data[i*numDOFs+j] = floatIt->GetValue(i*numComps+j);
                    //cout<<"Wert: "<<floatIt->GetValue(i*numComps+j)<<endl;
                    }
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
        std::string fn = dir_itr->path().filename().string(); 
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
        std::string fn = dir_itr->path().filename().string();
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
        std::string fn = dir_itr->path().filename().string();
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
        std::string fn = dir_itr->path().filename().string();
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
        std::string fn = dir_itr->path().filename().string();
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

  void FileReader_OPENFOAM::InitElemNodeMapping()
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
        unstrucElemNodeMapping_[et][0] = 0;
        unstrucElemNodeMapping_[et][1] = 1;
        unstrucElemNodeMapping_[et][2] = 2;
        unstrucElemNodeMapping_[et][3] = 3;
        unstrucElemNodeMapping_[et][4] = 4;
        unstrucElemNodeMapping_[et][5] = 5;
        unstrucElemNodeMapping_[et][6] = 6;
        unstrucElemNodeMapping_[et][7] = 7;

        uniformElemNodeMapping_[et][0] = 0;
        uniformElemNodeMapping_[et][1] = 1;
        uniformElemNodeMapping_[et][2] = 2;
        uniformElemNodeMapping_[et][3] = 3;
        uniformElemNodeMapping_[et][4] = 4;
        uniformElemNodeMapping_[et][5] = 5;
        uniformElemNodeMapping_[et][6] = 6;
        uniformElemNodeMapping_[et][7] = 7;
        break;

      default:
        break;
      }
    }
    
  }
} // end of namespace
