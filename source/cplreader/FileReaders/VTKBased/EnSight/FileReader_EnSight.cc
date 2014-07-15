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
#include "FileReader_EnSight.hh"


// This is due to the fucking OLAS New operator!!!
#undef New

#include "vtkEnSight6Reader.h"
#include "vtkEnSight6BinaryReader.h"
#include "vtkEnSightGoldReader.h"
#include "vtkEnSightGoldBinaryReader.h"

#include "vtkDataArrayCollection.h"
#include "vtkDataArrayCollectionIterator.h"
#include "vtkInformation.h"
#include "vtkCellData.h"

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

  FileReader_EnSight::FileReader_EnSight(const std::string& name,
                                           const UInt dim,
                                           const UInt numFiles) :
    FileReader_VTKMultiBlock(name, dim, numFiles)
  {
  }

  FileReader_EnSight::~FileReader_EnSight()
  {
  }

  void FileReader_EnSight::CreateReader()
  {
    Settings& settings = Settings::Instance();
    std::stringstream caseFileName;

    /* TODO: check if file even exists, otherwise vtkEnSightGoldReader will hang */
    caseFileName << settings.GetString("basedir") << "/"
                    << name_.c_str() << "/" << name_.c_str() << ".case";

    // Create new EnSight reader
    vtkGenericEnSightReader* reader = vtkEnSightReader::New();
    reader->SetCaseFileName(caseFileName.str().c_str());
    reader_ = reader;

    // Set names for velocity and pressure.
    vx_ = settings.GetString("vx");
    vy_ = settings.GetString("vy");
    vz_ = settings.GetString("vz");
    pres_ = settings.GetString("pres");
    presD2_ = settings.GetString("presD2");

    // Either velocity  is a  single vector  or consists of  2 or  3 component
    // fields.  In the first case the names  for vx, vz and vz are required to
    // match. In the second case we require 2 or three names for the component
    // fields.
    if(vx_ == vy_) {
      if(dim_ == 3 && vz_ != vx_) 
      {
        EXCEPTION("Array names for vx, vy and vz do not match for vector "
                  << "valued velocity field and dimension is 3!" << std::endl
                    << "vx: " << vx_ << std::endl
                    << "vy: " << vy_ << std::endl
                    << "vz: " << vz_ << std::endl);
      }      
    } else {
      if(dim_ == 3)
      {
        if( vz_ == vx_ || vz_ == vy_ ) 
        {
          EXCEPTION("Array names for vx, vy and vz should not match for a "
                    << "velocity field which consists of scalar components "
                    << "and dimension is 3!" << std::endl
                    << "vx: " << vx_ << std::endl
                    << "vy: " << vy_ << std::endl
                    << "vz: " << vz_ << std::endl);
        }        
      }      
    }
    
    reader->ReadAllVariablesOff();
    if(requiredResults_[FLUIDMECH_VELOCITY]) 
    {
      reader->SetPointArrayStatus(vx_.c_str(), 1);
      reader->SetCellArrayStatus(vx_.c_str(), 1);
      reader->SetPointArrayStatus(vy_.c_str(), 1);
      reader->SetCellArrayStatus(vy_.c_str(), 1);
      reader->SetPointArrayStatus(vz_.c_str(), 1);
      reader->SetCellArrayStatus(vz_.c_str(), 1);

      //      std::cout << "vx: " << vx_ << " vy: " << vy_ << " vz: " << vz_ << " pres: " << pres_ << std::endl;
    } else 
    {
      vx_ = ""; vy_ = ""; vz_ = "";
    }
    
    reader->SetPointArrayStatus(pres_.c_str(), 1);
    reader->SetCellArrayStatus(pres_.c_str(), 1);

    //if(requiredResults_[FLUIDMECH_PRESSURE_DERIV_2]){
      reader->SetPointArrayStatus(presD2_.c_str(), 1);
      reader->SetCellArrayStatus(presD2_.c_str(), 1);
    //}
  }
  
  void FileReader_EnSight::EnableRegions() 
  {
  }

  void FileReader_EnSight::GetTimeValues() 
  {
    // Get the time sets from the reader.
    std::set<Double> timeSteps;
    std::set<Double>::const_iterator tsIt, tsEnd;
    vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

    vtkDataArrayCollection* timeSets = reader->GetTimeSets();

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
  }

  void FileReader_EnSight::SetTimeValue(Double val) 
  {
    vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

    if(val < 0.0) 
    {
      reader->SetTimeValue(0.0);
    } else
    {
      reader->SetTimeValue(val);
    }    
  }

  void FileReader_EnSight::ReadElemValues(std::vector<FlowDataType>& nodalFlowData,
                                          const std::vector<bool>& activeParts,
                                          const UInt timeStepIdx){
  }

  
  /* get nodal values from the corresponding fluid datafile the new way */
  void FileReader_EnSight::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const std::vector<bool>& activeParts,
                                            const UInt timeStepIdx)
  {

    Settings& settings = Settings::Instance();
    //std::cout << "Entering FileReader_EnSight::ReadNodalValues..." << std::endl;

    SetTimeValue(timeValues_[timeStepIdx]);
    reader_->Modified();
    reader_->Update();

    if(settings.GetInt("verbose"))
      reader_->PrintSelf(std::cerr, vtkIndent());

    vtkCompositeDataIterator* iter = reader_->GetOutput()->NewIterator();
    vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();

    iter->GoToFirstItem();

    UInt node = 0;

    for (UInt actRegion=0; actRegion < numRegions_; ++actRegion)
    {
      vtkDataSet* ds = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      //      vtkInformation* info = iter->GetCurrentMetaData();
      //      std::cout << "##### Region: " << info->Get(vtkCompositeDataSet::NAME())<< std::endl; 
      std::vector<vtkDataSet*> datasets;
      
      vtkPointData* pointData = ds->GetPointData();
      UInt numArrays = pointData->GetNumberOfArrays();
      if(numArrays) datasets.push_back(ds);
      
      vtkCellData* cellData = ds->GetCellData();
      UInt numCellArrays = cellData->GetNumberOfArrays();
      
      if(numCellArrays) {
        c2p->SetInput(ds);
        c2p->Update();

        datasets.push_back( c2p->GetOutput() );
      }

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

      for(UInt  dset=0; dset < datasets.size(); dset++) {
         pointData = datasets[dset]->GetPointData();
         numArrays = pointData->GetNumberOfArrays();

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

        /* Get access to the fluid velocity data */
        /* check if the first array is fluid velocity data */
        if ((dsName == vx_ || dsName == vy_ || dsName == vz_) &&
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
        if (dsName == pres_ &&
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

        /* check if the array is fluid pressure data */
        if (dsName == presD2_ &&
            (requiredResults_[FLUIDMECH_PRESSURE_DERIV_2] ||
             requiredResults_[NO_SOLUTION_TYPE]))
        {
          fdps = &fd[FLUIDMECH_PRESSURE_DERIV_2];
          fdps->isActive = 1; // all partitions have results
          if (fdps->dofNames.empty())
          {
            fdps->definedOn = ResultInfo::NODE; // nodes
            fdps->entryType = ResultInfo::SCALAR;
            fdps->dofNames.push_back("-");
            fdps->unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE_DERIV_2);
            fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE_DERIV_2);
            fdps->data.resize(numDOFs * numTuples);
          }
        }

        if (((dsName == vx_ || dsName == vy_ || dsName == vz_) &&
             (requiredResults_[FLUIDMECH_VELOCITY] ||
              requiredResults_[NO_SOLUTION_TYPE])) ||
            (dsName == pres_ &&
             (requiredResults_[FLUIDMECH_PRESSURE] ||
              requiredResults_[NO_SOLUTION_TYPE])) ||
              (dsName == presD2_ &&
               (requiredResults_[FLUIDMECH_PRESSURE_DERIV_2] ||
                requiredResults_[NO_SOLUTION_TYPE])))
        {
          numDOFs = fdps->dofNames.size();
          fdps->data.resize(numDOFs * numTuples);

          //          std::cout << dsName << " " << numComps << " " << numTuples << " numNodes " << numNodes << std::endl;

          UInt offset = 0;
          if(numComps == 1) 
          {
            if (dsName == vx_)
              offset = 0;
            if (dsName == vy_)
              offset = 1;
            if (dsName == vz_)
              offset = 2;
          }

          /* copy the fluid velocity values */
          for(UInt i=0; i<numTuples; i++)
          {
            //            UInt node = actualRegionNodes_[actRegion][i];
            //            UInt node = 0;
            //            if( numTuples != numNodes )
            //            {
            //              node = actRegion == 0 ? 0 : actualRegionNodes_[actRegion-1].size();
            //            }
            UInt in = actualRegionNodes_[actRegion][i];
            UInt out = actualRegionNodesToNewIdx_[actRegion][in];
            
            //            std::cout << "in: " << in << " -> " << out << std::endl;
            
#define EXTRACT_DATA_ARRAY(TYPE) \
            { \
              vtkArrayIteratorTemplate<TYPE>* TYPEIt = NULL; \
              TYPEIt = static_cast< vtkArrayIteratorTemplate<TYPE>* >(dataIt); \
              for(UInt j=0; j<numComps; j++) { \
                TYPE val = TYPEIt->GetValue(in*numComps+j); \
                fdps->data[out*numDOFs+offset+j] = val; \
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
      }

      if( numTuples != numNodes )
      {
        node += numTuples;
      }

      iter->GoToNextItem();
    }

    iter->Delete();
    c2p->Delete();
  }

  //! get user data from file reader
  void FileReader_EnSight::GetUserData(std::map<std::string, std::string>& userData)
  {
    std::ifstream fin;
    std::stringstream sstr;
    vtkGenericEnSightReader* reader = dynamic_cast<vtkGenericEnSightReader*>(reader_);

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
  }

  void FileReader_EnSight::InitElemNodeMapping()
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
        unstrucElemNodeMapping_[et][0] = 0;
        unstrucElemNodeMapping_[et][1] = 1;
        unstrucElemNodeMapping_[et][2] = 2;
        unstrucElemNodeMapping_[et][3] = 3;
        unstrucElemNodeMapping_[et][4] = 4;
        unstrucElemNodeMapping_[et][5] = 5;
        break;
      case Elem::HEXA8:
        //FLUENT  & KMT ensight
        //unstrucElemNodeMapping_[et][0] = 4;
        //unstrucElemNodeMapping_[et][1] = 5;
        //unstrucElemNodeMapping_[et][2] = 6;
        //unstrucElemNodeMapping_[et][3] = 7;
        //unstrucElemNodeMapping_[et][4] = 0;
        //unstrucElemNodeMapping_[et][5] = 1;
        //unstrucElemNodeMapping_[et][6] = 2;
        //unstrucElemNodeMapping_[et][7] = 3;

        //uniformElemNodeMapping_[et][0] = 0;
        //uniformElemNodeMapping_[et][1] = 1;
        //uniformElemNodeMapping_[et][2] = 3;
        //uniformElemNodeMapping_[et][3] = 2;
        //uniformElemNodeMapping_[et][4] = 4;
        //uniformElemNodeMapping_[et][5] = 5;
        //uniformElemNodeMapping_[et][6] = 7;
        //uniformElemNodeMapping_[et][7] = 6;

        //FOAM
         uniformElemNodeMapping_[et][0] = 0;
         uniformElemNodeMapping_[et][1] = 1;
         uniformElemNodeMapping_[et][2] = 3;
         uniformElemNodeMapping_[et][3] = 2;
         uniformElemNodeMapping_[et][4] = 4;
         uniformElemNodeMapping_[et][5] = 5;
         uniformElemNodeMapping_[et][6] = 7;
         uniformElemNodeMapping_[et][7] = 6;
         
         unstrucElemNodeMapping_[et][0] = 0;
         unstrucElemNodeMapping_[et][1] = 1;
         unstrucElemNodeMapping_[et][2] = 2;
         unstrucElemNodeMapping_[et][3] = 3;
         unstrucElemNodeMapping_[et][4] = 4;
         unstrucElemNodeMapping_[et][5] = 5;
         unstrucElemNodeMapping_[et][6] = 6;
         unstrucElemNodeMapping_[et][7] = 7;
        break;

      default:
        break;
      }
    }
    
  }

} // end of namespace
