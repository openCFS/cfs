
#include "vtkCFSReader.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkInformationStringKey.h"
#include "vtkObjectFactory.h"
#include "vtkErrorCode.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"

#include "vtkStreamingDemandDrivenPipeline.h"


#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkQuadraticEdge.h>
#include <vtkTriangle.h>
#include <vtkQuadraticTriangle.h>
#include <vtkQuad.h>
#include <vtkQuadraticQuad.h>
#include <vtkBiQuadraticQuad.h>
#include <vtkTetra.h>
#include <vtkQuadraticTetra.h>
#include <vtkHexahedron.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkTriQuadraticHexahedron.h>
#include <vtkPyramid.h>
#include <vtkQuadraticPyramid.h>
#include <vtkWedge.h>
#include <vtkQuadraticWedge.h>
#include <vtkBiQuadraticQuadraticWedge.h>


#include <sstream>
#include <cmath>
#include <boost/format.hpp>


// **************************************************
//  For debugging, please define macro in next line
// **************************************************
#define DEBUG 0 /* set to 1 for debugging */


#if DEBUG >  0 
#if _WIN32
// Note: Under windows we normally do not get access to the 
// standard streams. Thus we use the vtk / paraview output window.
#define DBG_OUT(arg) {                                                        \
vtkOStreamWrapper::EndlType endl;                                             \
    vtkOStreamWrapper::UseEndl(endl);                                         \
    vtkOStrStreamWrapper vtkmsg;                                              \
    vtkmsg << "Debug: In " __FILE__ ", line " << __LINE__ << "\n"             \
           << this->GetClassName() << " (" << this << "): " arg  << "\n\n";   \
    vtkOutputWindowDisplayDebugText(vtkmsg.str());                            \
	vtkmsg.rdbuf()->freeze(0);   }
#else
#define DBG_OUT(arg) std::cerr << "line " << __LINE__ << ": " << arg << std::endl;
#endif
#else
#define DBG_OUT(arg)
#endif

vtkStandardNewMacro(vtkCFSReader)

const double PI = std::acos( -1.0 );

void vtkCFSReader::SetMultiSequenceStep( int step  ) {
  DBG_OUT("Setting multisequence step " << step << " for a valid range of "
           << this->MultiSequenceRange[0]  << ", " << this->MultiSequenceRange[1] );
  
  if( step == this->MultiSequenceStep )
    return;
  if( step > this->MultiSequenceRange[1] || 
      step < this->MultiSequenceRange[0] ) {
    DBG_OUT("Please enter a valid multisequence step between "
        << this->MultiSequenceRange[0] << " and " 
        << this->MultiSequenceRange[1] << "!");
    vtkErrorMacro(<< "Please enter a valid multisequence step between "
                  << this->MultiSequenceRange[0] << " and " 
                  << this->MultiSequenceRange[1] << "!");
  }
 
  this->MultiSequenceStep = step;
  this->msStepChanged = true;
  
  // Important: Set state to "modified", so that the result information
  // gets updated
  this->Modified();
  DBG_OUT("Step has changed");
}

void vtkCFSReader::SetTimeStep (unsigned int step) {
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting TimeStep to " << step);
  DBG_OUT("=> Set Time Step to " << step);
  DBG_OUT("\tPrevious step was " << this->TimeStep)
  DBG_OUT("\tSize of stepVals is" << this->stepVals.size())
  
  // security check: immediately leave, if no time / frequency steps are available
  if( this->stepVals.size() == 0 )
    return;

  if (this->TimeStep != step && this->stepVals.size() >= step ) 
  { 
    this->TimeStep = step; 
    this->timeFreqVal = stepVals[step-1];
    if( harmonicData )
      PrettyPrintFreq(this->timeFreqVal, this->timeFreqValStr );
    else 
      PrettyPrintTime(this->timeFreqVal, this->timeFreqValStr );
  
    
    // update pipeline
    this->Modified();
  }
}

const char* vtkCFSReader::GetTimeFreqValStr() {
  return timeFreqValStr.c_str();
}

void vtkCFSReader::SetAddDimensionsToArrayNames(int flag) {
  this->AddDimensionsToArrayNames = flag;

  // In addition trigger resetting the data value arrays 
  this->resetDataArrays = true;
  
  // update pipeline
  this->Modified();
}

void vtkCFSReader::SetHarmonicDataAsModeShape(int flag) {
  this->HarmonicDataAsModeShape = flag;
  // In addition trigger resetting the data value arrays 
  this->resetDataArrays = true;
  
  // update pipeline
  this->Modified();
}

void vtkCFSReader::SetFillMissingResults(int flag) {
  this->FillMissingResults = flag;
  
  // In addition trigger resetting the data value arrays 
  this->resetDataArrays = true;  
  
  // update pipeline
  this->Modified();
}

int vtkCFSReader::GetNumberOfRegionArrays() {
  return regionNames_.size();
}

int vtkCFSReader::GetRegionArrayStatus( const char * name ){
  std::map<std::string,int>::const_iterator it = regionSwitch.find(name);
  if( it == regionSwitch.end() ) {
    vtkErrorMacro(<< "Region '" << name << "' not found.");
  }
  return it->second;
  
}

void vtkCFSReader::SetRegionArrayStatus( const char * name, int status ){
  regionSwitch[name] = status;
  
  // vtkDebugMacro("Region '" << name << "' status " << status);
  
  // IMPORTANT: Let VTK know, that the pipeline needs an update
  // (otherwise nothing would happen)
  this->Modified();
  
  // In addition, set flag to false, that regions have to be updated
  this->activeRegionsChanged = true;
}
const char * vtkCFSReader::GetRegionArrayName( int index) {
  return regionNames_[index].c_str();
}

int vtkCFSReader::GetNumberOfNamedNodeArrays() {
  return namedNodeNames_.size();
}

int vtkCFSReader::GetNamedNodeArrayStatus( const char * name ){
  std::map<std::string,int>::const_iterator it = namedNodeSwitch.find(name);
  if( it == namedNodeSwitch.end() ) {
    vtkErrorMacro(<< "Named nodes '" << name << "' not found.");
  }
  return it->second;
}

void vtkCFSReader::SetNamedNodeArrayStatus( const char * name, int status ){
  namedNodeSwitch[name] = status;
  this->Modified();
  // In addition, set falg to false, that regions have to be updated
  this->activeRegionsChanged = true;
}

const char * vtkCFSReader::GetNamedNodeArrayName( int index) {
  return namedNodeNames_[index].c_str();
}

int vtkCFSReader::GetNumberOfNamedElemArrays() {
  return namedElemNames_.size();
}

int vtkCFSReader::GetNamedElemArrayStatus( const char * name ){
  std::map<std::string,int>::const_iterator it = namedElemSwitch.find(name);
  if( it == namedElemSwitch.end() ) {
    vtkErrorMacro(<< "Named elems '" << name << "' not found.");
  }
  return it->second;
}

void vtkCFSReader::SetNamedElemArrayStatus( const char * name, int status ){
  namedElemSwitch[name] = status;
  this->Modified();
  // In addition, set falg to false, that regions have to be updated
   this->activeRegionsChanged = true;
}

const char * vtkCFSReader::GetNamedElemArrayName( int index) {
  return namedElemNames_[index].c_str();
}

//----------------------------------------------------------------------------
vtkCFSReader::vtkCFSReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  
  // Define Debugging
  this->SetDebug(0);
  this->mbDataSet = NULL;

  this->isInitialized = false;
  this->hdf5InfoRead = false;
  this->activeRegionsChanged = true;
  this->resetDataArrays = false;
  this->timeFreqVal = 0.0;
  this->timeFreqValStr = "0.0";
  this->TimeStep = 1;
  this->MultiSequenceStep = 1;
  this->MultiSequenceRange[0] = 1;
  this->MultiSequenceRange[1] = 1;
  this->TimeStepNumRange[0] = 1;
  this->TimeStepNumRange[1] = 1;
  this->TimeStepValuesRange[0] = 0.0;
  this->TimeStepValuesRange[1] = 1.0;
  this->msStepChanged = true;
  this->HarmonicDataAsModeShape = false;
  this->harmonicData = false;
  this->FillMissingResults = 0;
}

//----------------------------------------------------------------------------
vtkCFSReader::~vtkCFSReader()
{
  if (this->FileName)
    {
    delete [] this->FileName;
    }
}

//----------------------------------------------------------------------------

#if 0
int vtkCFSReader::CanReadFile(const char* fname) {
  
  // Try to "load" the file (which internally just opens the mesh group)
  // to see, if this is a true CFS HDF5 file. 
  std::string fileName = std::string(fname);
  H5CFS::Hdf5Reader temp; 
  try {
    temp.LoadFile( fileName );
    // If this succeeds, we assume that this is a true CFS HDF5 file.
    return true;
  } catch (...) {
    return false;
  }
}
#endif

//----------------------------------------------------------------------------
int vtkCFSReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{

  // Load hdf file (only done first time
  if( !hdf5InfoRead ) {
    try {
      in_.LoadFile(this->FileName);
      in_.GetAllRegionNames( regionNames_ );
      in_.GetNodeNames( namedNodeNames_ );
      in_.GetElemNames( namedElemNames_ );
      
      for( unsigned int i = 0, n = regionNames_.size(); i < n; i++ ) {
        std::map<std::string,int>::iterator it = regionSwitch.find(regionNames_[i]);
        if( it == regionSwitch.end() ) {
          regionSwitch[regionNames_[i]] = 1;
        }

        // vtkDebugMacro("**** " << regionNames_[i] << " " << regionSwitch[regionNames_[i]] << " ****");
      }

      for( unsigned int i = 0, n = namedNodeNames_.size(); i < n; i++ ) {
        std::map<std::string,int>::iterator it = namedNodeSwitch.find(namedNodeNames_[i]);
        if( it == namedNodeSwitch.end() ) {
          namedNodeSwitch[namedNodeNames_[i]] = 0;
        }

        // vtkDebugMacro("+++++ " << namedNodeNames_[i] << " " << namedNodeSwitch[namedNodeNames_[i]] << " ++++");
      }

      for( unsigned int i = 0, n = namedElemNames_.size(); i < n; i++ ) {
        std::map<std::string,int>::iterator it = namedElemSwitch.find(namedElemNames_[i]);
        if( it == namedElemSwitch.end() ) {
          namedElemSwitch[namedElemNames_[i]] = 0;
        }

        // vtkDebugMacro("---- " << namedElemNames_[i] << " " << namedElemSwitch[namedElemNames_[i]] << " ----");
      }
    } catch( std::exception& ex ) {
      vtkErrorMacro(<< "Caught exception when trying to read file '"
                    << this->FileName << "' : '" 
                    << ex.what() << "'");
    } catch(std::string& ex ) {
      vtkErrorMacro(<< "Caught exception when trying to read file '"
                    << this->FileName << "' : '" 
                    << ex  << "'");
    }
    
    hdf5InfoRead = true;
    isInitialized = false;
  } // info read
  
 DBG_OUT(  "============================\n"
            << " Calling RequestInformation\n"
            << "============================\n")
 try {
   
   
   // ----------------------------------
   //  Check for new multisequence step
   // ----------------------------------
   if( this->msStepChanged ) {
     DBG_OUT("===> New MS Step: " << this->MultiSequenceStep );
     // Request number of multiSequenceSteps
     std::map<unsigned int,H5CFS::AnalysisType> analysis;
     std::map<unsigned int, unsigned int> numSteps;
     in_.GetNumMultiSequenceSteps( analysis, numSteps, false );

     // adjust range for multisequence steps
     if( analysis.size() ) {
       this->MultiSequenceRange[0] = analysis.begin()->first;
       this->MultiSequenceRange[1] = analysis.rbegin()->first;
       
       // ensure, that the current multisequence step is within the
       // valid range. Otherwise, set it to the smallest one defined
       if( this->MultiSequenceStep < this->MultiSequenceRange[0] ) {
         this->MultiSequenceStep = this->MultiSequenceRange[0];
       }
     }

     // check, if any result at all are present
     // Run over all resulttypes
     std::set<double> globStepVals;
     std::set<unsigned int> globStepNums;
     std::vector<shared_ptr<H5CFS::ResultInfo> >infos;

     if( analysis.size() > 0 ) {

       // store current analsis type
       this->analysisType = analysis[MultiSequenceStep]; 
       if( this->analysisType == H5CFS::HARMONIC || 
           this->analysisType == H5CFS::EIGENFREQUENCY) {
         harmonicData = true;
       } else {
         harmonicData = false;
       }

       in_.GetResultTypes( MultiSequenceStep, infos );
       for( unsigned int i = 0; i < infos.size(); i++ ) {
         std::map<unsigned int, double> actStepVals;
         in_.GetStepValues( MultiSequenceStep, infos[i],
                          actStepVals, false );
         std::map<unsigned int, double>::const_iterator it;
         for( it = actStepVals.begin(); it != actStepVals.end(); it++ ) {
           globStepNums.insert( it->first );
           globStepVals.insert( it->second );
           DBG_OUT("Timestep: " << it->first << " = " << it->second);
         }
       } // loop over resulttypes

       this->stepNums.resize( globStepNums.size() );
       std::copy( globStepNums.begin(), globStepNums.end(),
                  this->stepNums.begin() );
       this->stepVals.resize( globStepVals.size() );
       std::copy( globStepVals.begin(), globStepVals.end(), this->stepVals.begin() );
     } // if analysis size   
     this->NumberOfTimeSteps = this->stepVals.size();
     DBG_OUT( "Simulation has " << NumberOfTimeSteps << " time / frequency steps");
     this->TimeStepNumRange[1] = this->NumberOfTimeSteps;
     if(this->NumberOfTimeSteps) 
     {
       this->TimeStepValuesRange[0] = this->stepVals[0];
       this->TimeStepValuesRange[1] = this->stepVals[NumberOfTimeSteps-1];
     }
     else
     {
       this->TimeStepValuesRange[0] = 0;
       this->TimeStepValuesRange[1] = 0;
     }
   } // if multisequence step has changed
   
 } catch(std::string& ex ) {
   vtkErrorMacro(<< "Caught exception when trying to read file '"
                 << this->FileName << "' : '" 
                 << ex  << "'");
 } catch(const char *  ex ) {
    vtkErrorMacro(<< "Caught exception when trying to read file '"
                  << this->FileName << "' : '" 
                  << ex  << "'");
 }
  
  if( this->NumberOfTimeSteps > 0 ) {
    if(harmonicData && HarmonicDataAsModeShape) {
       // From harmonic data we can produce a continuous time data set which ranges from 0 to one
       // All other kinds of data are discrete in time/freq.

      this->TimeStepValuesRange[0] = 0;
      this->TimeStepValuesRange[1] = 1.0;
      
      outputVector->GetInformationObject(0)->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    } else {    
      outputVector->GetInformationObject(0)->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
                                                 &(this->stepVals[0]), this->NumberOfTimeSteps);
    }
    
    outputVector->GetInformationObject(0)->
    Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), TimeStepValuesRange, 2);
  }

  return 1;
}
//----------------------------------------------------------------------------
int vtkCFSReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{


  DBG_OUT( "============================\n"
           << " Calling RequestData\n"
           << "============================\n" );

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut
  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(
     outInfo->Get(vtkMultiBlockDataSet::DATA_OBJECT()));


  // ===================
  // Time Stepping Stuff
  // (taken from OpenFAOMReader)
  // ===================
  DBG_OUT( "Current time step before is " << TimeStep );
  
  // Attention: Only change value of pipeline / step values, if 
  // there are any results present at all.
  if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
      && stepVals.size() > 0 )
     {
     requestedTimeValue
       = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
     int nSteps
       = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
     double* rsteps
       = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
     int timeI = 0;
     while(timeI < nSteps - 1 && rsteps[timeI] < requestedTimeValue)
       {
       timeI++;
       }
     DBG_OUT( "Current time step is " << TimeStep );
     
     if(harmonicData && HarmonicDataAsModeShape) {
       outInfo->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
     } else {
       if(NumberOfTimeSteps) 
       {
         outInfo->Set(vtkDataObject::DATA_TIME_STEP(), stepVals[timeI]);
         this->TimeStep = timeI+1;
         this->timeFreqVal = stepVals[timeI];
         if( harmonicData )
           PrettyPrintFreq(this->timeFreqVal, this->timeFreqValStr );
         else 
           PrettyPrintTime(this->timeFreqVal, this->timeFreqValStr );
       }
       else 
       {
         outInfo->Set(vtkDataObject::DATA_TIME_STEP(), 0.0);
       }       
     }
   }

  // Read in data from file
  this->ReadFile(output);

  return 1;
}

//----------------------------------------------------------------------------
void vtkCFSReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " 
     << (this->FileName ? this->FileName : "(none)") << "\n";
}

//----------------------------------------------------------------------------
void vtkCFSReader::ReadFile(vtkMultiBlockDataSet *output)
{
  

  DBG_OUT( "==========================================\n"
           << " R E A D I N G   I N   F I L E \n"
           << "==========================================\n" );
  // ================================
  // The first time, the reader is initialized, we have to
  // fill the  multiblock data set
  // ================================
  if( !this->isInitialized ) {
    DBG_OUT( "==== Reading Grid for First Time ===");
    mbDataSet = vtkMultiBlockDataSet::New();
    
    mbDataSet->ShallowCopy(output);

    for( unsigned  int i = 0, n = regionNames_.size(); i < n; i++ ) {
      vtkUnstructuredGrid * newGrid = vtkUnstructuredGrid::New();
      mbDataSet->SetBlock( i, newGrid );
      newGrid->Delete();
    }

    // initialize mapping from 
    // (region,globalNodeNumber)->regionLocalNodeNumber
    nodeMap_.resize( regionNames_.size() );

    // read nodal connectivity
    this->ReadNodes( mbDataSet );

    // read element definition  
    this->ReadCells( mbDataSet );

    // Now the grid is finalized and we store it initially also 
    // to the activeSet
    this->mbActiveDataSet = vtkMultiBlockDataSet::New();
    this->mbActiveDataSet->ShallowCopy(this->mbDataSet);
    
    // now it should be definitly initialized
    this->isInitialized = true;
    
  } 

  // If region status has changed, we have to perform an update
  // of the active regions
  if( this->activeRegionsChanged || this->resetDataArrays ) {
    this->UpdateActiveRegions();
  }

  // only read data, if time/frequency steps are defined
  if( this->NumberOfTimeSteps > 0 ) {
    // read nodal data
    DBG_OUT( "Reading nodal data");
    this->ReadNodeCellData( mbActiveDataSet, true);

    // read cell data
    DBG_OUT( "Reading cell data");
    this->ReadNodeCellData( mbActiveDataSet, false);
  }

  // Use the active set as base for the current set
  output->ShallowCopy(mbActiveDataSet);
}

//----------------------------------------------------------------------------
void vtkCFSReader::UpdateActiveRegions() {
  
  // Idea:
  // Create new multiblock data set, based on the internal stored one and copy only those
  // grid-blocs, which are active.
  
  vtkMultiBlockDataSet *newSet  = vtkMultiBlockDataSet::New();
  
  // run over all regionNames_
  unsigned int index = 0;
  for(unsigned int i = 0, n = regionNames_.size(); i < n; i++ ) {

    // vtkDebugMacro("#### " << regionNames_[i] << " " << regionSwitch[regionNames_[i]] << " ####")
    // check, if region is active
    if( regionSwitch[regionNames_[i]] != 0 ) {

      // create new grid with output and make shallow copy
      vtkUnstructuredGrid * newGrid = vtkUnstructuredGrid::New();
      newGrid->ShallowCopy(mbDataSet->GetBlock(i));
      
      // try to find the corresponding block of this region in the current
      // active set in order to copy its data to this block.
      unsigned int numActBlocks = this->mbActiveDataSet->GetNumberOfBlocks();
      unsigned int activeIndex = 0;
      bool found = false;
      for( unsigned int j = 0; j < numActBlocks; ++j ) {
        const char * blockName = 
            this->mbActiveDataSet->GetMetaData(j)->Get(vtkCompositeDataSet::NAME());
        if( blockName != NULL ) {
          if( std::string(blockName) == regionNames_[i]) {
            found = true;
            activeIndex = j;
          }
        }
      }
      
      if( found && !this->resetDataArrays) {
        vtkUnstructuredGrid * actGrid = 
            vtkUnstructuredGrid::SafeDownCast(this->mbActiveDataSet->GetBlock(activeIndex));
        newGrid->GetPointData()->ShallowCopy( actGrid->GetPointData());
        newGrid->GetCellData()->ShallowCopy( actGrid->GetCellData());
      }
      
      // Get the array containing the actual node numbers for the current region
      newGrid->GetPointData()->SetActiveScalars("origNodeNums");
      vtkUnsignedIntArray *origNodeNums;
      origNodeNums = vtkUnsignedIntArray::SafeDownCast(newGrid->GetPointData()->GetScalars());
      newGrid->GetPointData()->SetActiveScalars("");
      
      // Iterate over all named node arrays and find out the named nodes which
      // are contained in the current region. For these nodes a value of 1 will
      // be set in the nodeNumsArray. For all other nodes 0 will be set. This makes
      // it possible to apply a glyph filter to the named node array and set the
      // scaling  mode to scalar or to create threshold based selections.
      for(unsigned int j = 0, m = namedNodeNames_.size(); j < m; j++ ) {
        if( namedNodeSwitch[namedNodeNames_[j]] != 0 ) {
          // Read named node array from HDF5 file
          std::vector<unsigned int> nodes;
          in_.GetNamedNodes(namedNodeNames_[j], nodes);

          // Create new attribute by copying the origNodeNums array
          vtkUnsignedIntArray *nodeNumsArray = vtkUnsignedIntArray::New();
          nodeNumsArray->DeepCopy(origNodeNums);
          std::string arrayName = namedNodeNames_[j];
          arrayName += " (named nodes)";
          nodeNumsArray->SetName( arrayName.c_str() );
          newGrid->GetPointData()->AddArray(nodeNumsArray);

          // Iterate over origNodeNums and find out which named nodes are contained in it.
          for(unsigned int k = 0, o = origNodeNums->GetNumberOfTuples(); k < o; k++ ) {
            std::vector<unsigned int>::const_iterator it;
            it = std::find(nodes.begin(), nodes.end(), origNodeNums->GetValue(k));
            nodeNumsArray->SetValue(k, (it != nodes.end()));
          }
        }
      }

      // Get the array containing the actual node numbers for the current region
      newGrid->GetCellData()->SetActiveScalars("origElemNums");
      vtkUnsignedIntArray *origElemNums;
      origElemNums = vtkUnsignedIntArray::SafeDownCast(newGrid->GetCellData()->GetScalars());
      newGrid->GetCellData()->SetActiveScalars("");
      
      // Iterate over all named element arrays and find out the named elements which
      // are contained in the current region. For these elems a value of 1 will
      // be set in the elemNumsArray. For all other nodes 0 will be set. This makes
      // it possible to apply to create threshold based selections.
      for(unsigned int j = 0, m = namedElemNames_.size(); j < m; j++ ) {
        if( namedElemSwitch[namedElemNames_[j]] != 0 ) {
          // Read named node array from HDF5 file
          std::vector<unsigned int> elems;
          in_.GetNamedElems(namedElemNames_[j], elems);

          // Create new attribute by copying the origNodeNums array
          vtkUnsignedIntArray *elemNumsArray = vtkUnsignedIntArray::New();
          elemNumsArray->DeepCopy(origElemNums);
          std::string arrayName = namedElemNames_[j];
          arrayName += " (named elems)";
          elemNumsArray->SetName( arrayName.c_str() );
          newGrid->GetCellData()->AddArray(elemNumsArray);

          // Iterate over origNodeNums and find out which named nodes are contained in it.
          for(unsigned int k = 0, o = origElemNums->GetNumberOfTuples(); k < o; k++ ) {
            std::vector<unsigned int>::const_iterator it;
            it = std::find(elems.begin(), elems.end(), origElemNums->GetValue(k));
            elemNumsArray->SetValue(k, (it != elems.end()));
          }
        }
      }

      newSet->SetBlock(index, newGrid);
      
      // Tag new block with region name
      newSet->GetMetaData(index)->Set(vtkCompositeDataSet::NAME(), regionNames_[i].c_str());
      index++;

      newGrid->Delete();
    }
  }
  
  // in the end, replace the active set by the current one
  this->mbActiveDataSet->ShallowCopy(newSet);
  newSet->Delete();

  // in the end, re-set the status of the update flag
  this->activeRegionsChanged = false;
}

//----------------------------------------------------------------------------
void vtkCFSReader::ReadNodeCellData(vtkMultiBlockDataSet *output, bool isNode)
{
  DBG_OUT("=========================");
  DBG_OUT( " Reading Node/Cell Data" );
  DBG_OUT("=========================");
  H5CFS::EntityType entType = H5CFS::ELEMENT;
  static const double H180DEG_OVER_PI = 180.0 / PI;

  if( isNode ) {
    DBG_OUT("EntityType is NODE");
    entType = H5CFS::NODE;
  } else {
    DBG_OUT("EntityType is ELEMENT");
  }

  try {
    // get result infos for 1ms step and step 1
    std::vector<shared_ptr<H5CFS::ResultInfo> > infos;
    in_.GetResultTypes( MultiSequenceStep, infos );
    DBG_OUT("Reading " << infos.size() << " result types");
    
    // Read all result types occuring in the current multisequence step
    std::set<std::string> resultNames;
    std::map<std::string, shared_ptr<H5CFS::ResultInfo> > nameInfo;
    for( unsigned int i = 0; i < infos.size(); ++i ) {
      std::string & resultName = infos[i]->name;
      resultNames.insert( resultName );
      nameInfo[resultName] = infos[i];
    }

    std::vector<std::vector<shared_ptr<H5CFS::ResultInfo> > > regionInfos;
    shared_ptr<H5CFS::Result> result;
    regionInfos.resize( regionNames_.size() );

    // ----------------------------
    //  Loop over all result names
    // -----------------------------
    std::set<std::string>::const_iterator it = resultNames.begin();
    for( ; it != resultNames.end(); ++it ) {
      std::string  resultName = *it;

      // We distinguish between volume and surface elements, but VTK does not.
      H5CFS::EntityType actualEntType = nameInfo[resultName]->listType;
      switch(nameInfo[resultName]->listType) 
      {
      case H5CFS::SURF_ELEM:
        actualEntType = H5CFS::ELEMENT;
        break;
      default:
        break;
      }        

      // check if requested entity type (node, cell) matches result one
      if( ( isNode && (actualEntType != H5CFS::NODE ) )
          || (!isNode && (actualEntType != H5CFS::ELEMENT ) ) )
        continue;
      
      DBG_OUT( "\tReading Result : '" << resultName << "'" );
      
      // ------------------------------
      //  Loop over all active regions
      // ------------------------------
      unsigned int blockIndex = 0;
      for( unsigned int iRegion = 0; iRegion < regionNames_.size(); iRegion++ ) {

        // variable used to write the result 
        std::string  writeName = resultName;
        
        // check, if region is active. If not, just leave
        if( regionSwitch[regionNames_[iRegion]] == 0)
          continue;

        DBG_OUT( "\t\tReading Region: '" << regionNames_[iRegion] << "'");

        // get grip of grid of current region
        vtkUnstructuredGrid *actGrid = 
            vtkUnstructuredGrid::SafeDownCast( output->GetBlock(blockIndex++) );
        
        
        // -------------------------------------------------
        //  Find resultinfo for this result and this region
        // -------------------------------------------------
        bool found = false;
        shared_ptr<H5CFS::ResultInfo> actInfo;
        for( unsigned int i = 0; i < infos.size(); ++i ) {
          switch(infos[i]->listType) 
          {
          case H5CFS::SURF_ELEM:
            actualEntType = H5CFS::ELEMENT;
            break;
          default:
            break;
          }  

          if( infos[i]->name == resultName &&
              infos[i]->listName == regionNames_[iRegion] &&
              actualEntType == entType ) {
            actInfo = infos[i];
            found = true;
            DBG_OUT("\t\t=> Found");
            break;
          }
        }
        
        // if result for this region was not found and we should not generate
        // "empty" datasets for missing result, we just leave
        if( found == false ) {
          DBG_OUT("\t\t=> Not Found!");
          if (!FillMissingResults ) {
            DBG_OUT("Leaving");
            continue;
          } else {
            actInfo = nameInfo[resultName];
            DBG_OUT("\t\t => Filling missing results with zeros");
          }
        } 
        

        std::vector<unsigned int> entities;
        in_.GetEntities( entType, regionNames_[iRegion], entities );
        unsigned int numEntities = entities.size();
        unsigned int numDofs = actInfo->dofNames.size();
        std::vector<double> * ptRealVals = NULL;
        std::vector<double> * ptImagVals = NULL;
        bool isComplex = false;
        std::vector<double> dummyVals;
        H5CFS::EntryType entryType = actInfo->entryType;
        if( found ) {
          result = shared_ptr<H5CFS::Result>(new H5CFS::Result());
          result->resultInfo = actInfo;
          try {
            in_.GetResult( MultiSequenceStep, stepNums[TimeStep-1], result);
          } catch (...) {
            // In case the result is not defined in every step, we simply
            // keep the result from the last step.
            break;
          }
          isComplex = result->isComplex;
          numDofs = actInfo->dofNames.size();
          ptRealVals = &(result->realVals);
          ptImagVals = &(result->imagVals);
        } else {
          dummyVals.resize( numEntities * numDofs);
          ptRealVals = &dummyVals;
          ptImagVals = &dummyVals;
          // NOTE: we have to find a different strategy how to 
          // find out about complex valued results
          isComplex = (this->analysisType == H5CFS::HARMONIC);
        }
       
        std::vector<double> & realVals = *ptRealVals;

        vtkDoubleArray *vals1 = NULL, *vals2 = NULL;

        // REAL-VALUED
        DBG_OUT( "\t\tCopying data to vtk array");
        if( !isComplex ) {
          if( ! (harmonicData && HarmonicDataAsModeShape) ) {
            DBG_OUT( "\t\tResult is REAL-VALUED");

            if(AddDimensionsToArrayNames) {
              writeName += " [";
              writeName += actInfo->unit;
              writeName += "]";
            }
            vals1 = SaveToArray( realVals, numDofs, numEntities,
                                 entryType, writeName );
          } else { // We obviously have a real-valued eigenvalue solution
            //std::string freqNormed;
            unsigned int numEntries = realVals.size();

            std::string writeName = actInfo->name;
            writeName += "-mode";
            if(AddDimensionsToArrayNames) {
              writeName += " [";
              writeName += actInfo->unit;
              writeName += "]";
            }

              double reRot = cos(-2*PI*requestedTimeValue);

            for(unsigned int i = 0; i < numEntries; i++ ) {
              realVals[i] *= reRot;
            }

            vals1 = SaveToArray( realVals, numDofs, numEntities,
                                 entryType, writeName );

          }
        } else {
          std::vector<double> & imagVals = *ptImagVals;
          // REAL-IMAG
          if( !HarmonicDataAsModeShape ) {
            DBG_OUT( "\t\tResult is COMPLEX");
            if( ComplexMode == 0 ) {
              std::string realName = actInfo->name + "-real";
              std::string imagName = actInfo->name + "-imag";
              if(AddDimensionsToArrayNames) {
                realName += " [";
                realName += actInfo->unit;
                realName += "]";
                imagName += " [";
                imagName += actInfo->unit;
                imagName += "]";
              }
              vals1 = SaveToArray( realVals, numDofs, numEntities,
                                   entryType, realName );
              vals2 = SaveToArray( imagVals, numDofs, numEntities,
                                   entryType, imagName );
            } else {
              // AMPL-PHASE
              unsigned int numEntries = realVals.size();
              std::string amplName = actInfo->name + "-ampl";
              std::string phaseName = actInfo->name + "-phase";
              if(AddDimensionsToArrayNames) {
                amplName += " [";
                amplName += actInfo->unit;
                amplName += "]";
                phaseName += " [deg]";
              }
              std::vector<double> ampl(numEntries), phase(numEntries);
              for(unsigned int i = 0; i < numEntries; i++ ) {
                ampl[i] = hypot( realVals[i], imagVals[i] );
                phase[i] = (std::abs(imagVals[i]) > 1e-16) ?                   
                    std::atan2( imagVals[i], realVals[i] ) * H180DEG_OVER_PI : 
                    ( realVals[i] < 0.0 ) ? 180 : 0 ; 
              }
              vals1 = SaveToArray( ampl, numDofs, numEntities, entryType, amplName );
              vals2 = SaveToArray( phase, numDofs, numEntities, entryType, phaseName );
            }
          } else { // !HarmonicDataAsModeShape

            // Compute time continuous field from harmonic data by multiplying complex
            // result field with rotating phasor

            unsigned int numEntries = realVals.size();
            std::stringstream sstr;

            //sstr << actInfo->name << "-mode at " << freqNormed;
            sstr << actInfo->name << "-mode";
            std::string writeName = sstr.str();
            sstr.clear(); sstr.str("");
            //sstr << actInfo->name << "-ampl at " << freqNormed;
            sstr << actInfo->name << "-ampl";
            std::string amplName = sstr.str();
            if(AddDimensionsToArrayNames) {
              writeName += " [";
              writeName += actInfo->unit;
              writeName += "]";
              amplName += " [";
              amplName += actInfo->unit;
              amplName += "]";
            }
            std::vector<double> modeShape(numEntries);
            std::vector<double> ampl(numEntries);
            // Compute the rotation pointer from the requestedTimeValue which is in the
            // interval [0,1]
            double reRot = cos(2*PI*requestedTimeValue);
            double imRot = sin(2*PI*requestedTimeValue);
            for(unsigned int i = 0; i < numEntries; i++ ) {
              // Just compute the real part of the field
              modeShape[i] = realVals[i] * reRot - imagVals[i] * imRot;

              // Also compute the amplitude field
              ampl[i] = hypot( realVals[i], imagVals[i] );
            }
            vals1 = SaveToArray( modeShape, numDofs, numEntities, entryType, writeName );
            vals2 = SaveToArray( ampl, numDofs, numEntities, entryType, amplName );
          } // harmonic data as mode shape
        } // if isComplex

        // save values to grid
        if( entType == H5CFS::NODE) { 
          DBG_OUT( "\t\tPassing nodal result to grid");
          actGrid->GetPointData()->AddArray(vals1);
          if( vals2)
            actGrid->GetPointData()->AddArray(vals2);
        } else {
          DBG_OUT( "\t\tPassing element result to grid");
          actGrid->GetCellData()->AddArray(vals1);
          if( vals2 )
            actGrid->GetCellData()->AddArray(vals2);
        }
        // delete values (de-crement reference counter)
        DBG_OUT("\t\tCleaning up pointer data");
        vals1->Delete();
        if( vals2 )
          vals2->Delete();

      } // loop over active regions
    } // loop over result types
  } catch(std::exception& ex ) {
    vtkErrorMacro(<< "Caught exception when reading results: '" 
                  << ex.what() << "'");
  } catch(std::string& ex ) {
    vtkErrorMacro(<< "Caught exception when reading results: '" 
                  << ex << "'");
  }

}

vtkDoubleArray* vtkCFSReader::SaveToArray( const std::vector<double>& vals,
                                           unsigned int numDofs,
                                           unsigned int numEntities,
                                           H5CFS::EntryType entryType,
                                           const std::string& name ) {
  
  vtkDoubleArray * ret = vtkDoubleArray::New();
  
  // Case 1: numDofs is 1 (scalar) or >= 3 (real vector/tensor)
  if( numDofs == 1 || numDofs >= 3) {
    unsigned int nDofs = numDofs;
    bool isTensor = (entryType == H5CFS::TENSOR);

    if(isTensor) {
      nDofs = 6;
    }
    
    ret->SetNumberOfComponents( nDofs );
    ret->SetNumberOfTuples( numEntities);
    ret->SetName( name.c_str() );
    double * retPtr = ret->GetPointer(0);
    unsigned int numEntries = nDofs * numEntities;
    DBG_OUT("numEntries should be " << numEntries);
    DBG_OUT("length of vals array is " << vals.size());
    
    if(isTensor) {
      // Since VTK and ParaView have a different ordering of Voigt vectors,
      // we fill our components accordingly.
      switch(numDofs) {
        case 3:
          {
            unsigned idx;
            for( unsigned int iEnt = 0; iEnt < numEntities; iEnt++ ) { 
              idx = iEnt*nDofs;
              retPtr[idx+0] = vals[iEnt*3+0];
              retPtr[idx+1] = vals[iEnt*3+1];
              retPtr[idx+2] = 0.0;
              retPtr[idx+3] = vals[iEnt*3+2];
              retPtr[idx+4] = 0.0;
              retPtr[idx+5] = 0.0;
            }
          }
          break;
        case 6:
          {
            unsigned idx;
            for( unsigned int iEnt = 0; iEnt < numEntities; iEnt++ ) { 
              idx = iEnt*nDofs;
              retPtr[idx+0] = vals[idx+0];
              retPtr[idx+1] = vals[idx+1];
              retPtr[idx+2] = vals[idx+2];
              retPtr[idx+3] = vals[idx+5];
              retPtr[idx+4] = vals[idx+3];
              retPtr[idx+5] = vals[idx+4];
            }
          }
          break;
      }
    } else {
      for( unsigned int j = 0; j < numEntries; j++ ) { 
        retPtr[j] = vals[j];
      }
    }
  } else {
    // Case 2: numDofs (vector field in 2D)
    // -> add artificial 3rd component, which is 0
    ret->SetNumberOfComponents( 3 );
    ret->SetNumberOfTuples( numEntities );
    ret->SetName( name.c_str() );
    double * retPtr = ret->GetPointer(0);
    unsigned int index = 0;
    for( unsigned int iEnt = 0; iEnt < numEntities; iEnt++ ) { 
      retPtr[index++] = vals[iEnt*2+0];
      retPtr[index++] = vals[iEnt*2+1];
      retPtr[index++] = 0.0;
    }
  }
  return ret;
  
  // Alternative, slower way of filling a double array
  // =================================================
  //      std::vector<double> tuple (numDofs);
  //      for( unsigned int e = 0; e < entities.size(); e++ ) {
  //        for( unsigned int d = 0; d < numDofs; d++ ) {
  //          tuple[d] = result->realVals[numDofs*e+d];
  //        }
  //vals->SetTuple(e, &tuple[0]);
  //vals->SetArray(&(result->realVals[0]), 
  //               result->realVals.size(), 1);
  //}
}





//----------------------------------------------------------------------------
void vtkCFSReader::ReadCells(vtkMultiBlockDataSet *output)
{

  std::vector<H5CFS::ElemType> types;
  std::vector<std::vector<unsigned int> > connect;
  
  try {
    in_.GetElems( types, connect );
    vtkDebugMacro(<< "Read in the element definitions");

    // unsigned int numElems = connect.size();

    VTKCellType elemType = VTK_EMPTY_CELL;

    // loop over all regions
    for( unsigned int iRegion = 0, n = regionNames_.size(); 
    iRegion < n; iRegion++ ) {

      std::vector<unsigned int> regionElems;
      in_.GetElemsOfRegion( regionNames_[iRegion], regionElems);
      unsigned int numRegionElems = regionElems.size();

      //vtkDebugMacro(<< "Reading region '" << regionNames_[iRegion] << "'");

      vtkUnstructuredGrid *actGrid = 
        vtkUnstructuredGrid::SafeDownCast( output->GetBlock(iRegion) );
      actGrid->Allocate( numRegionElems );

      // Add a result vector containing the original element numbers
      // to the current unstructured grid.
      vtkUnsignedIntArray * origElemNums = vtkUnsignedIntArray::New();
      origElemNums->SetNumberOfValues( numRegionElems );
      origElemNums->SetName( "origElemNums" );
      actGrid->GetCellData()->AddArray(origElemNums);

      // Add a result vector containing the CFS++ element types
      // to the current unstructured grid.
      vtkUnsignedIntArray * elemTypes = vtkUnsignedIntArray::New();
      elemTypes->SetNumberOfValues( numRegionElems );
      elemTypes->SetName( "elemTypes" );
      actGrid->GetCellData()->AddArray(elemTypes);
      
      // loop over all elements
      for( unsigned int i = 0; i < numRegionElems; i++ ) {

        unsigned int actElem = regionElems[i] - 1;
        // vtkDebugMacro(<< "Treating element"
        //              << actElem );
        origElemNums->SetValue(i, actElem+1);
        elemTypes->SetValue(i, types[actElem]);
        std::vector<unsigned int>& conn = connect[actElem];
        
        // create correct vtk element
        elemType = GetCellIdType( types[actElem] );
        if( elemType == VTK_EMPTY_CELL ) {
          vtkErrorMacro(<<"FE type " << types[actElem] << " not implemented yet");
        }
        // set connectivity
        vtkIdType nodeIds[27];

	for(unsigned int j = 0, n = 27; j < n; j++) {
	  nodeIds[j] = 0;
	}
	
        for(unsigned int j = 0, n = conn.size(); j < n && conn[j] > 0; j++)
        {
          // vtkDebugMacro(<< "addding nodeNum" << conn[j] );
	  unsigned int connJ = conn[j];
	  unsigned int nm = nodeMap_[iRegion][connJ-1];
          nodeIds[j] = nm - 1;
        }

        if(elemType == VTK_TRIQUADRATIC_HEXAHEDRON) 
        {
          // top 
          //   7--14--6
          //   |      |
          //  15  25  13
          //   |      |
          //   4--12--5

          //   middle
          //  19--23--18
          //   |      |
          //  20  26  21
          //   |      |
          //  16--22--17

          //  bottom
          //   3--10--2
          //   |      |
          //  11  24  9 
          //   |      |
          //   0-- 8--1

          nodeIds[20] = nodeMap_[iRegion][conn[23]-1]-1;
          nodeIds[21] = nodeMap_[iRegion][conn[21]-1]-1;
          nodeIds[22] = nodeMap_[iRegion][conn[20]-1]-1;
          nodeIds[23] = nodeMap_[iRegion][conn[22]-1]-1;
        }

        // unsigned int id = 
        actGrid->InsertNextCell( elemType, conn.size(),nodeIds );
      }
    }
  } catch(std::exception& ex ) {
    vtkErrorMacro(<< "Caught exception when reading cells: '" 
                  << ex.what() << "'");
  } catch(std::string& ex ) {
    vtkErrorMacro(<< "Caught exception when reading cells: '" 
                  << ex << "'");
  }
}


//----------------------------------------------------------------------------
void vtkCFSReader::ReadNodes(vtkMultiBlockDataSet *output)
{

  try {
    // get all nodal coordinates
    std::vector<std::vector<double> > coords;
    in_.GetNodeCoords( coords );


    // loop over all regions
    for( unsigned int iRegion = 0, n = regionNames_.size(); 
    iRegion < n; iRegion++ ) {
      nodeMap_[iRegion].resize( coords.size() );
      std::vector<unsigned int> regionNodes;
      in_.GetNodesOfRegion( regionNames_[iRegion], regionNodes);
      unsigned int numRegionNodes = regionNodes.size();

      vtkUnstructuredGrid *actGrid = 
        vtkUnstructuredGrid::SafeDownCast( output->GetBlock(iRegion) );

      vtkPoints* points = vtkPoints::New();
      points->SetNumberOfPoints( numRegionNodes );
      
      // Add a result vector containing the original node numbers
      // to the current unstructured grid.
      vtkUnsignedIntArray * origNodeNums = vtkUnsignedIntArray::New();
      origNodeNums->SetNumberOfValues( numRegionNodes );
      origNodeNums->SetName( "origNodeNums" );
      actGrid->GetPointData()->AddArray(origNodeNums);

      for( unsigned int i = 0; i < numRegionNodes; i++ ) {
        const std::vector<double> & p = coords[regionNodes[i]-1];
        origNodeNums->SetValue(i, regionNodes[i]);
        nodeMap_[iRegion][regionNodes[i]-1] = i+1;
        points->SetPoint(i, p[0], p[1], p[2]);
      }
      actGrid->SetPoints( points );
      points->Delete();
    }
  } catch(std::exception& ex ) {
    vtkErrorMacro(<< "Caught exception when reading nodes: '" 
                   << ex.what() << "'");
  } catch(std::string& ex ) {
    vtkErrorMacro(<< "Caught exception when reading nodes: '" 
                   << ex << "'");
  }

  vtkDebugMacro(<< "Finished reading nodes");

}

VTKCellType  vtkCFSReader::GetCellIdType( H5CFS::ElemType cfsType) {
  VTKCellType type =  VTK_EMPTY_CELL;
  switch( cfsType ) {
    case H5CFS::ET_UNDEF:
      break;
    case H5CFS::ET_POINT:
      type = VTK_VERTEX;
      break;
    case H5CFS::ET_LINE2:
      type = VTK_LINE;
      break;
    case H5CFS::ET_LINE3:
      type = VTK_QUADRATIC_EDGE;
      break;
    case H5CFS::ET_TRIA3:
      type = VTK_TRIANGLE;
      break;
    case H5CFS::ET_TRIA6:
      type = VTK_QUADRATIC_TRIANGLE;
      break;
    case H5CFS::ET_QUAD4:
      type = VTK_QUAD;
      break;
    case H5CFS::ET_QUAD8:
      type = VTK_QUADRATIC_QUAD;
      break;
    case H5CFS::ET_QUAD9:
      type = VTK_BIQUADRATIC_QUAD;
      break;
    case H5CFS::ET_TET4:
      type = VTK_TETRA;
      break;
    case H5CFS::ET_TET10:
      type = VTK_QUADRATIC_TETRA;
      break;
    case H5CFS::ET_HEXA8:
      type = VTK_HEXAHEDRON;
      break;
    case H5CFS::ET_HEXA20:
      type = VTK_QUADRATIC_HEXAHEDRON;
      break;
    case H5CFS::ET_HEXA27:
      type = VTK_TRIQUADRATIC_HEXAHEDRON;
      break;
    case H5CFS::ET_PYRA5:
      type = VTK_PYRAMID;
      break;
    case H5CFS::ET_PYRA13:
      type = VTK_QUADRATIC_PYRAMID;
      break;
    case H5CFS::ET_PYRA14:
      type = VTK_QUADRATIC_PYRAMID;
      break;
    case H5CFS::ET_WEDGE6:
      type = VTK_WEDGE;
      break;
    case H5CFS::ET_WEDGE15:
      type = VTK_QUADRATIC_WEDGE;
      break;
    case H5CFS::ET_WEDGE18:
      type = VTK_BIQUADRATIC_QUADRATIC_WEDGE;
      break;
  }
  return type;
      
 }

const char* freqSuffixes[] = 
{"fHz", "nHz", "muHz", "mHz", "Hz", "kHz", "MHz", "GHz", "THz", "PHz", "EHz", "Zhz", "YHz"};
 
void vtkCFSReader::PrettyPrintFreq(double freq, std::string& val) {
  
  // explicitly check for freq = 0
  if( std::abs(freq) < 1e-16) {
    val = "0.0 s";
    return;
  }
  
  int expo3 = floor(log10(std::abs(freq))/3);
  double reduced = freq/pow(double(10.0),expo3*3.0);
  boost::format fmt ("%3.3f %s");
  fmt % reduced % freqSuffixes[expo3+4];
  val = fmt.str();
}

const char* timeSuffixes[] = {"fs", "ns", "mus", "ms", "s"};

void vtkCFSReader::PrettyPrintTime(double time, std::string& val) {

  // explicitly check for time = 0
  if( std::abs(time) < 1e-16) {
    val = "0.0 s";
    return;
  }
    
  int expo3 = floor(log10(std::abs(time))/3);
  double reduced = time/pow(double(10.0),expo3*3);
  boost::format fmt ("%3.3f %s");
  fmt % reduced % timeSuffixes[expo3+4];
  val = fmt.str();
}
