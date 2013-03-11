// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Settings.hh"

#include "FileReader_fvuns.hh"
#include "fv_reader_tags.h"

//Filling Start

#include <stdlib.h> // For exit()
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <string>


//Manipulation Start
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>


#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
namespace fs=boost::filesystem;

#include "Domain/resultInfo.hh"

#include "cplreader/Settings.hh"

#include <def_cplreader.hh>
#include "cplreader/FileReader.hh"

//Manipulation End

int VERBOSE;
#define TRACE 0
#define DEBUG 0
#define COORDDEBUG 0
#define VARDEBUG 0

namespace CoupledField {

// constructor
FileReader_fvuns::FileReader_fvuns(
    const std::string& name,
    const UInt dim,
    const UInt numFiles,
    const UInt startIndex) :
    FileReader(name, dim, numFiles, startIndex) {

  // get the settings
  Settings& settings = Settings::Instance();

  if( settings.GetInt("verbose") == 1 ) {
    VERBOSE = 1;
  } else {
    VERBOSE = 0;
  }

  // TODO: only dim == 3 supported right now
  if( dim != 3 ) {
    std::cerr << "CPLREADER_FIELDVIEW currently only supports \"--dim 3\"" << std::endl;
    exit(EXIT_FAILURE);
  }

  // output information
  if( VERBOSE ) {
    std::cout << "Creating FileReader_fvuns: " << name << " " << dim << " " << numFiles << " " << startIndex << std::endl;
  }

  // redefine basename; .fvuns files are expected to be named baseName.fvuns or baseName_timestep.fvuns
  baseName_ = settings.GetString("basedir");
  baseName_ += "/";
  baseName_ += name_;
  baseName_ += "/";
  baseName_ += name_;

  if( VERBOSE ) {
    std::cout << "Setting baseName to '" << baseName_ << "'" << std::endl;
  }

}

FileReader_fvuns::~FileReader_fvuns() {

}

void FileReader_fvuns::Init() {

  if( TRACE ) std::cout << "TRACE: FileReader_fvuns::Init()" << std::endl;

  // filename for initialization is {baseName_}.fvuns
  std::string filename = baseName_ + ".fvuns";

  if( VERBOSE ) std::cout << "Initializing with '" << filename << "'" << std::endl;

  // read data from the .fvuns file
  fvunsFile_.openGrid(filename);

  // set some values
  numRegions_ = (UInt)fvunsFile_.getNumGrids();

  // output some information
  std::cout << "File contains " << numRegions_ << " region(s)." << std::endl;

  // nodes and elements per region
  numNodesPerRegion_.resize(numRegions_);
  numElemsPerRegion_.resize(numRegions_);

  for( UInt i = 0; i < numRegions_; i++ ) {
    numNodesPerRegion_[i] = (UInt)fvunsFile_.getNumNodes(i);
    numElemsPerRegion_[i] = (UInt)fvunsFile_.getNumElements(i);
    std::cout << "Region " << (i+1) << " has " << numNodesPerRegion_[i] << " nodes / " << numElemsPerRegion_[i] << " elements." << std::endl;

  }

}

void FileReader_fvuns::ReadNodalCoords(std::vector<Double> & NODECOORD) {

  /*
   * NODECOORD is to be filled with the coordinates of the nodes from ALL regions
   * format is: x,y,z, x,y,z, x,y,z ...
   *
   * serialize nodes from all grids in fvuns file AND change order (x,x,x,...y,y,y,...z,z,z,... in fvuns)
   */

  // resize NODECOORD
  UInt numNodes = 0;
  for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {
    numNodes += numNodesPerRegion_[regionIdx];
  }
  NODECOORD.resize(3 * numNodes);

  // fill NODECOORD
  UInt nodeOffset = 0;

  for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {

    float *nodeCoords = fvunsFile_.getNodeCoordinates(regionIdx);

    // check if nodeCoords are available from meshFile
    if( nodeCoords == NULL ) {

      std::cerr << "ERROR: " << fvunsFile_.getFileName() << " did not contain node coordinates!" << std::endl;
      exit(EXIT_FAILURE);

    }

    for( UInt nodeIdx = 0; nodeIdx < numNodesPerRegion_[regionIdx]; nodeIdx++ ) {

      NODECOORD[(nodeOffset * 3) + (3 * nodeIdx)] = nodeCoords[nodeIdx];
      NODECOORD[(nodeOffset * 3) + (3 * nodeIdx) + 1] = nodeCoords[numNodesPerRegion_[regionIdx] + nodeIdx];
      NODECOORD[(nodeOffset * 3) + (3 * nodeIdx) + 2] = nodeCoords[(numNodesPerRegion_[regionIdx] * 2) + nodeIdx];

      if( COORDDEBUG ) {
        std::cout << "NODECOORD[" << (nodeOffset * 3) + (3 * nodeIdx) << "] = " << nodeCoords[nodeIdx] << std::endl;
        std::cout << "NODECOORD[" << (nodeOffset * 3) + (3 * nodeIdx) + 1<< "] = " << nodeCoords[numNodesPerRegion_[regionIdx] + nodeIdx] << std::endl;
        std::cout << "NODECOORD[" << (nodeOffset * 3) + (3 * nodeIdx) + 2 << "] = " << nodeCoords[(numNodesPerRegion_[regionIdx] * 2) + nodeIdx] << std::endl;
      }

    }

    nodeOffset += numNodesPerRegion_[regionIdx];

  }

}

void FileReader_fvuns::GetRegionElements(std::vector<UInt> & regionElements, const UInt regionIdx) {

  // it is important to note that the elements in the fvuns file are NOT numbered
  // their number/id is defined HERE

  // the easiest approach: first element in first grid is element 1, last element in first grid
  // is element numElemsPerRegion_[0], first element in second grid is numElemsPerRegion_[0]+1 etc.

  // offset i.e. total number of elements in regions of lower index
  UInt elemOffset = 0;

  // add up
  for( UInt i = 0; i < regionIdx; i++ ) {
    elemOffset += numElemsPerRegion_[i];
  }

  // resize vector
  regionElements.resize(numElemsPerRegion_[regionIdx]);

  // fill with offseted numbers
  for ( UInt i = 0; i < numElemsPerRegion_[regionIdx]; i++ ) {
    // +1 since we start at 0
    regionElements[i] = elemOffset + i + 1;
  }

}

void FileReader_fvuns::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
    const std::vector<bool>& activeParts,
    const UInt timeStepIdx) {

  Settings& settings = Settings::Instance();

  // filename for timestep is {baseName_}_{startIndex+timeStepIdx}.fvuns
  UInt currentIdx = startIndex_ + timeStepIdx;
  std::stringstream filename;
  filename << baseName_ << "_" << currentIdx << ".fvuns";
  std::cout << "Processing " << filename.str() << std::endl;

  // read data from the .fvuns file
  fvunsFile_.openVars(filename.str().c_str());

  // only first timestep: find out which data is contained in the .fvuns file
  if( timeStepIdx == 0 ) {

    std::vector<std::string> nodalVarNames = fvunsFile_.getNodalVarNames();

    for( UInt i = 0; i < nodalVarNames.size(); i++ ) {

      if( ( strcmp(nodalVarNames.at(i).c_str(), "Velocity_0; Velocity") == 0 ) ||
          ( strcmp(nodalVarNames.at(i).c_str(), "MappedVertexVelocity_0; MappedVertexVelocity") == 0 ) ) {

        if( VERBOSE ) std::cout << "Found velocity data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(FLUIDMECH_VELOCITY);
        dofIndices_.push_back(0);

      } else if( ( strcmp(nodalVarNames.at(i).c_str(), "Velocity_1") == 0 ) ||
          ( strcmp(nodalVarNames.at(i).c_str(), "MappedVertexVelocity_1") == 0 ) ) {

        if( VERBOSE ) std::cout << "Found velocity data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(FLUIDMECH_VELOCITY);
        dofIndices_.push_back(1);

      } else if( ( strcmp(nodalVarNames.at(i).c_str(), "Velocity_2") == 0 ) ||
          ( strcmp(nodalVarNames.at(i).c_str(), "MappedVertexVelocity_2") == 0 ) ) {

        if( VERBOSE ) std::cout << "Found velocity data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(FLUIDMECH_VELOCITY);
        dofIndices_.push_back(2);

      } else if( ( strcmp(nodalVarNames.at(i).c_str(), "Pressure") == 0 ) ||
          ( strcmp(nodalVarNames.at(i).c_str(), "MappedVertexPressure") == 0 ) ) {

        if( VERBOSE ) std::cout << "Found pressure data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(FLUIDMECH_PRESSURE);
        dofIndices_.push_back(0);

      } else if( strcmp(nodalVarNames.at(i).c_str(), "MappedVertexdLHT_0; MappedVertexdLHT") == 0 ) {

        if( VERBOSE ) std::cout << "Found divergence of lighthill tensor data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(ACOU_DIV_LH_TENSOR_NODAL);
        dofIndices_.push_back(0);

      } else if( strcmp(nodalVarNames.at(i).c_str(), "MappedVertexdLHT_1") == 0 ) {

        if( VERBOSE ) std::cout << "Found divergence of lighthill tensor data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(ACOU_DIV_LH_TENSOR_NODAL);
        dofIndices_.push_back(1);

      } else if( strcmp(nodalVarNames.at(i).c_str(), "MappedVertexdLHT_2") == 0 ) {

        if( VERBOSE ) std::cout << "Found divergence of lighthill tensor data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(ACOU_DIV_LH_TENSOR_NODAL);
        dofIndices_.push_back(2);

      } else {

        std::cout << "Found unrecognized data (variable " << (i+1) << ", '" << nodalVarNames.at(i) << "')" << std::endl;
        solutionTypes_.push_back(NO_SOLUTION_TYPE);
        dofIndices_.push_back(0);

      }

    }

  }

  // only first timestep: initialize nodalFlowData
  if( timeStepIdx == 0 ) {

    if( VERBOSE ) std::cout << "Initializing nodalFlowData object" << std::endl;

    // loop over regions
    for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {

      if( VERBOSE ) std::cout << "Now processing region " << (regionIdx + 1) << std::endl;

      // get the FlowDataType for this region
      FlowDataType& fd = nodalFlowData[regionIdx];
      FlowDataPartStruct* fdps = NULL;

      // set all the flow data types where data is available
      for( UInt i = 0; i < solutionTypes_.size(); i++ ) {

        switch( solutionTypes_[i] ) {

        case FLUIDMECH_VELOCITY:

          // velocity
          fdps = &fd[FLUIDMECH_VELOCITY];

          // if velocity was already processed (since it comes up three times), skip
          if( fdps->isActive == true ) {
            break;
          }

          // first time velocity comes up, set everything
          std::cout << "Processing velocity data...";
          fdps->isActive = true;
          fdps->definedOn = ResultInfo::NODE;
          fdps->entryType = ResultInfo::VECTOR;
          fdps->dofNames.push_back("x");
          fdps->dofNames.push_back("y");
          fdps->dofNames.push_back("z");
          fdps->unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
          fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
          fdps->data.resize(numNodesPerRegion_[regionIdx] * fdps->dofNames.size());
          std::cout << "done." << std::endl;
          break;

        case ACOU_DIV_LH_TENSOR_NODAL:

          // divergence of lighthill tensor
          if( settings.GetInt("useDivLHT") ) {

            fdps = &fd[ACOU_DIV_LH_TENSOR_NODAL];

            // if data was already processed (since it comes up three times), skip
            if( fdps->isActive == true ) {
              break;
            }

            // first time data comes up, set everything
            std::cout << "Processing div of Lighthill-Tensor data...";
            fdps->isActive = true;
            fdps->definedOn = ResultInfo::NODE;
            fdps->entryType = ResultInfo::VECTOR;
            fdps->dofNames.push_back("x");
            fdps->dofNames.push_back("y");
            fdps->dofNames.push_back("z");
            fdps->unit = MapSolTypeToUnit(ACOU_DIV_LH_TENSOR_NODAL);
            fdps->resultName = SolutionTypeEnum.ToString(ACOU_DIV_LH_TENSOR_NODAL);
            fdps->data.resize(numNodesPerRegion_[regionIdx] * fdps->dofNames.size());
            std::cout << "done." << std::endl;

          } else {

            if( VERBOSE ) std::cout << "div of Lighthill-Tensor data found, but useDivLH = false, ignoring data." << std::endl;

          }

          break;

        default:

          if( VERBOSE ) std::cout << "Ignoring data for " << SolutionTypeEnum.ToString(solutionTypes_[i]) << std::endl;

          break;

        }

      }

    }

  }

  // every timestep: copy the nodal values
  if( VERBOSE ) std::cout << "Filling nodalFlowData object" << std::endl;

  // loop over regions
  for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {

    if( VERBOSE ) std::cout << "Now processing region " << (regionIdx + 1) << "...";

    // get the FlowDataType for this region
    FlowDataType& fd = nodalFlowData[regionIdx];
    FlowDataPartStruct* fdps = NULL;

    // set all the flow data types where data is available
    for( UInt i = 0; i < solutionTypes_.size(); i++ ) {

      switch( solutionTypes_[i] ) {

      case FLUIDMECH_VELOCITY:

        if( VERBOSE ) std::cout << "Reading velocity data...";

        // velocity
        fdps = &fd[FLUIDMECH_VELOCITY];

        // copy nodal values
        for( UInt nodeIdx = 0; nodeIdx < numNodesPerRegion_[regionIdx]; nodeIdx++ ) {

          // get nodal vars
          float* nodalVars = fvunsFile_.getNodalVars(regionIdx);
          fdps->data[(nodeIdx * 3) + 0] = nodalVars[((i+0) * numNodesPerRegion_[regionIdx]) + nodeIdx];
          fdps->data[(nodeIdx * 3) + 1] = nodalVars[((i+1) * numNodesPerRegion_[regionIdx]) + nodeIdx];
          fdps->data[(nodeIdx * 3) + 2] = nodalVars[((i+2) * numNodesPerRegion_[regionIdx]) + nodeIdx];

          if( VARDEBUG ) {
            std::cout << "fdps->data[" << (nodeIdx * 3) + 0 << "] = " << nodalVars[((i+0) * numNodesPerRegion_[regionIdx]) + nodeIdx] << std::endl;
            std::cout << "fdps->data[" << (nodeIdx * 3) + 1 << "] = " << nodalVars[((i+1) * numNodesPerRegion_[regionIdx]) + nodeIdx] << std::endl;
            std::cout << "fdps->data[" << (nodeIdx * 3) + 2 << "] = " << nodalVars[((i+2) * numNodesPerRegion_[regionIdx]) + nodeIdx] << std::endl;
          }

        }

        // advance i by additional 2 (3 in total) for 3 components
        i += 2;

        if( VERBOSE ) std::cout << "done." << std::endl;

        break;

      case ACOU_DIV_LH_TENSOR_NODAL:

        // divergence of lighthill tensor
        if( settings.GetInt("useDivLHT") ) {

          if( VERBOSE ) std::cout << "Reading divLHT data...";

          fdps = &fd[ACOU_DIV_LH_TENSOR_NODAL];

          // copy nodal values
          for( UInt nodeIdx = 0; nodeIdx < numNodesPerRegion_[regionIdx]; nodeIdx++ ) {

            // get nodal vars
            float* nodalVars = fvunsFile_.getNodalVars(regionIdx);
            fdps->data[(nodeIdx * 3) + 0] = nodalVars[((i+0) * numNodesPerRegion_[regionIdx]) + nodeIdx];
            fdps->data[(nodeIdx * 3) + 1] = nodalVars[((i+1) * numNodesPerRegion_[regionIdx]) + nodeIdx];
            fdps->data[(nodeIdx * 3) + 2] = nodalVars[((i+2) * numNodesPerRegion_[regionIdx]) + nodeIdx];

            if( VARDEBUG ) {
              std::cout << "fdps->data[" << (nodeIdx * 3) + 0 << "] = " << nodalVars[((i+0) * numNodesPerRegion_[regionIdx]) + nodeIdx] << std::endl;
              std::cout << "fdps->data[" << (nodeIdx * 3) + 1 << "] = " << nodalVars[((i+1) * numNodesPerRegion_[regionIdx]) + nodeIdx] << std::endl;
              std::cout << "fdps->data[" << (nodeIdx * 3) + 2 << "] = " << nodalVars[((i+2) * numNodesPerRegion_[regionIdx]) + nodeIdx] << std::endl;
            }

          }

        } else {

          std::cout << "div of Lighthill-Tensor data found, but useDivLH = false, ignoring data." << std::endl;

        }

        // advance i by additional 2 (3 in total) for 3 components
        i += 2;

        if( VERBOSE ) std::cout << "done." << std::endl;

        break;

      default:

        std::cout << "Ignoring data for " << SolutionTypeEnum.ToString(solutionTypes_[i]) << std::endl;

        break;

      }

    }

    //if( VERBOSE ) std::out << "done." << std::endl;

  }

}

void FileReader_fvuns::ReadTopology(std::vector<UInt> & TOPOLOGYDATA, std::vector<UInt> & elemTypes) {

  // find out maximum number of element nodes per region
  UInt *maxNumElemNodes = new UInt[numRegions_];

  for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {

    if( fvunsFile_.getNumHexElems(regionIdx) > 0 ) {
      maxNumElemNodes[regionIdx] = 8;
      // 8 is maximum, in this case we know maxNumElemNodes_ already
      maxNumElemNodes_ = 8;
      break;
    } else if( fvunsFile_.getNumPriElems(regionIdx) > 0 ) {
      maxNumElemNodes[regionIdx] = 6;
      continue;
    } else if( fvunsFile_.getNumPyrElems(regionIdx) > 0 ) {
      maxNumElemNodes[regionIdx] = 5;
      continue;
    } else if( fvunsFile_.getNumTetElems(regionIdx) > 0 ) {
      maxNumElemNodes[regionIdx] = 4;
    }

  }

  if( maxNumElemNodes_ != 8 ) {

    for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {

      if( maxNumElemNodes[regionIdx] > maxNumElemNodes_ ) maxNumElemNodes_ = maxNumElemNodes[regionIdx];

    }

  }

  std::cout << "Maximum number of element nodes is " << maxNumElemNodes_ << std::endl;
  delete[] maxNumElemNodes;

  // find out global number of elems
  UInt numElems = 0;

  for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {
    numElems += numElemsPerRegion_[regionIdx];
  }

  // resize arrays
  TOPOLOGYDATA.resize(maxNumElemNodes_ * numElems);
  elemTypes.resize(numElems);

  // fill arrays

  // global element index
  int elemIdx = 0;
  // global node offset; how many nodes in the regions before?
  int nodeOffset = 0;


  // walk regions
  for( UInt regionIdx = 0; regionIdx < numRegions_; regionIdx++ ) {

    // tet
    for( int i = 0; i < fvunsFile_.getNumTetElems(regionIdx); i++ ) {

      int* elemNodes = fvunsFile_.getTetElems(regionIdx);

      TOPOLOGYDATA[(elemIdx * 8) + 0] = elemNodes[(i * 4) + 0] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 1] = elemNodes[(i * 4) + 1] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 2] = elemNodes[(i * 4) + 2] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 3] = elemNodes[(i * 4) + 3] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 4] = 0;
      TOPOLOGYDATA[(elemIdx * 8) + 5] = 0;
      TOPOLOGYDATA[(elemIdx * 8) + 6] = 0;
      TOPOLOGYDATA[(elemIdx * 8) + 7] = 0;

      elemTypes[elemIdx] = Elem::TET4;

      elemIdx++;

    }

    // pyr
    for( int i = 0; i < fvunsFile_.getNumPyrElems(regionIdx); i++ ) {

      int* elemNodes = fvunsFile_.getPyrElems(regionIdx);

      TOPOLOGYDATA[(elemIdx * 8) + 0] = elemNodes[(i * 5) + 0] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 1] = elemNodes[(i * 5) + 1] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 2] = elemNodes[(i * 5) + 2] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 3] = elemNodes[(i * 5) + 3] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 4] = elemNodes[(i * 5) + 4] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 5] = 0;
      TOPOLOGYDATA[(elemIdx * 8) + 6] = 0;
      TOPOLOGYDATA[(elemIdx * 8) + 7] = 0;

      elemTypes[elemIdx] = Elem::PYRA5;

      elemIdx++;

    }

    // pri
    for( int i = 0; i < fvunsFile_.getNumPriElems(regionIdx); i++ ) {

      int* elemNodes = fvunsFile_.getPriElems(regionIdx);

      TOPOLOGYDATA[(elemIdx * 8) + 0] = elemNodes[(i * 6) + 0] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 1] = elemNodes[(i * 6) + 3] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 2] = elemNodes[(i * 6) + 5] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 3] = elemNodes[(i * 6) + 1] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 4] = elemNodes[(i * 6) + 2] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 5] = elemNodes[(i * 6) + 4] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 6] = 0;
      TOPOLOGYDATA[(elemIdx * 8) + 7] = 0;

      elemTypes[elemIdx] = Elem::WEDGE6;

      elemIdx++;

    }

    // hex
    for( int i = 0; i < fvunsFile_.getNumHexElems(regionIdx); i++ ) {

      int* elemNodes = fvunsFile_.getHexElems(regionIdx);
      /*
      TOPOLOGYDATA[(elemIdx * 8) + 0] = elemNodes[(i * 8) + 2] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 1] = elemNodes[(i * 8) + 6] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 2] = elemNodes[(i * 8) + 4] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 3] = elemNodes[(i * 8) + 0] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 4] = elemNodes[(i * 8) + 3] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 5] = elemNodes[(i * 8) + 7] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 6] = elemNodes[(i * 8) + 5] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 7] = elemNodes[(i * 8) + 1] + nodeOffset;
      */
      TOPOLOGYDATA[(elemIdx * 8) + 0] = elemNodes[(i * 8) + 0] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 1] = elemNodes[(i * 8) + 1] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 2] = elemNodes[(i * 8) + 3] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 3] = elemNodes[(i * 8) + 2] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 4] = elemNodes[(i * 8) + 4] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 5] = elemNodes[(i * 8) + 5] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 6] = elemNodes[(i * 8) + 7] + nodeOffset;
      TOPOLOGYDATA[(elemIdx * 8) + 7] = elemNodes[(i * 8) + 6] + nodeOffset;
      elemTypes[elemIdx] = Elem::HEXA8;

      elemIdx++;

    }

    nodeOffset += fvunsFile_.getNumNodes(regionIdx);

  }

}

}



