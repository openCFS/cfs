/* 
 * File:   FlowDataLoader.cpp
 * Author: Mace
 * 
 * Created on 4. Juni 2013, 23:05
 */

#include "FlowDataLoader.hh"

#include <vector>
#include <string>

#include "Utilities.hh"
#include "DataFile.hh"
#include "PreAllocatedDataAcceptor.hh"


#include "General/environment.hh"
#include "Domain/resultInfo.hh"
#include "cplreader/FlowDataTypes.hh"

#include "cplreader/Settings.hh"
#include <def_cplreader.hh>
#include "cplreader/FlowDataTypes.hh"
#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "integlib/elemIntegr.hh"

#include "FlowDataTypes.hh"
#include "OutputWriter.hh"
#include "CouplingHandler.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace CoupledField;

namespace CCM {
  
  FlowDataLoader::FlowDataLoader(uint cellCount, std::vector<ConsecutiveMap>& maps)
      :cellCount_(cellCount),acceptor_(maps) {
  }
  
  FlowDataLoader::FlowDataLoader(const FlowDataLoader& orig):cellCount_(orig.cellCount_),acceptor_(acceptor_) {
  }
  
  FlowDataLoader::~FlowDataLoader() {
  }
  
  void FlowDataLoader::LoadScalar(FlowDataPartStruct* fdps, DataFile& file, std::string name, SolutionType solType) {
    fdps->data.resize(cellCount_);
    acceptor_.data = &fdps->data[0];
    acceptor_.SetFieldName(name);
    file.ReadFields(NULL, &acceptor_);
    
    fdps->isActive = true;
    fdps->definedOn = ResultInfo::NODE;
    fdps->entryType = ResultInfo::SCALAR;
    fdps->dofNames.push_back("-");
    fdps->unit = MapSolTypeToUnit(solType);
    fdps->resultName = SolutionTypeEnum.ToString(solType);
  }
  
  void FlowDataLoader::LoadVector(FlowDataPartStruct* fdps, DataFile& file, std::string name, SolutionType solType) {
    fdps->data.resize(cellCount_*3);
    double* cache = new double[cellCount_];
    acceptor_.data = cache;
    
    std::string nameC = name + "_0";
    acceptor_.SetFieldName(nameC);
    file.ReadFields(NULL, &acceptor_);
#pragma omp parallel for
    for (uint i=0; i < cellCount_; i++) {
      fdps->data[i*3] = cache[i];
    }
    
    nameC = name + "_1";
    acceptor_.SetFieldName(nameC);
    file.ReadFields(NULL, &acceptor_);
#pragma omp parallel for
    for (uint i=0; i < cellCount_; i++) {
      fdps->data[i*3+1] = cache[i];
    }

    nameC = name + "_2";
    acceptor_.SetFieldName(nameC);
    file.ReadFields(NULL, &acceptor_);
#pragma omp parallel for
    for (uint i=0; i < cellCount_; i++) {
      fdps->data[i*3+2] = cache[i];
    }
    
    delete[] cache;
    
    
    fdps->isActive = true;
    fdps->definedOn = ResultInfo::NODE;
    fdps->entryType = ResultInfo::VECTOR;
    fdps->dofNames.push_back("x");
    fdps->dofNames.push_back("y");
    fdps->dofNames.push_back("z");
    fdps->unit = MapSolTypeToUnit(solType);
    fdps->resultName = SolutionTypeEnum.ToString(solType);
  }

}

