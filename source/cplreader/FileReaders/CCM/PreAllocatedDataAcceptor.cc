/* 
 * File:   PreAllocatedDataAcceptor.cpp
 * Author: Matthias Tautz
 * 
 * Created on 3. Juni 2013, 01:29
 */

#include <vector>
#include <string>
#include <algorithm>

#include "PreAllocatedDataAcceptor.hh"
#include "FieldDataAcceptor.hh"
#include "Utilities.hh"
#include "Meshes.hh"

namespace CCM {
  
  PreAllocatedDataAcceptor::PreAllocatedDataAcceptor(std::vector<ConsecutiveMap> maps):maps_(maps),numMaps_(maps.size()) {
    loadedDataForMaps_.resize(maps.size(), false);
  }
  
  PreAllocatedDataAcceptor::PreAllocatedDataAcceptor(const PreAllocatedDataAcceptor& orig) {
  }
  
  PreAllocatedDataAcceptor::~PreAllocatedDataAcceptor() {
  }
  
  std::string PreAllocatedDataAcceptor::GetFieldName() { 
    return fieldName_;
  }
  
  void PreAllocatedDataAcceptor::SetFieldName(std::string fieldName) {
    fieldName_ = fieldName;
  }
  
  bool PreAllocatedDataAcceptor::Accept(FieldData& fielddata, bool& alreadyAllocated) {
    if (fielddata.fieldLabel.compare(fieldName_) == 0) {
      for (uint i = 0; i < maps_.size(); i++) {
        std::string l = maps_[i].label;
        if (l.compare(fielddata.mapLabel) == 0) {
          fielddata.doubleData = &data[maps_[i].startValue - 1];
	  loadedDataForMaps_[i] = true;
          alreadyAllocated = true;
          return true;
        }
      }
    }
    return false;
  }
  
  void PreAllocatedDataAcceptor::ResetUnloadedData() {
    std::fill (loadedDataForMaps_.begin(),loadedDataForMaps_.end(),false);
  }
  
  bool PreAllocatedDataAcceptor::HasUnloadedData() {
    for (uint i = 0; i < numMaps_; i++) {
      if (!loadedDataForMaps_[i]) {
	return true;
      }
    }
    return false;
  }
  
  void PreAllocatedDataAcceptor::GetUnloadedDataMaps(std::vector<ConsecutiveMap>& unloadedMaps) {
    unloadedMaps.clear();
    for (uint i = 0; i < numMaps_; i++) {
      if (!loadedDataForMaps_[i]) {
	unloadedMaps.push_back(maps_[i]);
      }
    }
  }
  
}

