/* 
 * File:   PreAllocatedDataAcceptor.cpp
 * Author: Mace
 * 
 * Created on 3. Juni 2013, 01:29
 */

#include <vector>
#include <string>

#include "PreAllocatedDataAcceptor.hh"
#include "FieldDataAcceptor.hh"
#include "Utilities.hh"
#include "Meshes.hh"

namespace CCM {
  
  PreAllocatedDataAcceptor::PreAllocatedDataAcceptor(std::vector<ConsecutiveMap> maps):maps_(maps) {
  }
  
  PreAllocatedDataAcceptor::PreAllocatedDataAcceptor(const PreAllocatedDataAcceptor& orig) {
  }
  
  PreAllocatedDataAcceptor::~PreAllocatedDataAcceptor() {
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
          alreadyAllocated = true;
          return true;
        }
      }
    }
    return false;
  }
  
}

