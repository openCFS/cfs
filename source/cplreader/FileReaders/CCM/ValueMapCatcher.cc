/* 
 * File:   ValueMapCatcher.cpp
 * Author: Mace
 * 
 * Created on 4. Juni 2013, 22:31
 */

#include "ValueMapCatcher.hh"

#include <vector>

#include "Utilities.hh"
#include "FieldDataAcceptor.hh"
#include <iostream>

namespace CCM {
  
  ValueMapCatcher::ValueMapCatcher(std::vector<Map>& maps, std::vector<Map>& cellMaps, 
            std::vector<Map>& boundaryFaceMaps)
        :maps_(maps),cellMaps_(cellMaps),boundaryFaceMaps_(boundaryFaceMaps) {
  }
  
  ValueMapCatcher::ValueMapCatcher(const ValueMapCatcher& orig)
    :maps_(orig.maps_),cellMaps_(orig.cellMaps_),boundaryFaceMaps_(orig.boundaryFaceMaps_){
  }
  
  ValueMapCatcher::~ValueMapCatcher() {
  }
  
  bool ValueMapCatcher::Accept(FieldData& data, bool& alreadyAllocated) {
    for (uint i = 0; i < maps_.size(); i++) {
      if (data.mapLabel.compare(maps_[i].label) == 0) {
        std::vector<Map>& list = data.location == CellDataLoc ? cellMaps_ : boundaryFaceMaps_;
        bool found = false;
        for (uint j=0; j<list.size() && !found;j++) {
          found = list[j].label.compare(maps_[i].label) == 0;
        }
        if (!found) {
          list.push_back(maps_[i]);
        }
      }
    }
    return false;
  }
  
}

