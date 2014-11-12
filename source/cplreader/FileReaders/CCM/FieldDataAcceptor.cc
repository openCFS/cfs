/* 
 * File:   CCMFieldDataAcceptor.cpp
 * Author: Mace
 * 
 * Created on 1. Juni 2013, 23:55
 */

#include "FieldDataAcceptor.hh"

namespace CCM {
  
  FieldDataAcceptor::FieldDataAcceptor() {
  }
  
  FieldDataAcceptor::FieldDataAcceptor(const FieldDataAcceptor& orig) {
  }
  
  FieldDataAcceptor::~FieldDataAcceptor() {
  }
  
  bool FieldDataAcceptor::Accept(FieldData& data, bool& alreadyAllocated) {
    alreadyAllocated = false;
    return true;
  }
  
}

