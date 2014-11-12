/* 
 * File:   NameFieldDataAcceptor.cpp
 * Author: Mace
 * 
 * Created on 2. Juni 2013, 00:36
 */

#include "NameFieldDataAcceptor.hh"

namespace CCM {
  
  NameFieldDataAcceptor::NameFieldDataAcceptor(std::string name):name_(name) {
  }
  
  NameFieldDataAcceptor::NameFieldDataAcceptor(const NameFieldDataAcceptor& orig) {
  }
  
  NameFieldDataAcceptor::~NameFieldDataAcceptor() {
  }
  
  bool NameFieldDataAcceptor::Accept(FieldData& data, bool& alreadyAllocated) {
    alreadyAllocated = false;
    return data.fieldLabel.compare(name_) == 0;
  }
  
}

