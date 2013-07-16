/* 
 * File:   PreAllocatedCCMDoubleDataAcceptor.h
 * Author: Mace
 *
 * Created on 3. Juni 2013, 01:29
 */

#ifndef PREALLOCATEDDATAACCEPTOR_H
#define	PREALLOCATEDDATAACCEPTOR_H

#include <vector>
#include <string>

#include "FieldDataAcceptor.hh"
#include "Utilities.hh"
#include "Meshes.hh"

namespace CCM {

  class PreAllocatedDataAcceptor : public FieldDataAcceptor {
  public:
    PreAllocatedDataAcceptor(std::vector<ConsecutiveMap> maps);
    PreAllocatedDataAcceptor(const PreAllocatedDataAcceptor& orig);
    virtual ~PreAllocatedDataAcceptor();
    
    double* data;
    void SetFieldName(std::string fieldName);
    bool Accept(FieldData& fielddata, bool& alreadyAllocated);
  private:
    std::string fieldName_;
    std::vector<ConsecutiveMap> maps_;
    
  
  };
  
}

#endif	/* PREALLOCATEDDATAACCEPTOR_H */

