/* 
 * File:   ValueMapCatcher.h
 * Author: Mace
 *
 * Created on 4. Juni 2013, 22:31
 */

#ifndef VALUEMAPCATCHER_H
#define	VALUEMAPCATCHER_H

#include <vector>

#include "Utilities.hh"
#include "FieldDataAcceptor.hh"

namespace CCM {
  
  class ValueMapCatcher : public FieldDataAcceptor {
  public:
    ValueMapCatcher(std::vector<Map>& maps, std::vector<Map>& cellMaps, 
            std::vector<Map>& boundaryFaceMaps);
    ValueMapCatcher(const ValueMapCatcher& orig);
    virtual ~ValueMapCatcher();
    
    virtual bool Accept(FieldData& data, bool& alreadyAllocated);
  private:
    std::vector<Map>& maps_;
    std::vector<Map>& cellMaps_;
    std::vector<Map>& boundaryFaceMaps_;
  };
  
}

#endif	/* VALUEMAPCATCHER_H */

