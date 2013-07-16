/* 
 * File:   FlowDataLoader.h
 * Author: Mace
 *
 * Created on 4. Juni 2013, 23:05
 */

#ifndef FLOWDATALOADER_H
#define	FLOWDATALOADER_H

#include <vector>
#include <string>

#include "Utilities.hh"
#include "DataFile.hh"
#include "PreAllocatedDataAcceptor.hh"

#include "cplreader/Settings.hh"
#include <def_cplreader.hh>
#include "cplreader/FlowDataTypes.hh"
#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "integlib/elemIntegr.hh"

namespace CCM {

  class FlowDataLoader {
  public:
    FlowDataLoader(uint cellCount, std::vector<ConsecutiveMap>& maps);
    FlowDataLoader(const FlowDataLoader& orig);
    virtual ~FlowDataLoader();
    void LoadScalar(FlowDataPartStruct* fdps, DataFile& file, std::string name, SolutionType solType);
    void LoadVector(FlowDataPartStruct* fdps, DataFile& file, std::string name, SolutionType solType);
  private:
    uint cellCount_;
    PreAllocatedDataAcceptor acceptor_;
  };
  
}

#endif	/* FLOWDATALOADER_H */

