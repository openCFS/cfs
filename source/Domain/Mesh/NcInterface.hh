/*
 * NcInterface.hh
 *
 *  Created on: 19.02.2013
 *      Author: jens
 */

#ifndef _NCINTERFACE_HH_
#define _NCINTERFACE_HH_

#include <string>

#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  // forward declarations
  class Grid;
  
  class NcInterface {
    
  public:

    NcInterface(Grid* grid, PtrParamNode nciNode);
    
    virtual ~NcInterface();
    
    std::string& GetName() { return name_; }
    
    RegionIdType GetRegionId() { return region_; }
    
    RegionIdType GetMasterSurfRegion() { return masterSurfRegion_; }
    
    RegionIdType GetSlaveSurfRegion() { return masterSurfRegion_; }

    RegionIdType GetMasterVolRegion() { return masterVolRegion_; }
    
    RegionIdType GetSlaveVolRegion() { return masterVolRegion_; }
    
    bool IsPlanar() { return coplanar_; }
    
    bool IsMoving() { return needsUpdate_; }
    
    std::string& GetCoordSys() { return coordSysId_; }

  protected:

    enum NcIntersectAlgo {
      NCI_INTERSECT_LINE,
      NCI_INTERSECT_RECT,
      NCI_INTERSECT_POLYGON,
    };

    Grid* ptGrid_;
    std::string name_;
    RegionIdType region_;
    RegionIdType masterSurfRegion_;
    RegionIdType slaveSurfRegion_;
    RegionIdType masterVolRegion_;
    RegionIdType slaveVolRegion_;
    bool coplanar_;
    bool needsUpdate_;
    NcIntersectAlgo intersectAlgo_;
    Double tolAbs_;
    Double tolRel_;
    Double rpm_;
    std::string coordSysId_;

  }; // class NcInterface
  
} // namespace CoupledField

#endif /* _NCINTERFACE_HH_ */
