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
    
    enum NcMotionType {
      NCI_NOT_MOVING,
      NCI_LINEAR_MOTION,
      NCI_ROTATION,
    };
    
    std::string& GetName() const { return name_; }
    
    RegionIdType GetRegionId() const { return region_; }
    
    RegionIdType GetMasterSurfRegion() const { return masterSurfRegion_; }
    
    RegionIdType GetSlaveSurfRegion() const { return masterSurfRegion_; }

    RegionIdType GetMasterVolRegion() const { return masterVolRegion_; }
    
    RegionIdType GetSlaveVolRegion() const { return masterVolRegion_; }
    
    bool IsPlanar() const { return coplanar_; }
    
    NcMotionType GetMotionType() const { return motionType_; }
    
    std::string& GetCoordSys() const { return coordSysId_; }
    
    void MoveInterface();

  protected:

    void UpdateIntersection();
    
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
    NcMotionType motionType_;
    NcIntersectAlgo intersectAlgo_;
    Double tolAbs_;
    Double tolRel_;
    Double rpm_;
    std::string coordSysId_;

  }; // class NcInterface
  
} // namespace CoupledField

#endif /* _NCINTERFACE_HH_ */
