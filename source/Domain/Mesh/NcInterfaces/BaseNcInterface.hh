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
#include "Domain/Mesh/Grid.hh"

namespace CoupledField {

class BaseNcInterface {

  public:

    BaseNcInterface(Grid* grid) : ptGrid_(grid), region_(NO_REGION_ID) {};

    virtual ~BaseNcInterface() { ptGrid_ = NULL; };

    const std::string& GetName() const { return name_; }
    
    void SetName(const std::string & name) {
      this->name_ = name;
      region_ = ptGrid_->AddSurfaceRegion(name_);
    }
    
    RegionIdType GetRegionId() const { return region_; }

    virtual bool NeedsUpdate() const = 0;
    
    virtual void UpdateInterface() = 0;

  protected:

    // =======================================================================
    // Class variables
    // =======================================================================

    Grid* ptGrid_;
    std::string name_;
    RegionIdType region_;

}; // class BaseNcInterface

} // namespace CoupledField

#endif /* _NCINTERFACE_HH_ */
