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

namespace CoupledField {

// forward declarations
class Grid;

class BaseNcInterface {

  public:

    BaseNcInterface(Grid* grid, std::string &name);

    virtual ~BaseNcInterface();

    std::string& GetName() const { return name_; }
    
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
