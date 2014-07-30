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

class BiLinFormContext;

class BaseNcInterface {

  public:

    BaseNcInterface(Grid* grid)
      : ptGrid_(grid),
        elemList_( new NcSurfElemList(grid) )
    {};

    virtual ~BaseNcInterface() { ptGrid_ = NULL; };

    const std::string& GetName() const { return name_; }
    
    void SetName(const std::string &name) {
      name_ = name;
    }

    const shared_ptr<NcSurfElemList> GetElemList() const {
      return elemList_;
    }
    
    virtual bool NeedsUpdate() const = 0;
    
    virtual void ResetInterface() = 0;

    virtual void UpdateInterface() = 0;

    void RegisterIntegrator(BiLinFormContext* context) {
      integrators_.Push_back(context);
    }
    
  protected:

    void UpdateIntegrators();
    
    // =======================================================================
    // Class variables
    // =======================================================================

    Grid* ptGrid_;
    std::string name_;
    shared_ptr<NcSurfElemList> elemList_;
    StdVector<BiLinFormContext*> integrators_;

}; // class BaseNcInterface

} // namespace CoupledField

#endif /* _NCINTERFACE_HH_ */
