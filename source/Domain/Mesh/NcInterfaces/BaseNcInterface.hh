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
      // base class constructor
      BaseNcInterface(Grid* grid): ptGrid_(grid), elemList_(new NcSurfElemList(grid)) {};
      // base class destructor
      virtual ~BaseNcInterface() { ptGrid_ = nullptr; elemList_ = nullptr;};
      // set the name of the interface
      void SetName(const std::string &name) { name_ = name; }
      // return the list of elements
      const shared_ptr<NcSurfElemList> GetElemList() const { return elemList_; }
      // return the name of the interface
      const std::string& GetName() const { return name_; }
      // return if the interface is a moving interface
      virtual bool IsMoving() const = 0;
      // reset the interface
      virtual void ResetInterface() = 0;
      // update the interface
      virtual void UpdateInterface() = 0;
      // add passed integrator to the stored integrators_
      void RegisterIntegrator(BiLinFormContext* context) { integrators_.Push_back(context); }
    protected:
      // update the integrators by assigning the cached elements in elemList_
      void UpdateIntegrators();
      // =======================================================================
      // Class variables
      // =======================================================================
      //! pointer to the grid
      Grid* ptGrid_;
      //! name of the interface
      std::string name_;
      //! pointer to the entity list containing all intersection grid elements
      shared_ptr<NcSurfElemList> elemList_;
      //! pointer to the information for defining the surface bilinear integrators
      StdVector<BiLinFormContext*> integrators_;
  }; // class BaseNcInterface
} // namespace CoupledField

#endif /* _NCINTERFACE_HH_ */
