/*
 * BaseNcInterface.cc
 *
 *  Created on: 22.09.2013
 *      Author: jens
 */
#include "BaseNcInterface.hh"
#include "Driver/FormsContexts.hh"

namespace CoupledField {

void BaseNcInterface::UpdateIntegrators() {
  // iterate over integrators and assign the interface elemList_
  StdVector<BiLinFormContext*>::iterator it = integrators_.Begin(), itEnd = integrators_.End();
  for ( ; it != itEnd; ++it ) {
    (*it)->SetEntities(elemList_, elemList_);
  }
}

} // end of namespace CoupledField
