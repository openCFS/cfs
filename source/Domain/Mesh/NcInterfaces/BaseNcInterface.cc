/*
 * NcInterface.cc
 *
 *  Created on: 20.02.2013
 *      Author: jens
 */

#include "BaseNcInterface.hh"
#include "Domain/Mesh/Grid.hh"

namespace CoupledField {

BaseNcInterface::BaseNcInterface(Grid* grid, std::string &name) :
  ptGrid_(grid),
  name_(name)
{
  region_ = ptGrid_->AddSurfaceRegion(name_);
}

BaseNcInterface::~BaseNcInterface() {
  ptGrid_ = NULL;
}

} // namespace CoupledField
