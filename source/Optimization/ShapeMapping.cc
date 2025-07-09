/*
 * ShapeMapping.cc
 *
 *  Created on: Mar 14, 2016
 *      Author: fwein
 */

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/ShapeMapping.hh"

namespace CoupledField {

DEFINE_LOG(ShapeMap, "shapeMap")

ShapeMapping::ShapeMapping()
 : AcouSIMP()
{

}

ShapeMapping::~ShapeMapping()
{

}


void ShapeMapping::PostInit()
{
  AcouSIMP::PostInit();
}

} // end of namespace

