/*
 * ShapeMapping.cc
 *
 *  Created on: Mar 14, 2016
 *      Author: fwein
 */

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Optimization/ShapeMapping.hh"


namespace CoupledField {

DECLARE_LOG(ShapeMap)
DEFINE_LOG(ShapeMap, "shapeMap")

ShapeMapping::ShapeMapping()
 : SIMP()
{

}

ShapeMapping::~ShapeMapping()
{

}


void ShapeMapping::PostInit()
{
  SIMP::PostInit();
}

} // end of namespace

