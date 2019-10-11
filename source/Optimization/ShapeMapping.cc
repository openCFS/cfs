/*
 * ShapeMapping.cc
 *
 *  Created on: Mar 14, 2016
 *      Author: fwein
 */

#include "Optimization/ShapeMapping.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


namespace CoupledField {

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

} // enf of namespace

