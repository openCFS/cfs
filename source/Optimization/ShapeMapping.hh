/*
 * ShapeMapping.hh
 *
 *  Created on: Mar 11, 2016
 *      Author: fwein
 */
#ifndef OPTIMIZATION_SHAPEMAPPING_HH_
#define OPTIMIZATION_SHAPEMAPPING_HH_

#include "Optimization/AcouSIMP.hh"

namespace CoupledField {

/** Perform parametric shape mapping. The parametric shape maps to the ersatz material. That means a standard SIMP like
 * parameterization and gradients which is a function of the mapping. Only the shape parameters are given to the optimizer but the
 * DesignElements are used as in standard ErsatzMaterial */
class ShapeMapping : public AcouSIMP
{
public:
  ShapeMapping();

  virtual ~ShapeMapping();

  virtual void PostInit();

};

} // end of namespace

#endif /* OPTIMIZATION_SHAPEMAPPING_HH_ */
