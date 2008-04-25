#ifndef SHAPEGRAD_HH_
#define SHAPEGRAD_HH_

#include "Optimization/ErsatzMaterial.hh"

namespace CoupledField
{

/** Optimization via ShapeGradient and Level-Set method, not by Parametrization */
class ShapeGrad : public ErsatzMaterial
{
public:
  /** Up to now w/o parameters */
  ShapeGrad();
  
  virtual ~ShapeGrad()
  {
  }


  /** ... */
  void CalcObjectiveGradient(double* grad_out);

private:
  
  /** stores in DesignElement */
  void CalcTopGrad();
  
  /** evaluate */
  bool topgrad;
};


} // namespace


#endif /*SHAPEGRAD_HH_*/
