#ifndef HEATSIMP_HH_
#define HEATSIMP_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField
{
class HeatSIMP: public SIMP
{
public:
  HeatSIMP() {};

  virtual ~HeatSIMP() {};

private:

  void InitSecondMaterialCache() override
  {
    AddSecondMaterialCache(MaterialClass::THERMIC, MaterialType::HEAT_CAPACITY);
    AddSecondMaterialCache(MaterialClass::THERMIC, MaterialType::DENSITY);
  }
};

} // namespace

#endif /*HEATSIMP_HH_*/