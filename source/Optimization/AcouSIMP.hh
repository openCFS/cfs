#ifndef ACOUSIMP_HH_
#define ACOUSIMP_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField
{
class AcouSIMP: public SIMP
{
public:
  AcouSIMP();

  virtual ~AcouSIMP();

  void PostInit() override;

private:
  const Complex GetExcitationPressure(Function* f) override;

  const RegionIdType GetExcitationRegion(Function* f, const std::string& attr = "surfRegion") override;

};

} // namespace

#endif /*ACOUSIMP_HH_*/