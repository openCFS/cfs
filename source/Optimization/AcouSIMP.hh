#pragma once
#include "Optimization/SIMP.hh"

namespace CoupledField
{
class AcouSIMP: public SIMP
{
public:
  AcouSIMP() {};

  virtual ~AcouSIMP() {};

private:

  void InitSecondMaterialCache() override;

  StdVector<Complex> GetExcitationPressureVector(Excitation& excite, Function* f) override;

  const RegionIdType GetExcitationRegion(Function* f, const std::string& attr = "surfRegion") override;

};

} // namespace
