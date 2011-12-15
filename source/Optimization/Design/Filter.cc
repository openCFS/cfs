#include <assert.h>
#include <cmath>
#include <limits>
#include <ostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/Filter.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

DECLARE_LOG(ds)

using namespace CoupledField;

Filter::Filter()
{
  type_          = NO_FILTERING;
  sensitivity_   = SIGMUND;
  density_       = STANDARD;
  beta_          = -1.0;
  heaviside_corr = -1.0;
  eta            = -1.0;
  explicit_lower_bound_ = std::numeric_limits<double>::min();
}

double Filter::GetLowerBound(const DesignElement* de) const
{
  if(explicit_lower_bound_ != std::numeric_limits<double>::min())
    return explicit_lower_bound_;
  else
    return de->GetLowerBound();
}


void Filter::SetBeta(double val, const DesignSpace* space)
{
  beta_ = val;
  // is there a need to set heavside_corr?
  if(type_ != DENSITY || density_ != HEAVISIDE) return;

  // in case we do bimaterial optimization disable the correction by setting it to 0 which makes it transparent
  assert(space->regions.GetSize() == 1 && space->regions[0].GetSize() == 1);
  if(space->regions[0][space->GetRegionId()].HasBiMaterial())
  {
    this->heaviside_corr = 0.0;
  }
  else
  {
    double lb = GetLowerBound(&(space->data[0]));
    // we assume a constant lower bound, if this is not he case we need another implementation
    for(unsigned int i = 0, n = space->data.GetSize(); i < n; i++)
      assert(GetLowerBound(&(space->data[i])) == lb);

    DesignElement de = DesignElement();
    de.simp = new SIMPElement(&de);
    de.simp->filter = *this;
    de.elem = space->data[0].elem;
    explicit_lower_bound_ = lb;
    de.SetDesign(lb);

    // find haviside_corr by a simple bisection (results in 0.2 ... 0.99999)
    double lower = 1e-6;
    double upper = 2.0;
    while(upper - lower > 1e-6) // so close that it doesn't matter in the end which value [lower:upper] is best
    {
      // this is nonsense!! TODO
      // why not like
      // double mid = 0.5 * (ub + lb);
      // double test = tf->Transform(NULL, DesignElement::PLAIN, mid);
      // if(test < physical)

      double mid = 0.5 * (upper + lower);
      de.simp->filter.heaviside_corr = 0.5 * (upper + mid);
      double err_u = lb - de.simp->CalcHeaviside(lb);
      de.simp->filter.heaviside_corr = 0.5 * (mid + lower);
      double err_l = lb - de.simp->CalcHeaviside(lb);
      if(std::abs(err_u) < std::abs(err_l))
        lower = mid;
      else
        upper = mid;
      LOG_DBG2(ds) << "Filter:SB val=" << val << " hc=" << de.simp->filter.heaviside_corr << " lower=" << lower
                   << " beta=" << de.simp->filter.GetBeta() << " upper=" << upper << " err_l=" << err_l << " err_u=" << err_u
                   << " H(" << lb << ")=" << de.simp->CalcHeaviside(lb);
    }
    this->heaviside_corr = de.simp->filter.heaviside_corr;
  }
}
