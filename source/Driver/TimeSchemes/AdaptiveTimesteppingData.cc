#include "AdaptiveTimesteppingData.hh"

namespace CoupledField {

void AdaptiveTimesteppingData::InitFromXml(PtrParamNode node) {
    enabled_ = true;
    dtMin_   = node->Get("deltaTmin")->MathParse<Double>();
    dtMax_  = node->Get("deltaTmax")->MathParse<Double>();
    sigma_   = node->Get("cutbackfactor")->MathParse<Double>();

    std::string s = node->Get("Stepsizesmoothing")->As<std::string>();
    smoothing_    = (s == "ON");
    minStepFactor_ = smoothing_ ? 0.2 : 0.25;

    PtrParamNode minSFNode = node->Get("minStepFactor", ParamNode::PASS);
    if (minSFNode) minStepFactor_ = minSFNode->MathParse<Double>();

    PtrParamNode tolNode = node->Get("directTol", ParamNode::PASS);
    tol_ = tolNode ? tolNode->MathParse<Double>() : 1e-6;

    PtrParamNode rtolNode = node->Get("relativeTol", ParamNode::PASS);
    if (rtolNode) {
        rtol_ = rtolNode->MathParse<Double>();
        PtrParamNode atolNode = node->Get("absoluteTol", ParamNode::PASS);
        atol_ = atolNode ? atolNode->MathParse<Double>() : 1e-10;
        tol_  = 1.0;  // dimensionless threshold when rtol/atol are used
    }

    std::string scheme = node->Get("scheme")->As<std::string>();
    errorScheme_ = (scheme == "normalizedError") ? 2 : 1;
}

bool AdaptiveTimesteppingData::is_error_finite(Double Error)
{
    return std::isfinite(Error);
}

void AdaptiveTimesteppingData::mark_saturated()
{
    postSaturationMode_ = true;
    healthyStepCount_   = 0;
}

Double AdaptiveTimesteppingData::apply_post_saturation_cap(
        Double h_next, Double h, Double local_error, bool accepted)
{
    if (!postSaturationMode_)
        return h_next;

    const Double cap = 1.5;
    h_next = std::min(h_next, h * cap);

    if (accepted && local_error <= tol_) {
        if (++healthyStepCount_ >= 3)
            postSaturationMode_ = false;
    } else {
        healthyStepCount_ = 0;
    }
    return h_next;
}

 Double AdaptiveTimesteppingData::standardStepsize(bool* accepted, Double local_error_, Double dtCurrent_)
  {
    Double Rtol  = tol_;
    Double est   = local_error_;
    Double h     = dtCurrent_;

    const Double z_U = 0.1;
    const Double z_S = sigma_ + 1;  // zs has to be between 1 and 2
    const Double F_U      = 1.5;
    Double h_next;

    if (est == 0.0) {
      h_next   = F_U * h;
      *accepted = true;
    } else {
      Double z = z_S * std::pow((est / Rtol), (1.0 / 3.0));
      if (z <= z_U) {
        h_next = F_U * h;  *accepted = true;
      } else if (z <= z_S) {
        h_next = h / z;    *accepted = true;
      } else {
        h_next = h / z;    *accepted = false;
      }
    }

    return h_next;
  }

  Double AdaptiveTimesteppingData::smoothStepsize(bool* accepted, Double local_error_,Double dtCurrent_ )
  {
    Double Rtol  = tol_;
    Double est   = local_error_;
    Double h     = dtCurrent_;
    const Double F_U = 1.5;
    Double h_next;
    Double sigma       = sigma_;
    Double prev_error_ = prevError_;
    double k = 3; // For BDF2
    Double alpha = 0.3/k ;
    Double beta  = 0.6/k;

    if (est == 0.0) {
      h_next = F_U * h;
      *accepted = true;
    } else {
      Double ratio_I = std::pow(sigma * Rtol / est, alpha);
      Double ratio_P = 1.0;

      if (prev_error_ > 0.0) {
        ratio_P = std::pow(prev_error_ / est, beta);
      }
      Double factor = ratio_I * ratio_P;

      if (est > Rtol && factor > 0.95)
        {
            *accepted = true;
            h_next = h*0.9;   // keep current dt — don't reduce to dtMin
            toleranceNotReachable_ = true;
            return h_next;
        }
      h_next = h * factor;

      *accepted = (est <= Rtol);
    }
  
    return h_next;
  }

} // namespace CoupledField
