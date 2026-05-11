#include "AdaptiveTimesteppingData.hh"

namespace CoupledField {

void AdaptiveTimesteppingData::InitFromXml(PtrParamNode node) {
    enabled_ = true;
    dtMin_   = node->Get("deltaTmin")->MathParse<Double>();
    dtMax_  = node->Get("deltaTmax")->MathParse<Double>();
    PtrParamNode sigmaNode = node->Get("cutbackfactor", ParamNode::PASS);
    if (sigmaNode) sigma_ = sigmaNode->MathParse<Double>();

    std::string ctrl = node->Get("controller")->As<std::string>();
    if      (ctrl == "PI")  controllerType_ = 1;
    else if (ctrl == "PID") controllerType_ = 2;
    else                    controllerType_ = 0;  // "I"
    minStepFactor_ = (controllerType_ > 0) ? 0.2 : 0.25;
  
    if(controllerType_ == 2)
    {
      minStepFactor_ = 0.10;
    }

    PtrParamNode startNode = node->Get("StartAtmin", ParamNode::PASS);
    startFromDtMin_ = startNode && (startNode->As<std::string>() == "ON");

    PtrParamNode warmUpNode = node->Get("warmUpLTE", ParamNode::PASS);
    if (warmUpNode) {
        if (startFromDtMin_)
            EXCEPTION("warmUpLTE and StartAtmin=\"ON\" are mutually exclusive.");
        warmUpEnabled_  = true;
        warmUpLTETarget_ = warmUpNode->MathParse<Double>();
        inWarmUpPhase_  = true;
    }

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

void AdaptiveTimesteppingData::registerFieldLTE(double lte)
{
    fieldLocalErrors_.push_back(lte);
    localError_ = getControllingError();
    lteCollected_ = true;
}

double AdaptiveTimesteppingData::getControllingError() const
{
    if (fieldLocalErrors_.empty()) return localError_;
    return *std::max_element(fieldLocalErrors_.begin(), fieldLocalErrors_.end());
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

 Double AdaptiveTimesteppingData::iController(bool* accepted, Double local_error_, Double dtCurrent_)
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

  Double AdaptiveTimesteppingData::piController(bool* accepted, Double local_error_,Double dtCurrent_ )
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

Double AdaptiveTimesteppingData::pidController(bool* accepted,
        Double local_error_, Double dtCurrent_)
{
    // H312 PID (Söderlind 2005, eq. 38): kk_I=2/9, k=3 → kI=2/27
    // h_{n+1} = h × (σ·tol/r̂_n)^e1 × (σ·tol/r̂_{n-1})^e2 × (σ·tol/r̂_{n-2})^e1
    // Fall back to PI when not yet enough history.
    if (prevPrevError_ <= 0.0)
        return piController(accepted, local_error_, dtCurrent_);

    const Double kI  = (2.0/9.0) / 3.0;  // kk_I/k = 2/27
    const Double e1  = kI / 4.0;
    const Double e2  = kI / 2.0;

    Double h   = dtCurrent_;
    Double est = local_error_;

    if (est == 0.0) { *accepted = true; return 1.5 * h; }

    Double base   = sigma_ * tol_;
    Double h_next = h
        * std::pow(base / est,            e1)
        * std::pow(base / prevError_,     e2)
        * std::pow(base / prevPrevError_, e1);

    *accepted = (est <= tol_);
    return h_next;
}

} // namespace CoupledField
