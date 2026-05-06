#pragma once
#include <algorithm>
#include <vector>
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

class AdaptiveTimesteppingData {
public:
    //! AdaptiveTimesteppingData contains all Data neccesery for the Adaptive Schemes, is part of domain. Includes xml parsing.
    void InitFromXml(PtrParamNode node);

    // config (read-only after init)
    bool   enabled_         = false;
    int    errorScheme_     = 1;       // 1 = maxlocalError, 2 = normalizedError
    double tol_             = 1e-6;
    double dtMin_           = 0.0;
    double dtMax_           = 1e10;
    double sigma_           = 0.9;
    double rtol_            = 0.0;     // 0 = disabled
    double atol_            = 0.0;
    double minStepFactor_   = 0.2;
    bool   smoothing_       = false;
    bool   startFromDtMin_  = false;

    // per-step state (written by TimeSchemeGLM, read by TransientDriver)
    double localError_            = 0.0;
    bool   stepRejected_          = false;
    bool   toleranceNotReachable_ = false;
    double prevRetryError_        = 0.0;  // LTE from previous retry; 0 = first attempt

    // Multi-field LTE collection: each field registers its LTE via FinishStepLTE()
    // before FinishStep() makes a single consistent step-size decision.
    // fieldLocalErrors_ holds one entry per field; getControllingError() returns
    // the worst (currently max) — future: per-domain field selection.
    std::vector<double> fieldLocalErrors_;
    bool lteCollected_     = false;  // true once at least one field has registered
    bool stepDecisionMade_ = false;  // true once the first field called ComputeAdaptiveStepSize

    void   registerFieldLTE(double lte);
    double getControllingError() const;

    // per-step feedback (written by TransientDriver, read by TimeSchemeGLM)
    double prevError_ = 0.0;

    bool is_error_finite(Double Error);

    //! Post-saturation growth limiter: activated by mark_saturated() whenever
    // a force-accept (toleranceNotReachable) occurs.  
    void   mark_saturated();

    //! apply_post_saturation_cap()
    // restricts h_next to at most 1.5×h until 3 consecutive healthy steps have
    // been seen.  Call mark_saturated() at every toleranceNotReachable site;
    // call apply_post_saturation_cap() after the BDF2 stability cap, before
    // the dtMin/dtMax clamp.
    Double apply_post_saturation_cap(Double h_next, Double h,
                                     Double local_error, bool accepted);

    bool postSaturationMode_ = false;
    int  healthyStepCount_   = 0;

    // warm-up phase: run at fixed deltaT until LTE/tol drops to warmUpLTETarget_
    bool   warmUpEnabled_   = false;
    double warmUpLTETarget_ = 2.0;
    bool   inWarmUpPhase_   = false;
    int    warmUpRampSteps_ = 0;   // counts accepted adaptive steps during post-warm-up ramp (0-3)

    //! Calculates next Step (Standart fomular)
    Double standardStepsize(bool* accepted, Double local_error_, Double dtCurrent_);
    
    //! Calculates next Step (Söderlin PI Controller)
    Double smoothStepsize(bool* accepted, Double local_error_,Double dtCurrent_ );
};

} // namespace CoupledField
