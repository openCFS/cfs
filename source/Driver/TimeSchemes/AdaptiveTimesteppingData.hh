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
    int    controllerType_   = 0;    // 0 = I, 1 = PI, 2 = PID
    bool   startFromDtMin_  = false;

    // per-step state (written by TimeSchemeGLM, read by TransientDriver)
    double localError_            = 0.0;
    bool   stepRejected_          = false;
    bool   toleranceNotReachable_ = false;
    double prevRetryError_        = 0.0;  // LTE from previous retry; 0 = first attempt
    int    consecutiveNaN_       = 0;     // consecutive NaN-solution steps; abort after threshold

    // Multi-field LTE: each field registers via FinishStepLTE(); FinishStep() makes one consistent decision.
    // fieldLocalErrors_ per field; getControllingError() returns max.
    std::vector<double> fieldLocalErrors_;
    bool lteCollected_     = false;  // true once at least one field has registered
    bool stepDecisionMade_ = false;  // true once the first field called ComputeAdaptiveStepSize

    void   registerFieldLTE(double lte);
    double getControllingError() const;

    // per-step feedback (written by TransientDriver, read by TimeSchemeGLM)
    double prevError_     = 0.0;
    double prevPrevError_ = 0.0;

    bool is_error_finite(Double Error);

    //! Post-saturation growth limiter: caps h_next to 1.5×h for 3 steps after a force-accept.
    void   mark_saturated();

    //! Call mark_saturated() at every toleranceNotReachable site; apply_post_saturation_cap() after BDF2 cap.
    Double apply_post_saturation_cap(Double h_next, Double h,
                                     Double local_error, bool accepted);

    bool postSaturationMode_ = false;
    int  healthyStepCount_   = 0;

    struct StepResult {
        double h_next;
        bool   accepted;
    };

    StepResult computeNextStep(double h, double h_prev, double est,
                               double dtMin, double dtMax);

    // bound-saturation counters (reset each sequence step by TransientDriver)
    int stepsAtDtMin_      = 0;
    int stepsAtDtMax_      = 0;
    int totalAcceptedSteps_ = 0;

    // warm-up phase: run at fixed deltaT until LTE/tol drops to warmUpLTETarget_
    bool   warmUpEnabled_   = false;
    double warmUpLTETarget_ = 2.0;
    bool   inWarmUpPhase_   = false;

    //! I-controller: simple single-step power-law formula
    Double iController(bool* accepted, Double local_error_, Double dtCurrent_);

    //! PI.3.4 controller (Söderlind 2002) (includes error rejektion)
    Double piController(bool* accepted, Double local_error_, Double dtCurrent_);

    //! H312 PID controller (Söderlind 2005, eq. 38) — for non-smooth/noisy problems
    Double pidController(bool* accepted, Double local_error_, Double dtCurrent_);

    //! 3-step LTE stability factor: 1=stable (full controller freedom), <1=damp step change.
    //! Detects monotone LTE growth and large oscillations; applies to all controllers.
    double lteStabilityFactor(double est) const;
};

} // namespace CoupledField
