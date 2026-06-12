#pragma once
#include <algorithm>
#include <cmath>
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

    //! Per-DOF tolerance scaling: lte / (atol + maxAbsY·rtol) in rtol/atol mode, raw lte otherwise.
    double scaledLTE(double lte, double maxAbsY) const {
        return (rtol_ > 0.0) ? lte / (atol_ + maxAbsY * rtol_) : lte;
    }

    //! Accumulates per-DOF (scaled) LTE values into the configured error norm.
    struct ErrorNorm {
        int    scheme;        // 1 = maxlocalError, 2 = normalizedError (RMS)
        double maxv = 0.0;
        double sum  = 0.0;
        std::size_t n = 0;
        void   add(double lte) { if (scheme == 2) sum += lte*lte; if (lte > maxv) maxv = lte; n++; }
        double result() const  { return (scheme == 2 && n > 0) ? std::sqrt(sum/n) : maxv; }
    };
    ErrorNorm newErrorNorm() const { return ErrorNorm{errorScheme_}; }

    //! Warm-up hold: true while running at fixed deltaT until LTE/tol <= warmUpLTETarget_.
    //! Resets localError_ and the caller's scheme-local error if LTE is non-finite.
    bool warmUpHold(double& schemeLocalError);

    //! Scheme-dependent stability cap on step growth h_next/h; 0 = unset (falls back to 1+sqrt(2)).
    //! Set at adaptive activation; with mixed schemes the most restrictive bound wins.
    double maxGrowthRatio_ = 0.0;
    void restrictMaxGrowthRatio(double r) {
        maxGrowthRatio_ = (maxGrowthRatio_ > 0.0) ? std::min(maxGrowthRatio_, r) : r;
    }

    //! LTE-trend damping (lteStabilityFactor) on/off. Newmark's ZX estimate alternates with
    //! 2*dt (ü update has a -ü_n term for beta=1/4), which falsely triggers the damper.
    bool lteDampingEnabled_ = true;

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
