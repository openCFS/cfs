#pragma once
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

class AdaptiveTimesteppingData {
public:
    void InitFromXml(PtrParamNode node);

    // config (read-only after init)
    bool   enabled         = false;
    int    errorScheme     = 1;       // 1 = maxlocalError, 2 = normalizedError
    double tol             = 1e-6;
    double dtMin           = 0.0;
    double dtMax           = 1e10;
    double sigma           = 0.9;
    double rtol            = 0.0;     // 0 = disabled
    double atol            = 0.0;
    double minStepFactor   = 0.25;
    bool   smoothing       = false;

    // per-step state (written by TimeSchemeGLM, read by TransientDriver)
    double localError            = 0.0;
    bool   stepRejected          = false;
    bool   toleranceNotReachable = false;

    // per-step feedback (written by TransientDriver, read by TimeSchemeGLM)
    double prevError = 0.0;
};

} // namespace CoupledField
