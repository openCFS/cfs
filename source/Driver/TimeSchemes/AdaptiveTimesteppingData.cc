#include "AdaptiveTimesteppingData.hh"

namespace CoupledField {

void AdaptiveTimesteppingData::InitFromXml(PtrParamNode node) {
    enabled = true;
    dtMin   = node->Get("deltaTmin")->MathParse<Double>();
    dtMax   = node->Get("deltaTmax")->MathParse<Double>();
    sigma   = node->Get("cutbackfactor")->MathParse<Double>();

    std::string s = node->Get("Stepsizesmoothing")->As<std::string>();
    smoothing     = (s == "ON");
    minStepFactor = smoothing ? 0.2 : 0.25;

    PtrParamNode minSFNode = node->Get("minStepFactor", ParamNode::PASS);
    if (minSFNode) minStepFactor = minSFNode->MathParse<Double>();

    PtrParamNode tolNode = node->Get("directTol", ParamNode::PASS);
    tol = tolNode ? tolNode->MathParse<Double>() : 1e-6;

    PtrParamNode rtolNode = node->Get("relativeTol", ParamNode::PASS);
    if (rtolNode) {
        rtol = rtolNode->MathParse<Double>();
        PtrParamNode atolNode = node->Get("absoluteTol", ParamNode::PASS);
        atol = atolNode ? atolNode->MathParse<Double>() : 1e-10;
        tol  = 1.0;  // dimensionless threshold when rtol/atol are used
    }

    std::string scheme = node->Get("scheme")->As<std::string>();
    errorScheme = (scheme == "normalizedError") ? 2 : 1;
}

} // namespace CoupledField
