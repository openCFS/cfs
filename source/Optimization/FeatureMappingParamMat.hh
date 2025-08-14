#ifndef OPTIMIZATION_FEATURE_MAPPING_PARAM_MAT_HH_
#define OPTIMIZATION_FEATURE_MAPPING_PARAM_MAT_HH_

#include "Optimization/ParamMat.hh"

namespace CoupledField 
{

class FeatureMappingDesign;

/** Anisotropic feature mapping. 
 * FeatureMappingDesign (derived from DesignSpace) creates the material
 * tensor (and it's derivative) itself. 
 * FeatureMappingParamMat (indirectly derived from SIMP) uses some ParamMat features 
 * to provide SetElementK() for the complex sensitivity analysis of aniostropic feature mapping */  
class FeatureMappingParamMat : public ParamMat
{
public:
  FeatureMappingParamMat();

  virtual ~FeatureMappingParamMat() { }

  void PostInit() override;

  /** In the anisotropic case we are ParamMat use here FeatureMappingDesign::CalcAnisoTensor() based on the aggregated oriented features in the FeatureToGradient function.
   * However, when the optimizer evaluates the gradient of state functions like the gradient to compute d_J/d_rho, we have here a dummy implementation
   * as what actually is essential, is FeatureMappingDesign::FeatureToGradient() 
   * @see aniso_shape_derivative */
  void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, ErsatzMaterial::CalcMode calcMode = STANDARD, double ev = -1.0) override;

private:

  FeatureMappingDesign* fmd = nullptr; // set by PostInit()
};

} // end of namespace

#endif /* OPTIMIZATION_FEATURE_MAPPING_PARAM_MAT_HH_ */
