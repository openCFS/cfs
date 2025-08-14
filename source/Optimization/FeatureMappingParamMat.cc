#include "Optimization/FeatureMappingParamMat.hh"
#include "Optimization/Design/FeatureMappingDesign.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(fmpm, "featureMappingPM")

using namespace CoupledField;

FeatureMappingParamMat::FeatureMappingParamMat()
{
  design->SetDesignMaterial(DesignMaterial::FEATURE_MAPPING_ANISO); // non param mat design variables constructor
}

void FeatureMappingParamMat::PostInit()
{
  ParamMat::PostInit();
  fmd = dynamic_cast<FeatureMappingDesign*>(domain->GetDesign());
  assert(fmd != nullptr);


}

void FeatureMappingParamMat::SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev)
{
  // set FeatureMappingDesign::AnisotropicGradientHelper()

  const FeatureMappingDesign::AnisoCurrentVariable& cur = fmd->aniso_current_variable;
  
  switch(cur.var_type)
  {
    case BaseDesignElement::NO_DERIVATIVE:
      SIMP::SetElementK(f, de, tf, app, mat_out, derivative, calcMode, ev);
      break;

    case DesignElement::FEATURE_MAPPING_PX:
    case DesignElement::FEATURE_MAPPING_PY:
    case DesignElement::FEATURE_MAPPING_QX:
    case DesignElement::FEATURE_MAPPING_QY:
    case DesignElement::FEATURE_MAPPING_P:  
    {
      assert(cur.var_idx == cur.var_type - DesignElement::FEATURE_MAPPING_PX); // does not hold for NO_DERIVATIVE
      OptimizationMaterial* mat = Optimization::context->mat;
      Matrix<double>& out = *dynamic_cast<Matrix<double>*>(mat_out);
      // in the end calls FeatureMappingDesign::GetAnisoMechTensor()
      mat->ComputeElementMatrix(out, mat->stiff.integrator, de->elem, false, cur.var_type );
      LOG_DBG(fmpm) << "SEK: f=" << f->ToString() << " asp=" << cur.var_type << " s=" << cur.var_idx 
                    << " afi=" << cur.feature_idx << " de=" << de->elem->elemNum << " -> " << mat_out->ToString();
      break;
    }

    default:
      assert(false);
  }
}
