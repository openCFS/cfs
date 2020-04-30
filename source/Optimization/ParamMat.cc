#include <assert.h>
#include <stddef.h>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/Exception.hh"
#include "MatVec/DenseMatrix.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/ParamMat.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/tools.hh"

namespace CoupledField {
class TransferFunction;
}  // namespace CoupledField

using namespace CoupledField;

EXTERN_LOG(em)

ParamMat::ParamMat() : ErsatzMaterial()
{
  // Note: this constructor is also called from constructor of ShapeOpt even when no ParamMat is used, in this case, nothing may be done
  
  if((method_ == PARAM_MAT || method_ == SHAPE_PARAM_MAT) && pn->Has("paramMat")){ 
    design->SetDesignMaterial(pn->Get("paramMat/designMaterial"), OptimizationMaterial::system.Parse(pn->Get("material")->As<std::string>()), this);
  }
}

void ParamMat::PostInit()
{
  ErsatzMaterial::PostInit();
}


template <class T1, class T2>
void ParamMat::SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode mode, double ev)
{
  // this is only called from CalcU1KU2 which is only used in derivative calculation (compliance, tracking, volume)
  // therefore we always return a derivative, de indicating which
  // for transient problems, this does also need to return the derivative of the mass matrix
  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);
  int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;
  MechMat* mech_mat = dynamic_cast<MechMat*>(context->mat); // don't cache

  switch(app)
  {
  case App::MECH:
  case App::BUCKLING:
  {
    const Matrix<T2>& tmp = dynamic_cast<const Matrix<T2>& >(mech_mat->Stiffness(de->elem, false, mm, derivative ? de->GetType() : DesignElement::NO_DERIVATIVE));
    Assign(out, tmp, 1.0);

    if(context->DoBuckling())
    {
//      out *= -ev;

      AddGeometricStiffnessToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, false, mode, ev);
      // LOG_DBG3(simp) << "SetElementK: m_factor " << m_factor << " -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddGeometricStiffnessToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, true, mode, ev); // bimaterial
        // LOG_DBG3(simp) << "SetElementK: m_bi_factor " << m_factor << " -> " << out.ToString();
      }
    }

    if(context->IsComplex())
    {
      AddMassToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, false, mode, ev); // no bimaterial
      // LOG_DBG3(simp) << "SetElementK: m_factor " << m_factor << " -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddMassToStiffness(f->ctxt, tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, true, mode, ev); // bimaterial
        // LOG_DBG3(simp) << "SetElementK: m_bi_factor " << m_factor << " -> " << out.ToString();
      }
    }
    break;
  }
  case App::MASS:
  {
    const Matrix<T2>& tmp = dynamic_cast<const Matrix<T2>& >(mech_mat->Mass(de->elem, false, mm, derivative ? de->GetType() : DesignElement::NO_DERIVATIVE));
    Assign(out, tmp, 1.0);
    break;
  }
  default:
    Exception("Only mech and mass matrix are available for paramMat");
    break;
  }
  LOG_DBG3(em) << "PM:SEK de=" << de->ToString() << " app=" << application.ToString(app) << " d=" << derivative << " out=" << mat_out->ToString(0, false);
}


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template void ParamMat::SetElementK<double, double>( Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev);
template void ParamMat::SetElementK<Complex, Complex>(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev);
template void ParamMat::SetElementK<Complex, double>(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev);
#endif

