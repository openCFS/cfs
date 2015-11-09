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
#include "DataInOut/Logging/log.hpp"
#include "Utils/tools.hh"

namespace CoupledField {
class TransferFunction;
}  // namespace CoupledField

using namespace CoupledField;

DECLARE_LOG(em)

ParamMat::ParamMat() : ErsatzMaterial()
{
  // Note: this constructor is also called from constructor of ShapeOpt even when no ParamMat is used, in this case, nothing may be done
  
  if((method_ == PARAM_MAT || method_ == SHAPE_PARAM_MAT) && pn->Has("paramMat")){ 
    design->SetDesignMaterial(pn->Get("paramMat/designMaterial"), OptimizationMaterial::system.Parse(pn->Get("material")->As<std::string>()), this);
  }
  
  mech_mat_ = NULL; // set in PostInit()
}

void ParamMat::PostInit()
{
  if(pn->Has("filters"))
  {
    ParamNodeList list = pn->Get("filters")->GetList("filter");
    // this is save for design=polarization
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      if(structure_ == NULL)
        structure_ = new DesignStructure(this);
      structure_->SetFilter(list[i], this->optInfoNode);
    }
  }

  ErsatzMaterial::PostInit();
  
  mech_mat_ = dynamic_cast<MechMat*>(material); // just set in EM:PostInit()
  assert(mech_mat_ != NULL);
}


template <class T1, class T2>
void ParamMat::SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, bool derivative, CalcMode mode, double ev)
{
  // this is only called from CalcU1KU2 which is only used in derivative calculation (compliance, tracking, volume)
  // therefore we always return a derivative, de indicating which
  // for transient problems, this does also need to return the derivative of the mass matrix
  Matrix<T1>& out = dynamic_cast<Matrix<T1>& >(*mat_out);
  int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;
  switch(app)
  {
  case MECH:
  {
    const Matrix<T2>& tmp = dynamic_cast<const Matrix<T2>& >(mech_mat_->MechStiffness(de->elem, false, mm, derivative ? de->GetType() : DesignElement::NO_DERIVATIVE));
    Assign(out, tmp, 1.0);
    if(context.IsComplex())
    {
      AddMassToStiffness(tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, false, mode, ev); // no bimaterial

      // LOG_DBG3(simp) << "SetElementK: m_factor " << m_factor << " -> " << out.ToString();

      if(design->GetRegion(de->elem->regionId)->HasBiMaterial())
      {
        // rho^3 * E1 + (1-rho^3) * E2, in the derivative case 3*rho^2 * E1 - 3*rho^2 * E2
        AddMassToStiffness(tf, de, dynamic_cast<Matrix<Complex>& >(out), derivative, true, mode, ev); // bimaterial

        // LOG_DBG3(simp) << "SetElementK: m_bi_factor " << m_factor << " -> " << out.ToString();
      }
    }
    break;
  }
  case MASS:
  {
    const Matrix<T2>& tmp = dynamic_cast<const Matrix<T2>& >(mech_mat_->MechMass(de->elem, false, mm, derivative ? de->GetType() : DesignElement::NO_DERIVATIVE));
    Assign(out, tmp, 1.0);
    break;
  }
  default:
    Exception("Only mech and mass matrix are available for paramMat");
    break;
  }
  LOG_DBG3(em) << "PM:SEK de=" << de->ToString() << " app=" << application.ToString(app) << " d=" << derivative << " out=" << mat_out->ToString(0, false);
}


void ParamMat::SetElementKMapping(DesignElement* de, BaseDesignElement::Type type, const TransferFunction* tf, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  // this is only called from CalcU1KU2 which is only used in derivative calculation (compliance, tracking, volume)
  // therefore we always return a derivative, de indicating which
  // for transient problems, this does also need to return the derivative of the mass matrix
  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);
  int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;

  DesignElement::Type t = derivative ? type : DesignElement::NO_DERIVATIVE;

  switch(app)
  {
  case MECH:
    out = dynamic_cast<Matrix<double> &>(mech_mat_->MechStiffness(de->elem, false, mm, t));
    break;
  case MASS:
    out = dynamic_cast<Matrix<double> &>(mech_mat_->MechMass(de->elem, false, mm, t));
    break;
  default:
    Exception("Only mech and mass matrix are available for paramMat");
    break;
  }
  LOG_DBG3(em) << "PM:SEK de=" << de->ToString() << " d=" << derivative << " out=" << mat_out->ToString(0, false);
}

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template void ParamMat::SetElementK<double, double>(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev);
template void ParamMat::SetElementK<Complex, Complex>(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev);
template void ParamMat::SetElementK<Complex, double>(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev);
#endif

