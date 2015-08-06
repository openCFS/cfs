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
    if(complex_)
    {
      // in SIMP we would call AddMassToStiffness(). As we assume to have no damping here and no bimaterial it is simpleer here!
      assert(mode != EIGENFREQ || (derivative == true && ev > 0)); // EIGENVALUE only for derivative
      assert(derivative == true);
      assert(pde->GetDamping(de->elem->regionId) == NONE);
      assert(de->multimaterial == NULL);
      double m_factor = mode == EIGENFREQ ? ev : 1.0;
      double omega = mode != EIGENFREQ ? 2.0 * M_PI * pde->GetSolveStep()->GetActFreq() : 1.0 ;  // todo: check with multiple excitation frequencies!
      assert(mode != EIGENFREQ || (omega == 1.0 && m_factor != 0)); // note that we might have very_small negative eigenvalues!
      T1 damp_mass = -1.0 *omega*omega*m_factor;
      Matrix<double> M;
      if(material->ComplexElementMatrix(de->elem->regionId)) // true for bloch
      {
        const Matrix<Complex>& T = dynamic_cast<const Matrix<Complex>&>(mech_mat_->MechMass(de->elem, false, -1, de->GetType())); // no bimat, multimatindx -1
                assert(IsNoise(T.GetPart(Global::IMAG).NormL2())); // up to now we have nothing interesting in the imaginary part of a mass matrix
        M = T.GetPart(Global::REAL);
      }
      else
        M = dynamic_cast<const Matrix<double>& >(mech_mat_->MechMass(de->elem, false, -1, de->GetType()));
      Add<T1, double>(out, damp_mass, M);
      LOG_DBG(em) << "PM:SEK de=" << de->ToString() << " app=" << application.ToString(app) << " d=" << derivative << " m_factor=" << m_factor << " omega=" << omega << " -> " << damp_mass << " K=" << tmp.NormL2() << " M=" << M.NormL2();
      LOG_DBG3(em) << "PM:SEK M=" << M.ToString();
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

