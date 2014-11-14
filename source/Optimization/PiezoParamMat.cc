#include <assert.h>
#include <stddef.h>
#include <cmath>

#include "DataInOut/Logging/cfslog.hh"
#include "Domain/elem.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/elecPDE.hh"
#include "PDE/eqnMap.hh"
#include "Utils/StdVector.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/PiezoParamMat.hh"
#include "Optimization/TransferFunction.hh"
#include "boost/lexical_cast.hpp"

using namespace CoupledField;

using std::complex;

DECLARE_LOG(ppm)
DEFINE_LOG(ppm, "piezo_para_mat")


PiezoParamMat::PiezoParamMat() : PiezoSIMP()
{
  design->SetDesignMaterial(pn->Get("paramMat/designMaterial"), OptimizationMaterial::system.Parse(pn->Get("material")->As<std::string>()));
}

PiezoParamMat::~PiezoParamMat()
{
}


void PiezoParamMat::SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  // we assume to have no interpolation
  assert(tf->GetType() == TransferFunction::IDENTITY);
  assert(calcMode == STANDARD);

  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);

  DesignElement::Type dt = derivative ? de->GetType() : DesignElement::NO_DERIVATIVE;

  switch(app)
  {
  case MECH:
    out = piezo_mat_->MechStiffness(de->elem, false, de->multimaterial != NULL ? de->multimaterial->index : -1, dt);
    break;

  case ELEC:
    out = piezo_mat_->ElecStiffnessNeg(de, dt); // we need the -K_pp matrix
    break;

  case PIEZO_COUPLING:
    // out needs to be defined
    assert(out.GetNumCols() != out.GetNumRows());

    if(out.GetNumCols() > out.GetNumRows())
      out = piezo_mat_->CoupledStiffnessTransposed(de, dt);
    else
      out = piezo_mat_->CoupledStiffness(de, dt);
    break;

  default:
    assert(false);
    break;
  }

  LOG_DBG2(ppm) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " d=" << derivative << " dt=" << dt;
}

void PiezoParamMat::SetElementKMapping(DesignElement* de, BaseDesignElement::Type type, const TransferFunction* tf, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  // we assume to have no interpolation
  assert(tf->GetType() == TransferFunction::IDENTITY);
  assert(calcMode == STANDARD);

  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);

  DesignElement::Type dt = derivative ? type : DesignElement::NO_DERIVATIVE;

  switch(app)
  {
  case MECH:
    out = piezo_mat_->MechStiffness(de->elem, false, de->multimaterial != NULL ? de->multimaterial->index : -1, dt);
    break;

  case ELEC:
    out = piezo_mat_->ElecStiffnessNeg(de, dt); // we need the -K_pp matrix
    break;

  case PIEZO_COUPLING:
    // out needs to be defined
    assert(out.GetNumCols() != out.GetNumRows());

    if(out.GetNumCols() > out.GetNumRows())
      out = piezo_mat_->CoupledStiffnessTransposed(de, dt);
    else
      out = piezo_mat_->CoupledStiffness(de, dt);
    break;

  default:
    assert(false);
    break;
  }

  LOG_DBG2(ppm) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " d=" << derivative << " dt=" << dt;
}
