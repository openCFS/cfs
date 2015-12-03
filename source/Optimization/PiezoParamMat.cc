#include <assert.h>
#include <stddef.h>
#include <cmath>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/ElecPDE.hh"
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
  design->SetDesignMaterial(pn->Get("paramMat/designMaterial"), OptimizationMaterial::system.Parse(pn->Get("material")->As<std::string>()), this);
}

PiezoParamMat::~PiezoParamMat()
{
}


void PiezoParamMat::SetElementK(DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode, double ev)
{
  // we assume to have no interpolation
  assert(tf->GetType() == TransferFunction::IDENTITY);
  assert(calcMode == STANDARD);

  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);
  PiezoElecMat* pem = dynamic_cast<PiezoElecMat*>(context->mat); // don't cache!

  DesignElement::Type dt = derivative ? de->GetType() : DesignElement::NO_DERIVATIVE;

  switch(app)
  {
  case App::MECH:
    out = dynamic_cast<const Matrix<double>& >(pem->MechStiffness(de->elem, false, de->multimaterial != NULL ? de->multimaterial->index : -1, dt));
    break;

  case App::ELEC:
    out = pem->ElecStiffnessNeg(de, dt); // we need the -K_pp matrix
    break;

  case App::PIEZO_COUPLING:
    // out needs to be defined
    assert(out.GetNumCols() != out.GetNumRows());

    if(out.GetNumCols() > out.GetNumRows())
      out = pem->CoupledStiffnessTransposed(de, dt);
    else
      out = pem->CoupledStiffness(de, dt);
    break;

  default:
    assert(false);
    break;
  }

  LOG_DBG2(ppm) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " d=" << derivative << " dt=" << dt;
}

void PiezoParamMat::SetElementKMapping(DesignElement* de, BaseDesignElement::Type type, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  // we assume to have no interpolation
  assert(tf->GetType() == TransferFunction::IDENTITY);
  assert(calcMode == STANDARD);
  PiezoElecMat* pem = dynamic_cast<PiezoElecMat*>(context->mat); // don't cache outside the function because of multiple seqeuence issues

  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);

  DesignElement::Type dt = derivative ? type : DesignElement::NO_DERIVATIVE;

  switch(app)
  {
  case App::MECH:
    out = dynamic_cast<Matrix<double> &>(pem->MechStiffness(de->elem, false, de->multimaterial != NULL ? de->multimaterial->index : -1, dt));
    break;

  case App::ELEC:
    out = pem->ElecStiffnessNeg(de, dt); // we need the -K_pp matrix
    break;

  case App::PIEZO_COUPLING:
    // out needs to be defined
    assert(out.GetNumCols() != out.GetNumRows());

    if(out.GetNumCols() > out.GetNumRows())
      out = pem->CoupledStiffnessTransposed(de, dt);
    else
      out = pem->CoupledStiffness(de, dt);
    break;

  default:
    assert(false);
    break;
  }

  LOG_DBG2(ppm) << "PiezoSIMP::SetElementK elem: " << de->elem->elemNum << " app: " << application.ToString(app) << " d=" << derivative << " dt=" << dt;
}
