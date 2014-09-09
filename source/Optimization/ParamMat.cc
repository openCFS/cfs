#include <assert.h>
#include <stddef.h>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/Exception.hh"
#include "MatVec/denseMatrix.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/ParamMat.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

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
  PtrParamNode simp_pn = pn->Get("SIMP", ParamNode::PASS);

  // There might be a filter regularization based on the design element.
  if(simp_pn)
  {
    if(simp_pn->HasByVal("regularization", "type", "filter"))
    {
      ParamNodeList list = simp_pn->Get("regularization")->GetList("filter");
      // this is save for design=polarization
      for(unsigned int i = 0; i < list.GetSize(); i++)
      {
        if(structure_ == NULL)
          structure_ = new DesignStructure(this);
        structure_->SetFilters(list[i], this->optInfoNode);
      }
    }
    else
    {
      if(simp_pn->Has("regularization"))
        throw Exception("regularization not implemented");
    }
  }

  ErsatzMaterial::PostInit();
  
  mech_mat_ = dynamic_cast<MechMat*>(material); // just set in EM:PostInit()
  assert(mech_mat_ != NULL);
}

void ParamMat::SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  // this is only called from CalcU1KU2 which is only used in derivative calculation (compliance, tracking, volume)
  // therefore we always return a derivative, de indicating which
  // for transient problems, this does also need to return the derivative of the mass matrix
  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);
  int mm = de->multimaterial != NULL ? de->multimaterial->index : -1;

  switch(app)
  {
  case MECH:
    out = dynamic_cast<Matrix<double>& >(mech_mat_->MechStiffness(de->elem, false, mm, derivative ? de->GetType() : DesignElement::NO_DERIVATIVE));
    break;
  case MASS:
    out = mech_mat_->MechMass(de->elem, false, mm, derivative ? de->GetType() : DesignElement::NO_DERIVATIVE);
    break;
  default:
    Exception("Only mech and mass matrix are available for paramMat");
    break;
  }
  LOG_DBG3(em) << "PM:SEK de=" << de->ToString() << " d=" << derivative << " out=" << mat_out->ToString(0, false);
}
