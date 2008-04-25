#include "Optimization/ShapeGrad.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include <string>

using namespace CoupledField;


DECLARE_LOG(shapeGrad)
DEFINE_LOG(shapeGrad, "shapeGrad")

ShapeGrad::ShapeGrad() : ErsatzMaterial()
{
  ParamNode* sg_pn = pn->Get("shapeGrad");
  
  topgrad = sg_pn != NULL && sg_pn->Has("topGrad") && sg_pn->Get("topGrad")->Get("enabled")->AsBool();
}


void ShapeGrad::CalcObjectiveGradient(double* grad_out)
{
  CalcTopGrad();
  
  
}
