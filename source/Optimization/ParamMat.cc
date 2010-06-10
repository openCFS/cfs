#include "Optimization/ParamMat.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/OptimizationMaterial.hh"

using namespace CoupledField;

ParamMat::ParamMat() : ErsatzMaterial()
{
  // Note: this constructor is also called from constructor of ShapeOpt even when no ParamMat is used, in this case, nothing may be done
  if(pn->Has("paramMat")){ 
    design->SetDesignMaterial(pn->Get("paramMat")->Get("designMaterial"));  
  }
  
  mech_mat_ = NULL; // set in PostInit()
}

void ParamMat::PostInit()
{
  ErsatzMaterial::PostInit();
  
  mech_mat_ = dynamic_cast<OptMechMat*>(material); // just set in EM:PostInit()
  assert(mech_mat_ != NULL);
}

void ParamMat::SetElementK(DesignElement* de, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative)
{
  // this is only called from CalcU1KU2 which is only used in derivative calculation (compliance, tracking, volume)
  // therefore we always return a derivative, de indicating which
  // for transient problems, this does also need to return the derivative of the mass matrix
  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);
  switch(app){
  case MECH:
    out = mech_mat_->MechStiffness(de->elem, de->GetType());
    break;
  case MASS:
    out = mech_mat_->MechMass(de->elem, de->GetType());
    break;
  default:
    Exception("Only mech and mass matrix are available for paramMat");
  }
}
