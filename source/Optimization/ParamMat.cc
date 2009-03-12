#include "Optimization/ParamMat.hh"
#include "Optimization/DesignSpace.hh"

using namespace CoupledField;

ParamMat::ParamMat() : ErsatzMaterial()
{
  ParamNode* pmpn = pn->Get("paramMat");

  design->SetDesignMaterial(pmpn->Get("designMaterial"));  
}

void ParamMat::SetElementK(DesignElement* de, Application app, DenseMatrix* mat_out)
{
  // this is only called from CalcU1KU2 which is only used in derivative calculation (compliance, tracking, volume)
  // therefore we always return a derivative, de indicating which
  Matrix<double>& out = dynamic_cast<Matrix<double>& >(*mat_out);
  const Matrix<double> mechStiffness = MechStiffness(de->elem, de->GetType());
  out.Resize(mechStiffness.GetNumRows(), mechStiffness.GetNumCols());
  Assign(out, mechStiffness, 1);
}
