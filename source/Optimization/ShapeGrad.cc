#include "Optimization/ShapeGrad.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Domain/domain.hh"
#include "PDE/mechPDE.hh"
#include "Forms/mechStressStrain.hh"
#include "Forms/linElastInt.hh"
#include "Utils/nodestoresol.hh"
#include <string>
#include <memory>

namespace CoupledField
{

DECLARE_LOG(shapeGrad)
DEFINE_LOG(shapeGrad, "shapeGrad")

ShapeGrad::ShapeGrad() : ErsatzMaterial(), mech_mat_(NULL)
{
}

void ShapeGrad::GetSubTensorType(SubTensorType &stt) const
{ 
  String2Enum(mech->GetSubType(), stt);
}

void ShapeGrad::GetMaterialParameters(double &lambda, double &mu) const
{
  // TODO: extend for multi-region-optimization if necessary
  BaseMaterial* material = mech->getPDEMaterialData()[regionIds[0]];
  material->GetScalar(lambda, MECH_LAME_LAMBDA, Global::REAL);
  material->GetScalar(mu, MECH_LAME_MU, Global::REAL);
  LOG_DBG3(shapeGrad) << "lame parameters:  lambda = " << lambda << ", mu = " << mu;
}

void ShapeGrad::GetStrainsOnElement(Vector<double> &forward_strain, Vector<double> &adjoint_strain, 
                                    const unsigned int e, const SubTensorType type) const
{
  Vector<double> intPoint;
  Matrix<double> elem_sol_forward_matrix;
  Matrix<double> elem_sol_adjoint_matrix;

  // we need a NodeStoreSol to get a element*matrix*. We use
  // the current one and explicitly set the raw solution vector from the forward problem
  NodeStoreSol<double>* node_store_sol = dynamic_cast<NodeStoreSol<double>* >(ErsatzMaterial::mech->getPDESolution());
  
  // for now, the application we consider is only mech
  // Application app1 = MECH;
  // Application app2 = MECH;
  
  // get current design element
  DesignElement* de = &design->data[e];
  
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  elemList.SetElement(de->elem);
  EntityIterator it = elemList.GetIterator();
  
  it.GetElem()->ptElem->GetCoordMidPoint(intPoint);

  const int idx(0);
  
  // set element solution to matrix
	node_store_sol->GetElemSolutionAsMatrix(elem_sol_forward_matrix, it, 0, forward.data[idx]->raw[MECH]);
	node_store_sol->GetElemSolutionAsMatrix(elem_sol_adjoint_matrix, it, 0, adjoint.data[idx]->raw[MECH]);

  BaseMaterial* material = mech->getPDEMaterialData()[regionIds[idx]]; // TODO: extend for multi-region-optimization if necessary
  
  // will contain the calculated strains = partial derivatives of displacement u
  shared_ptr<MechStressStrain<double> > strain(new MechStressStrain<double>(material,type));

  strain->SetIntPoint(intPoint);

  // calculate strain for forward problem
  strain->SetActElemSol(elem_sol_forward_matrix);
  strain->CalcStrainVec(forward_strain, 1, it);

  // FIXME adjoint problem
  //strain->SetActElemSol(elem_sol_adjoint_matrix);
  //strain->CalcStrainVec(adjoint_strain, 1, it);
  adjoint_strain = forward_strain;

  LOG_DBG3(shapeGrad) << "elem " << de->elem->elemNum << " forward strain: " << forward_strain.ToString();
  LOG_DBG3(shapeGrad) << "elem " << de->elem->elemNum << " adjoint strain: " << adjoint_strain.ToString();
}

linElastInt* ShapeGrad::getBDBForm()
{
  return dynamic_cast<linElastInt*>(GetForm(regionIds[0], mech, mech, "linElastInt"));
}


void ShapeGrad::PostInit()
{
  ErsatzMaterial::PostInit();

  mech_mat_ = dynamic_cast<OptMechMat*>(material); // just created in PostInit()
  assert(material != NULL);
}

} //namespace
