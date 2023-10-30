#include <assert.h>
#include <stddef.h>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "FeBasis/BaseFE.hh"
#include "General/defs.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Materials/BaseMaterial.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/ShapeGrad.hh"
#include "PDE/SinglePDE.hh"

namespace CoupledField
{

DEFINE_LOG(shapeGrad, "shapeGrad")

ShapeGrad::ShapeGrad() : ErsatzMaterial()
{
  PtrParamNode pncf = domain->GetParamRoot()->Get("optimization")->Get("costFunction");
  if(pncf->Has("multipleExcitation"))
  {
    std::cout << "shapegrad has multiple excitation" << std::endl;
  }
}

void ShapeGrad::GetMaterialParameters(double &lambda, double &mu) const
{
  if(context->pde->GetName() != "mechanic")
    return;

  // TODO: extend for multi-region-optimization if necessary
  assert(false);
  // const BaseMaterial* bm = NULL; // FIXME = pde->getPDEMaterialData()[design->GetRegionIds()[0]];
  // bm->GetScalar(lambda, MECH_LAME_LAMBDA, Global::REAL);
  // bm->GetScalar(mu, MECH_LAME_MU, Global::REAL);
  // LOG_DBG3(shapeGrad) << "lame parameters:  lambda = " << lambda << ", mu = " << mu;
}

void ShapeGrad::GetElementSolution(Vector<double> &vecforward, Vector<double> &vecadjoint, 
                                    const unsigned int e, const SubTensorType type,
                                    App::Type app)
{ 
  assert(false);
  /* FIXME
  Vector<double> intPoint;
  Matrix<double> elem_sol_forward_matrix;
  Matrix<double> elem_sol_adjoint_matrix;

  // we need a NodeStoreSol to get a element*matrix*. We use
  // the current one and explicitly set the raw solution vector from the forward problem
  assert(false);
  NodeStoreSol<double>* node_store_sol = NULL; // FIXME = dynamic_cast<NodeStoreSol<double>* >(ErsatzMaterial::pde->getPDESolution());
  
  // get current design element
  const DesignElement* de = &design->data[e];
  
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  elemList.SetElement(de->elem);
  EntityIterator it = elemList.GetIterator();
  
  it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
  
  const int idx(0); // region index; TODO: extend for multi-region-optimization if necessary

  switch(app)
  {
  case App::MECH:
    {
      // set element solution to matrix
      // use forward as adjoint because it is selfadjoint anyway (and make sign mistake!)
      node_store_sol->GetElemSolutionAsMatrix(elem_sol_forward_matrix, it, 0, forward.Get(idx)->GetVector(Solution::RAW_VECTOR));
      node_store_sol->GetElemSolutionAsMatrix(elem_sol_adjoint_matrix, it, 0, forward.Get(idx)->GetVector(Solution::RAW_VECTOR));

      BaseMaterial* material = pde->getPDEMaterialData()[design->GetRegionIds()[idx]];
  
      // will contain the calculated strains = partial derivatives of displacement u
      shared_ptr<MechStressStrain<double> > strain(new MechStressStrain<double>(material,type));
  
      strain->SetIntPoint(intPoint);
  
      // calculate strains
      strain->SetActElemSol(elem_sol_forward_matrix);
      strain->CalcStrainVec(vecforward, 1, it);
      strain->SetActElemSol(elem_sol_adjoint_matrix);
      strain->CalcStrainVec(vecadjoint, 1, it);
    }
    break;
  case App::HEAT:
    node_store_sol->GetElemSolution(vecforward, it, 0);
    node_store_sol->GetElemSolution(vecadjoint, it, 0);
    break;
  default:
    assert(false);
  }
  
  // debug output
  LOG_DBG3(shapeGrad) << "elem " << de->elem->elemNum << " forward strain: " << vecforward.ToString();
  LOG_DBG3(shapeGrad) << "elem " << de->elem->elemNum << " adjoint strain: " << vecadjoint.ToString();
  */
}

linElastInt* ShapeGrad::getBDBForm()
{
  if(context->pde->GetName() != "mechanic")
    return NULL;
  assert(false);
  return NULL; // FIXME dynamic_cast<linElastInt*>(GetForm(design->GetRegionIds()[0], pde, pde, "LinElastInt"));
}


void ShapeGrad::PostInit()
{
  ErsatzMaterial::PostInit();

  if(context->pde->GetName() != "mechanic")
    return;
}

} //namespace
