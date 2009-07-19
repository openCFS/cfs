#include "Forms/baseForm.hh"
#include "Forms/linElecInt.hh"
#include "Domain/domain.hh"
#include "PDE/mechPDE.hh"
#include "PDE/elecPDE.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/DesignSpace.hh"
#include "DataInOut/Logging/cfslog.hh"

using namespace CoupledField;

DECLARE_LOG(om)
DEFINE_LOG(om, "optimizationMaterial")


Enum<OptimizationMaterial::System> OptimizationMaterial::system;


OptimizationMaterial::OptimizationMaterial(ErsatzMaterial* em) : regionIds(em->regionIds)
{
  this->opt = em;
}

OptimizationMaterial::~OptimizationMaterial()
{
}


void OptimizationMaterial::GetElementMatrix(BaseForm* form, Matrix<double>& out, Elem* elem, DesignElement::Type direction, double factor)
{
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());

  // when there is no element given, we use the first from our design element.
  // For this process we disable all transfer functions. - they shoul be enabled!
  // This works only for isotropic material.

  // temporarily disables the transfer functions
  opt->GetDesign()->DisableTransferFunctions();

  // We have to do this because CalcElementMatrix asks the density via domain and
  // multiplies it with it's matrix!

  if(elem == NULL) elemList.SetElement(opt->GetDesign()->data[0].elem);
              else elemList.SetElement(elem);

  const EntityIterator& it = elemList.GetIterator();

  // get our element stiffness matrix -> it could be saved over the iteration loops!
  form->CalcElementMatrix(out,const_cast<EntityIterator&>(it), const_cast<EntityIterator&>(it), direction);

  // in piezoelectricity K_pp is -1.0* BDB
  out *= factor;

  LOG_DBG3(om) << "CalcElemMatrix for " << form->GetName() << " factor=" << factor << " -> " << out.ToString();

  // enable again our transfer functions
  opt->GetDesign()->EnableTransferFunctions();
}



OptMechMat::OptMechMat(ErsatzMaterial* em) : OptimizationMaterial(em)
{
  system_ = MECHANIC;
  mech = dynamic_cast<MechPDE*>(opt->ToPDE(Optimization::MECH));
  assert(mech != NULL);

  for(unsigned int r=0; r < regionIds.GetSize(); r++)
  {
    GetElementMatrix(opt->GetForm(regionIds[r], mech, mech, "linElastInt"), mechStiffness_map[regionIds[r]]);
    if(opt->IsHarmonic())
      GetElementMatrix(opt->GetForm(regionIds[r], mech, mech, "MassInt"), mechMass_map[regionIds[r]]);
  }
}


const Matrix<double>& OptMechMat::MechStiffness(Elem* elem, DesignElement::Type direction)
{
  if(!opt->IsDomainStructured())
    GetElementMatrix(opt->GetForm(elem->regionId, mech, mech, "linElastInt"), mechStiffness_map[elem->regionId], elem, direction);

  return mechStiffness_map[elem->regionId];
}

const Matrix<double>& OptMechMat::MechMass(Elem* elem)
{
  if(!opt->IsDomainStructured())
    GetElementMatrix(opt->GetForm(elem->regionId, mech, mech, "MassInt"), mechMass_map[elem->regionId], elem);

  return mechMass_map[elem->regionId];
}


OptPiezoMat::OptPiezoMat(ErsatzMaterial* em) : OptMechMat(em)
{
  system_ = PIEZO;
  elec = dynamic_cast<ElecPDE*>(opt->ToPDE(Optimization::ELEC));
  assert(elec != NULL);

  for(unsigned int r = 0; r < regionIds.GetSize(); r++)
  {
    // in the piezoelectric case K_pp is set to -K_pp. As we use the linElecInt from a piezoelectric problem this is done!
    GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linElecInt"), elecStiffness_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, -1.0); // see above
    GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linElecInt"), elecStiffness_neg_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, 1.0); // see above

    // wrong!!!
    //GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linElecInt"), elecStiffness_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, 1.0); // see above
    //GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linElecInt"), elecStiffness_neg_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, -1.0); // see above


    GetElementMatrix(opt->GetForm(regionIds[r], mech, elec, "linPiezoCoupling"), coupledStiffness_map[regionIds[r]]);
    coupledStiffness_map[regionIds[r]].Transpose(coupledStiffnessTransposed_map[regionIds[r]]);
  }

  // validate that we indeed have the real form
#ifndef NDEBUG
  // create a real linElecInt object

  SubTensorType tensor = domain->GetGrid()->GetDim() == 2 ? PLANE_STRAIN : FULL; // assume AXI not possible
  BaseForm* elec_int =  new linElecInt(elec->getPDEMaterialData()[regionIds[0]], tensor);

  ElemList list(domain->GetGrid());
  list.SetRegion(regionIds[0]);
  EntityIterator elemIt = list.GetIterator();
  elemIt.Begin();

  Matrix<Double> mat;
  GetElementMatrix(elec_int, mat, const_cast<Elem*>(elemIt.GetElem())); // checks the pseudo densitiy disable stuff
  LOG_DBG3(om) << "OptPiezoMat: true K_pp=" << mat.ToString();

  // we just need to check the sign
  Matrix<Double>& our = elecStiffness_map[regionIds[0]];
  LOG_DBG3(om) << "OptPiezoMat: our K_pp=" << our.ToString();
  assert(mat[0][0] == our[0][0]);

  Matrix<Double> our_neg = elecStiffness_neg_map[regionIds[0]];
  assert(mat[0][0] == -1.0 * our_neg[0][0]);
  LOG_DBG3(om) << "OptPiezoMat: our -K_pp=" << our_neg.ToString();

  delete elec_int;

#endif

}

const Matrix<double>& OptPiezoMat::ElecStiffness(Elem* elem, int factor)
{
  assert(factor == 1 || factor == -1);

  std::map<RegionIdType, Matrix<double> >& map = factor == 1 ? elecStiffness_map : elecStiffness_neg_map;
  // overwrite the element
  if(!opt->IsDomainStructured())
    GetElementMatrix(opt->GetForm(elem->regionId, elec, elec, "linElecInt"), map[elem->regionId], NULL, DesignElement::NO_DERIVATIVE, (double) factor);

  return map[elem->regionId];
}

const Matrix<double>& OptPiezoMat::CoupledStiffness(Elem* elem)
{
  if(!opt->IsDomainStructured())
    GetElementMatrix(opt->GetForm(elem->regionId, mech, elec, "linPiezoCoupling"), coupledStiffness_map[elem->regionId]);

  return coupledStiffness_map[elem->regionId];
}


const Matrix<double>& OptPiezoMat::CoupledStiffnessTransposed(Elem* elem)
{
  if(!opt->IsDomainStructured())
  {
    GetElementMatrix(opt->GetForm(elem->regionId, mech, elec, "linPiezoCoupling"), coupledStiffness_map[elem->regionId]);
    coupledStiffness_map[elem->regionId].Transpose(coupledStiffnessTransposed_map[elem->regionId]);
  }

  return coupledStiffnessTransposed_map[elem->regionId];
}

