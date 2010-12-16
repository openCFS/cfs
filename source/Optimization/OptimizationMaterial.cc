#include "Driver/assemble.hh"
#include "Forms/baseForm.hh"
#include "Forms/linearForm.hh"
#include "Forms/linGradBDBInt.hh"
#include "Domain/domain.hh"
#include "PDE/mechPDE.hh"
#include "PDE/elecPDE.hh"
#include "PDE/heatCondPDE.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/cfslog.hh"

using namespace CoupledField;

DECLARE_LOG(om)
DEFINE_LOG(om, "optimizationMaterial")


Enum<OptimizationMaterial::System> OptimizationMaterial::system;


OptimizationMaterial::OptimizationMaterial(ErsatzMaterial* em) :
  regionIds(em->regionIds),
  opt(em)
{
}

OptimizationMaterial::~OptimizationMaterial()
{
}


void OptimizationMaterial::GetElementMatrix(BaseForm* form, Matrix<double>& out, const Elem* elem, DesignElement::Type direction, double factor)
{
  GetElementEntity(form, &out, NULL, elem, direction);
  // in piezoelectricity K_pp is -1.0* BDB
  if(factor != 1.0){
    out *= factor;
  }

  LOG_DBG3(om) << "CalcElemMatrix for " << form->GetName() << " factor=" << factor << " -> " << out.ToString();
}

void OptimizationMaterial::GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem, const Vector<double>* ts)
{
  GetElementEntity(form, NULL, &out, elem,DesignElement::NO_DERIVATIVE, ts);

  LOG_DBG3(om) << "CalcElemVector for " << form->GetName() << " -> " << out.ToString();
}


void OptimizationMaterial::GetElementEntity(BaseForm* form, Matrix<double>* mat_out, Vector<double>* vec_out, const Elem* elem, DesignElement::Type direction, const Vector<double>* ts)
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

  assert(mat_out == NULL || vec_out == NULL);
  if(mat_out != NULL)
  {
    // get our element stiffness matrix -> it could be saved over the iteration loops!
    form->CalcElementMatrix(*mat_out, const_cast<EntityIterator&>(it), const_cast<EntityIterator&>(it), direction);
  }
  else
  {
    if(form->GetName() == "AddStrainRHSInt")
    {
      AddStrainRHSInt* lf = dynamic_cast<AddStrainRHSInt*>(form);
      lf->CalcElemVector(*vec_out, const_cast<EntityIterator&>(it), ts);
    }
    else
    {
      LinearForm* lf = dynamic_cast<LinearForm*>(form);
      lf->CalcElemVector(*vec_out, const_cast<EntityIterator&>(it));
    }
  }

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
    LOG_DBG(om) << "OptMechMat MechStiffness region=" << domain->GetGrid()->GetRegion().ToString(regionIds[r])
                << std::endl << mechStiffness_map[regionIds[r]].ToString(0,true);
    if(opt->IsHarmonic())
    {
      GetElementMatrix(opt->GetForm(regionIds[r], mech, mech, "MassInt"), mechMass_map[regionIds[r]]);
      LOG_DBG(om) << "OptMechMat MechMaxx region=" << domain->GetGrid()->GetRegion().ToString(regionIds[r])
                  << std::endl << mechMass_map[regionIds[r]].ToString(0,true);
    }
  }
}


const Matrix<double>& OptMechMat::MechStiffness(const Elem* elem, DesignElement::Type direction)
{
  if(!opt->IsDomainStructured() || direction != DesignElement::NO_DERIVATIVE)
    GetElementMatrix(opt->GetForm(elem->regionId, mech, mech, "linElastInt"), mechStiffness_map[elem->regionId], elem, direction);

  return mechStiffness_map[elem->regionId];
}


const Vector<double>& OptMechMat::MechStrainRHS(const Elem* elem, MechPDE::TestStrain testStrain)
{
  // in homogenization we always set/replace the actual AddStrainRHSInt which contains the current test strain,
  // therefore we do not cache!
  SinglePDE* pde = opt->ToPDE(Optimization::MECH);
  LinearForm* lf = pde->getPDE_assemble()->GetLinearForm(opt->GetDesign()->GetRegionId(), pde, "AddStrainRHSInt")->GetIntegrator();
  // this is really inefficient -> but won't cost to much! the vector is created too often!
  if(testStrain != MechPDE::NOT_SET)
  {
    Vector<double> ts = mech->CalcTestStrainVector(testStrain, true);
    GetElementVector(lf, mechStrainRHS, elem, &ts);
  }
  else
  {
    GetElementVector(lf, mechStrainRHS, elem);
  }

  return mechStrainRHS;
}

const Matrix<double>& OptMechMat::MechMass(const Elem* elem, DesignElement::Type direction)
{
  if(!opt->IsDomainStructured() || direction != DesignElement::NO_DERIVATIVE)
    GetElementMatrix(opt->GetForm(elem->regionId, mech, mech, "MassInt"), mechMass_map[elem->regionId], elem, direction);

  return mechMass_map[elem->regionId];
}


OptPiezoMat::OptPiezoMat(ErsatzMaterial* em) :
  OptMechMat(em)
{
  system_ = PIEZO;
  elec = dynamic_cast<ElecPDE*>(opt->ToPDE(Optimization::ELEC));
  assert(elec != NULL);

  for(unsigned int r = 0; r < regionIds.GetSize(); r++)
  {
    // in the piezoelectric case K_pp is set to -K_pp. As we use the linGradBDBInt from a piezoelectric problem this is done!
    GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linGradBDBInt"), elecStiffness_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, -1.0); // see above
    GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linGradBDBInt"), elecStiffness_neg_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, 1.0); // see above

    // wrong!!!
    //GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linGradBDBInt"), elecStiffness_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, 1.0); // see above
    //GetElementMatrix(opt->GetForm(regionIds[r], elec, elec, "linGradBDBInt"), elecStiffness_neg_map[regionIds[r]], NULL, DesignElement::NO_DERIVATIVE, -1.0); // see above


    GetElementMatrix(opt->GetForm(regionIds[r], mech, elec, "linPiezoCoupling"), coupledStiffness_map[regionIds[r]]);
    coupledStiffness_map[regionIds[r]].Transpose(coupledStiffnessTransposed_map[regionIds[r]]);
  }

  // validate that we indeed have the real form
#ifndef NDEBUG
  // create a real linGradBDBInt object

  SubTensorType tensor = domain->GetGrid()->GetDim() == 2 ? PLANE_STRAIN : FULL; // assume AXI not possible
  BaseForm* elec_int =  new linGradBDBInt(elec->getPDEMaterialData()[regionIds[0]], ELEC_PERMITTIVITY, tensor);

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
    GetElementMatrix(opt->GetForm(elem->regionId, elec, elec, "linGradBDBInt"), map[elem->regionId], NULL, DesignElement::NO_DERIVATIVE, (double) factor);

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

HeatMat::HeatMat(ErsatzMaterial* em) :
  OptimizationMaterial(em)
{
  system_ = HEAT;
  heat = dynamic_cast<HeatCondPDE*>(opt->ToPDE(Optimization::HEAT));
  assert(heat != NULL);
}
