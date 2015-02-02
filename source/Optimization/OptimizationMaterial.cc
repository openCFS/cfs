#include <assert.h>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/FormsContexts.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "General/defs.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/AcousticPDE.hh"
#include "PDE/BasePDE.hh"
#include "PDE/ElecPDE.hh"
#include "PDE/HeatPDE.hh"
#include "PDE/MechPDE.hh"
// #include "PDE/LatticeBoltzmannPDE.hh"

namespace CoupledField {
class BaseMaterial;
}  // namespace CoupledField

using namespace CoupledField;

DECLARE_LOG(om)
DEFINE_LOG(om, "optimizationMaterial")


Enum<OptimizationMaterial::System> OptimizationMaterial::system;


OptimizationMaterial::OptimizationMaterial(ErsatzMaterial* em, DesignSpace* space)
{
  assert((em == NULL && space != NULL) || (em != NULL && (space == NULL || space == em->GetDesign())));

  this->opt = em;
  this->space = space != NULL ? space : em->GetDesign();

  regionIds = this->space->GetRegionIds();

  harmonic_ =  domain->GetDriver()->IsComplex();
  transient_ = em->IsTransient();
  structured_ = em != NULL ? em->IsDomainStructured() : false; // we would have a problem with the L-Shape!
}

OptimizationMaterial::~OptimizationMaterial()
{
}

OptimizationMaterial* OptimizationMaterial::CreateInstance(System sys, ErsatzMaterial* em)
{
  switch(sys)
  {
  case PIEZOCOUPLING:
    return new PiezoElecMat(em);

  case MECH:
    return new MechMat(em);

  case ELEC:
    return new ElecMat(em);

  case HEAT:
    return new HeatMat(em);

  case ACOUSTIC:
    return new AcouMat(em);

  case LBM:
    return new LBMMat(em);

  default:
    assert(false);
    return NULL;
  }
}


shared_ptr<CoefFunctionOpt> OptimizationMaterial::GetMatCoef(const std::string& integrator, BiLinFormContext* context, RegionIdType reg_id)
{
  if(context != NULL)
  {
    BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->GetIntegrator());
    if(bdb != NULL)
    {
      PtrCoefFct coef = bdb->GetCoef();
      return dynamic_pointer_cast<CoefFunctionOpt>(bdb->GetCoef());
    }
  }
  EXCEPTION("not implemented");
}

void OptimizationMaterial::GetElementMatrix(BiLinearForm* form, Matrix<double>& out, const Elem* elem, BaseMaterial* bimaterial, DesignElement::Type direction, double factor)
{
  GetElementEntity(form, &out, NULL, elem, bimaterial, direction);

  // FIXME LOG_DBG3(om) << "CalcElemMatrix for " << form->GetName() << " factor=" << factor << " -> " << out.ToString();
}

void OptimizationMaterial::GetElementMatrix(Matrix<double>& out, const std::string& integrator, const Elem* elem, BaseMaterial* bimaterial, DesignElement::Type direction,  Global::ComplexPart entryType)
{
  assert(bimaterial == NULL);
  assert(direction == DesignElement::NO_DERIVATIVE || direction == DesignElement::NO_MULTIMATERIAL);
  assert(entryType != Global::IMAG);

  if(elem == NULL)
    elem = space->data[0].elem;

  BiLinFormContext* c = GetPDE()->GetAssemble()->GetBiLinForm(integrator, elem->regionId, GetPDE(), GetPDE(), false);

  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  elemList.SetElement(elem);
  EntityIterator it = elemList.GetIterator();

  shared_ptr<CoefFunctionOpt> coef = GetMatCoef(integrator, c, elem->regionId);
  coef->SetEnable(false);

  c->GetIntegrator()->CalcElementMatrix(out, it, it);

  coef->SetEnable(true);
}


BiLinearForm* OptimizationMaterial::GetForm(RegionIdType regionId, const std::string& integrator, Global::ComplexPart entryType)
{
   Assemble* ass = GetPDE()->GetAssemble();
   BiLinFormContext* c = ass->GetBiLinForm(integrator, regionId, GetPDE(), GetPDE(), false);
   // TODO
   // ignore the condition -> in legacy cfs this was part of the search routine in assemble!
   assert(c->GetEntryType() == entryType);

   return c->GetIntegrator();
}

void OptimizationMaterial::GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem, BaseMaterial* bimaterial, const Vector<double>* ts)
{
  assert(false);
  // FIXME GetElementEntity(form, NULL, &out, elem, bimaterial, DesignElement::NO_DERIVATIVE, ts);

  LOG_DBG3(om) << "CalcElemVector for " << form->GetName() << " -> " << out.ToString();
}


void OptimizationMaterial::GetElementEntity(BiLinearForm* form, Matrix<double>* mat_out, Vector<double>* vec_out, const Elem* elem, BaseMaterial* bimaterial,
                                            DesignElement::Type direction, const Vector<double>* ts)
{
  assert(false);
  /*

  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());

  // when there is no element given, we use the first from our design element.
  // For this process we disable all transfer functions. - they shoul be enabled!
  // This works only for isotropic material.

  // temporarily disables the transfer functions
  space->DisableTransferFunctions();

  // We have to do this because CalcElementMatrix asks the density via domain and
  // multiplies it with it's matrix!

  if(elem == NULL) elemList.SetElement(space->data[0].elem);
              else elemList.SetElement(elem);

  // form needs to be the right one for the region!
  assert(false);
  BaseMaterial* org_mat = NULL; // FIXME form->GetMaterial(); // for bimaterial

  if(bimaterial != NULL)
  {
    assert(false);
    // FIXME form->SetMaterial(bimaterial);
  }

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
  space->EnableTransferFunctions();

  if(bimaterial != NULL)
    form->SetMaterial(org_mat);

  */
}

Matrix<double>& OptimizationMaterial::GeneralStiffness(std::map<RegionIdType, StdVector<Matrix<double> > >& map, const DesignElement* de, MaterialClass mc,
                                                       StdPDE* pde1, StdPDE* pde2, DesignElement::Type direction, double factor, bool transposed)
{
  unsigned int index = de->multimaterial != NULL ? de->multimaterial->index : 0;

  StdVector<Matrix<double> >& mat_vec = map[de->elem->regionId];

  // generate entries on the fly
  for(unsigned int i = mat_vec.GetSize(); i <= index; i++)
    mat_vec.Push_back(Matrix<double>());

  Matrix<double>& mat = mat_vec[index];

  if(space->HasMultiMaterial())
    direction = DesignElement::NO_MULTIMATERIAL;

  std::string form;
  switch(mc)
  {
  case MECHANIC:
    form = "LinElastInt";
    break;
  case ELECTROSTATIC:
    form = "linGradBDBInt";
    break;
  case PIEZO:
    form = "linPiezoCoupling";
    break;
  default:
    assert(false);
    break;
  }

  if(mat.GetNumRows() == 0 || !structured_ || direction != DesignElement::NO_DERIVATIVE)
  {
    BaseMaterial* bm = NULL;
    if(de != NULL && de->multimaterial != NULL) // TODO is actually nonsense but due to non-region dependency of DesignSpace::multimaterial
      bm = space->GetRegion(de->elem->regionId, &space->GetMultiMaterials()[index])->multimaterial->GetMultiMaterial(mc);

    GetElementMatrix(GetForm(de->elem->regionId, form), mat, NULL, bm, direction, factor);
    if(transposed)
    {
      Matrix<double> tmp = mat;
      tmp.Transpose(mat);
    }
  }

  return mat;
}




MechMat::MechMat(ErsatzMaterial* em) : OptimizationMaterial(em)
{
  Init();
}

MechMat::MechMat(DesignSpace* space) : OptimizationMaterial(NULL, space)
{
  Init();
}

SinglePDE* MechMat::GetPDE() {
  return dynamic_cast<SinglePDE*>(mech);
}

void MechMat::Init()
{
  system_ = MECH;
  mech = dynamic_cast<MechPDE*>(opt != NULL ? opt->ToPDE(Optimization::MECH) : domain->GetSinglePDE("mechanic"));
  assert(mech != NULL);
  StdVector<MultiMaterial>& mm = space->GetMultiMaterials();

  for(unsigned int r=0; r < regionIds.GetSize(); r++)
  {
    RegionIdType reg_id = regionIds[r];
    StdVector<Matrix<double> >& K = mechStiffness_map[reg_id];
    StdVector<Matrix<Complex> >& KC = mechStiffness_mapC[reg_id];

    // three cases: multimaterial, (legacy) bimaterial, normal. Run once for the last two
    for(unsigned int m = 0; m < (space->HasMultiMaterial() ? mm.GetSize() : 1); m++)
    {
      DesignSpace::DesignRegion* dr = space->HasMultiMaterial() ? space->GetRegion(reg_id, DesignElement::MULTIMATERIAL, mm[m].index) : space->GetRegion(reg_id);

      // either system material or multimaterial
      BaseMaterial* bm = space->HasMultiMaterial() ? mm[m].GetMultiMaterial(MECHANIC) : NULL;
      // prepare target

      //GetElementMatrix(ErsatzMaterial::GetForm(reg_id, mech, mech, "LinElastInt"), K.Last(), NULL, bm, DesignElement::NO_MULTIMATERIAL);
      Matrix<double> tmp_mat;
      Matrix<Complex> tmp_mat2;

      GetElementMatrix(tmp_mat, "LinElastInt", NULL, bm, DesignElement::NO_MULTIMATERIAL);
      if (mech->HasComplexMatData(reg_id))
      {
        tmp_mat2.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
        tmp_mat2.SetPart(Global::REAL, tmp_mat);
        GetElementMatrix(GetForm(reg_id, "LinElastInt", Global::IMAG), tmp_mat, NULL, bm, DesignElement::NO_MULTIMATERIAL);
        tmp_mat2.SetPart(Global::IMAG, tmp_mat);
        KC.Push_back(tmp_mat2);
        LOG_DBG(om) << "OptMechMat: MechStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << std::endl << KC.Last().ToString(0,true);
      }
      else
      {
        K.Push_back(tmp_mat);
        LOG_DBG(om) << "OptMechMat: MechStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << std::endl << K.Last().ToString(0,true);
      }

      if(dr->HasBiMaterial())
      {
//        K.Push_back(Matrix<double>());
//
//        GetElementMatrix(ErsatzMaterial::GetForm(reg_id, mech, mech, "LinElastInt"), K.Last(), NULL, dr->GetBiMaterial(MECHANIC));
        GetElementMatrix(GetForm(reg_id, "LinElastInt", Global::REAL), tmp_mat, NULL, dr->GetBiMaterial(MECHANIC));
        if (mech->HasComplexMatData(reg_id)){
          tmp_mat2.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
          tmp_mat2.SetPart(Global::REAL, tmp_mat);
          GetElementMatrix(GetForm(reg_id, "LinElastInt", Global::IMAG), tmp_mat, NULL, dr->GetBiMaterial(MECHANIC));
          tmp_mat2.SetPart(Global::IMAG, tmp_mat);
          KC.Push_back(tmp_mat2);
        }
        else K.Push_back(tmp_mat);
        LOG_DBG(om) << "OptMechMat: MechStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << " bimaterial" << std::endl << K.Last().ToString(0,true);
      }

      if(harmonic_ || transient_)
      {
        mechMass_map[reg_id].Push_back(Matrix<double>());

        GetElementMatrix(mechMass_map[reg_id].Last(), "MassInt");
        LOG_DBG(om) << "OptMechMat MechMass region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << std::endl << mechMass_map[reg_id].Last().ToString(0,true);

        if(dr->HasBiMaterial())
        {
          mechMass_map[reg_id].Push_back(Matrix<double>());

          GetElementMatrix(GetForm(reg_id, "MassInt"), mechMass_map[reg_id].Last(), NULL, dr->GetBiMaterial(MECHANIC));
          LOG_DBG(om) << "OptMechMat MechMass region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << " bimaterial" << std::endl << mechMass_map[reg_id].Last().ToString(0,true);
        }
      }
    }
  }
}




DenseMatrix& MechMat::MechStiffness(const Elem* elem, bool bimaterial, int multimaterial, DesignElement::Type direction)
{
  unsigned int index = multimaterial < 0 ? 0 : (unsigned int) multimaterial;
  bool complexMaterial = mech->HasComplexMatData(elem->regionId);

  // in the multimaterial case we do not want the form to obtain the multimaterial tensor but only the one of the currently material
  if(space->HasMultiMaterial())
    direction = DesignElement::NO_MULTIMATERIAL;

  LOG_DBG3(om) << "MS: el=" << elem->elemNum << " bi=" << bimaterial << " mm=" << multimaterial << " d=" << direction << " index=" << index;

  // we need to get the data on request
  if(!structured_ || direction != DesignElement::NO_DERIVATIVE)
  {
    BaseMaterial* bm = NULL;
    if(multimaterial >= 0) // TODO is actually nonsense but due to non-region dependency of DesignSpace::multimaterial
      bm = space->GetRegion(elem->regionId, &space->GetMultiMaterials()[multimaterial])->multimaterial->GetMultiMaterial(MECHANIC);

    Matrix<double> tmp_mat;
    Matrix<Complex> tmp_mat2;
    if(!bimaterial)
    {
      GetElementMatrix(tmp_mat, "LinElastInt", elem, bm, direction);

      if (complexMaterial)
      {
        tmp_mat2.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
        tmp_mat2.SetPart(Global::REAL, tmp_mat);
        GetElementMatrix(GetForm(elem->regionId, "LinElastInt", Global::IMAG), tmp_mat, elem, bm, direction);
        tmp_mat2.SetPart(Global::IMAG, tmp_mat);
        mechStiffness_mapC[elem->regionId][index] = tmp_mat2;
      }
      else mechStiffness_map[elem->regionId][index] = tmp_mat;
    }
    else
    {
      BaseMaterial* bm = space->GetRegion(elem->regionId)->GetBiMaterial(MECHANIC);
      //GetElementMatrix(ErsatzMaterial::GetForm(elem->regionId, mech, mech, "LinElastInt"), mechStiffness_map[elem->regionId][1], elem, bm, direction);
      GetElementMatrix(GetForm(elem->regionId, "LinElastInt", Global::REAL), tmp_mat, elem, bm, direction);
      if (complexMaterial){
        tmp_mat2.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
        tmp_mat2.SetPart(Global::REAL, tmp_mat);
        GetElementMatrix(GetForm(elem->regionId, "LinElastInt", Global::IMAG), tmp_mat, elem, bm, direction);
        tmp_mat2.SetPart(Global::IMAG, tmp_mat);
        mechStiffness_mapC[elem->regionId][index] = tmp_mat2;
      }
      else mechStiffness_map[elem->regionId][1] = tmp_mat;
    }
  }

  if(bimaterial)
  {
    if (complexMaterial)
      return dynamic_cast<DenseMatrix& >(mechStiffness_mapC[elem->regionId][1]);
    else
      return dynamic_cast<DenseMatrix& >(mechStiffness_map[elem->regionId][1]);
  }

  if (complexMaterial)
  {
    LOG_DBG3(om) << "MS: el=" << elem->elemNum << " -> " << mechStiffness_mapC[elem->regionId][index].ToString();
    return dynamic_cast<DenseMatrix& >(mechStiffness_mapC[elem->regionId][index]); // no multimaterial is frst material
  }
  else
    LOG_DBG3(om) << "MS: el=" << elem->elemNum << " -> " << mechStiffness_map[elem->regionId][index].ToString();
    return dynamic_cast<DenseMatrix& >(mechStiffness_map[elem->regionId][index]);
}

const Matrix<double>& MechMat::MechMass(const Elem* elem, bool bimaterial, int multimaterial, DesignElement::Type direction)
{
   // FIXME assert(ErsatzMaterial::GetForm(elem->regionId, mech, mech, "MassInt")->GetMaterialDescriptor().Enabled());
  unsigned int index = multimaterial < 0 ? 0 : (unsigned int) multimaterial;

  if(!structured_ || direction != DesignElement::NO_DERIVATIVE)
  {
    BaseMaterial* bm = NULL;
    if(multimaterial >= 0) // TODO is actually nonsense but due to non-region dependency of DesignSpace::multimaterial
      bm = space->GetRegion(elem->regionId, &space->GetMultiMaterials()[multimaterial])->multimaterial->GetMultiMaterial(MECHANIC);

    if(!bimaterial)
      GetElementMatrix(GetForm(elem->regionId, "MassInt"), mechMass_map[elem->regionId][index], elem, bm, direction);
    else
    {
      BaseMaterial* bm = space->GetRegion(elem->regionId)->GetBiMaterial(MECHANIC);
      GetElementMatrix(GetForm(elem->regionId, "MassInt"), mechMass_map[elem->regionId][1], elem, bm, direction);
    }
  }

  if(bimaterial)
    return mechMass_map[elem->regionId][1];

  return mechMass_map[elem->regionId][index];
}

const Vector<double>& MechMat::MechStrainRHS(const Elem* elem, MechPDE::TestStrain testStrain)
{
  // in homogenization we always set/replace the actual AddStrainRHSInt which contains the current test strain,
  // therefore we do not cache!
  assert(false);
  LinearForm* lf = NULL; // FIXME mech->GetAssemble()->GetLinearForm(space->GetRegionId(), mech, "AddStrainRHSInt")->GetIntegrator();
  // this is really inefficient -> but won't cost to much! the vector is created too often!
  if(testStrain != MechPDE::NOT_SET)
  {
    assert(false);
    Vector<double> ts; // FIXME = mech->CalcTestStrainVector(testStrain, true);
    GetElementVector(lf, mechStrainRHS, elem, NULL, &ts);
  }
  else
  {
    GetElementVector(lf, mechStrainRHS, elem);
  }

  return mechStrainRHS;
}


AcouMat::AcouMat(ErsatzMaterial* em) : OptimizationMaterial(em)
{
  system_ = ACOUSTIC;
  acou = dynamic_cast<AcousticPDE*>(opt != NULL ? opt->ToPDE(Optimization::ACOUSTIC) : domain->GetSinglePDE("acoustic"));
  assert(acou != NULL);

  for(unsigned int r=0; r < regionIds.GetSize(); r++)
  {
    RegionIdType reg_id = regionIds[r];
    DesignSpace::DesignRegion* dr = space->GetRegion(reg_id);
    BaseMaterial* bm = dr->HasBiMaterial() ? dr->GetBiMaterial(FLUID) : NULL;

    GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].first, NULL, NULL);

    if(dr->HasBiMaterial())
      GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].second, NULL, bm);

    if(harmonic_)
    {
      GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].first, NULL, NULL);

      if(dr->HasBiMaterial())
        GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].second, NULL, bm);
    }
  }
}


SinglePDE* AcouMat::GetPDE() {
  return dynamic_cast<SinglePDE*>(acou);
}


Matrix<double>& AcouMat::AcouStiffness(const Elem* elem, bool bimaterial)
{
  RegionIdType reg_id = elem->regionId;

  if(!structured_)
  {
    if(!bimaterial)
      GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].first, elem, NULL);
    else
    {
      BaseMaterial* bm = space->GetRegion(reg_id)->GetBiMaterial(FLUID);
      GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].second, elem, bm);
    }
  }

  return !bimaterial ? acouStiffness_map[reg_id].first : acouStiffness_map[reg_id].second;
}

const Matrix<double>& AcouMat::AcouMass(const Elem* elem, bool bimaterial)
{
  RegionIdType reg_id = elem->regionId;

  // FIXME assert(ErsatzMaterial::GetForm(reg_id, acou, acou, "MassInt")->GetMaterialDescriptor().Enabled());

  if(!structured_)
  {
    if(!bimaterial)
      GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].first, elem, NULL);
    else
    {
      BaseMaterial* bm = space->GetRegion(reg_id)->GetBiMaterial(FLUID);
      GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].second, elem, bm);
    }
  }

  return !bimaterial ? acouMass_map[reg_id].first : acouMass_map[reg_id].second;
}



PiezoElecMat::PiezoElecMat(ErsatzMaterial* em) :
  MechMat(em)
{
  system_ = PIEZOCOUPLING;
  elec = dynamic_cast<ElecPDE*>(opt != NULL ? opt->ToPDE(Optimization::ELEC) : domain->GetSinglePDE("electrostatic"));
  assert(elec != NULL);

  for(unsigned int r = 0; r < regionIds.GetSize(); r++)
  {
    Elem elem;
    DesignElement de;
    de.elem = &elem;
    de.elem->regionId = regionIds[r];

    // three cases: multimaterial, (legacy) bimagterial, normal. Run once for the last two
    for(int m = -1; m < (int) space->GetMultiMaterials().GetSize(); m++)
    {
      de.multimaterial = (m == -1) ? NULL : &(space->GetMultiMaterials()[m]);

      ElecStiffnessPos(&de, DesignElement::NO_DERIVATIVE); // has assert
      ElecStiffnessNeg(&de, DesignElement::NO_DERIVATIVE);
      CoupledStiffness(&de, DesignElement::NO_DERIVATIVE);
      CoupledStiffnessTransposed(&de, DesignElement::NO_DERIVATIVE);
    }
  }
}

SinglePDE* PiezoElecMat::GetPDE() {
  return dynamic_cast<SinglePDE*>(elec);
}



const Matrix<double>& PiezoElecMat::ElecStiffnessPos(const DesignElement* de, DesignElement::Type direction)
{
  Matrix<double>& mat = GeneralStiffness(elecStiffness_map, de, ELECTROSTATIC, elec, elec, direction, -1.0 , false); // piezo form already negates
  mat *= -1.0;
  assert(direction != DesignElement::NO_DERIVATIVE || mat.Trace() >= 0.0); // allow zero in case of piezo-fmo for non-DIELEC direction
  return mat;
}

const Matrix<double>& PiezoElecMat::ElecStiffnessNeg(const DesignElement* de, DesignElement::Type direction)
{
  Matrix<double>& mat = GeneralStiffness(elecStiffness_neg_map, de, ELECTROSTATIC, elec, elec, direction, 1.0, false); // piezo form already negates
  //   type.Add(NO_MULTIMATERIAL, "no_multimaterial");LOG_DBG3(om) << "PEM:ESN de=" << de->ToString() << " dt=" << DesignElement::type.ToString(direction) << " -> " << mat.ToString();
  assert(direction != DesignElement::NO_DERIVATIVE || mat.Trace() <= 0.0); // allow zero in case of piezo-fmo for non-DIELEC direction
  return mat;
}


const Matrix<double>& PiezoElecMat::CoupledStiffness(const DesignElement* de, DesignElement::Type direction)
{
  return GeneralStiffness(coupledStiffness_map, de, PIEZO, mech, elec, direction, 1.0, false);
}


const Matrix<double>& PiezoElecMat::CoupledStiffnessTransposed(const DesignElement* de, DesignElement::Type direction)
{
  return GeneralStiffness(coupledStiffnessTransposed_map, de, PIEZO, mech, elec, direction, 1.0, true);
}

HeatMat::HeatMat(ErsatzMaterial* em) :
  OptimizationMaterial(em)
{
  system_ = HEAT;
  assert(false);
  heat = NULL; // FIXME dynamic_cast<HeatCondPDE*>(opt != NULL ? opt->ToPDE(Optimization::HEAT) : domain->GetSinglePDE("heatConduction"));
  assert(heat != NULL);
}

SinglePDE* HeatMat::GetPDE() {
  assert(false);
  return NULL;
  // FIXME return dynamic_cast<SinglePDE*>(heat);
}



ElecMat::ElecMat(ErsatzMaterial* em) : OptimizationMaterial(em)
{
  system_ = ELEC;
  elec = dynamic_cast<ElecPDE*>(opt->ToPDE(Optimization::ELEC));
  assert(elec != NULL);

  for(unsigned int r=0; r < regionIds.GetSize(); r++)
  {
    RegionIdType reg_id = regionIds[r];
    DesignSpace::DesignRegion* dr = opt->GetDesign()->GetRegion(reg_id);

    Matrix<double> tmp_mat;
    GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat);
    elecStiffness_map[reg_id].first.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
    elecStiffness_map[reg_id].first.SetPart(Global::REAL, tmp_mat);
    if (elec->HasComplexMatData(reg_id)){
      GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat);
      elecStiffness_map[reg_id].first.SetPart(Global::IMAG, tmp_mat);
    }
    LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id)
                        << std::endl << elecStiffness_map[reg_id].first.ToString(0,true);

    if(dr->HasBiMaterial())
    {
      elecStiffness_map[reg_id].second.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
      elecStiffness_map[reg_id].second.SetPart(Global::REAL, tmp_mat);
      if (elec->HasComplexMatData(reg_id)){
        GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
        elecStiffness_map[reg_id].second.SetPart(Global::IMAG, tmp_mat);
      }
      LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << " bimaterial"
          << std::endl << elecStiffness_map[reg_id].second.ToString(0,true);
    }

  }
}


SinglePDE* ElecMat::GetPDE() {
  return dynamic_cast<SinglePDE*>(elec);
}

void ElecMat::ReInit(){
  assert(elec != NULL);

  for(unsigned int r=0; r < regionIds.GetSize(); r++)
  {
    RegionIdType reg_id = regionIds[r];
    DesignSpace::DesignRegion* dr = opt->GetDesign()->GetRegion(reg_id);

    Matrix<double> tmp_mat;
    GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat);
    elecStiffness_map[reg_id].first.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
    elecStiffness_map[reg_id].first.SetPart(Global::REAL, tmp_mat);
    if (elec->HasComplexMatData(reg_id)){
      GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat);
      elecStiffness_map[reg_id].first.SetPart(Global::IMAG, tmp_mat);
    }
    LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id)
                        << std::endl << elecStiffness_map[reg_id].first.ToString(0,true);

    if(dr->HasBiMaterial())
    {
      elecStiffness_map[reg_id].second.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
      elecStiffness_map[reg_id].second.SetPart(Global::REAL, tmp_mat);
      if (elec->HasComplexMatData(reg_id)){
        GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
        elecStiffness_map[reg_id].second.SetPart(Global::IMAG, tmp_mat);
      }
      LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << " bimaterial"
          << std::endl << elecStiffness_map[reg_id].second.ToString(0,true);
    }

  }
}


Matrix<std::complex<double> >& ElecMat::ElecStiffness(const Elem* elem, bool bimaterial, DesignElement::Type direction)
{
  if(!opt->IsDomainStructured() || direction != DesignElement::NO_DERIVATIVE)
  {
    Matrix<double> tmp_mat;
    if(!bimaterial){
      GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::REAL), tmp_mat, elem, NULL, direction);
      elecStiffness_map[elem->regionId].first.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      elecStiffness_map[elem->regionId].first.SetPart(Global::REAL, tmp_mat);
      if (opt->IsHarmonic()){
        GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::IMAG), tmp_mat, elem, NULL, direction);
        elecStiffness_map[elem->regionId].first.SetPart(Global::IMAG, tmp_mat);
        //  GetElementMatrix(opt->GetForm(elem->regionId, elec, elec, "linGradBDBInt"), elecStiffness_map[elem->regionId].first, elem, NULL, direction);
      }
    }
    else
    {
      BaseMaterial* bm = opt->GetDesign()->GetRegion(elem->regionId)->GetBiMaterial(ELECTROSTATIC);
      GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::REAL), tmp_mat, elem, bm, direction);
      elecStiffness_map[elem->regionId].second.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      elecStiffness_map[elem->regionId].second.SetPart(Global::REAL, tmp_mat);
      if (opt->IsHarmonic()){
        GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::IMAG), tmp_mat, elem, bm, direction);
        elecStiffness_map[elem->regionId].second.SetPart(Global::IMAG, tmp_mat);
        //  GetElementMatrix(opt->GetForm(elem->regionId, elec, elec, "linGradBDBInt"), elecStiffness_map[elem->regionId].second, elem, bm, direction);
      }
    }
  }

  return !bimaterial ? elecStiffness_map[elem->regionId].first : elecStiffness_map[elem->regionId].second;
}

LBMMat::LBMMat(ErsatzMaterial* em) :
  OptimizationMaterial(em)
{
  system_ = LBM;
  assert(false);
  lbm = NULL; // FIXME dynamic_cast<LatticeBoltzmannPDE*>(opt != NULL ? opt->ToPDE(Optimization::LBM) : domain->GetSinglePDE("LatticeBoltzmann"));
  assert(lbm != NULL);
}

