#include <assert.h>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Driver/Assemble.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Driver/FormsContexts.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "General/defs.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/LocalElementCache.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "PDE/SinglePDE.hh"

namespace CoupledField {
class BaseMaterial;
}  // namespace CoupledField

using std::complex;

using namespace CoupledField;

DECLARE_LOG(om)
DEFINE_LOG(om, "optimizationMaterial")


Enum<OptimizationMaterial::System> OptimizationMaterial::system;


OptimizationMaterial::OptimizationMaterial(ErsatzMaterial* em, Context* ctxt, DesignSpace* space)
{
  assert((em == NULL && space != NULL) || (em != NULL && (space == NULL || space == em->GetDesign())));

  this->ctxt_ = ctxt;
  this->opt = em;
  this->space = space != NULL ? space : em->GetDesign();
  this->system_ = NO_SYSTEM;

  regionIds = this->space->GetRegionIds();

  assert(em->context == ctxt); // shall hold for initialization but is not really important

  needs_mass_ = ctxt->IsComplex() || em->IsTransient() || ctxt->IsEigenvalue();

  transient_ = em->IsTransient();
}

OptimizationMaterial::~OptimizationMaterial()
{
}



OptimizationMaterial* OptimizationMaterial::CreateInstance(System sys, ErsatzMaterial* em, Context* ctxt)
{
  switch(sys)
  {
  case PIEZOCOUPLING:
    return new PiezoElecMat(em, ctxt);

  case MECH:
    return new MechMat(em, ctxt);

  case ELEC:
    return new ElecMat(em, ctxt);

  case HEAT:
    return new HeatMat(em, ctxt);

  case ACOUSTIC:
    return new AcouMat(em, ctxt);

  case LBM:
    return new LBMMat(em, ctxt);

  default:
    assert(false);
    return NULL;
  }
}

shared_ptr<CoefFunctionOpt> OptimizationMaterial::GetMatCoef(const std::string& integrator, BiLinFormContext* context, RegionIdType reg_id)
{
  assert(context != NULL);
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->GetIntegrator());
  assert(bdb != NULL);
  assert(bdb->GetCoef());
  // LOG_DBG3(om) << "GMC int=" << integrator << " coef=" << bdb->GetCoef()->ToString();
  return dynamic_pointer_cast<CoefFunctionOpt>(bdb->GetCoef());
}


template <class T>
void OptimizationMaterial::GetElementMatrix(Matrix<T>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType)
{
  LOG_DBG3(om) << "GEM int=" << integrator << " elem=" << (elem != NULL ? elem->elemNum : 4711) << " lb=" << lower_bimat << " d=" << direction << " et=" << entryType;

  assert(entryType != Global::IMAG);
  assert(elem != NULL);

  SinglePDE* pde = ctxt_->pde;
  assert(pde != NULL);

  BiLinFormContext* c = pde->GetAssemble()->GetBiLinForm(integrator, elem->regionId, pde, pde, false);

  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  elemList.SetElement(elem);
  EntityIterator it = elemList.GetIterator();

  shared_ptr<CoefFunctionOpt> coef = GetMatCoef(integrator, c, elem->regionId);
  assert(!(coef == NULL));

  // we temporarily switch the coef to one of three states and after evaluating the element matrix switch it back to optimization
  assert(!(lower_bimat && (direction != DesignElement::NO_DERIVATIVE && direction != DesignElement::NO_MULTIMATERIAL)));

  if(lower_bimat)
  {
    DesignSpace::DesignRegion* dr = space->GetRegion(elem->regionId, DesignElement::NO_TYPE, -1, false); // tolerant for off-design stress constraints
    if(dr == NULL)
      throw Exception("cannot find bimaterial for region " + domain->GetGrid()->GetRegion().ToString(elem->regionId));

    assert(dr->HasBiMaterial());

    MaterialClass mc = NO_CLASS;
    MaterialType  mt = NO_MATERIAL;
    switch(system_)
    {
    case MECH:
      mc = MECHANIC;
      assert(integrator == "LinElastInt" || integrator == "MassInt");
      mt = integrator == "LinElastInt" ? MECH_STIFFNESS_TENSOR : DENSITY;
      break;
    case ELEC:
      mc = ELECTROSTATIC;
      break;
    case PIEZOCOUPLING:
      mc = PIEZO;
      break;
    default:
      assert(false);
    }

    PtrCoefFct bimat = dr->GetBiMaterial(mc, mt);

    LOG_DBG3(om) << "GEM: shadow=" << bimat->ToString();

    coef->SetToShadow(bimat);
  }
  if(!lower_bimat)
    coef->SetToOrgMaterial();
  if(direction != DesignElement::NO_DERIVATIVE && direction != DesignElement::NO_MULTIMATERIAL)
    coef->SetToMaterialDerivative(direction);

  c->GetIntegrator()->CalcElementMatrix(out, it, it);
  coef->SetToOptimization(); // removes the shadow material and direction
}


template <class T>
const Matrix<T>& OptimizationMaterial::GetCachedElementMatrix(const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction)
{
  assert(!(lower_bimat && direction != DesignElement::NO_DERIVATIVE)); // what is with NO_MULTIMATERIAL?!

  LocalElementCache::Type state = lower_bimat ? LocalElementCache::SHADOW : LocalElementCache::ORG;
  if(direction != DesignElement::NO_DERIVATIVE && direction != DesignElement::NO_MULTIMATERIAL)
    state = LocalElementCache::DIRECTION;

  return space->elementCache->CachedElement<T>(integrator, state, elem, direction);
}

const DenseMatrix& OptimizationMaterial::GetCachedElementMatrixDM(const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction)
{
  if(ComplexElementMatrix(elem->regionId))
    return GetCachedElementMatrix<complex<double> >(integrator, elem, lower_bimat, direction);
  else
    return GetCachedElementMatrix<double>(integrator, elem, lower_bimat, direction);
}

void OptimizationMaterial::GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem, BaseMaterial* bimaterial, const Vector<double>* ts)
{
  assert(false);
  // FIXME GetElementEntity(form, NULL, &out, elem, bimaterial, DesignElement::NO_DERIVATIVE, ts);

  LOG_DBG3(om) << "CalcElemVector for " << form->GetName() << " -> " << out.ToString();
}


/*
const Matrix<double>& OptimizationMaterial::GeneralStiffness(std::map<RegionIdType, StdVector<Matrix<double> > >& map, const DesignElement* de, MaterialClass mc,
                                                       DesignElement::Type direction, double factor, bool transposed)
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
    //BaseMaterial* bm = NULL;
    // if(de != NULL && de->multimaterial != NULL) // TODO is actually nonsense but due to non-region dependency of DesignSpace::multimaterial
    //  bm = space->GetRegion(de->elem->regionId, &space->GetMultiMaterials()[index])->multimaterial->GetMultiMaterial(mc);

    assert(false);
    // GetElementMatrix(GetForm(de->elem->regionId, form), mat, NULL, bm, direction, factor);
    if(transposed)
    {
      Matrix<double> tmp = mat;
      tmp.Transpose(mat);
    }
  }

  return mat;
}

*/


bool OptimizationMaterial::ComplexElementMatrix(RegionIdType reg) const
{
  if(ctxt_->DoBloch())
    return true;

  assert(reg != NO_REGION_ID);

  if(opt != NULL)
    return ctxt_->pde->HasComplexMatData(reg);

  assert(false); // why should opt be NULL??
  return false;
}

MechMat::MechMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  Init();
}

MechMat::MechMat(DesignSpace* space) : OptimizationMaterial(NULL, Optimization::context, space)
{
  Init();
}

void MechMat::Init()
{
  system_ = MECH;
}


const DenseMatrix& MechMat::MechStiffness(const Elem* elem, bool bimaterial, int multimaterial, DesignElement::Type direction)
{

  // index is 0 for standard material, 1 for bimaterial and anything >= 0 for a real multimaterial index
  //unsigned int index = multimaterial < 0 ? (bimaterial ? 1 : 0) : (unsigned int) multimaterial;

  // in the multimaterial case we do not want the form to obtain the multimaterial tensor but only the one of the currently material
  if(space->HasMultiMaterial()) // ??? and why not in mass??
    direction = DesignElement::NO_MULTIMATERIAL;

  // in the bloch case a change of the wave vector requires to calculate new stiffness matrices
  //bool new_wave_vector = ctxt_->DoBloch() && ctxt_->GetEigenFrequencyDriver()->GetCurrentWaveVector().NormL2() != current_wave_vector_[reg_id][index];

  LOG_DBG3(om) << "MS: el=" << elem->elemNum << " bi=" << bimaterial << " mm=" << multimaterial << " d=" << direction; // << " index=" << index;

  return GetCachedElementMatrixDM("LinElastInt", elem, bimaterial, direction);
}

const DenseMatrix& MechMat::MechMass(const Elem* elem, bool bimaterial, int multimaterial, DesignElement::Type direction)
{
  // see MechStiffness, almost copy & paste :(
  // unsigned int index = multimaterial < 0 ? (bimaterial ? 1 : 0) : (unsigned int) multimaterial;

  return GetCachedElementMatrixDM("MassInt", elem, bimaterial, direction);
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


AcouMat::AcouMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  system_ = ACOUSTIC;
}



const Matrix<double>& AcouMat::AcouStiffness(const Elem* elem, bool bimaterial)
{
  return GetCachedElementMatrix<double>("LaplaceInt", elem, bimaterial);
}

const Matrix<double>& AcouMat::AcouMass(const Elem* elem, bool bimaterial)
{
  return GetCachedElementMatrix<double>("MassInt", elem, bimaterial);
}



PiezoElecMat::PiezoElecMat(ErsatzMaterial* em, Context* ctxt) : MechMat(em, ctxt)
{
  system_ = PIEZOCOUPLING;
/*
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
  */
}



const Matrix<double>& PiezoElecMat::ElecStiffnessPos(const DesignElement* de, DesignElement::Type direction)
{
  /*
  Matrix<double>& mat = GeneralStiffness(elecStiffness_map, de, ELECTROSTATIC, direction, -1.0 , false); // piezo form already negates
  mat *= -1.0;
  assert(direction != DesignElement::NO_DERIVATIVE || mat.Trace() >= 0.0); // allow zero in case of piezo-fmo for non-DIELEC direction
  return mat;
  */
  assert(false);
  return GetCachedElementMatrix<double>("killme", NULL, false);
}

const Matrix<double>& PiezoElecMat::ElecStiffnessNeg(const DesignElement* de, DesignElement::Type direction)
{
/*
  Matrix<double>& mat = GeneralStiffness(elecStiffness_neg_map, de, ELECTROSTATIC, direction, 1.0, false); // piezo form already negates
  //   type.Add(NO_MULTIMATERIAL, "no_multimaterial");LOG_DBG3(om) << "PEM:ESN de=" << de->ToString() << " dt=" << DesignElement::type.ToString(direction) << " -> " << mat.ToString();
  assert(direction != DesignElement::NO_DERIVATIVE || mat.Trace() <= 0.0); // allow zero in case of piezo-fmo for non-DIELEC direction
  return mat;
  */
  assert(false);
  return GetCachedElementMatrix<double>("killme", NULL, false);
}


const Matrix<double>& PiezoElecMat::CoupledStiffness(const DesignElement* de, DesignElement::Type direction)
{
  //return GeneralStiffness(coupledStiffness_map, de, PIEZO, direction, 1.0, false);
  assert(false);
  return GetCachedElementMatrix<double>("killme", NULL, false);
}


const Matrix<double>& PiezoElecMat::CoupledStiffnessTransposed(const DesignElement* de, DesignElement::Type direction)
{
  //return GeneralStiffness(coupledStiffnessTransposed_map, de, PIEZO, direction, 1.0, true);
  assert(false);
  return GetCachedElementMatrix<double>("killme", NULL, false);
}

HeatMat::HeatMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  system_ = HEAT;
}

const DenseMatrix& HeatMat::HeatStiffness(const Elem* elem, bool bimaterial, int multimaterial)
{
  return GetCachedElementMatrixDM("HeatConductivity", elem, bimaterial);
}

ElecMat::ElecMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  system_ = ELEC;
/*
  SinglePDE* elec = ctxt_->ToPDE(App::ELEC);
  assert(elec != NULL);

  for(unsigned int r=0; r < regionIds.GetSize(); r++)
  {
    RegionIdType reg_id = regionIds[r];
    DesignSpace::DesignRegion* dr = opt->GetDesign()->GetRegion(reg_id);

    Matrix<double> tmp_mat;
    assert(false);
    // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat);
    elecStiffness_map[reg_id].first.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
    elecStiffness_map[reg_id].first.SetPart(Global::REAL, tmp_mat);
    if (elec->HasComplexMatData(reg_id)){
      assert(false);
      // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat);
      elecStiffness_map[reg_id].first.SetPart(Global::IMAG, tmp_mat);
    }
    LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id)
                        << std::endl << elecStiffness_map[reg_id].first.ToString(0,true);

    if(dr->HasBiMaterial())
    {
      elecStiffness_map[reg_id].second.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      assert(false);
      // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
      elecStiffness_map[reg_id].second.SetPart(Global::REAL, tmp_mat);
      if (elec->HasComplexMatData(reg_id)){
        assert(false);
        // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
        elecStiffness_map[reg_id].second.SetPart(Global::IMAG, tmp_mat);
      }
      LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << " bimaterial"
          << std::endl << elecStiffness_map[reg_id].second.ToString(0,true);
    }

  }
  */
}



void ElecMat::ReInit()
{
  /*
  SinglePDE* elec = ctxt_->ToPDE(App::ELEC);
  assert(elec != NULL);

  for(unsigned int r=0; r < regionIds.GetSize(); r++)
  {
    RegionIdType reg_id = regionIds[r];
    DesignSpace::DesignRegion* dr = opt->GetDesign()->GetRegion(reg_id);

    Matrix<double> tmp_mat;
    assert(false);
    // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat);
    elecStiffness_map[reg_id].first.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
    elecStiffness_map[reg_id].first.SetPart(Global::REAL, tmp_mat);
    if (elec->HasComplexMatData(reg_id)){
      assert(false);
      // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat);
      elecStiffness_map[reg_id].first.SetPart(Global::IMAG, tmp_mat);
    }
    LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id)
                        << std::endl << elecStiffness_map[reg_id].first.ToString(0,true);

    if(dr->HasBiMaterial())
    {
      elecStiffness_map[reg_id].second.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      assert(false);
      // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::REAL), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
      elecStiffness_map[reg_id].second.SetPart(Global::REAL, tmp_mat);
      if (elec->HasComplexMatData(reg_id)){
        assert(false);
        // GetElementMatrix(GetForm(reg_id, "linGradBDBInt", Global::IMAG), tmp_mat, NULL, dr->GetBiMaterial(ELECTROSTATIC));
        elecStiffness_map[reg_id].second.SetPart(Global::IMAG, tmp_mat);
      }
      LOG_DBG(om) << "OptElecMat: ElecStiffness region=" << domain->GetGrid()->GetRegion().ToString(reg_id) << " bimaterial"
          << std::endl << elecStiffness_map[reg_id].second.ToString(0,true);
    }
  }
  */
}


const Matrix<complex<double> >& ElecMat::ElecStiffness(const Elem* elem, bool bimaterial, DesignElement::Type direction)
{
  /*
  if(!opt->IsDomainStructured() || direction != DesignElement::NO_DERIVATIVE)
  {
    Matrix<double> tmp_mat;
    if(!bimaterial){
      assert(false);
      // GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::REAL), tmp_mat, elem, NULL, direction);
      elecStiffness_map[elem->regionId].first.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      elecStiffness_map[elem->regionId].first.SetPart(Global::REAL, tmp_mat);
      if (ctxt_->IsComplex()){
        assert(false);
        // GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::IMAG), tmp_mat, elem, NULL, direction);
        elecStiffness_map[elem->regionId].first.SetPart(Global::IMAG, tmp_mat);
        //  GetElementMatrix(opt->GetForm(elem->regionId, elec, elec, "linGradBDBInt"), elecStiffness_map[elem->regionId].first, elem, NULL, direction);
      }
    }
    else
    {
      //BaseMaterial* bm = opt->GetDesign()->GetRegion(elem->regionId)->GetBiMaterial(ELECTROSTATIC);
      assert(false);
      // GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::REAL), tmp_mat, elem, bm, direction);
      elecStiffness_map[elem->regionId].second.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      elecStiffness_map[elem->regionId].second.SetPart(Global::REAL, tmp_mat);
      if (ctxt_->IsComplex()){
        assert(false);
        // GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::IMAG), tmp_mat, elem, bm, direction);
        elecStiffness_map[elem->regionId].second.SetPart(Global::IMAG, tmp_mat);
        //  GetElementMatrix(kate GetForm(elem->regionId, elec, elec, "linGradBDBInt"), elecStiffness_map[elem->regionId].second, elem, bm, direction);
      }
    }
  }

  return !bimaterial ? elecStiffness_map[elem->regionId].first : elecStiffness_map[elem->regionId].second;
  */
  assert(false);
  return GetCachedElementMatrix<complex<double> >("killme", NULL, false);
}

LBMMat::LBMMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  system_ = LBM;
}


// explicit template instantiation for GCC compiler
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template void OptimizationMaterial::GetElementMatrix<double>(Matrix<double>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType);
template void OptimizationMaterial::GetElementMatrix<Complex>(Matrix<Complex>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType);
#endif


