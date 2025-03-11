#include <assert.h>
#include <ostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
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
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "General/defs.hh"
#include "Optimization/Context.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/LocalElementCache.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/H1Elems.hh"

namespace CoupledField {
class BaseMaterial;
}  // namespace CoupledField

using std::complex;

using namespace CoupledField;

DEFINE_LOG(om, "optimizationMaterial")


Enum<OptimizationMaterial::System> OptimizationMaterial::system;


OptimizationMaterial::OptimizationMaterial(ErsatzMaterial* em, Context* ctxt, DesignSpace* space)
{
  assert((em == NULL && space != NULL) || (em != NULL && (space == NULL || space == em->GetDesign())));

  this->ctxt_ = ctxt;
  this->opt = em;
  this->space = space != NULL ? space : em->GetDesign();

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

  case MAG:
    return new MagMat(em, ctxt);

  case ACOUSTIC:
    return new AcouMat(em, ctxt);

  case LBM:
    return new LBMMat(em, ctxt);

  default:
    assert(false);
    return NULL;
  }
}

CoefFunction* OptimizationMaterial::GetOrgMatCoef(BaseBDBInt* bdb)
{
  assert(bdb != NULL);
  assert(bdb->GetCoef());

  CoefFunctionOpt* opt = dynamic_cast<CoefFunctionOpt*>(bdb->GetCoef().get());
  if(opt != NULL)
    return opt->orgMat.get();
  else
    return bdb->GetCoef().get();
}

CoefFunctionOpt* OptimizationMaterial::GetMatCoef(BiLinFormContext* context)
{
  assert(context != NULL);
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->GetIntegrator());
  assert(bdb != NULL);
  assert(bdb->GetCoef());
  return dynamic_cast<CoefFunctionOpt*>(bdb->GetCoef().get());
}

inline CoefFunctionOpt* OptimizationMaterial::GetMatCoef(LinearFormContext* context)
{
  assert(context != NULL);
  App::Type app =  ctxt_->ToApp();
  PtrCoefFct coef;
  if (app == App::MAG) {
    BUIntegrator<>* bu = dynamic_cast<BUIntegrator<>*>(context->GetIntegrator());
    assert(bu != NULL);
    assert(bu->GetCoef());
    coef = bu->GetCoef();
  } else if (app == App::HEAT) {
    // FIXME Nasty stuff
    if (domain->GetGrid()->GetDim() == 2) {
      BDUIntegrator<GradientOperator<FeH1,2> >* bdu = dynamic_cast<BDUIntegrator<GradientOperator<FeH1,2> >*>(context->GetIntegrator());
      assert(bdu != NULL);
      assert(bdu->GetDCoef());
      coef = bdu->GetDCoef();
    } else {
      assert(domain->GetGrid()->GetDim() == 3);
      BDUIntegrator<GradientOperator<FeH1,3> >* bdu = dynamic_cast<BDUIntegrator<GradientOperator<FeH1,3> >*>(context->GetIntegrator());
      assert(bdu != NULL);
      assert(bdu->GetDCoef());
      coef = bdu->GetDCoef();
    }
  } else
    EXCEPTION("OM: GetMatCoef only implemented for application MAG or HEAT.");

  return dynamic_cast<CoefFunctionOpt*>(coef.get());
}

CoefFunctionOpt* OptimizationMaterial::GetMatCoef(const string& integrator, RegionIdType reg_id)
{
  SinglePDE* pde = ctxt_->pde;
  assert(pde != NULL);

  BiLinFormContext* c = pde->GetAssemble()->GetBiLinForm(integrator, reg_id, pde, pde, false);
  return GetMatCoef(c);
}

template <class T>
const Matrix<T>& OptimizationMaterial::ComputeElementMatrix(Matrix<T>& out, const FormID& form_id, const Elem* elem, shared_ptr<CoefFunction> shadow)
{
  LOG_DBG3(om) << "GEM form=" << form_id.integrator << " elem=" << (elem != NULL ? elem->elemNum : 4711) << " shadow=" << shadow->ToString();

    assert(elem != NULL);

  SinglePDE* pde = ctxt_->pde;
  assert(pde != NULL);

  BiLinFormContext* c = pde->GetAssemble()->GetBiLinForm(form_id.integrator, elem->regionId, pde, pde, false);

  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  elemList.SetElement(elem);
  EntityIterator it = elemList.GetIterator();

  CoefFunctionOpt* coef = GetMatCoef(c);
  assert(!(coef == NULL));

  CoefFunctionOpt::State old_state = coef->GetState();
  assert(old_state == CoefFunctionOpt::ORG || old_state == CoefFunctionOpt::OPT); // otherwise we would need to store also shadow/direction
  coef->SetToShadow(shadow);

  // let CFS do the hard stuff
  // TODO! Check that the local element cache does not do any harm!!
  c->GetIntegrator()->CalcElementMatrix(out, it, it);

  if(old_state == CoefFunctionOpt::ORG)
    coef->SetToOrgMaterial();
  else
    coef->SetToOptimization();

  return out;
}


template <class T>
const Matrix<T>& OptimizationMaterial::ComputeElementMatrix(Matrix<T>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType)
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

  CoefFunctionOpt* coef = GetMatCoef(c);
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
      // TODO: add buckling with integrator "PreStressInt"
      assert(integrator == "LinElastInt" || integrator == "MassInt" || integrator == "PreStressInt");
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

    PtrCoefFct bimat = dr->GetScndMaterial(mc, mt);

    LOG_DBG3(om) << "GEM: shadow=" << bimat->ToString();

    coef->SetToShadow(bimat);
  }
  if(!lower_bimat)
    coef->SetToOrgMaterial();
  if(direction != DesignElement::NO_DERIVATIVE && direction != DesignElement::NO_MULTIMATERIAL)
    coef->SetToMaterialDerivative(direction);

  // let CFS do the hard stuff
  c->GetIntegrator()->CalcElementMatrix(out, it, it);

  coef->SetToOptimization(); // removes the shadow material and direction

  return out;
}

const DenseMatrix& OptimizationMaterial::GetElementMatrix(OptimizationMaterial::FormID& id, const Elem* elem, bool bimaterial, int multimaterial, DesignElement::Type mat_deriv)
{

  // index is 0 for standard material, 1 for bimaterial and anything >= 0 for a real multimaterial index
  //unsigned int index = multimaterial < 0 ? (bimaterial ? 1 : 0) : (unsigned int) multimaterial;


  // in the bloch case a change of the wave vector requires to calculate new stiffness matrices
  //bool new_wave_vector = ctxt_->DoBloch() && ctxt_->GetEigenFrequencyDriver()->GetCurrentWaveVector().NormL2() != current_wave_vector_[reg_id][index];

  LOG_DBG(om) << "OM:S: sys=" << system_ << " el=" << elem->elemNum << " bi=" << bimaterial << " mm=" << multimaterial; // << " index=" << index;

  LocalElementCache* lec = space->elementCache;
  bool is_complex = ComplexElementMatrix(elem->regionId);
  int design_id = ctxt_->DoBuckling() ? space->GetCurrentDesignId() : -1;

  if(bimaterial)
  {
    DesignSpace::DesignRegion* dr = space->GetRegion(elem->regionId, id.mc, id.mt);
    assert(dr != NULL);

    if(lec == NULL || !lec->HasCachedData(id.integrator, lec->SHADOW, DesignElement::NO_DERIVATIVE, dr->GetScndMaterial(id.mc, id.mt)))
    {
      if(is_complex)
        return ComputeElementMatrix(id.calc_cplx.Mine(), id.integrator, elem, true);
      else
        return ComputeElementMatrix(id.calc_real.Mine(), id.integrator, elem, true);
    }
    else
    {
      if(is_complex)
        return lec->CachedShadowElement<Complex>(id.integrator, elem, dr->GetScndMaterial(id.mc, id.mt));
      else
        return lec->CachedShadowElement<double>(id.integrator, elem, dr->GetScndMaterial(id.mc, id.mt));
    }
  }

  if(mat_deriv != DesignElement::NO_DERIVATIVE && mat_deriv != DesignElement::NO_MULTIMATERIAL)
  {
    if(ctxt_->DoBuckling())
    {
      // One (macroscopic) buckling constraint depends on the derivative of the
      // stiffness matrix twice and the derivative of the stress stiffness
      // matrix once (see ErsatzMaterial::CalcEigenvalueDerivativeBuckling). If
      // we have multiple constraints, these derivatives don't change. Thus, we
      // use the LocalElementCache to store dK_e and dG_e for an iteration.
      if(lec != NULL && !lec->HasCachedData(id.integrator, lec->DIRECTION, mat_deriv, PtrCoefFct(), design_id))
      {
        // remove the previously cached data to save memory
        if(space->GetCurrentDesignId() > 0)
          lec->ClearMatDeriv(id.integrator, elem->regionId, mat_deriv, design_id - 1);

        // compute element matrix and store it in cache
        lec->SetMatDeriv(id.integrator, elem->regionId, mat_deriv, design_id);

        // call GetElementMatrix again. this time, lec->HasCachedData will be True
        return GetElementMatrix(id, elem, bimaterial, multimaterial, mat_deriv);
      }
    }

    if(lec == NULL || !lec->HasCachedData(id.integrator, lec->DIRECTION, mat_deriv, PtrCoefFct(), design_id))
    {
      if(is_complex)
        return ComputeElementMatrix(id.calc_cplx.Mine(), id.integrator, elem, false, mat_deriv);
      else
        return ComputeElementMatrix(id.calc_real.Mine(), id.integrator, elem, false, mat_deriv);
    }
    else
    {
      if(is_complex)
        return lec->CachedMatDerivElement<Complex>(id.integrator, elem, mat_deriv, design_id);
      else
        return lec->CachedMatDerivElement<double>(id.integrator, elem, mat_deriv, design_id);
    }
  }

  // the standard simp case with org
  if(lec == NULL || !lec->HasCachedData(id.integrator, lec->ORG))
  {
    if(is_complex)
      return ComputeElementMatrix(id.calc_cplx.Mine(), id.integrator, elem);
    else
      return ComputeElementMatrix(id.calc_real.Mine(), id.integrator, elem);
  }

  if(is_complex)
    return lec->CachedOrgElement<Complex>(id.integrator, elem);
  else
    return lec->CachedOrgElement<double>(id.integrator, elem);
}


void OptimizationMaterial::GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem, BaseMaterial* bimaterial, const Vector<double>* ts)
{
  assert(false);
  // FIXME GetElementEntity(form, NULL, &out, elem, bimaterial, DesignElement::NO_DERIVATIVE, ts);

  LOG_DBG3(om) << "CalcElemVector for " << form->GetName() << " -> " << out.ToString();
}

void OptimizationMaterial::GetOrgElementVector(LinearFormContext* lc, const Elem* elem, Vector<double>& elemVec){
  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  elemList.SetElement(elem);
  EntityIterator it = elemList.GetIterator();

  CoefFunctionOpt* coef = GetMatCoef(lc);
  assert(!(coef == NULL));

  CoefFunctionOpt::State old_state = coef->GetState();
  assert(old_state == CoefFunctionOpt::ORG || old_state == CoefFunctionOpt::OPT); // otherwise we would need to store also shadow/direction

  coef->SetToOrgMaterial();

  lc->GetIntegrator()->CalcElemVector(elemVec, it);

  if(old_state == CoefFunctionOpt::ORG)
    coef->SetToOrgMaterial();
  else
    coef->SetToOptimization();
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
  system_   = MECH;
  stiff.integrator = "LinElastInt";
  stiff.mc = MECHANIC;
  stiff.mt = MECH_STIFFNESS_TENSOR;

  mass.integrator = "MassInt";
  mass.mc = MECHANIC;
  mass.mt = DENSITY;

  gstiff.integrator = "PreStressInt";
  gstiff.mc = MECHANIC;
  gstiff.mt = MECH_STIFFNESS_TENSOR;
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
  system_   = ACOUSTIC;
  assert(false); // mc_!
  stiff.integrator = "LaplaceInt";
  stiff.mc = MECHANIC; // What could match better?!
  stiff.mt = MECH_STIFFNESS_TENSOR;

  mass.integrator = "MassInt";
  mass.mc = stiff.mc;
  mass.mt = DENSITY;
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
  return ComputeElementMatrix<double>(stiff.calc_real.Mine(), "killme", de->elem);
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
  return ComputeElementMatrix<double>(stiff.calc_real.Mine(), "killme", de->elem);
}


const Matrix<double>& PiezoElecMat::CoupledStiffness(const DesignElement* de, DesignElement::Type direction)
{
  //return GeneralStiffness(coupledStiffness_map, de, PIEZO, direction, 1.0, false);
  assert(false);
  return ComputeElementMatrix<double>(stiff.calc_real.Mine(), "killme", de->elem);
}


const Matrix<double>& PiezoElecMat::CoupledStiffnessTransposed(const DesignElement* de, DesignElement::Type direction)
{
  //return GeneralStiffness(coupledStiffnessTransposed_map, de, PIEZO, direction, 1.0, true);
  assert(false);
  return ComputeElementMatrix<double>(stiff.calc_real.Mine(), "killme", de->elem);
}

HeatMat::HeatMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  system_ = HEAT;
  stiff.integrator = "HeatConductivity";
  stiff.mc = THERMIC;
  stiff.mt = HEAT_CONDUCTIVITY_TENSOR;
  // mass does not apply yet
}

Vector<double> HeatMat::CalcElementTemperature(Context* ctxt, const Elem* elem, HeatPDE::TestStrain testStrain)
{
  if (!domain->GetGrid()->IsGridRegular())
    EXCEPTION("CalcElementTemperature is only implemented for regular meshes!");

  RegionIdType regionId = elem->regionId;
  // each region has a different stiffness matrix and material
  const Matrix<double>& Ke = dynamic_cast<const Matrix<double>&>(Stiffness(elem));

  LOG_DBG3(om) << "\nCET: Get elem temp for region : " << domain->GetGrid()->GetRegionName(regionId) << " testStrain:" << testStrain;
  LOG_DBG3(om) << "region= " << domain->GetGrid()->GetRegionName(regionId) << "  ORG Element stiffness matrix: \n" << Ke.ToString();

  // don't compute multiple times
  if (elem_temp.size() > 0 && elem_temp[regionId].GetSize() > 0){
    assert(elem_temp[regionId].GetSize() >= testStrain);
    return elem_temp[regionId][testStrain];
  }

  assert(elem_temp[regionId].GetSize() == 0);

  UInt num_ex = opt->me->GetNumberHomogenization(App::HEAT);
  UInt meta = ctxt->GetExcitation()->meta_index;
  // each column j stores respective solution for test temperature gradient \eps^j
  StdVector<Vector<double> >& temp_on_reg = elem_temp[regionId];
  temp_on_reg.Resize(num_ex);

  for (UInt ex = 0; ex < num_ex; ex++)
  {
    LinearFormContext* lc = NULL;

    // find correct linear form depending on region
    for (LinearFormContext* it_lfc: ctxt->GetExcitation(ex,meta)->forms) {
      if (it_lfc->GetEntities()->GetRegion() == regionId)
        lc = it_lfc;
    }

    assert(lc != NULL);
    assert(lc->GetIntegrator()->GetName() == "TestTempGradInt");

    LOG_DBG3(om) << "for excitation=" << testStrain << " on region " << domain->GetGrid()->GetRegionName(lc->GetEntities()->GetRegion()) << " found form: " << lc->ToString();

    Vector<double> rhs;
    GetOrgElementVector(lc,elem,rhs);
    assert(rhs.GetSize() > 0);


    LOG_DBG3(om) << "elem " << elem->elemNum << " test strain:" << testStrain << " rhs: " << rhs.ToString();
    UInt dim_sol = rhs.GetSize();

    // as Ke is singular, we fix last row and column and solve only for reduced system
    Matrix<double> K_tmp(Ke.GetNumRows()-1,Ke.GetNumCols()-1);
    Ke.GetSubMatrix(K_tmp,1,1);
    //    LOG_DBG3(om) << "submatrix of K_e:\n" << K_tmp.ToString(0);

    Vector<double> f_tmp;
    for (UInt i = 1; i < dim_sol; i++)
      f_tmp.Push_back(rhs[i]);

    //    LOG_DBG3(om) << "sub rhs:\n" << f_tmp.ToString();

    Vector<double> sol_tmp(dim_sol-1);
    K_tmp.DirectSolve(sol_tmp,f_tmp);
    //    LOG_DBG3(om) << "sub of chi_0:\n" << sol_tmp.ToString();

    temp_on_reg[ex].Resize(dim_sol);
    temp_on_reg[ex].Init();
    for (UInt i = 1; i < dim_sol; i++)
      temp_on_reg[ex][i] = sol_tmp[i-1];

    LOG_DBG3(om) << "Test strain:" << testStrain << " chi_0:" << temp_on_reg[ex].ToString();
  }

  return elem_temp[elem->regionId][testStrain];
}

MagMat::MagMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  system_ = MAG;

  stiff.integrator = "CurlCurlIntegrator"; // or -NL!!
  stiff.mc = ELECTROMAGNETIC;
  stiff.mt = MAG_RELUCTIVITY_TENSOR;

  // mass does not apply yet
}

const Vector<double>& MagMat::MagExcitationRHS(const std::string& integrator, const Elem* elem)
{
  LOG_DBG3(om) << "MER int=" << integrator << " elem=" << (elem != NULL ? elem->elemNum : 4711);

  assert(elem != NULL);

  SinglePDE* pde = ctxt_->pde;
  assert(pde != NULL);
  LinearFormContext* lc = pde->GetAssemble()->GetLinForm(integrator, elem->regionId, pde, false);

  GetOrgElementVector(lc,elem,mag_excitation);

  return mag_excitation;
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
  return ComputeElementMatrix<Complex>(stiff.calc_cplx.Mine(), "killme", elem, bimaterial, direction);
}

LBMMat::LBMMat(ErsatzMaterial* em, Context* ctxt) : OptimizationMaterial(em, ctxt)
{
  system_ = LBM;
  // this is a just a fake system. LBM ist not FE-based, hence no element matrices are required
}


// explicit template instantiation
template const Matrix<double>&  OptimizationMaterial::ComputeElementMatrix<double>(Matrix<double>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType);
template const Matrix<Complex>& OptimizationMaterial::ComputeElementMatrix<Complex>(Matrix<Complex>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType);

template const Matrix<double>& OptimizationMaterial::ComputeElementMatrix(Matrix<double>& out, const FormID& form_id, const Elem* elem, shared_ptr<CoefFunction> shadow);
template const Matrix<Complex>& OptimizationMaterial::ComputeElementMatrix(Matrix<Complex>& out, const FormID& form_id, const Elem* elem, shared_ptr<CoefFunction> shadow);


