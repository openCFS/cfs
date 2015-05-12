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

using std::complex;

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

  needs_mass_ = em->IsComplex() || em->IsTransient() || em->IsEigenvalue();
  complex_ =  em->IsComplex();
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


template <class T>
void OptimizationMaterial::GetElementMatrix(Matrix<T>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType)
{
  LOG_DBG3(om) << "GEM int=" << integrator << " elem=" << (elem != NULL ? elem->elemNum : 4711) << " lb=" << lower_bimat << " d=" << direction << " et=" << entryType;

  assert(entryType != Global::IMAG);
  assert(elem != NULL);

  BiLinFormContext* c = GetPDE()->GetAssemble()->GetBiLinForm(integrator, elem->regionId, GetPDE(), GetPDE(), false);

  // create an element list to gain the iterator in the loop
  ElemList elemList(domain->GetGrid());
  elemList.SetElement(elem);
  EntityIterator it = elemList.GetIterator();

  shared_ptr<CoefFunctionOpt> coef = GetMatCoef(integrator, c, elem->regionId);

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
    coef->SetToTensorDerivative(direction);

  c->GetIntegrator()->CalcElementMatrix(out, it, it);

  coef->SetToOptimization(); // removes the shadow material and direction
}


void OptimizationMaterial::GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem, BaseMaterial* bimaterial, const Vector<double>* ts)
{
  assert(false);
  // FIXME GetElementEntity(form, NULL, &out, elem, bimaterial, DesignElement::NO_DERIVATIVE, ts);

  LOG_DBG3(om) << "CalcElemVector for " << form->GetName() << " -> " << out.ToString();
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

bool OptimizationMaterial::ComplexElementMatrix(RegionIdType reg) const
{
  if(domain->GetDriver()->DoBlochModeEigenfrequency())
    return true;

  assert(reg != NO_REGION_ID);

  if(opt != NULL)
    return opt->pde->HasComplexMatData(reg);

  assert(false); // why should opt be NULL??
  return false;
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
    bool bimat = space->GetRegion(reg_id)->HasBiMaterial();
    Elem* elem = dynamic_cast<GridCFS*>(domain->GetGrid())->SearchFistRegionElement(reg_id);

    // 1 for normal, 2 for bimaterial, n for multimaterial
    unsigned int max = space->HasMultiMaterial() ? mm.GetSize() : (bimat ? 2 : 1);

    // first the bloch case
    if(domain->GetDriver()->DoBlochModeEigenfrequency())
      current_wave_vector_[reg_id].Resize(max, -1.0);

    if(ComplexElementMatrix(reg_id)) {
      mechStiffness_mapC[reg_id].Resize(max);
      if(needs_mass_)
        mechMass_mapC[reg_id].Resize(max);
    }
    else {
      mechStiffness_map[reg_id].Resize(max);
      if(needs_mass_)
        mechMass_map[reg_id].Resize(max);
    }

    // the normal and multimaterial case, not the bimaterial case
    for(unsigned int m = 0; m < (space->HasMultiMaterial() ? mm.GetSize() : 1); m++) {
      MechStiffness(elem, false, m, DesignElement::NO_DERIVATIVE, true); // no bimaterial and enforce_unstructured
      if(needs_mass_)
        MechMass(elem, false, m, DesignElement::NO_DERIVATIVE, true);
    }

    // now the bimaterial case
    if(bimat) {
      MechStiffness(elem, true, -1, DesignElement::NO_DERIVATIVE, true);
      if(needs_mass_)
        MechMass(elem, true, -1, DesignElement::NO_DERIVATIVE, true);
    }
  }
}


const DenseMatrix& MechMat::MechStiffness(const Elem* elem, bool bimaterial, int multimaterial, DesignElement::Type direction, bool enforce_unstructured)
{
  RegionIdType reg_id = elem->regionId;

  // index is 0 for standard material, 1 for bimaterial and anything >= 0 for a real multimaterial index
  unsigned int index = multimaterial < 0 ? (bimaterial ? 1 : 0) : (unsigned int) multimaterial;
  bool complexMaterial = ComplexElementMatrix(reg_id);

  // in the multimaterial case we do not want the form to obtain the multimaterial tensor but only the one of the currently material
  if(space->HasMultiMaterial()) // ??? and why not in mass??
    direction = DesignElement::NO_MULTIMATERIAL;

  // in the bloch case a change of the wave vector requires to calculate new stiffness matrices
  bool new_wave_vector = domain->GetDriver()->DoBlochModeEigenfrequency() && dynamic_cast<EigenFrequencyDriver*>(domain->GetDriver())->GetCurrentWaveVector().NormL2() != current_wave_vector_[reg_id][index];

  LOG_DBG3(om) << "MS: el=" << elem->elemNum << " bi=" << bimaterial << " mm=" << multimaterial << " d=" << direction << " index=" << index << " es=" << enforce_unstructured << " s=" << structured_ << " nwv=" << new_wave_vector;

  StdVector<Matrix<double> >&  K =  mechStiffness_map[reg_id];
  StdVector<Matrix<Complex> >& KC = mechStiffness_mapC[reg_id];

  // we need to get the data on request
  if(new_wave_vector || enforce_unstructured || !structured_ || direction != DesignElement::NO_DERIVATIVE)
  {
   if(ComplexElementMatrix(reg_id))
      GetElementMatrix<complex<double> >(KC[index], "LinElastInt", elem, bimaterial, direction); // in the bimaterial case the standard (upper) material
    else
      GetElementMatrix<double>(K[index], "LinElastInt", elem, bimaterial, direction); // in the bimaterial case the standard (upper) material

   if(domain->GetDriver()->DoBlochModeEigenfrequency())
     current_wave_vector_[reg_id][index] = dynamic_cast<EigenFrequencyDriver*>(domain->GetDriver())->GetCurrentWaveVector().NormL2();
  }

  // e.g. structured case

  if(complexMaterial)  {
    LOG_DBG3(om) << "MS: el=" << elem->elemNum << " -> " << KC[index].ToString();
    return dynamic_cast<DenseMatrix& >(KC[index]); // no multimaterial is frst material
  }
  else {
    LOG_DBG3(om) << "MS: el=" << elem->elemNum << " -> " << K[index].ToString();
    return dynamic_cast<DenseMatrix& >(K[index]);
  }
}

const DenseMatrix& MechMat::MechMass(const Elem* elem, bool bimaterial, int multimaterial, DesignElement::Type direction, bool enforce_unstructured)
{
  // see MechStiffness, almost copy & paste :(
  unsigned int index = multimaterial < 0 ? (bimaterial ? 1 : 0) : (unsigned int) multimaterial;
  bool complexMaterial = ComplexElementMatrix(elem->regionId);

  StdVector<Matrix<double> >&  M =  mechMass_map[elem->regionId];
  StdVector<Matrix<Complex> >& MC = mechMass_mapC[elem->regionId];

  if(enforce_unstructured || !structured_ || direction != DesignElement::NO_DERIVATIVE)
  {
    if(complexMaterial)
      GetElementMatrix<complex<double> >(MC[index], "MassInt", elem, bimaterial, direction); // in the bimaterial case the standard (upper) material
    else
      GetElementMatrix<double>(M[index], "MassInt", elem, bimaterial, direction);
  }

  if(complexMaterial) {
    LOG_DBG3(om) << "MM el=" << elem->elemNum << " c=" << complexMaterial << " es=" << enforce_unstructured << " s=" << structured_<< " bm=" << bimaterial << " d=" << direction << " i=" << index << " -> " << MC[index].ToString();
    return dynamic_cast<DenseMatrix& >(MC[index]);
  }
  else {
    LOG_DBG3(om) << "MM el=" << elem->elemNum << " c=" << complexMaterial << " es=" << enforce_unstructured << " s=" << structured_ << " bm=" << bimaterial << " d=" << direction << " i=" << index << " -> " << M[index].ToString();
    return dynamic_cast<DenseMatrix& >(M[index]);
  }
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
    // BaseMaterial* bm = dr->HasBiMaterial() ? dr->GetBiMaterial(FLUID) : NULL;

    assert(false);
    // GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].first, NULL, NULL);

    if(dr->HasBiMaterial())
    {
      assert(false);
      // GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].second, NULL, bm);
    }

    if(complex_)
    {
      assert(false);
      //GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].first, NULL, NULL);

      if(dr->HasBiMaterial())
      {
        assert(false);
        // GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].second, NULL, bm);
      }
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
    {
      assert(false);
      //GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].first, elem, NULL);
    }
    else
    {
      assert(false);
      // BaseMaterial* bm = space->GetRegion(reg_id)->GetBiMaterial(FLUID);
      // GetElementMatrix(GetForm(reg_id, "LaplaceInt"), acouStiffness_map[reg_id].second, elem, bm);
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
    {
      assert(false);
      // GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].first, elem, NULL);
    }
    else
    {
      assert(false);
      // BaseMaterial* bm = space->GetRegion(reg_id)->GetBiMaterial(FLUID);
      // GetElementMatrix(GetForm(reg_id, "MassInt"), acouMass_map[reg_id].second, elem, bm);
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
}


Matrix<std::complex<double> >& ElecMat::ElecStiffness(const Elem* elem, bool bimaterial, DesignElement::Type direction)
{
  if(!opt->IsDomainStructured() || direction != DesignElement::NO_DERIVATIVE)
  {
    Matrix<double> tmp_mat;
    if(!bimaterial){
      assert(false);
      // GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::REAL), tmp_mat, elem, NULL, direction);
      elecStiffness_map[elem->regionId].first.Resize(tmp_mat.GetNumRows(), tmp_mat.GetNumCols());
      elecStiffness_map[elem->regionId].first.SetPart(Global::REAL, tmp_mat);
      if (opt->IsComplex()){
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
      if (opt->IsComplex()){
        assert(false);
        // GetElementMatrix(GetForm(elem->regionId, "linGradBDBInt", Global::IMAG), tmp_mat, elem, bm, direction);
        elecStiffness_map[elem->regionId].second.SetPart(Global::IMAG, tmp_mat);
        //  GetElementMatrix(opt->GetForm(elem->regionId, elec, elec, "linGradBDBInt"), elecStiffness_map[elem->regionId].second, elem, bm, direction);
      }
    }
  }

  return !bimaterial ? elecStiffness_map[elem->regionId].first : elecStiffness_map[elem->regionId].second;
}

LBMMat::LBMMat(ErsatzMaterial* em) : OptimizationMaterial(em)
{
  system_ = LBM;
  assert(false);
  lbm = NULL; // FIXME dynamic_cast<LatticeBoltzmannPDE*>(opt != NULL ? opt->ToPDE(Optimization::LBM) : domain->GetSinglePDE("LatticeBoltzmann"));
  assert(lbm != NULL);
}


// explicit template instantiation for GCC compiler
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template void OptimizationMaterial::GetElementMatrix<double>(Matrix<double>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType);
template void OptimizationMaterial::GetElementMatrix<Complex>(Matrix<Complex>& out, const std::string& integrator, const Elem* elem, bool lower_bimat, DesignElement::Type direction, Global::ComplexPart entryType);
#endif


