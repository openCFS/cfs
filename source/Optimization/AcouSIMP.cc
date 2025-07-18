
#include "Optimization/AcouSIMP.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"

namespace CoupledField {
  class DenseMatrix;
  class AcousticMaterial;
} // namespace CoupledField

using namespace CoupledField;
using std::complex;

DEFINE_LOG(acs, "acouSimp");

AcouSIMP::AcouSIMP()
{
}

AcouSIMP::~AcouSIMP()
{
}

void AcouSIMP::PostInit(){
  SIMP::PostInit();
  assert(design != nullptr);
  assert(context->pde != nullptr);
  // we get the second material once to have it cached for later use which might be in parallel
  // and causes exception in MathParserOMP::GetNewHandle()
  for (StdVector<DesignSpace::DesignRegion>& drv: design->regions)
    for (DesignSpace::DesignRegion& dr: drv)
      if (dr.HasBiMaterial())
      {
        dr.GetScndMaterial(MaterialClass::ACOUSTIC, MaterialType::DENSITY, context->pde);
        dr.GetScndMaterial(MaterialClass::ACOUSTIC, MaterialType::ACOU_BULK_MODULUS, context->pde);
      }
}

const Complex AcouSIMP::GetExcitationPressure(Function* f) { 
  PtrParamNode info = optInfoNode->Get(ParamNode::HEADER)->Get("getExcitationPressure");
  // get surfRegion
  const RegionIdType reg = AcouSIMP::GetExcitationRegion(f);
  info->Get("surfRegionId")->SetValue(reg);

  // get normalVelocity (ofc this only works if it is a number)
  // todo: expand this for coefFunction
  double vel = 0.0;
  PtrParamNode bcs = f->ctxt->pde->GetParamNode()->Get("bcsAndLoads");
  PtrParamNode velNode = bcs->GetByVal("normalVelocity", "name", grid->GetRegionName(reg));
  // If value is not a double, the As() method throws an Exception.
  vel = velNode->Get("value")->As<double>();
  LOG_DBG2(acs) << "AS: vel=" << vel;
  info->Get("vel")->SetValue(vel);

  // get volumeRegion
  const RegionIdType volreg = AcouSIMP::GetExcitationRegion(f, "volumeNeighbour");
  info->Get("volumeRegionId")->SetValue(volreg);

  // get c0 and K coef functions
  std::map<RegionIdType, BaseMaterial*> mat = context->pde->GetMaterialData();
  auto it = mat.find(volreg);
  if (it == mat.end())
    throw Exception("Given surfRegion not assigned to a material.");
  PtrCoefFct dens = it->second->GetScalCoefFnc(DENSITY, Global::REAL);
  PtrCoefFct blk = it->second->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
  const CoefFunctionConst<Double>* cfcdens = dynamic_cast<const CoefFunctionConst<Double>*>(dens.get());
  const CoefFunctionConst<Double>* cfcblk = dynamic_cast<const CoefFunctionConst<Double>*>(blk.get());
  if(cfcdens == nullptr || cfcblk == nullptr)
    throw Exception("reflectedWave only works for constant real valued material paramters.");

  // get c0 and K
  double rho = cfcdens->GetScalar();
  double K = cfcblk->GetScalar();
  LOG_DBG2(acs) << "AS: rho=" << rho << " K=" << K;
  info->Get("rho")->SetValue(rho);
  info->Get("K")->SetValue(K);

  // compute ep
  Complex ep = - rho * std::sqrt(K/rho) * vel;
  info->Get("excitationPressure_ideal")->SetValue(ep);

  // if ABC exists divide by 2
  if (bcs->HasByVal("absorbingBCs", "name", grid->GetRegionName(reg))) {
    ep /= 2.0;
    LOG_DBG2(acs) << "AS: Found absorbingBC";
  }
  info->Get("excitationPressure_abc")->SetValue(ep);
  return ep;
}

const RegionIdType AcouSIMP::GetExcitationRegion(Function* f, const std::string& attr) {
  // the paramnode of the function should have an output
  if(!f->pn->Has("output/acoustic")) 
    throw Exception("Output must have an acoustic child element.");
  PtrParamNode acoustic = f->pn->Get("output/acoustic");

  // Return the attr Id
  if(!acoustic->Has(attr))
    throw Exception("For reflectedWave output needs to specify the " + attr +" of the normalVelocity b.c.!");
  return grid->GetRegionId(acoustic->Get(attr)->As<std::string>());
}
