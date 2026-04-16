
#include "Optimization/AcouSIMP.hh"
#include "Materials/BaseMaterial.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Optimization/Excitation.hh"
#include "PDE/AcousticPDE.hh"

DEFINE_LOG(acs, "acouSimp");

void AcouSIMP::InitSecondMaterialCache()
{
  // we need only stiffness and density for now
  AddSecondMaterialCache(MaterialClass::ACOUSTIC, MaterialType::DENSITY);
  AddSecondMaterialCache(MaterialClass::ACOUSTIC, MaterialType::ACOU_BULK_MODULUS);
}

StdVector<Complex> AcouSIMP::GetExcitationPressureVector(Excitation& excite, Function* f)
{
  assert(context->sequence == excite.sequence);
  assert(f->ctxt == context);
  assert(f->GetType() == Function::REFLECTED_WAVE);

  // parse regions
  const RegionIdType surfReg = AcouSIMP::GetExcitationRegion(f);
  const RegionIdType volreg = AcouSIMP::GetExcitationRegion(f, "volumeNeighbour");
  // get pde
  AcousticPDE* acoupde = dynamic_cast<AcousticPDE*>(context->pde);
  
 // extract the CoefFunction of the normalVelocity boundary condition
  BUIntegrator<Complex, true>* integ = nullptr;
  assert(excite.forms.GetSize() > 0);
  for(LinearFormContext* lfctx : excite.forms) {
    if(lfctx->GetIntegrator()->GetName() == "NormalVelocityIntegrator") {
      if(integ == nullptr)
        integ = dynamic_cast<BUIntegrator<Complex, true>*>(lfctx->GetIntegrator());
      else {
        // maybe also multiple regions?
        throw Exception("Multiple normalVelocities for one excitation not handled");
      }
    }
  }
  assert(integ != nullptr);
  PtrCoefFct coef = integ->GetCoef();
  

  // get c0 and K coef functions
  std::map<RegionIdType, BaseMaterial*> mat = acoupde->GetMaterialData();
  auto it = mat.find(volreg);
  if (it == mat.end())
    throw Exception("Given region not assigned to a material.");
  PtrCoefFct dens = it->second->GetScalCoefFnc(DENSITY, Global::REAL);
  PtrCoefFct blk = it->second->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);

  // check if ABC xor PML exist
  PtrParamNode bcs = acoupde->GetParamNode()->Get("bcsAndLoads");
  PtrParamNode dmp = acoupde->GetParamNode()->Get("dampingList", ParamNode::ActionType::PASS);
  if (bcs->HasByVal("absorbingBCs", "name", grid->GetRegionName(surfReg)) ==
      (dmp != nullptr && dmp->Has("pml")))
    throw Exception("For accurate reflectedWave results, please place an ABC or PML behind the exciation surfRegion (not both)!");

  // get the mapping function
  const shared_ptr<FeSpace>& feSpace = acoupde->GetFeFunction(acoupde->GetFormulation())->GetFeSpace();
  FeSpace::SingleEqnMap& map = feSpace->GetNodeMap();

  // init our output vector with the same size as u
  Vector<Complex>& l = dynamic_cast<Vector<Complex>&>(*(adjoint.Get(excite, f)->GetVector(StateSolution::SEL_VECTOR)));
  StdVector<Complex> zvec(l.GetSize());  // does this init to {0, 0}?

  // init scaling factor to convert coef to vn
  Complex coef2vn{0, 1/excite.GetOmega()};

  StdVector<unsigned> nodeList;
  LocPointMapped lpm;
  for(RegionIdType regId : acoupde->GetRegions()){
    grid->GetNodesByRegion(nodeList, regId);
    for (unsigned nnr : nodeList) {
      // we use the node number to set the location
      lpm.lp.number = nnr;  

      // map node number to euqation number to align it with RHS storage type
      StdVector<int>& eqnrs = map.eqns[(int)nnr];
      assert(eqnrs.GetSize() == 1);
      int eqnr = eqnrs[0] - 1;  // switch to zero based

      if(eqnr < 1 ||  // clamp eqnr, e.g. for nodes outside of the opt region
        eqnr >= (int)(l.GetSize()) ||
        l[eqnr] == 0.0)  // if this DOF is anyways not selected, just skip it
        continue;
      
      // get result
      Complex val;
      coef->GetScalar(val, lpm);
      val *= coef2vn;  // scale val to vn
      
      LOG_DBG2(acs) << "GEPV: nnr=" << nnr << "; eqnr=" << eqnr << "; w=" << excite.GetOmega() << "; vn=" << val
        << "; exidx=" << excite.index << "; form=" << coef->GetName();

      // scale val to pI
      double rho, K;
      dens->GetScalar(rho, lpm);
      blk->GetScalar(K, lpm);
      val *= -rho * std::sqrt(K/rho) / 2;

      zvec[eqnr] = val;
    }
  }

  return zvec;
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
