#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Optimization/Design/FeatureMappingDesign.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "Optimization/Design/MaterialTensor.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "General/Environment.hh" // CFS_NUM_THREADS for the parallel mapping
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include <map>

using std::string;
using std::pair;
using std::make_pair;
using std::to_string;
using DE = DesignElement; // shortcut

namespace CoupledField {

DEFINE_LOG(fm, "featureMapping")

Matrix<double> FeatureMappingDesign::aniso_base_tensor; // static member  

FeatureMappingDesign::FeatureMappingDesign(StdVector<RegionIdType>& regionIds, PtrParamNode empn, ErsatzMaterial::Method method)
: FeaturedDesign(regionIds, empn, method)
{
  setup_timer_->Start();

  assert(method == ErsatzMaterial::FEATURE_MAPPING || method == ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT);

  anisotropic = method == ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT;
  if(anisotropic && dim_ == 3)
    throw Exception("anisotropic feature mapping is only implemented in 2D");
  PtrParamNode pn = empn->Get(ErsatzMaterial::method.ToString(ErsatzMaterial::FEATURE_MAPPING)); // featureMappingAniso has also a featureMapping element

  combine_ = combine.Parse(pn->Get("combine")->As<string>());
  if(combine_ != MAX)
  {
    if(!pn->Has("p"))
      throw Exception("featureMapping attribute 'p' is mandatory for combine " + combine.ToString(combine_));
    cmb_param = pn->Get("p")->As<double>();
  } 
  boundary_ = boundary.Parse(pn->Get("boundary")->As<string>());

  order = pn->Get("order")->As<unsigned int>();
  transition = pn->Get("transition")->As<double>(); // make optional
  extension = pn->Get("extension")->As<double>(); // asymmetric outer zone for bezier, 0 is symmetric

  OpenGradPlot(pn); // check "gradplot" and create file
  hessexport_ = pn->Get("hessexport")->As<bool>();

  // the optional geometry variable alpha (Norato's size variable) scales each feature density as alpha^q * rho
  if(pn->Has("alpha"))
  {
    has_alpha_ = true;
    alpha_q_ = pn->Get("alpha")->Get("q")->As<double>();
    num_var_by_feature = ALPHA_VAR + 1; // nodes, profile, alpha - needs to be set before SetupDesign()
    if(anisotropic)
      throw Exception("the geometry variable 'alpha' is not yet implemented for anisotropic feature mapping");
  }

  // map
  SetupMapping();

  // setup pills with pills and shape_param_ and opt_shape_param_
  SetupDesign(pn);

  // the framework for efficient distance field determination - we need to know the number of features 
  SetupFixedIntegrationPoints();

  // get rhomin from first density element and check transfer function
  assert(!data.IsEmpty());
  const DesignElement& de = data.First();
  rhomin = de.GetLowerBound(); // at least in python use for all
  rhomax = de.GetUpperBound();

  switch(boundary_)
  {
  case POLY:
    bnd_fnc = std::make_unique<Poly>(transition, rhomin);
    break;
  case QUINTIC:
    bnd_fnc = std::make_unique<Quintic>(transition, rhomin);
    break;
  case BEZIER:
    bnd_fnc = std::make_unique<Bezier>(transition, rhomin, extension);
    break;
  default:
    throw Exception("featureMapping boundary '" + boundary.ToString(boundary_) + "' not supported, use poly, quintic or bezier");
  }
  if(extension != 0.0 && boundary_ != BEZIER)
    throw Exception("featureMapping attribute 'extension' is only supported for boundary 'bezier'");

  switch(combine_)
  {
  case MAX:
    cmb_fnc = std::make_unique<Max>();
    if(features_.GetSize() > 1)
      info_->Get("features/combine")->SetWarning("Using non differentiable combine=max with multiple features is not recommended, use p-norm, KS or softmax");
    break;
  case P_NORM:
    cmb_fnc = std::make_unique<P_Norm>(cmb_param);
    break;
  case KS:
    cmb_fnc = std::make_unique<KreisselmeierSteinhauser>(cmb_param);
    break;
  case SOFTMAX:
    cmb_fnc = std::make_unique<SoftMax>(cmb_param);
    break;
  default:
    assert(false);
    break;
  }

  assert(GetNumberOfVariables() > 0);

  if(anisotropic) 
    aniso_current_variable.Reset(); 


  setup_timer_->Stop(); // python and map are subtimers, so count them here
}


void FeatureMappingDesign::PostInit(int objectives, int constraints)
{
  FeaturedDesign::PostInit(objectives, constraints);

  setup_timer_->Start();
  assert(mapped_design_ != design_id);
  
  if (Optimization::context->mat->GetSystem() == OptimizationMaterial::MECH)
  {
    CoefFunctionOpt* coef = Optimization::context->mat->GetMatCoef("LinElastInt", domain->GetDesign()->GetRegionId());
    assert(dynamic_cast<CoefFunctionConst<double>*>(coef->orgMat.get()) != nullptr);
    aniso_base_tensor = dynamic_cast<CoefFunctionConst<double>*>(coef->orgMat.get())->GetTensor();
  } 
  else if(domain->GetDesign()->GetMethod() == ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT)
    throw Exception("Anisotropic material only supported for MechPDE.");
  
  // this calls Pill::Update() and in the anisotropic case we need DesignSpace set up
  MapFeatureToDensity(); // only so late because of python -> calls PythonUpdateSpaghetti()
  LOG_DBG(fm) << "PI data=" << data.GetSize() << " aux=" << aux_design_.GetSize() << " shape=" << opt_shape_param_.GetSize() << " -> N=" << GetNumberOfVariables();
  setup_timer_->Stop(); // python and map are subtimers, so count them here
}


void FeatureMappingDesign::ToInfo(ErsatzMaterial* em)
{
  FeaturedDesign::ToInfo(em);

  PtrParamNode in = info_->Get("features");
  assert(combine_ != NO_COMBINE);
  if(combine_ != MAX)
    in->Get("combine/param")->SetValue(cmb_param);
  in->Get("boundary/transition")->SetValue(transition);
  in->Get("boundary/extension")->SetValue(extension);
}


void FeatureMappingDesign::SetupMapping()
{
  FeaturedDesign::SetupMapping();
  assert(map.GetSize() <= nx_ * ny_ * nz_); // nz_ is 1 for 2D. Not all cells of the virtual box need to be design

  for(Item& item : map)
  {
    ItemIP* item_ip = new ItemIP(); // we use ItemIP as extension
    if(order > 1)
      item_ip->corner.Reserve(dim_ == 2 ? 4 : 8); // max 4 corners for 2D, 8 for 3D
    item.extension = item_ip; // ~Item() will delete this
  }

  // inverse of Item::lexicographic_pos: design element by virtual cell, -1 where the domain is not
  // design (e.g. a fixed region). The element/node numbering is assumed lexicographic (regular mesh)
  lex_to_item_.Resize(n_.Product());
  lex_to_item_.Init(-1);
  for(unsigned int i = 0; i < map.GetSize(); i++)
  {
    assert(map[i].lexicographic_pos >= 0 && map[i].lexicographic_pos < (int) n_.Product());
    lex_to_item_[map[i].lexicographic_pos] = i;
  }
}

const FeatureMappingDesign::Pill& FeatureMappingDesign::GetFeatureForResult(const ResultDescription& rd) const
{
  string value = DE::valueSpecifier.ToString(rd.value);
  
  assert(rd.detail == DE::NONE || rd.detail == DE::GRAD_DISTANCE || rd.detail == DE::FEATURE_PART);

  if(rd.detail != DE::GRAD_DISTANCE && rd.design != DE::FEATURE)
    throw Exception("result '" + value + "'/'none' needs 'design' 'allFeatures' or 'feature'/'generic'");
  if(rd.generic == "")
    throw Exception("result '" + value + "'/'feature' needs 0-based feature index in 'generic'");
  if(!IsType<unsigned int>(rd.generic))  
    throw Exception("result '" + value + "'/'feature' has no number 'generic' but '" + rd.generic + "'");
  
  int fi = std::stoi(rd.generic); 

  if(fi < 0 || fi >= (int) pills.GetSize())
    throw Exception("result '" + value + "'/'feature' needs 0-based feature index smaller " + std::to_string(pills.GetSize()) + " in 'generic'");

  return pills[fi];
}

int FeatureMappingDesign::GetVariableForResult(const ResultDescription& rd) const
{
  // the per-feature variable slot in GetAllVariables()/GradDistance() order:
  // 2D [Px Py Qx Qy p], 3D [Px Py Pz Qx Qy Qz p]. The enum order differs as PZ/QZ were appended
  if((rd.design == DE::FEATURE_MAPPING_PZ || rd.design == DE::FEATURE_MAPPING_QZ) && dim_ == 2)
    throw Exception("result 'design' 'feature_var_Pz'/'feature_var_Qz' is 3D only");
  switch(rd.design)
  {
  case DE::FEATURE_MAPPING_PX: return 0;
  case DE::FEATURE_MAPPING_PY: return 1;
  case DE::FEATURE_MAPPING_PZ: return 2;
  case DE::FEATURE_MAPPING_QX: return dim_;
  case DE::FEATURE_MAPPING_QY: return dim_ + 1;
  case DE::FEATURE_MAPPING_QZ: return 5;
  case DE::FEATURE_MAPPING_P:  return 2 * dim_;
  case DE::FEATURE_MAPPING_ALPHA:
    if(!has_alpha_)
      throw Exception("result 'design' 'feature_var_alpha' needs the 'alpha' element in 'featureMapping'");
    return ALPHA_VAR;
  default:
    throw Exception("expect result 'design' from 'feature_var_Px', ..., 'feature_var_alpha'");
  }
}


StdVector<int> FeatureMappingDesign::GetSpecialResultIndices(const Function* f, DesignElement::Detail detail, bool d_s) const
{
  assert((f != nullptr && detail == DE::NONE) || (f == nullptr && (detail == DE::COMBINE || detail == DE::FEATURE_RHO)));

  assert(shape_param_.GetSize() == pills.GetSize() * num_var_by_feature);

  // we have the case d_J/d_s, d_rho/d_s, d_mrho/d_s and d_mrho/d_f 
  // where m_rho is the combined rho 
  // and J is compliance,volume,... In case you need to add the function string to CFS_Optimization.xsd and the enum.
  
  StdVector<int> res;

  // test for d_mrho/d_feature
  if(d_s == false)
  {
    assert(f == nullptr && detail == DE::COMBINE);
    res.Resize(pills.GetSize(), -1);

    for(const Pill& pill : pills)
    {
      for(const ResultDescription& rd : resultDescriptions)
      {
        bool detail_match = rd.detail == DE::COMBINE || rd.detail == DE::FEATURE_RHO;
        if(rd.value == DE::FEATURE_GRAD && detail_match && rd.design == DE::FEATURE && std::stoi(rd.generic) == pill.idx) 
        {
          res[pill.idx] = rd.solutionType - OPT_RESULT_1;
          break; 
        }
      }
    }
  }
  else
  {
    assert(d_s == true);
    // we have the case d_J/d_s or d_mrho/d_s
    res.Resize(shape_param_.GetSize(), -1);

    DE::Detail test = detail;
    if(f != nullptr)
    {
      std::string name = Function::type.ToString(f->GetType()); // 'compliance', 'volume', ...
      // python (and other non-detail) cost functions have no featureGrad special result to map here;
      // their objective gradient is assembled elsewhere, so just return the empty (-1) index list.
      if(!DE::detail.IsValid(name))
        return res;
      test = DE::detail.Parse(name);
    }
    // we loop over all possible combinations
    for(unsigned int si = 0; si < shape_param_.GetSize(); si++)
    {
      const FeatureVariable* var = dynamic_cast<FeatureVariable*>(shape_param_[si]);
      assert(var != nullptr);
      assert(var->feature >= 0 && var->feature < (int) pills.GetSize());
      assert(var->var != DE::NO_TYPE);
      // check every result (if we have any)
      for(const ResultDescription& rd : resultDescriptions)
      {
        // <result value="featureGrad" detail="compliance" design="feature_var_Px" generic="1" id="optResult_9" />
        LOG_DBG2(fm) << "GSRI si=" << si << " rd.v=" << DE::valueSpecifier.ToString(rd.value) << " var=" << var->var 
                     << " rd.d=" << DE::type.ToString(rd.design) << " f=" << var->feature << " rd.g=" << rd.generic;
        if(rd.value == DE::FEATURE_GRAD && rd.detail == test && rd.design == var->var && std::stoi(rd.generic) == var->feature) 
        {
          res[si] = rd.solutionType - OPT_RESULT_1;
          LOG_DBG2(fm) << "GSRI res[" << si << "] <- " << res[si];
          break; 
        }
      }
    }
  }
  LOG_DBG2(fm) << "GSRI f=" << (f ? f->ToString() : "-") << " detail=" << detail << " d_s=" << d_s << " -> " << res.ToString();
  return res;
}


double FeatureMappingDesign::GetNodalSpecialResult(unsigned int nodeNumber, const ResultDescription& rd) 
{
  // we assume non-parallel execution!
  if(node_ip_result_map.IsEmpty())
    SetupNodeIPMap();
 
  bool prj = rd.value == DE::FEATURE_PROJECTED;
  bool all = rd.design == DE::ALL_FEATURES;
  assert(all || rd.design == DE::FEATURE || (rd.design >= DE::FEATURE_MAPPING_PX && rd.design <= DE::FEATURE_MAPPING_P));

  assert(rd.value == DE::FEATURE_DISTANCE || prj);
  assert(nodeNumber <= node_ip_result_map.GetSize()); // we assume 1-base nodeNumber
  assert(nodeNumber > 0); 
  assert(node_ip_result_map[nodeNumber] != -1);
  LOG_DBG2(fm) << "GNSR nodeNumber=" << nodeNumber << " node_ip_result_map[nn]=" << node_ip_result_map[nodeNumber-1];
 
  const IntegrationPoint& ip = corners[node_ip_result_map[nodeNumber]];

  if(rd.detail != DE::NONE && rd.detail != DE::GRAD_DISTANCE && rd.detail != DE::FEATURE_PART)
    throw Exception("result 'featureDistance/Projected' allows 'detail' 'none'/'gradDistance'/'allows 'detail' 'none'/'gradDistance'");

  // dummy pill if not needed
  const Pill& pill = (rd.design == DE::FEATURE || FeatureVariable::IsFeatureVariable(rd.design))? GetFeatureForResult(rd) : pills[0]; 
  
   if(rd.detail == DE::FEATURE_PART)
   {
     if(all)
     {
       unsigned int idx = ip.dist.ArgMin();
       assert(idx < ip.part.GetSize());
       return (double) 3 * idx + ip.part[idx]; // P, inner, Q
     }
     else
       return (double) ip.part[pill.idx]; // technically -1,0,1,2
   }
    

  // common for detail 'none' (default) and 'gradDistance'
  double dist = -1.0;
  if(rd.design == DE::ALL_FEATURES)
    dist = ip.dist.Min();
  else // rd.design == DE::FEATURE
    dist = ip.dist[pill.idx];
  

  if(rd.detail == DE::NONE)
  {
    return prj ? bnd_fnc->Eval(dist) : dist;
  }
  else
  {
    assert(rd.detail == DE::GRAD_DISTANCE);

    Vector<double> ddist_ds(num_var_by_feature); // we have only the vector grad
    pill.GradDistance(ip.loc, ip.part[pill.idx], ddist_ds); // output is the vector ddist_ds (including fixed)

    double dH_ddist = prj ? bnd_fnc->Grad(dist) : 1.0; // we multiply below    
    int s = GetVariableForResult(rd); // 0 ... 4
    return dH_ddist * ddist_ds[s];
  }
}

void FeatureMappingDesign::SetupNodeIPMap()
{
  // don't call in parallel!
  // when we have nodal special results,  DesignSpace::GetNodalValue() calls GetNodalSpecialResult() for each node number
  // we do here a mapping from the node coordinates to the integration points loc in corners.
  // we need to assume that not the full mesh is design and then search by coordinate.
  if(order < 2)
    throw Exception("order at least 2 required for result featureDistance");

  assert(node_ip_result_map.IsEmpty());

  Grid* grid = domain->GetGrid();
  unsigned int grid_nodes = grid->GetNumNodes();
  node_ip_result_map.Resize(grid_nodes+1, -1); // we have 1-based nodeNumber, -1 means no mapping 

  // we compare by node coordinate and integration point loc.
  // we need this, if not all regions are design domain
  // note that the integration point loc is calculated. One could also go via grid Elem corner nodes

  // note that corners is almost ordered by coordinates but has a hick-up in the beginning

  int next_idx = 0; // index within integration points corner vector
  StdVector<unsigned int> nodeList; // nodes of current region
  Vector<double> coord; // grid coordinate of mesh node (2D)
  for(Grid::RegionData& reg : grid->regionData) 
  {
    grid->GetNodesByRegion(nodeList, reg.id);
    for(unsigned int nn : nodeList)
    {
        grid->GetNodeCoordinate(coord, nn, false); // not updated
        // search an hope for more or less consecutive but have hick-ups when a horizontal line starts (left side).
        next_idx = std::max(next_idx - 2, 0); 
        int cnt = 0;

        while(true)
        {
          IntegrationPoint& ip = corners[next_idx];
          //LOG_DBG3(fm) << "SNIM: nn=" << nn << " si=" << start_idx << " next_idx=" << next_idx << " c=" << coord.ToString() << " ip.l=" << ip.loc.ToString() << " d=" << ip.loc.Dist(coord);

          if(ip.loc.Close(coord))
          {
            node_ip_result_map[nn] = next_idx;
            LOG_DBG3(fm) << "SNIM: nn=" << nn << " idx=" << next_idx << " cnt=" << cnt << " c=" << coord.ToString() << " ip.l=" << ip.loc.ToString();
            break;
          }
          else
          {
            next_idx = next_idx <  (int) corners.GetSize()-1 ? next_idx + 1 : 0; // increment with wrap
            if(++cnt > (int) corners.GetSize())
              throw Exception("cannot match node " + to_string(nn) + " with coord " + coord.ToString() + " to integration points. Wrong nodeResult region?");     
        }
      }
    } // region finished
  } // all regions finished
  LOG_DBG3(fm) << "SNIPM: node_ip_result_map=" << node_ip_result_map.ToString();
}    


bool FeatureMappingDesign::GetAnisoMechTensor(const Elem* elem, MaterialTensor<double>& mt, DesignElement::Type direction, const CoefFunctionOpt* coef, DesignMaterial* dm) const
{
  // be aware, needs to be thread safe!

  // single feature orientation: c_e = R(phi) rho_e c_0 R(phi)^T = rho_e R c_0 R^T = rho_e * rot_mat
  // aggregated features with mrho_e is (smooth) max of rho_e^f for all features f
  // c_e = (sum_f c_e^f / sum_f rho_e^f) * mrho_e
  assert(coef->GetState() == CoefFunctionOpt::OPT || coef->GetState() == CoefFunctionOpt::DIRECTION); 
  assert(dm != nullptr);
  assert(dim_ == 2);
  assert(Find(elem, false) != -1);
  const Item& item = map[Find(elem, false)];
  const ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
  assert(ext->rho.GetSize() == pills.GetSize());

  // sum_f c_e^f / sum_f rho_e^f -> not scaled yet
  // we use mt to store the sum -> for evaluation we need not other variable
  Matrix<double>& out = mt.GetMatrix(VOIGT); // output
  out.Resize(3,3); // 2D
  out.Init(); // set to zero
  // sum_f c_e^f / sum_f rho_e^f) -> sum_c_e / sum_rho where c_e = rho_e * rot_mat
  double sum_rho = ext->rho.Sum();
  assert(sum_rho >= 0);
  for(unsigned int i = 0; i < pills.GetSize(); i++)
    out.Add(ext->rho[i], pills[i].rot_mat); // rot_mat is precomputed in Pill::Update()
  
  // the (smooth) max rho_e to scale the avg tensor
  double mrho_e = cmb_fnc->Eval(ext->rho);

  LOG_DBG2(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " rho_e=" << ext->rho.ToString() << " mrho=" << mrho_e << " sum_c_e=" << out.ToString();

  if(coef->GetState() == CoefFunctionOpt::OPT)
  {
    assert(direction == DE::NO_DERIVATIVE);
    assert(sum_rho > 1e-12 || mrho_e <= 1e-12); // if sum_rho is zero, mrho_e must be zero too
    if(sum_rho > 1e-12)
      out.Mult(mrho_e * 1.0/sum_rho); // scale to get c_e
    else
      out.Init(); // we are in the ground material case, no feature is active -> zero  

    LOG_DBG2(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " m_rho_e=" << mrho_e << " sum_rho=" << sum_rho << " -> " << out.ToString();
  }
  else if(sum_rho > 1e-12) // we may divide by sum_rho
  {
    assert(coef->GetState() == CoefFunctionOpt::DIRECTION);
    // mc_e = (sum_f c_e / sum_f rho_e) * mrho_e
    // d_mc_e/d_s = (d_c_e/d_s * sum_rho - sum_c_e * d_rho_e/d_s) / (sum_f rho_e)^2 * mrho_e 
    //            + (sum_f c_e / sum_f rho_e) * d_mrho_e/d_rho_e * d_rho_e/d_s
    // in variables:
    // dce_ds = (dc_ds * sum_rho - sum_c_e * drho_ds) / (sum_rho)^2 * mrho_e
    //        + (sum_c_e / sum_rho) * dmrho_drho * drho_ds
    int s = direction - DE::FEATURE_MAPPING_PX;
    assert(direction != DE::NO_DERIVATIVE);
    assert(direction == aniso_current_variable.var_type);
    assert(s == aniso_current_variable.var_idx);
    unsigned int fi = aniso_current_variable.feature_idx;   
    
    assert(s >= 0 && s <= 4);
    assert(aniso_current_variable.feature_idx >= 0 && fi < pills.GetSize());  

    const Pill& pill = pills[fi];

    // d_rho_e/d_s is precomputed 
    double drho_ds = ext->drho_ds_full[fi * num_var_by_feature + s]; 

    // the single rho_e for the feature (before combining)
    double rho_e = ext->rho[fi];

    // sum_f c_e^f was stored in out above -> save it
    Matrix<double> sum_c_e(out); 

    // compute the derivate of a single material tensor by variable s of feature fi (no averaging) 
    // C_e = R^T rho_e C_0 R where R=R(phi(s)) and rho_e = rho_e(s)
    // d_C/d_s = 2 R^T rho_e C_0 d_R/d_s + R^T d_rho_e/d_s C_0 R
    //          = rho_e 2 R^T C_0 d_R/d_s + d_rho_e/d_s R^T C_0 R with d_R/d_s = d_R/d_phi d_phi/d_s
    //          = rho_e d_phi/d_s RCdR    + d_rho_e/d_s RCR  with dR = d_R/d_phi
    // where rho_e R^T C_0 d_R/d_s is significant everywhere, where rho_e != 0 (not rho_min)
    // and d_rho_e/d_s R^T C_0 R is significant only where rho_e is > rho_min and < 1
    out.Init(); // we use out as dc_ds
    if(rho_e != 0.0)
    { // we have it also in inner and rho_min
      double dphi_ds = pill.GradAngle(s);
      assert(aniso_current_variable.RCdR.NormL2() > 0);
      out.Add(rho_e * dphi_ds, aniso_current_variable.RCdR); // rho * d_phi/d_s * RCdR (constant for given angle and variable)
    }
    if(rho_e != rhomin && rho_e != 1)
      out.Add(drho_ds, pill.rot_mat);

    // includes angle stuff  
    if(aniso_current_variable.drho_ds_res_idx[fi * num_var_by_feature + s] != -1)
      item.elemval->specialResult[aniso_current_variable.drho_ds_res_idx[fi * num_var_by_feature + s]] = drho_ds;

    // the single derivate by feature and variable d_C_0/d_s is computed 
    LOG_DBG3(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " s=" << s << " fi=" << fi << " rho_e=" << rho_e << " dphi_ds=" << pill.GradAngle(s) << " -> out=dc_ds=" << out.ToString();

    double dmrho_drho = cmb_fnc->Grad(ext->rho, pill.idx);

    if(aniso_current_variable.dmrho_drho_res_idx[fi] != -1)
      item.elemval->specialResult[aniso_current_variable.dmrho_drho_res_idx[fi]] = dmrho_drho;
    if(aniso_current_variable.drho_ds_res_idx[fi * num_var_by_feature + s] != -1)
      item.elemval->specialResult[aniso_current_variable.drho_ds_res_idx[fi * num_var_by_feature + s]] = dmrho_drho * drho_ds;

    // dc_ds * sum_rho - sum_c_e * drho_ds)
    out.Mult(sum_rho);
    out.Add(-drho_ds, sum_c_e);

    LOG_DBG3(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " s=" << s << " fi=" << fi << " sum_rho=" << sum_rho 
                 << " sum_c_e=" << sum_c_e.ToString() << " drho_ds=" << drho_ds << " -> " << out.ToString();
    out.Mult(mrho_e / (sum_rho * sum_rho)); 
    LOG_DBG3(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " s=" << s << " mrho_e=" << mrho_e << " f=" << (mrho_e / (sum_rho * sum_rho)) << " -> " << out.ToString();
    out.Add(dmrho_drho * drho_ds / sum_rho, sum_c_e);
    LOG_DBG3(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " s=" << s << " dmrho_drho=" << dmrho_drho << " f=" << (dmrho_drho * drho_ds / sum_rho) << " -> " << out.ToString();
  }
  else
  {
    assert(coef->GetState() == CoefFunctionOpt::DIRECTION);
    assert(sum_rho <= 1e-12); // if sum_rho is zero, we must be in the ground material case -> no derivative
    out.Init(); 
    LOG_DBG3(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " sum_rho=" << sum_rho << " is zero -> " << out.ToString();
  }
  LOG_DBG3(fm) << "GAT e=" << elem->elemNum << " d=" << direction << " s=" << (direction < DE::FEATURE_MAPPING_PX ? -1 : (direction - DE::FEATURE_MAPPING_PX)) 
               << " fi=" << aniso_current_variable.feature_idx << " -> " << mt.GetMatrix(VOIGT).ToString();

  return true;   
}


inline FeatureMappingDesign::ItemIP* FeatureMappingDesign::GetItemIP(unsigned int item_idx)
{
  assert(map[item_idx].extension != nullptr);
  assert(dynamic_cast<ItemIP*>(map[item_idx].extension) != nullptr);
  return dynamic_cast<ItemIP*>(map[item_idx].extension);
}

inline const FeatureMappingDesign::ItemIP* FeatureMappingDesign::GetItemIP(unsigned int item_idx) const
{
  assert(map[item_idx].extension != nullptr);
  assert(dynamic_cast<const ItemIP*>(map[item_idx].extension) != nullptr);
  return dynamic_cast<const ItemIP*>(map[item_idx].extension);
}


void FeatureMappingDesign::SetupDesign(PtrParamNode pn)
{
  FeaturedDesign::SetupDesign(pn);
 
  // the order is nodes, profile, [alpha]
  for(Pill& pill : pills)
  {
    for(auto& p : pill.points) // P, [inner], Q
      for(auto& v : p)
        AddVariable(&v);

    AddVariable(&(pill.p));

    if(pill.has_alpha)
      AddVariable(&(pill.alpha));
  }

  // post process mapping
  for(Feature& pill : pills)
  {
    for(FeatureVariable* var : pill.GetAllVariables()) 
    {
      if(var->key != "")
        if(CountKey(var->key) > 1)
          throw Exception("feature key '" + var->key + "' occurs more than once. Use unique keys for each variable."); 
      if(var->map != "")
      {
        var->map_to = FindKey(var->map); // important to traverse a non-resizable vector
        if(var->map_to == nullptr)
          throw Exception("There is no feature key '" + var->map + "' to map to");
        var->SetDesign(var->map_to->GetPlainDesignValue()); // set the design value to the mapped variable
      }
    } 
  }

  // we did exact space "estimation"
  LOG_DBG(fm) << "SD: shape_param=" << shape_param_.GetSize() << " opt_shape_param=" << opt_shape_param_.GetSize() 
            << " opt_indices=" << opt_indices.GetSize() << " pills=" << pills.GetSize();
  assert(shape_param_.GetCapacity() == shape_param_.GetSize());
  assert(opt_shape_param_.GetCapacity() == opt_shape_param_.GetSize());
}

void FeatureMappingDesign::SetupFixedIntegrationPoints()
{
  assert(features_.GetSize() >= 1);
  Grid* grid = domain->GetGrid();
  Matrix<double> coords;
  Point center;

  assert(map.GetSize() == nx_ * ny_ * nz_); // nz_ is 1 for 2D

  // all integration points are anchored at the element barycenters.
  // we are not sure about node ordering for pathological cases and therefore recompute everything
  // from the element barycenter. The same time we assume that the map ordering is lexicographic!!
  for(Item& item : map)
  {
    assert(item.elemval != nullptr);
    grid->GetElemNodesCoord(coords, item.elemval->elem->connect, false); // obtain element nodal coords, no update
    LagrangeElemShapeMap::CalcBarycenter(center, coords, domain); // get the barycenter of the element
    dynamic_cast<ItemIP*>(item.extension)->center = center;
    LOG_DBG3(fm) << "SD e=" << item.elemval->elem->elemNum << " c=" << center.ToString();
  }

  // per-item data sizes, independent of the integration order
  unsigned int nf = features_.GetSize();
  for(Item& item : map)
  {
    ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
    ext->rho.Resize(nf);
    if(has_alpha_)
      ext->rho_org.Resize(nf);
    ext->drho_ds_full.Resize(nf * num_var_by_feature);
  }

  // origin (lower node) of the virtual lexicographic box from the first element's barycenter
  int fp = map.First().lexicographic_pos;
  unsigned int fx = fp % nx_, fy = (fp / nx_) % ny_, fz = fp / (nx_ * ny_);
  const Point& fc = GetItemIP(0)->center;
  origin_.Set(fc[0] - (fx + 0.5)*dx_, fc[1] - (fy + 0.5)*dx_, dim_ == 3 ? fc[2] - (fz + 0.5)*dx_ : 0.0);

  if(order == 1)
  {
    // one static barycenter ip per element
    inners_.Resize(map.GetSize());
    num_inner_ = map.GetSize();
    for(unsigned int i = 0; i < map.GetSize(); i++)
    {
      ItemIP* ext = GetItemIP(i);
      IntegrationPoint& ip = inners_[i];
      ip.loc = ext->center;
      ip.dist.Resize(nf);
      ip.part.Resize(nf, FeatureVariable::NO_TIP);
      ext->inner.Push_back(&ip);
    }
  }
  else // for order > 2 the corners are fixed and the inner ip are dynamically handled in MapFeatureToDensity()
  {
    // we traverse the virtual node grid: each node is a corner ip, shared by its up to 4/8
    // adjacent design elements. Nodes without any design element around are skipped (the domain
    // does not need to be fully design)
    unsigned int mz = dim_ == 3 ? nz_ + 1 : 1; // node counts per direction
    corners.Reserve((nx_+1) * (ny_+1) * mz); // upper bound - the ip pointers in the items must stay valid
    for(unsigned int z = 0; z < mz; z++)
    {
      for(unsigned int y = 0; y <= ny_; y++)
      {
        for(unsigned int x = 0; x <= nx_; x++)
        {
          // the adjacent design cells (clamped ranges, in 2D only z=0)
          ItemIP* adj[8];
          int n_adj = 0;
          for(unsigned int az = (z > 0 ? z-1 : 0); az <= std::min(z, nz_-1); az++)
            for(unsigned int ay = (y > 0 ? y-1 : 0); ay <= std::min(y, ny_-1); ay++)
              for(unsigned int ax = (x > 0 ? x-1 : 0); ax <= std::min(x, nx_-1); ax++)
              {
                int it = lex_to_item_[LexicographicPos(ax, ay, az)];
                if(it >= 0)
                  adj[n_adj++] = GetItemIP(it);
              }
          if(n_adj == 0)
            continue; // node of a non-design part of the box

          corners.Push_back(IntegrationPoint());
          IntegrationPoint& ip = corners.Last();
          ip.loc.Set(origin_[0] + x*dx_, origin_[1] + y*dx_, dim_ == 3 ? origin_[2] + z*dx_ : 0.0);
          ip.dist.Resize(nf);
          ip.part.Resize(nf, FeatureVariable::NO_TIP);
          for(int a = 0; a < n_adj; a++)
            adj[a]->corner.Push_back(&ip);
        }
      }
    }
  }

  // optional debug logging
  for(unsigned int i = 0; i < corners.GetSize(); i++)
    LOG_DBG3(fm) << "SFIP: corners[" << i << "].loc=" << corners[i].loc.ToString();

  for(Item& item : map)
    LOG_DBG3(fm) << "SFIP: map " << item.lexicographic_pos << " ext: " << item.extension->ToString();
}

std::string FeatureMappingDesign::ItemIP::ToString()
{
  std::stringstream ss;
  
  ss << "c=" << center.ToString();
  
  ss << " corner #" << corner.GetSize() << "/" << corner.GetCapacity() << "->[";
  for(IntegrationPoint* ip : corner)
    ss << ip->loc.ToString() << " d=" << ip->dist.ToString() << ", ";
  ss << "]";
  return ss.str();
}

void FeatureMappingDesign::ItemIP::FillDist(Vector<double>& data, unsigned int num) const
{
  data.Resize(corner.GetSize() + inner.GetSize()); // no data copied
  for(unsigned int i = 0; i < corner.GetSize(); i++)
    data[i]=corner[i]->dist[num];
  for(unsigned int i = 0; i < inner.GetSize(); i++)
    data[corner.GetSize()+i]=inner[i]->dist[num];
}

void FeatureMappingDesign::ItemIP::GetAllIP(StdVector<IntegrationPoint*>& out)
{
  out.ResizeNoCopy(corner.GetSize() + inner.GetSize());
  for(unsigned int i = 0; i < corner.GetSize(); i++)
    out[i]=corner[i];
  for(unsigned int i = 0; i < inner.GetSize(); i++)
    out[corner.GetSize()+i]=inner[i];
}

StdVector<FeatureMappingDesign::IntegrationPoint*> FeatureMappingDesign::ItemIP::GetAllIP()
{
  StdVector<IntegrationPoint*> out;
  GetAllIP(out);
  return out;
}

void FeatureMappingDesign::CalcIpDistances(StdVector<IntegrationPoint>& ips, unsigned int n) const
{
  assert(n <= ips.GetSize());
  #pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(int i = 0; i < (int) n; i++)
  {
    IntegrationPoint& ip = ips[i];
    for(unsigned int f = 0; f < pills.GetSize(); f++)
      ip.dist[f] = pills[f].Distance(ip.loc, &(ip.part[f]));
  }
}

void FeatureMappingDesign::MapFeatureToDensity()
{
  mapping_timer_->Start();
  LOG_DBG(fm) << "MFTD called";

  assert(!map.IsEmpty());

  // the DesignElement might be updated, update the Feature internal representation
  for(Feature* f : features_)
    f->Update();
  
  // order 1: static barycenter ip, order >= 2: static corner ip at the element nodes,
  // order > 2: additionally the dynamic inner lattice points for elements in a transition zone
  if(order == 1)
    CalcIpDistances(inners_, num_inner_);
  else
  {
    CalcIpDistances(corners, corners.GetSize());

    if(order > 2)
    {
      // Dynamic inner integration points for the elements crossed by a transition zone. The points
      // live on the global refined lattice with q = order-1 subdivisions per element, so shared
      // points (on element faces and edges) have a unique global lattice index and no ownership or
      // sequential creation is needed. Four passes, all but the cheap compaction parallel:
      // 1) mark the lattice points of all active elements (element corners excluded, they are the
      //    static corner ip), 2) compact: assign consecutive inners_ slots, 3) fill the slots -
      //    reused slots keep their dist/part allocations, so this is allocation free after the
      //    first mapping - and compute the distances, 4) wire ItemIP::inner of the active elements
      unsigned int q = order - 1;
      unsigned int mx = q*nx_ + 1; // lattice points per direction
      unsigned int my = q*ny_ + 1;
      unsigned int mz = dim_ == 3 ? q*nz_ + 1 : 1;
      double spacing = dx_ / q;
      lattice_slot_.Resize(mx * my * mz);
      lattice_slot_.Init(0);
      item_active_.Resize(map.GetSize());

      // 1) transition zone check by the corner distances, then mark the element's lattice points
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for(int e = 0; e < (int) map.GetSize(); e++)
      {
        ItemIP* ext = GetItemIP(e);
        ext->inner.Clear(); // keeps the capacity for the re-wiring in 4)
        item_active_[e] = !bnd_fnc->IsConstant(ext);
        if(!item_active_[e])
          continue;
        unsigned int pos = map[e].lexicographic_pos;
        unsigned int bx = (pos % nx_) * q; // lattice base of the cell
        unsigned int by = ((pos / nx_) % ny_) * q;
        unsigned int bz = (pos / (nx_ * ny_)) * q;
        for(unsigned int k = 0; k < (dim_ == 3 ? order : 1); k++)
          for(unsigned int j = 0; j < order; j++)
            for(unsigned int i = 0; i < order; i++)
            {
              if(i % q == 0 && j % q == 0 && (dim_ == 2 || k % q == 0))
                continue; // element corner
              lattice_slot_[((bz + k)*my + (by + j))*mx + (bx + i)] = 1; // write races store the same value
            }
      }

      // 2) compact: consecutive slots (1-based, 0 stays unused)
      unsigned int cnt = 0;
      for(unsigned int g = 0; g < lattice_slot_.GetSize(); g++)
        if(lattice_slot_[g] != 0)
          lattice_slot_[g] = ++cnt;
      num_inner_ = cnt;
      if(inners_.GetSize() < cnt)
        inners_.Resize(cnt); // grow only: the ItemIP::inner pointers are re-wired below in 4)

      // 3) fill the slots
      unsigned int nf = pills.GetSize();
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for(long g = 0; g < (long) lattice_slot_.GetSize(); g++)
      {
        if(lattice_slot_[g] == 0)
          continue;
        IntegrationPoint& ip = inners_[lattice_slot_[g] - 1];
        unsigned int gx = g % mx, gy = (g / mx) % my, gz = g / (mx * my);
        ip.loc.Set(origin_[0] + gx*spacing, origin_[1] + gy*spacing, dim_ == 3 ? origin_[2] + gz*spacing : 0.0);
        ip.dist.Resize(nf); // no-op for a reused slot
        ip.part.Resize(nf); // set by CalcIpDistances
      }
      CalcIpDistances(inners_, num_inner_);

      // 4) wire the active elements to their slots, same lattice enumeration as the marking
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for(int e = 0; e < (int) map.GetSize(); e++)
      {
        if(!item_active_[e])
          continue;
        ItemIP* ext = GetItemIP(e);
        unsigned int pos = map[e].lexicographic_pos;
        unsigned int bx = (pos % nx_) * q;
        unsigned int by = ((pos / nx_) % ny_) * q;
        unsigned int bz = (pos / (nx_ * ny_)) * q;
        for(unsigned int k = 0; k < (dim_ == 3 ? order : 1); k++)
          for(unsigned int j = 0; j < order; j++)
            for(unsigned int i = 0; i < order; i++)
            {
              if(i % q == 0 && j % q == 0 && (dim_ == 2 || k % q == 0))
                continue; // element corner
              unsigned int slot = lattice_slot_[((bz + k)*my + (by + j))*mx + (bx + i)];
              assert(slot > 0);
              ext->inner.Push_back(&inners_[slot - 1]);
            }
      }

      LOG_DBG(fm) << "MFTD order=" << order << " inner ip=" << num_inner_;
    } // end of order > 2
  } // end of order >= 2

  LOG_DBG(fm) << "MFTD test boundary func: " << bnd_fnc->DebugLog();

  // all integration points set and their distances calculated, we can now integrate the density.
  // the elements are independent, the working arrays are thread private
  assert(num_var_by_feature == features_.First()->GetAllVariables().GetSize());
  #pragma omp parallel num_threads(CFS_NUM_THREADS)
  {
    Vector<double> all_dist; // all dist from ItemIP::corner and inner for a given feature
    Vector<double> bnd;      // all_dist mapped by boundary function
    // working arrays for CalcGradRhoByFeature
    Vector<double> ddist_ds;
    StdVector<IntegrationPoint*> elem_ip;
    Vector<double> drho_ds_vec;

    #pragma omp for
    for(int e = 0; e < (int) map.GetSize(); e++)
    {
      Item& item = map[e];
      ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
      LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " center=" << ext->center.ToString() << " corner=" << ext->corner.GetSize() << " inner=" << ext->inner.GetSize() << " -> " << IntegrationPoint::ToString(ext->inner);
      assert(order <= 2 || (ext->GetAllIP().GetSize() == std::pow(2,dim_) || ext->GetAllIP().GetSize() == std::pow(order,dim_))); // we have at least corners, but for order > 2 we have also inner points
      for(unsigned int fi = 0; fi < pills.GetSize(); fi++)
      {
        // collect all computed distances from all element integration points by feature number
        ext->FillDist(all_dist, fi);
        LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " f=" << fi << " d=" << all_dist.ToString();
        // project these distances by the boundary function (rhomin ... 1)
        bnd_fnc->Eval(all_dist, bnd);
        // integrate these projections to a single rho by feature -> before combine
        ext->rho[fi] = bnd_fnc->Integrate(bnd);
        LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " c=" << ext->center.ToString() << " f=" << fi << " prj:#" << bnd.GetSize()
                     << " prj=" << bnd.ToString() << "->" << bnd.ToString() << " -> " << ext->rho[fi];
        assert(ext->rho[fi] >= 0);

        // calculate d_rho/d_s to be re-used in the gradient calculation.
        CalcGradRhoByFeature(item, pills[fi], drho_ds_vec, ddist_ds, elem_ip); // drho_ds_vec is output, ddist_ds and elem_ip are working arrays
        // sort drho_drho_ds_vec in the full d_rho_ds vector
        for(unsigned int si = 0; si < num_var_by_feature; si++)
          ext->drho_ds_full[fi * num_var_by_feature + si] = drho_ds_vec[si];

        // the geometry variable alpha scales the feature density as rho_hat = alpha^q * rho (Norato's size
        // variable). We store rho_hat and its chain (alpha^q * drho/ds and drho_hat/dalpha = q*alpha^(q-1)*rho)
        // so all consumers (combine, gradient, Jacobian) work on rho_hat. Note the alpha column lives on the
        // whole feature footprint, not only on the transition zone.
        if(has_alpha_)
        {
          double a = pills[fi].alpha.GetPlainDesignValue();
          double aq = std::pow(a, alpha_q_);
          ext->rho_org[fi] = ext->rho[fi]; // keep the unscaled rho for the alpha Hessian entries
          ext->drho_ds_full[fi * num_var_by_feature + ALPHA_VAR] = alpha_q_ * std::pow(a, alpha_q_ - 1.0) * ext->rho[fi];
          for(unsigned int si = 0; si < ALPHA_VAR; si++)
            ext->drho_ds_full[fi * num_var_by_feature + si] *= aq;
          ext->rho[fi] *= aq;
        }
      }
      // combine for the element the rho of all features to the remaining element rho
      double rho = cmb_fnc->Eval(ext->rho);
      LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " c=" << ext->center.ToString() << " proj=" << ext->rho.ToString() << " -> rho=" << rho;
      item.elemval->SetDesign(rho);
    }
  } // end omp parallel

  LOG_DBG(fm)  << "MFTD: data=" << data.GetSize() << " first=" << data[0].ToString() << " simp=" << data[0].simp;
  LOG_DBG(fm)  << "MFTD: old mid=" << mapped_design_ << " current did=" << design_id;

  mapped_design_ = design_id;
  mapping_timer_->Stop();
}


void FeatureMappingDesign::MapFeatureGradient(Function* func)
{
  // for many density based functions we could cache stuff, as all is same but d_J/d_rho_e where J is func

  assert(mapped_design_ == design_id);
  gradient_timer_->Start();

  LOG_DBG(fm) << "MFG md= " << mapped_design_;

  // The sensitivity for isotropic and anisotropic features reads quite different.
  // In the isotropic case almost all are products from chain rule 
  // d_J/d_s = sum_e d_J/d_rho_e * d_rho_e/d_s 
  //           with d_rho_e/d_s = 1/N_ip sum_i d_H/d_d * d_d/d_s
  // In the anisotropic case the state dependent functions are like
  // d_J/d_s = sum_e <lmbda_e, dK_e/d_s u_e> (ErsatzMaterial::CalcU1KU2())
  //           with everything in dK_e/d_s
  // in the non-state dependent but rho dependent case (volume) it is like isotropic

  // we contribute to d_J/d_s and in the end move it to the opt variables (e.g. considering linking)
  // fixed is 0

  // assume all features are equal - this is not opt variables. We skip later fixed optionally link later
  assert(num_var_by_feature == features_.First()->GetAllVariables().GetSize()); 
  Vector<double> dJ_ds(features_.GetSize() * num_var_by_feature, 0.0);
  StdVector<FeatureVariable*> f_si;     // all variables of a feature

  // the isotropic or non-state dependent (volume) is purely chain rule and only factor sum_e dJ/drho * ...
  if(!anisotropic || !func->IsStateDependent())
    IsotropicGradientHelper(func, dJ_ds);
  else  // anisotropic state has dJ/drho and dJ/dphi and is computed by sum_e <lmbd_e, dK/ds u_e> via ParamMat::SetElementK()
    AnisotropicGradientHelper(func, dJ_ds);

  // copy from dJ_ds to the opt variables - resolve mapping, skip constants (see FeatureVariable::EffectiveOptIndex)
  for(Feature& pill : pills)
  {
    pill.GetAllVariables(f_si);
    assert(f_si.GetSize() == num_var_by_feature);
    for(unsigned int i = 0; i < f_si.GetSize(); i++)
    {
      int oi = f_si[i]->EffectiveOptIndex();
      if(oi >= 0)
        opt_shape_param_[oi]->AddGradient(func, dJ_ds[pill.idx * num_var_by_feature + i]);
    }
  }

  assert(dJ_ds.GetSize() == features_.GetSize() * num_var_by_feature);
  assert(dJ_ds.GetSize() == shape_param_.GetSize()); 
  LOG_DBG(fm) << "MFG: func=" << func->ToString() << " dJ_ds=" << dJ_ds.GetSize() << " sp=" << shape_param_.GetSize()  << " osp=" << opt_shape_param_.GetSize();
  for(unsigned int i = 0; i < dJ_ds.GetSize(); i++)
    LOG_DBG2(fm) << "MFG: func=" << func->ToString() << " dJ_ds[" << i << "]=" << dJ_ds[i] << " var=" << dynamic_cast<FeatureVariable*>(shape_param_[i])->ToString(); 
  for(unsigned int i = 0; i < opt_shape_param_.GetSize(); i++)
    LOG_DBG2(fm) << "MFG: func=" << func->ToString() << " opt[" << i << "]=" << opt_shape_param_[i]->GetPlainGradient(func) << " var=" << dynamic_cast<FeatureVariable*>(shape_param_[i])->ToString(); 

  gradient_timer_->Stop();
}

void FeatureMappingDesign::IsotropicGradientHelper(const Function* func, Vector<double>& dJ_ds) const
{
  // for many density based functions we could cache stuff, as all is same but d_J/d_rho_e where J is func
  
  // the expensive part, the distance calculation is already done, therefore we parallelize only by feature, not by variable

  // this is the case for max_rho_e via p-norm: rho_e = (sum_f rho_f^p)^(1/p)
  // d_max_rho_r/d_s = (sum_f rho_f^p)^(1/p -1) * rho_f(s)^(p-1) * d_rho_f(s)/d_s 
  // d_J/d_s = sum_e d_J/d_rho_e * (sum_f rho_f^p)^(1/p -1) * rho_f(s)^(p-1) * 1/N_ip sum_i d_H/d_d * d_d/d_s

  // we contribute to d_J/d_s and in the end move it to the opt variables (e.g. considering linking)
  // fixed is 0

  // assume all features are equal - this is not opt variables. We skip later fixed optionally link later
  assert(num_var_by_feature == features_.First()->GetAllVariables().GetSize()); 
  // we use the precomputed ext->drho_ds_full

  StdVector<int> dJ_ds_res_idx = GetSpecialResultIndices(func);
  StdVector<int> dmrho_ds_res_idx = GetSpecialResultIndices(nullptr, DE::COMBINE, true); // d_mrho_d_s
  StdVector<int> drho_ds_res_idx = GetSpecialResultIndices(nullptr, DE::FEATURE_RHO, true); // d_rho_d_s
  assert(dJ_ds_res_idx.GetSize() == pills.GetSize() * num_var_by_feature && dJ_ds_res_idx.GetSize() == dmrho_ds_res_idx.GetSize() &&  dJ_ds_res_idx.GetSize() == drho_ds_res_idx.GetSize());
  
  StdVector<int> dmrho_drho_res_idx = GetSpecialResultIndices(nullptr, DE::COMBINE, false); // d_mrho_d_feature
  assert(dmrho_drho_res_idx.GetSize() == pills.GetSize()); 

  // for each element
  for(const Item& item : map)
  {
    const ItemIP* ext = dynamic_cast<const ItemIP*>(item.extension);
    // SIMP gradient: function by rho_e
    double dJ_drho_e = item.elemval->GetPlainGradient(func);

    // for each feature
    for(unsigned int f = 0; f < pills.GetSize(); f++)
    {
      // skip when the whole d_rho/d_z row of this feature vanishes (fully solid/void element, gradient is
      // zero) but stay if we have the special result. Without alpha this is equivalent to the former
      // rho == 1 || rho == rhomin check; with alpha the interior elements stay alive via the alpha column.
      bool zero_row = true;
      for(unsigned int s = 0; s < num_var_by_feature && zero_row; s++)
        zero_row = ext->drho_ds_full[f * num_var_by_feature + s] == 0.0;

      if(zero_row && dmrho_drho_res_idx[f] == -1)
        continue;

      // derivative of rho_f aggregation of all features (1.0 single feature)
      double dmax_drho_f = cmb_fnc->Grad(ext->rho, f);

      if(dmrho_drho_res_idx[f] != -1)
      {
        item.elemval->specialResult[dmrho_drho_res_idx[f]] = dmax_drho_f;
        LOG_DBG3(fm) << "IGH f=" << f << " dmax_drho_f=" << dmax_drho_f << " -> s_r[" << dmrho_drho_res_idx[f] << "]";
        if(zero_row) // in that case drho_ds is zero and all results contain this
          continue;
      }

      assert(dJ_ds.GetSize() == ext->drho_ds_full.GetSize());
      for(unsigned int s = 0; s < num_var_by_feature; s++) 
      {
        unsigned int idx = f * num_var_by_feature + s;
        double value = dJ_drho_e * dmax_drho_f * ext->drho_ds_full[idx]; 
        dJ_ds[idx] += value;

        if(drho_ds_res_idx[idx] != -1)
          item.elemval->specialResult[drho_ds_res_idx[idx]] = ext->drho_ds_full[idx]; // we automatically have the proper feature

        if(dmrho_ds_res_idx[idx] != -1)
          item.elemval->specialResult[dmrho_ds_res_idx[idx]] = dmax_drho_f * ext->drho_ds_full[idx];

        if(dJ_ds_res_idx[idx] != -1)
          item.elemval->specialResult[dJ_ds_res_idx[idx]] = value;

        LOG_DBG3(fm) << "IGH: item=" << item.lexicographic_pos << " f=" << f << " s=" << s 
                   << " dJ_drho_e=" << dJ_drho_e << " dmax_drho_f=" << dmax_drho_f << " drho_f_ds=" << ext->drho_ds_full[idx] << " dJ_ds[" << idx << "] += " << value <<  " -> " <<  dJ_ds[idx];
      } 
    }
  }
}

void FeatureMappingDesign::CalcShapeJacobian(Matrix<double>& D) const
{
  // Shape Jacobian D = d_mrho/d_s: derivative of each element's aggregated density mrho_e w.r.t.
  // all shape variables s (Px, Py, Qx, Qy, p per feature), full variable space (fixed vars as zero
  // columns). Chain rule per element e and feature f:
  //   D[e][f,s] = (d_mrho_e/d_rho_f) * (d_rho_f/d_s)
  //   - d_mrho_e/d_rho_f : aggregation sensitivity (cmb_fnc->Grad, couples overlapping features)
  //   - d_rho_f/d_s      : geometric, precomputed in ItemIP::drho_ds_full (with alpha incl. its column)
  // Rows with a vanishing chain (fully solid/void without alpha) are skipped (D sparse).
  // D carries the element coupling, so the objective Hessian term is D^T diag(curv) D.

  unsigned int N_full = pills.GetSize() * num_var_by_feature;
  D.Resize(map.GetSize(), N_full);
  D.Init();

  for(unsigned int e = 0; e < map.GetSize(); e++)
  {
    const ItemIP* ext = dynamic_cast<const ItemIP*>(map[e].extension);
    for(unsigned int f = 0; f < pills.GetSize(); f++)
    {
      // skip when the whole chain row vanishes (fully solid/void feature without alpha; with alpha the
      // alpha column keeps the footprint alive), as in IsotropicGradientHelper()
      bool zero_row = true;
      for(unsigned int s = 0; s < num_var_by_feature && zero_row; s++)
        zero_row = ext->drho_ds_full[f * num_var_by_feature + s] == 0.0;
      if(zero_row)
        continue;

      double dmrho_drho_f = cmb_fnc->Grad(ext->rho, f);
      for(unsigned int s = 0; s < num_var_by_feature; s++)
        D[e][f * num_var_by_feature + s] = dmrho_drho_f * ext->drho_ds_full[f * num_var_by_feature + s];
    }
  }
}

void FeatureMappingDesign::AssembleHessianTerms(const Function* func, Matrix<double>& H_agg, Matrix<double>& H_feat) const
{
  unsigned int nv = num_var_by_feature;
  unsigned int N_full = pills.GetSize() * nv;
  int N_ip = std::pow(order, dim_); // normalization as in CalcGradRhoByFeature()

  H_agg.Resize(N_full, N_full);
  H_agg.Init();
  H_feat.Resize(N_full, N_full);
  H_feat.Init();

  // working data
  Vector<double> ddist_ds(nv);
  Matrix<double> hd; // nv x nv distance Hessian per integration point
  Matrix<double> rh(nv, nv); // d^2 rho_e^f/(d_s_i d_s_j), eqn 'rho_hess' in the tracking paper
  Vector<double> rho_s(ALPHA_VAR); // unscaled d_rho_f/d_s for the mixed alpha-s Hessian entries
  StdVector<IntegrationPoint*> elem_ip;
  StdVector<unsigned int> active;
  active.Reserve(pills.GetSize());

  // per-feature alpha factors: rho_hat = a^q rho, first and second derivative of a^q
  Vector<double> aq(pills.GetSize(), 1.0), daq(pills.GetSize(), 0.0), ddaq(pills.GetSize(), 0.0);
  if(has_alpha_)
    for(unsigned int f = 0; f < pills.GetSize(); f++)
    {
      double a = pills[f].alpha.GetPlainDesignValue();
      aq[f]   = std::pow(a, alpha_q_);
      daq[f]  = alpha_q_ * std::pow(a, alpha_q_ - 1.0);
      // q(q-1) a^(q-2); exactly zero for q=1, and we take the (one-sided) limit 0 at a=0 for q < 2
      ddaq[f] = alpha_q_ == 1.0 || (a == 0.0 && alpha_q_ < 2.0) ? 0.0 : alpha_q_ * (alpha_q_ - 1.0) * std::pow(a, alpha_q_ - 2.0);
    }

  for(const Item& item : map)
  {
    double w_e = item.elemval->GetPlainGradient(func); // d_func/d_mrho_e
    if(w_e == 0.0)
      continue;

    ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);

    // the features contributing in this element: any non-zero chain entry (transition zone; with alpha
    // also the footprint via the alpha column) or a non-zero pure alpha curvature on the footprint
    active.Resize(0); // keeps capacity
    for(unsigned int f = 0; f < pills.GetSize(); f++)
    {
      bool zero_row = true;
      for(unsigned int s = 0; s < nv && zero_row; s++)
        zero_row = ext->drho_ds_full[f*nv + s] == 0.0;
      if(!zero_row || (has_alpha_ && ddaq[f] != 0.0 && ext->rho_org[f] != 0.0))
        active.Push_back(f);
    }
    if(active.IsEmpty())
      continue;

    // aggregation term: w_e * d^2 mrho/(d_rho_f d_rho_g) * d_rho_f/d_s_i * d_rho_g/d_s_j.
    // this is the only cross-feature coupling, present where features share elements
    for(unsigned int f : active)
      for(unsigned int g : active)
      {
        double agg = cmb_fnc->Hessian(ext->rho, f, g);
        for(unsigned int i = 0; i < nv; i++)
          for(unsigned int j = 0; j < nv; j++)
            H_agg[f*nv + i][g*nv + j] += w_e * agg * ext->drho_ds_full[f*nv + i] * ext->drho_ds_full[g*nv + j];
      }

    // feature term: w_e * d_mrho/d_rho_hat_f * d^2 rho_hat_f/(d_z_i d_z_j) with the geometric part
    // d^2 rho_f/(d_s_i d_s_j) = 1/N_ip * sum_ip (H''*dd_i*dd_j + H'*d^2d_ij)
    // and with alpha: rho_hat = a^q rho, so the s-s block scales by a^q, the mixed alpha-s entries are
    // (a^q)' * d_rho/d_s (first-order data) and the alpha-alpha entry is (a^q)'' * rho (eqn
    // 'geom_var_hess' in the tracking paper)
    ext->GetAllIP(elem_ip);
    for(unsigned int f : active)
    {
      rh.Init();
      rho_s.Init();
      for(const IntegrationPoint* ip : elem_ip)
      {
        double dist = ip->dist[f];
        if(bnd_fnc->IsConstant(dist))
          continue; // H' and H'' are zero outside the transition zone

        double Hp  = bnd_fnc->Grad(dist);
        double Hpp = bnd_fnc->Hessian(dist); // throws for poly which is only C^1
        pills[f].GradDistance(ip->loc, ip->part[f], ddist_ds);
        pills[f].HessDistance(ip->loc, ip->part[f], hd); // 5x5, the geometry variables only
        for(unsigned int i = 0; i < ALPHA_VAR; i++)
          for(unsigned int j = 0; j < ALPHA_VAR; j++)
            rh[i][j] += (Hpp * ddist_ds[i] * ddist_ds[j] + Hp * hd[i][j]) / N_ip;
        if(has_alpha_)
          for(unsigned int j = 0; j < ALPHA_VAR; j++)
            rho_s[j] += Hp * ddist_ds[j] / N_ip; // unscaled d_rho_f/d_s_j
      }

      if(has_alpha_)
      {
        for(unsigned int i = 0; i < ALPHA_VAR; i++)
        {
          for(unsigned int j = 0; j < ALPHA_VAR; j++)
            rh[i][j] *= aq[f];
          rh[ALPHA_VAR][i] = rh[i][ALPHA_VAR] = daq[f] * rho_s[i];
        }
        rh[ALPHA_VAR][ALPHA_VAR] = ddaq[f] * ext->rho_org[f];
      }

      double dmrho_drho_f = cmb_fnc->Grad(ext->rho, f);
      for(unsigned int i = 0; i < nv; i++)
        for(unsigned int j = 0; j < nv; j++)
          H_feat[f*nv + i][f*nv + j] += w_e * dmrho_drho_f * rh[i][j];
    }
  }
}

std::string FeatureMappingDesign::HessExportTermXML(const std::string& name, unsigned int nv,
    const Matrix<double>& H, bool diagonal_only) const
{
  unsigned int Nf = pills.GetSize();
  // slice the (f,g) 5x5 sub-block of the full-space term matrix, return true if it has content
  Matrix<double> sub(nv, nv);
  auto slice = [&](unsigned int f, unsigned int g) -> bool {
    bool nonzero = false;
    for(unsigned int i = 0; i < nv; i++)
      for(unsigned int j = 0; j < nv; j++)
      {
        sub[i][j] = H[f*nv + i][g*nv + j];
        nonzero |= sub[i][j] != 0.0;
      }
    return nonzero;
  };

  // term-major: <name> wraps the non-empty (f,g) blocks as children. The bulk content is verbatim and
  // only the first line (<name>) is auto-indented by the printer (to 6, under <function>); place the
  // blocks at 8 (data at 10, typeTag=false omits the <real> wrapper) and </name> at 6 by hand. Only
  // g >= f (lower triangle by symmetry); the feature/geometry term has diagonal (f==f) blocks only.
  std::stringstream ss;
  for(unsigned int f = 0; f < Nf; f++)
    for(unsigned int g = f; g < Nf; g++)
    {
      if(diagonal_only && g != f)
        continue;
      if(!slice(f, g))
        continue;
      std::string attr = "f=\"" + std::to_string(f) + "\" g=\"" + std::to_string(g) + "\"";
      ss << "        " << sub.ToXMLFormat("block", 8, false, attr) << "\n";
    }

  if(ss.str().empty())
    return "<" + name + "/>"; // term structurally zero (e.g. curvature of a function linear in rho)

  std::stringstream term;
  term << "<" << name << ">\n" << ss.str() << "      </" << name << ">";
  return term.str();
}

bool FeatureMappingDesign::CalcShapeHessianTerms(const Function* func, Matrix<double>& H_obj, Matrix<double>& H_agg, Matrix<double>& H_feat)
{
  unsigned int N_full = pills.GetSize() * num_var_by_feature;

  Vector<double> curv; // per-element diagonal d^2 func/d_rho_e^2
  if(!const_cast<Function*>(func)->CalcCurvature(curv))
    return false; // no curvature information for this function

  Matrix<double> D; // the Jacobian d_mrho/d_s
  CalcShapeJacobian(D);
  assert(curv.GetSize() == D.GetNumRows());

  // objective term D^T diag(curv) D, see eqn 'hessian_terms' in the tracking paper
  H_obj.Resize(N_full, N_full);
  H_obj.Init();
  for(unsigned int e = 0; e < D.GetNumRows(); e++)
  {
    if(curv[e] == 0.0)
      continue;
    for(unsigned int i = 0; i < N_full; i++)
    {
      if(D[e][i] == 0.0)
        continue; // D is sparse, only transition zone elements contribute
      double cD = curv[e] * D[e][i];
      for(unsigned int j = 0; j < N_full; j++)
        H_obj[i][j] += cD * D[e][j];
    }
  }

  AssembleHessianTerms(func, H_agg, H_feat);
  return true;
}

void FeatureMappingDesign::BuildOptIndexMap(StdVector<int>& opt_idx) const
{
  opt_idx.Resize(pills.GetSize() * num_var_by_feature);
  StdVector<FeatureVariable*> f_si;
  for(const Pill& pill : pills)
  {
    pill.GetAllVariables(f_si);
    assert(f_si.GetSize() == num_var_by_feature);
    for(unsigned int i = 0; i < f_si.GetSize(); i++)
      opt_idx[pill.idx * num_var_by_feature + i] = f_si[i]->EffectiveOptIndex();
  }
}

bool FeatureMappingDesign::CalcShapeHessian(const Function* func, Matrix<double>& H_opt)
{
  if(anisotropic)
    throw Exception("the exact shape Hessian is only implemented for the isotropic case");
  assert(mapped_design_ == design_id);

  Matrix<double> H_obj, H_agg, H_feat;
  if(!CalcShapeHessianTerms(func, H_obj, H_agg, H_feat))
    return false;

  // map the full feature space total (sum of the three terms) to optimization variable space:
  // linked variables accumulate on their target, fixed variables (-1) are dropped
  StdVector<int> opt_idx;
  BuildOptIndexMap(opt_idx);
  unsigned int N = opt_shape_param_.GetSize();
  unsigned int N_full = H_obj.GetNumRows();
  H_opt.Resize(N, N);
  H_opt.Init();
  for(unsigned int i = 0; i < N_full; i++)
  {
    if(opt_idx[i] < 0)
      continue;
    for(unsigned int j = 0; j < N_full; j++)
      if(opt_idx[j] >= 0)
        H_opt[opt_idx[i]][opt_idx[j]] += H_obj[i][j] + H_agg[i][j] + H_feat[i][j];
  }
  return true;
}

void FeatureMappingDesign::WriteHessExportFile()
{
  if(!hessexport_)
    return;

  if(anisotropic)
    throw Exception("the hessexport option is only implemented for the isotropic case");

  assert(mapped_design_ == design_id);

  unsigned int nv = num_var_by_feature;
  unsigned int Nf = pills.GetSize();

  // persistent root, created once: each call appends an <iteration> child so the file keeps the whole
  // history. Written (re-written with the grown tree) every iteration via ParamNode::ToFile(), with
  // self-written xml strings for the blocks.
  if(!hessexport_root_)
  {
    hessexport_root_ = PtrParamNode(new ParamNode(ParamNode::INSERT));
    hessexport_root_->SetName("hessian");
    hessexport_root_->Get("features")->SetValue(Nf);
    hessexport_root_->Get("varsPerFeature")->SetValue(nv);
    // the local variable order within a feature block, matches GradDistance()/HessDistance()
    std::string vars = dim_ == 2 ? "Px Py Qx Qy P" : "Px Py Pz Qx Qy Qz P";
    hessexport_root_->Get("variables")->SetValue(vars + (has_alpha_ ? " alpha" : ""));
  }

  PtrParamNode iteration = hessexport_root_->Get("iteration", ParamNode::APPEND);
  iteration->Get("nr")->SetValue(opt->GetCurrentIteration());

  // the design variables [Px Py Qx Qy P] per feature (one row per feature)
  Matrix<double> V(Nf, nv);
  StdVector<FeatureVariable*> fv;
  for(unsigned int f = 0; f < Nf; f++)
  {
    pills[f].GetAllVariables(fv);
    for(unsigned int i = 0; i < nv; i++)
      V[f][i] = fv[i]->GetPlainDesignValue();
  }
  iteration->Get("variables")->SetValue(V);

  // a function without a second derivative in the export gets a one-line <noHessian name=".."/> marker
  // so every function shows up (e.g. the C++ 'volume' constraint, or a local constraint not yet dispatched)
  auto noHessian = [&](const Function* f) {
    iteration->Get("noHessian", ParamNode::APPEND)->Get("name")->SetValue(f->ToString());
  };

  Matrix<double> H_obj, H_agg, H_feat; // the three Hessian terms from the paper. The sum is the real Heassian
  for(Function* func : opt->GetFunctions(false))
  {
    if(func->IsLocal())
    {
      if(func->GetType() != Function::DISTANCE)
      {
        noHessian(func);
        continue;
      }

      // local constraints up to now live only for the feature part and have a compact form
      // <function name="distance_(node)_lowerBound" local="0">
      //      <feature>
      //        <block f="0" g="0" dim1="5" dim2="5">
      //           2.490270e-01 -5.600441e-01 -2.490270e-01  5.600441e-01  0.000000e+00
      Matrix<double> hl;
      for(unsigned int f = 0; f < Nf; f++)
      {
        pills[f].HessLength(hl); 
        std::string attr = "f=\"" + std::to_string(f) + "\" g=\"" + std::to_string(f) + "\"";
        std::stringstream feature; // term-major: a local constraint has only the <feature> geometry term
        feature << "<feature>\n        " << hl.ToXMLFormat("block", 8, false, attr) << "\n      </feature>";
        PtrParamNode fn = iteration->Get("function", ParamNode::APPEND);
        fn->Get("name")->SetValue(func->ToString());
        fn->Get("local")->SetValue(f);
        fn->GetFastBulkBlock().Push_back(feature.str());
      }
      continue;
    }
    if(!CalcShapeHessianTerms(func, H_obj, H_agg, H_feat))
    {
      noHessian(func); // no curvature information for this function (e.g. C++ volume)
      continue;
    }

    // term-major: one element per Hessian term, each wrapping its (f,g) blocks as children. An empty
    // term becomes a self-closing marker, e.g. <curvature/> for a function linear in rho (reward).
    StdVector<std::string> terms;
    terms.Push_back(HessExportTermXML("curvature",   nv, H_obj,  false));
    terms.Push_back(HessExportTermXML("aggregation", nv, H_agg,  false));
    terms.Push_back(HessExportTermXML("feature",     nv, H_feat, true));

    PtrParamNode fn = iteration->Get("function", ParamNode::APPEND);
    fn->Get("name")->SetValue(func->ToString());
    fn->GetFastBulkBlock() = terms;
  }

  hessexport_root_->ToFile(progOpts->GetSimName() + ".hessian.xml", true); // grows every iteration
}

void FeatureMappingDesign::AnisotropicGradientHelper(Function* func, Vector<double>& dJ_ds)
{
  // for a single anisotropic feature we have
  // C_e= rho_e * R*C0*R^T with R=R(phi(s)) and rho_e = rho_e(H(d(s))) as in the isotropic case
  // for a state dependent function we have
  // d_J/d_s = sum_e <lmbda_e, dK_e/d_s u_e> 
  // where dK_e/d_s needs d_C_e/d_e = d(R*C0*R^T)/d_phi*d_phi/d_s rho_e + R^T C0 R d_rho_e/d_s
  // the implementation is quite different from the isotropic case as we are based on CalcU1KU2()
  // for every variable we execute the compliance, .... gradient calculation <lmbda_e, dK_e/d_s u_e> and sum it up
  // see FeatueMappingParamMat::SetElementK() which eventually calls GetAnisoMechTensor() 
  assert(num_var_by_feature == features_.First()->GetAllVariables().GetSize()); 
  assert(opt != nullptr);
  assert(func->IsStateDependent());

  // of size num_vars_by_feature * num_features
  aniso_current_variable.dJ_ds_res_idx = GetSpecialResultIndices(func);
  aniso_current_variable.dmrho_ds_res_idx = GetSpecialResultIndices(nullptr, DE::COMBINE, true); // d_mrho_d_s
  aniso_current_variable.drho_ds_res_idx = GetSpecialResultIndices(nullptr, DE::FEATURE_RHO, true); // d_rho_d_s
  // of size num_features
  aniso_current_variable.dmrho_drho_res_idx = GetSpecialResultIndices(nullptr, DE::COMBINE, false); // d_mrho_d_rho within feature

  dJ_ds.Init();
  assert(dJ_ds.GetSize() == features_.GetSize() * num_var_by_feature);

  for(unsigned int f = 0; f < pills.GetSize(); f++)
  {
    aniso_current_variable.feature_idx = f;
    for(unsigned int s = 0; s < num_var_by_feature; s++) 
    {
      assert(s <= (DE::FEATURE_MAPPING_P - DE::FEATURE_MAPPING_PX));
      aniso_current_variable.var_type = (DE::Type) (DE::FEATURE_MAPPING_PX + s);
      aniso_current_variable.var_idx = (int) s;
      
      // this speeds up the derivate for GetAnisoMechTensor() a lot
      aniso_current_variable.RCdR = aniso_base_tensor; // copy the unrotated anisotropic tensor
      MaterialTensor<double> RCdR(VOIGT, &(aniso_current_variable.RCdR), false); // wrap the MaterialTensor around the Matrix
      Optimization::context->dm->RotateTensor(RCdR, DE::ROTANGLE, true, pills[f].angle); // true means we provide the angle
      LOG_DBG2(fm) << "AGH: func=" << func->ToString() << " f=" << f << " s=" << s << " angle=" << pills[f].angle 
                   << " RCdR=" << RCdR.GetMatrix(VOIGT).ToString() << " acv.RCdR=" <<  aniso_current_variable.RCdR.ToString();
      assert(aniso_current_variable.RCdR.NormL2() > 0);

      for(Excitation* excite : func->ctxt->excitations)
      {
        // some objectives are only to be evaluated for the last excitation
        if(!func->DoEvaluate(excite))
          continue;
        excite->Apply(true); // set the correct context
        
        // shall only affect DesignSpace::data -> density elements
        ResetGradient(func);
          
        // calls FeatureMappingParamMat::SetElementK() for every element
        opt->CalcFunction(*excite, func, true); 

        // sum up all gradients -> they will later be handled for mapping
        for(const Item& item : map) 
        {
          dJ_ds[f * num_var_by_feature + s] += item.elemval->GetPlainGradient(func);
          LOG_DBG3(fm) << "AGH: func=" << func->ToString() << " f=" << f << " s=" << s << " + " << item.elemval->GetPlainGradient(func) << " -> dJ_s=" <<  dJ_ds[f * num_var_by_feature + s];
          if(aniso_current_variable.dJ_ds_res_idx[f * num_var_by_feature + s] != -1)
             item.elemval->specialResult[aniso_current_variable.dJ_ds_res_idx[f * num_var_by_feature + s]] = item.elemval->GetPlainGradient(func);
        }
        LOG_DBG2(fm) << "AGH: func=" << func->ToString() << " f=" << f << " s=" << s << " -> dJ_s=" <<  dJ_ds[f * num_var_by_feature + s];
      }
    } 
  }
  aniso_current_variable.Reset(); // The optimizer will compute the compliance gradient by rho -> set for FeatureMappingParamMat::SetElementK()
}

inline void FeatureMappingDesign::CalcGradRhoByFeature(const Item& item, const Feature& feature, Vector<double>& out, Vector<double>& ddist_ds, StdVector<IntegrationPoint*>& elem_ip) const
{
  int N_ip = std::pow(order,dim_); // number of integration points 
  
  // the expensive part, the distance calculation is already done

  // d_rho/d_s = 1/N_ip * sum_i d_H/d_d * d_d/d_s

  // assume all features are equal - this is not opt variables. We skip later fixed optionally link later
  assert(num_var_by_feature == features_.First()->GetAllVariables().GetSize()); 
  out.Resize(num_var_by_feature);
  out.Init(0.0); // we sum up later

  ddist_ds.Resize(num_var_by_feature);
  ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
  assert(feature.idx >= 0 && feature.idx <= (int) features_.GetSize());
  unsigned int f = (unsigned int) feature.idx;

  // for all integration points of the element
  ext->GetAllIP(elem_ip);
  LOG_DBG3(fm) << "CGRBF: item=" << item.lexicographic_pos << " f=" << f << " inner=" << ext->inner.GetSize() << " elem_ip=" << elem_ip.GetSize() << " num_var_by_feature=" << num_var_by_feature;
  assert(elem_ip.GetSize() == std::pow(2,dim_) || (int) elem_ip.GetSize() == N_ip); // corners only or full order lattice
  for(IntegrationPoint* ip : elem_ip)
  {
    // now we have the actual point X and distance for each IP know which part is closest

    // derivative of boundary function (poly) by distance
    double dH_ddist = bnd_fnc->Grad(ip->dist[f]);

    // for all variables of the specific feature (P_x, P_y, Q_x, Q_y and p at once)
    // this is the most tricky part in the derivatives
    feature.GradDistance(ip->loc, ip->part[f], ddist_ds); // output is the vector ddist_ds (including fixed)

    // vectorial operation. Contains also values for fixed variables!
    out.Add(1./N_ip * dH_ddist, ddist_ds);
  }
}


int FeatureMappingDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  // practically updates the optimization Variable objects within features. Set's design_id
  FeaturedDesign::ReadDesignFromExtern(space_in, setAndWriteCurrent); 

  // all variables are set, we can now perform the optional mapping to mapped variables which are not within opt_shape_param_
  for(ShapeParamElement* spe : shape_param_)
  {
    FeatureVariable* var = dynamic_cast<FeatureVariable*>(spe);
    assert(var != nullptr); // all shape_param_ are FeatureVariables  

    FeatureVariable* org = var->map_to;
    if(org) // yes, we map
    {
      assert(var->map != "" && org->key != "" && var->map == org->key);
      var->SetDesign(org->GetPlainDesignValue()); // here we could apply advanced mapping like -key, 1-key, key+5, etc.
    }
  }

  LOG_DBG(fm) << "RDFE mid=" << mapped_design_ << " did=" << design_id;

  if(mapped_design_ != design_id)
    MapFeatureToDensity();

  return design_id;
}

void FeatureMappingDesign::ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation)
{
  // we expect all nodes (start x, start y, end x, end y, ...) then all profiles then all normals (if any).
  // by the shape attribute we could alll collect from random order but for easy we expect this ordering
  assert(features_.GetSize() > 0);

  ParamNodeList list = set->GetList("shapeParamElement");
  if(list.GetSize() != shape_param_.GetSize())
    EXCEPTION("expect " << shape_param_.GetSize() << " 'shapeParamElement' in density.xml but got " << list.GetSize());

  for(unsigned int i = 0; i < shape_param_.GetSize(); i++)
  {
    FeatureVariable*    var = dynamic_cast<FeatureVariable*>(shape_param_[i]);
    PtrParamNode pn  = list[i];

    // do a lot of validation
    const std::string pnnr = pn->Get("nr")->As<string>();
    if(pn->Get("nr")->As<unsigned int>() != i)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'nr' value " << i);

    if(pn->Get("shape")->As<int>() != var->feature)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'shape' value " << var->feature);

    if(DesignElement::type.Parse(pn->Get("type")->As<string>()) != var->GetType())
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'type' value " << DE::type.ToString(var->GetType()));

    if(var->GetType() == FeatureVariable::NODE)
    {
      if(ShapeParamElement::dof.Parse(pn->Get("dof")->As<string>()) != var->dof_)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'dof' value " << ShapeParamElement::dof.ToString(var->dof_));
      if(FeatureVariable::tip_enum.Parse(pn->Get("tip")->As<string>()) != var->tip)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'tip' value " << FeatureVariable::tip_enum.ToString(var->tip));
    }

    var->SetDesign(pn->Get("design")->As<double>());
  }

  // it shall be save to update the design_id, because we changed stuff
  design_id++;

  MapFeatureToDensity();
}


void FeatureMappingDesign::AddToDensityHeader(PtrParamNode header)
{
  PtrParamNode pn = header->Get("featureMapping");
  pn->Get("transition")->SetValue(transition);
  pn->Get("extension")->SetValue(extension); // only relevant for bezier
  if(has_alpha_)
    pn->Get("alpha_q")->SetValue(alpha_q_); // featureviz.py scales the drawn pills by alpha^q
}


FeatureVariable* FeatureMappingDesign::FindKey(const std::string& key) const
{
  for(const Feature* f : features_)
    for(FeatureVariable* var : f->GetAllVariables()) // important that this variable is not subject to a resizable vector
      if(var->key == key)
        return var;
  return nullptr; // not found
}

int FeatureMappingDesign::CountKey(const std::string& key) const
{
  int count = 0;
  for(const Feature* f : features_)
    for(FeatureVariable* var : f->GetAllVariables())
      if(var->key == key)
        count++;
  return count;
}

void FeatureMappingDesign::SetupParsedFeatures(PtrParamNode base)
{
  // external feature files (root element 'features', validated against the same schema so the
  // defaults apply) and local pills - the schema fixes the order external*, pill*, so the external
  // features come first. E.g. an external base skeleton with a varying number of features plus
  // special local pills
  StdVector<PtrParamNode> pnl;
  for(PtrParamNode child : base->GetChildren())
  {
    if(child->GetName() == "pill")
      pnl.Push_back(child);
    if(child->GetName() == "external")
    {
      std::string file = child->Get("file")->As<std::string>();
      std::string schema = progOpts->GetSchemaPathStr() + "/CFS-Simulation/CFS.xsd";
      PtrParamNode root = XmlReader::ParseFile(file, schema, "http://www.cfs++.org/simulation");
      if(root->GetName() != "features")
        throw Exception("root element of external feature file '" + file + "' needs to be 'features'");
      for(PtrParamNode pill : root->GetList("pill"))
        pnl.Push_back(pill);
    }
  }
  if(pnl.IsEmpty()) // the schema cannot express at least one pill from either source
    throw Exception("no features given, neither as local pills nor via external files");
  pills.Resize(pnl.GetSize());
  features_.Resize(pnl.GetSize());
  for(unsigned int i = 0; i < pnl.GetSize(); i++)
  {
    pills[i].Parse(pnl[i], i);
    // the global <alpha> gives the defaults, an optional <alpha> in the pill overrides (e.g. fixed or map)
    if(has_alpha_)
      pills[i].ParseAlpha(pnl[i]->Has("alpha") ? pnl[i]->Get("alpha") : base->Get("alpha"), i);
    features_[i] = &pills[i];
  }
}

 std::string FeatureMappingDesign::IntegrationPoint::ToString(StdVector<IntegrationPoint*>& ips)
{
  std::stringstream ss;
  for(const IntegrationPoint* ip : ips)
    ss << ip->loc.ToString() + " ";
  return ss.str();
}

void FeatureMappingDesign::Pill::ParseAlpha(PtrParamNode pn, int idx)
{
  alpha.Parse(pn, idx, DesignElement::FEATURE_MAPPING_ALPHA);
  has_alpha = true;
  opt_variables_ += alpha.IsVariable() ? 1 : 0;
}

void FeatureMappingDesign::Pill::GetAllVariables(StdVector<FeatureVariable*>& out) const
{
  // as the base but with the optional alpha appended after the profile, reserve to avoid expansion
  out.Clear();
  out.Reserve(GetTotalVariables());
  for(const StdVector<FeatureVariable>& vec : points)
    for(const FeatureVariable& var : vec)
      out.Push_back(const_cast<FeatureVariable*>(&var));
  out.Push_back(const_cast<FeatureVariable*>(&p));
  if(has_alpha)
    out.Push_back(const_cast<FeatureVariable*>(&alpha));
}

void FeatureMappingDesign::Pill::Update()
{
  assert(points.GetSize() >= 2);
  P.Set(points.First()[0].GetPlainDesignValue(), points.First()[1].GetPlainDesignValue(),
        dim_ == 3 ? points.First()[2].GetPlainDesignValue() : 0.0);
  Q.Set(points.Last()[0].GetPlainDesignValue(), points.Last()[1].GetPlainDesignValue(),
        dim_ == 3 ? points.Last()[2].GetPlainDesignValue() : 0.0);

  length = P.Dist(Q); // || Q - P ||
  length2 = length*length;

  U = (Q-P) / length;

  if(dim_ == 2) // in 3D the perpendicular direction depends on the point, no helpers to precompute
  {
    V.Set(U[1], -U[0]);

    dx = Q[0] - P[0]; // x-component of the vector from start to end
    dy = Q[1] - P[1]; // y-component
    dp_norm = Q[0] * P[1] - Q[1] * P[0]; // end_x*start_y - end_y*start_x
  }

  // only required in the anisotropic case
  bool aniso = domain->GetDesign()->GetMethod() == ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT;
  if(aniso)
  {
    angle = std::atan2(dy, dx);

    // we assume, we do anisotropy only in elasticity and we assume that we have in all design regions the same material
    rot_mat = FeatureMappingDesign::aniso_base_tensor; // copy C_0
    // RotateTensor needs a MaterialTensor which wraps a Matrix and assignes it with a notation (VOIGT)
    MaterialTensor<double> RCR(VOIGT, &rot_mat, false); // false: do not copy 
    // actually rot_mat is rotated
    Optimization::context->dm->RotateTensor(RCR, DE::NO_DERIVATIVE, true, angle); // true because we provide the angle
  }

  LOG_DBG(fm) << "P:U: P=" << P.ToString() << " Q=" << Q.ToString() << " l=" << length << " U=" << U.ToString() << " V=" << V.ToString() 
              << " dx= " << dx << " dy=" << dy << " dp_nom=" << dp_norm << " angle=" << angle 
              << (aniso ? " rot_mat=" + rot_mat.ToString() : "");
}

double FeatureMappingDesign::Pill::Distance(const Point& X, FeatureVariable::Tip* part) const
{
  assert(points.GetSize() >= 2);

  double pval = this->p.GetPlainDesignValue();
  double dist = -1;

  // test if we are in the region of the straight bar
  // if np.dot((X - self.Q), self.U0) <= 0 and np.dot((X - self.P), self.U0) >= 0:
  //    projected_distance = abs(np.dot((X - self.Q), self.V0))
  // Point XQ = X-Q;
  // LOG_DBG(fm) << "P:CD U.Dot(X-Q)=" << U.Dot(X-Q) << " U.Dot(X-P)=" << U.Dot(X-P) << " X-Q=" << XQ.ToString();
  if(length > 0 && U.Dot(X-Q) <= 0 && U.Dot(X-P) >= 0)
  {
    if(dim_ == 2)
    {
      // distance to line segment https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
      dist = std::abs(V.Dot(X-Q)) - pval;
      // double dl = (std::abs(dy*X[0] - dx*X[0] + dp_norm) / length) - pval;
    }
    else
    {
      // 3D capsule: perpendicular distance to the axis ||w - (U.w) U|| via Pythagoras
      Point w = X - P;
      double t = U.Dot(w);
      dist = std::sqrt(std::max(w.Dot(w) - t*t, 0.0)) - pval; // max guards rounding for points on the axis
    }
    if(part)
      *part = FeatureVariable::INNER;
    LOG_DBG(fm) << "P:CD X=" << X.ToString() << " p=" << pval << " P=" << P.ToString() << " Q=" << Q.ToString() << " -> line=" << dist;
  }
  else
  {
    double dp =  X.Dist(P)-pval; // distance to start point
    double dq =  X.Dist(Q)-pval; 
  
    dist = dp < dq ? dp : dq;
    if(part)
      *part = dp < dq ? FeatureVariable::START : FeatureVariable::END;

    LOG_DBG(fm) << "P:CD from " << X.ToString() << " p=" << pval << " P=" << P.ToString() << "->"  << dp 
              << " Q=" << P.ToString() << "->" << dq << " tip min=" << dist << " part=" << (part ? *part : -1);
  }  
  return dist;
}

void FeatureMappingDesign::Pill::GradDistance(const Point& X, FeatureVariable::Tip part, Vector<double>& out) const
{
  assert(part != FeatureVariable::NO_TIP);
  assert(points.GetSize() == 2);
  assert(out.GetSize() == 2*dim_ + 1 + (has_alpha ? 1 : 0));
  if(has_alpha)
    out[2*dim_ + 1] = 0.0; // the distance does not depend on the geometry variable alpha

  switch(part)
  {
    // variable order [Px Py (Pz) Qx Qy (Qz) p] as in GetAllVariables()
  case FeatureVariable::START:
  case FeatureVariable::END:
  {
    // d = ||X-C|| - p with the cap center C = P or Q, the other node does not enter
    const Point& C = part == FeatureVariable::START ? P : Q;
    int o = part == FeatureVariable::START ? 0 : dim_; // index offset of the block
    out.Init(0.0);
    for(unsigned int d = 0; d < dim_; d++)
      out[o+d] = (C[d]-X[d])/X.Dist(C);
    out[2*dim_] = -1.0; // p
    break;
  }
  case FeatureVariable::INNER:
  {
    if(dim_ == 3)
    {
      // envelope theorem for the foot point c = P + beta*(Q-P) of the infinite axis: with the
      // unit normal n = (X-c)/||X-c|| we get d_d/d_P = -(1-beta)*n and d_d/d_Q = -beta*n
      Point w = X - P;
      double beta = U.Dot(w) / length;
      Point r = w - (Q-P)*beta; // = X - c, perpendicular to the axis
      double dn = std::sqrt(r.Dot(r)); // = dist + p
      // on the axis the direction n is undefined (cone tip of the distance function). Such a point
      // is strictly inside the pill and the H'(dist) factor is zero, we pick the subgradient 0
      out.Init(0.0);
      if(dn > 0.0)
      {
        for(unsigned int d = 0; d < 3; d++)
        {
          out[d]   = -(1.0-beta) * r[d]/dn;
          out[3+d] = -beta * r[d]/dn;
        }
      }
      out[6] = -1.0; // p
      break;
    }
    // wolframalpha
    // diff (V0*(X0-Q0)+V1*(X1-Q1)) by Q1
    // V0 = ((Q1-P1) / sqrt((Q0-P0)**2 + (Q1-P1)**2))
    // V1 = ((P0-Q0) / sqrt((Q0-P0)**2 + (Q1-P1)**2))
    // diff (((Q1-P1) / sqrt((Q0-P0)**2 + (Q1-P1)**2))*(X0-Q0)+((P0-Q0) / sqrt((Q0-P0)**2 + (Q1-P1)**2))*(X1-Q1)) by Q1    
    // 
    // d/dQ1(((Q1 - P1) (X0 - Q0))/sqrt((Q0 - P0)^2 + (Q1 - P1)^2) + ((P0 - Q0) (X1 - Q1))/sqrt((Q0 - P0)^2 + (Q1 - P1)^2)) 
    // = ((Q0 - P0) (P0^2 - P0 (Q0 + X0) + P1^2 - P1 (Q1 + X1) + Q0 X0 + Q1 X1))/((P0 - Q0)^2 + (P1 - Q1)^2)^(3/2)
    
    // d/dP0 =  ((Q1 - P1) (P0 (Q0 - X0) + (P1 - Q1) (Q1 - X1) - Q0^2 + Q0 X0)) / ||P-Q||^3
    // d/dP1 = -((Q0 - P0) (P0 (Q0 - X0) + (P1 - Q1) (Q1 - X1) - Q0^2 + Q0 X0)) / ||P-Q||^3
    // d/dQ0 = -((Q1 - P1) (P0^2 - P0 (Q0 + X0) + P1^2 - P1 (Q1 + X1) + Q0 X0 + Q1 X1)) / ||P-Q||^3
    // d/dQ1 =  ((Q0 - P0) (P0^2 - P0 (Q0 + X0) + P1^2 - P1 (Q1 + X1) + Q0 X0 + Q1 X1)) / ||P-Q||^3

    Point QmP = Q - P;
    Point QmX = Q - X;
    Point QpX = Q + X;

    double denom = length*length*length;
    double sg = V.Dot(X-Q) > 0 ? 1.0 : -1.0; // to differentiate the abs from std::abs(V.Dot(X-Q)) 

    double dP_common = (P[0]*QmX[0]-QmP[1]*QmX[1]-Q[0]*Q[0]+Q[0]*X[0]);
    double dQ_common = (P[0]*P[0]-P[0]*QpX[0]+P[1]*P[1]-P[1]*QpX[1]+Q.Dot(X));

    // we assume the compiler sees that 
    out[0] = sg *  (QmP[1] * dP_common) / denom; // dP_x
    out[1] = sg * -(QmP[0] * dP_common) / denom; // dP_y
    out[2] = sg * -(QmP[1] * dQ_common) / denom; // dQ_x
    out[3] = sg *  (QmP[0] * dQ_common) / denom; // dQ_y
    out[4] = -1.0; // p
    break;
  } 
  case FeatureVariable::NO_TIP:
    assert(false);
    break;
  }
}

void FeatureMappingDesign::Pill::HessDistance(const Point& X, FeatureVariable::Tip part, Matrix<double>& out) const
{
  assert(part != FeatureVariable::NO_TIP);
  assert(points.GetSize() == 2);
  out.Resize(2*dim_ + 1, 2*dim_ + 1);
  out.Init(); // the p row/column stays zero in all cases as the distance is linear in p

  switch(part)
  {
  case FeatureVariable::START:
  case FeatureVariable::END:
  {
    // d = ||X-C|| - p with C = P or Q. The Hessian w.r.t. C is the projector
    // (I - m*m^T)/R with R = ||X-C|| and m = (C-X)/R, all other entries are zero
    const Point& C = part == FeatureVariable::START ? P : Q;
    int o = part == FeatureVariable::START ? 0 : dim_; // index offset of the block
    double R = X.Dist(C);
    for(unsigned int i = 0; i < dim_; i++)
      for(unsigned int j = 0; j < dim_; j++)
        out[o+i][o+j] = ((i == j ? 1.0 : 0.0) - (C[i]-X[i])/R * (C[j]-X[j])/R)/R;
    break;
  }
  case FeatureVariable::INNER:
  {
    if(dim_ == 3)
    {
      // differentiate the envelope gradient from GradDistance: with the foot parameter
      // beta = U.(X-P)/length, r = X-P-beta*(Q-P), dn = ||r||, n = r/dn and the projector
      // Pi = I - n*n^T the blocks are (verified against finite differences)
      //   H_PP = -dn/length^2 * m1*m1^T + (1-beta)^2/dn * Pi   with m1 =  n + (1-beta)*(Q-P)/dn
      //   H_PQ = -dn/length^2 * m1*m2^T + beta*(1-beta)/dn * Pi with m2 = -n + beta*(Q-P)/dn
      //   H_QQ = -dn/length^2 * m2*m2^T + beta^2/dn * Pi
      Point v = Q - P;
      Point w = X - P;
      double beta = U.Dot(w) / length;
      Point r = w - v*beta;
      double dn = std::sqrt(r.Dot(r));
      if(dn == 0.0)
        break; // on the axis the distance is not differentiable, keep zeros as in GradDistance
      Point n = r / dn;
      Point m1 = n + v*((1.0-beta)/dn);
      Point m2 = v*(beta/dn) - n;
      double f = dn/length2;
      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
          double Pi = (i == j ? 1.0 : 0.0) - n[i]*n[j];
          out[i][j]     = -f*m1[i]*m1[j] + (1.0-beta)*(1.0-beta)/dn * Pi;
          out[i][3+j]   = -f*m1[i]*m2[j] + beta*(1.0-beta)/dn * Pi;
          out[3+i][j]   = -f*m2[i]*m1[j] + beta*(1.0-beta)/dn * Pi;
          out[3+i][3+j] = -f*m2[i]*m2[j] + beta*beta/dn * Pi;
        }
      break;
    }
    // d = sg*n/D - p with the numerator n = V.(X-Q)*D and D = ||Q-P|| = length, see GradDistance.
    // Product/chain rule in structured form (verified with sympy):
    //   Hess(n/D) = Hess(n)/D - (gn*gD^T + gD*gn^T)/D^2 + n*(2*gD*gD^T/D^3 - Hess(D)/D^2)
    // with derivatives w.r.t. (P_x, P_y, Q_x, Q_y):
    //   gn = (y-q_y, q_x-x, p_y-y, x-p_x), Hess(n) only (P_x,Q_y) = -1 and (P_y,Q_x) = +1
    //   gD = (-u_x, -u_y, u_x, u_y),       Hess(D) = blocks [[Hv,-Hv],[-Hv,Hv]], Hv = (I-u*u^T)/D
    double sg = V.Dot(X-Q) > 0 ? 1.0 : -1.0; // to differentiate the abs, as in GradDistance
    double n = V.Dot(X-Q) * length;
    const double gn[4] = {X[1]-Q[1], Q[0]-X[0], P[1]-X[1], X[0]-P[0]};
    const double gD[4] = {-U[0], -U[1], U[0], U[1]};
    const double Hv[2][2] = {{(1.0 - U[0]*U[0])/length, -U[0]*U[1]/length},
                             {-U[0]*U[1]/length, (1.0 - U[1]*U[1])/length}};
    for(int i = 0; i < 4; i++)
      for(int j = 0; j < 4; j++)
      {
        double HD = (i/2 == j/2 ? 1.0 : -1.0) * Hv[i%2][j%2];
        out[i][j] = sg * (-(gn[i]*gD[j] + gD[i]*gn[j])/length2 + n * (2.0*gD[i]*gD[j]/length - HD)/length2);
      }
    // the Hess(n)/D term: the bilinear numerator has only these two second derivatives
    out[0][3] += sg * -1.0/length;
    out[3][0] += sg * -1.0/length;
    out[1][2] += sg *  1.0/length;
    out[2][1] += sg *  1.0/length;
    break;
  }
  case FeatureVariable::NO_TIP:
    assert(false);
    break;
  }
}

void FeatureMappingDesign::Pill::HessLength(Matrix<double>& out) const
{
  // the native 'distance' constraint is c = ||Q-P|| = length. Its exact Hessian w.r.t.
  // (Px,Py,(Pz),Qx,Qy,(Qz)) is the [[Hv,-Hv],[-Hv,Hv]] block with Hv = (I - U U^T)/length and
  // U = (Q-P)/length - exactly the Hess(D) term documented in HessDistance()'s INNER case.
  // The profile p is linear (absent), so its row/column stays zero.
  int d = (int) dim_;
  out.Resize(2*d + 1, 2*d + 1);
  out.Init();
  assert(length > 0.0);
  for(int i = 0; i < 2*d; i++)
    for(int j = 0; j < 2*d; j++)
      out[i][j] = (i/d == j/d ? 1.0 : -1.0) * ((i%d == j%d ? 1.0 : 0.0) - U[i%d]*U[j%d])/length;
}

inline double FeatureMappingDesign::Pill::GradAngle(int s) const
{
  assert(s >= 0 && s <= 4);
  if(length < 1e-12)
    return 0.0;

  // angle = atan2(dy,dx)
  // dx = Q_x - P_x, dy = Q_y-P_y
  // d_atan2/d_y = x/(x^2+y^2) d_atan2/d_x = -y/(x^2+y^2)
  switch(s)
  {
  case 0:
    return  dy/length2; // P_x: -dy * -1
  case 1:
    return -dx/length2; // P_y: dx * -1
  case 2:
    return -dy/length2; // Q_x: -dy * 1
  case 3:
    return  dx/length2; // Q_y: dx * 1
  case 4:
    return  0.0;       // p: no angle impact
  default:
    assert(false);
  }
  return -1.0;
}

inline void FeatureMappingDesign::BoundaryFunction::Eval(const Vector<double>& dist, Vector<double>& eval) const
{
  eval.Resize(dist.GetSize());
  for(unsigned int i = 0; i < dist.GetSize(); i++)
    eval[i] = Eval(dist[i]);
}

inline bool FeatureMappingDesign::BoundaryFunction::IsConstant(const ItemIP* item_ip) const
{
  for(const IntegrationPoint* ip : item_ip->corner)
  {
    for(double d : ip->dist)
      if(!IsConstant(d))
        return false;
  }
  return true;
}

double FeatureMappingDesign::BoundaryFunction::Hessian(double dist) const
{
  throw Exception("second derivative not implemented for this boundary function, use quintic or bezier");
}

std::string FeatureMappingDesign::BoundaryFunction::DebugLog()
{
  std::stringstream ss;
  ss << "BF:DL t=" << 2*h << " h=" << h << " rhomin=" << rhomin << ": ";

  for(double d : {2*h, h, h/2, h/4, 0.0, -h/2, -h})
    ss << " f(" << d << ")=" << Eval(d) << ", ";
  return ss.str();
}

FeatureMappingDesign::Poly::Poly(double transition, double rhomin)
{
  this->h = transition/2; // half transition zone!
  this->rhomin = rhomin;
  this->fac = 3.0*(1.0-rhomin)/4.0;
  this->bias = ((1+rhomin)/2);
  this->h3 = h*h*h;
}


inline double FeatureMappingDesign::Poly::Eval(double d) const
{
  // note that the Feature-Mapping Review paper uses phi(x) which is -dist(x). We use here distance!
  if(IsOutside(d))
    return rhomin;
  if(IsInside(d))
    return 1;

  double res = fac * (((d*d*d)/(3*h3))-(d/h))+bias; // d = -phi!
  assert(res >= -1e-16);
  res = std::max(res,0.0); // remove numerical artefacts to be very slightly negative
  return res;
}

inline double FeatureMappingDesign::Poly::Grad(double d) const
{
  if(IsConstant(d))
    return 0.0;

  return fac * (d*d - h*h)/h3;
}

FeatureMappingDesign::Quintic::Quintic(double transition, double rhomin)
{
  this->h = transition/2; // half transition zone as for Poly
  this->rhomin = rhomin;
  this->fac = (1.0-rhomin)/16.0;
  this->bias = (1.0+rhomin)/2;
  this->grad_fac = 15.0*(1.0-rhomin)/(16.0*h);
  this->hess_fac = 15.0*(1.0-rhomin)/(4.0*h*h);
}

double FeatureMappingDesign::Quintic::Hessian(double d) const
{
  if(IsConstant(d)) // continuous transition as H''(-h) = H''(h) = 0
    return 0.0;

  double s = d/h;
  return hess_fac * s * (1.0 - s*s);
}

inline double FeatureMappingDesign::Quintic::Eval(double d) const
{
  if(IsOutside(d))
    return rhomin;
  if(IsInside(d))
    return 1;

  double s = d/h;
  double s2 = s*s;
  return bias - fac * s * (15.0 - s2*(10.0 - 3.0*s2));
}

inline double FeatureMappingDesign::Quintic::Grad(double d) const
{
  if(IsConstant(d))
    return 0.0;

  double t = 1.0 - (d*d)/(h*h);
  return -grad_fac * t*t;
}

FeatureMappingDesign::Bezier::Bezier(double transition, double rhomin, double extension)
{
  assert(transition > 0.0 && extension >= 0.0);
  this->a = transition/2; // as for Poly, transition is the full width of the symmetric transition zone
  this->b = transition/2 + extension;
  // the base class helpers IsConstant()/IsInside()/IsOutside() are symmetric, h = b is
  // conservative: dist in (-b, -a] is reported non-constant but Eval() is correctly 1 there
  this->h = b;
  this->rhomin = rhomin;

  gamma = FindGamma();
  double c = (a - b) / 30.0;
  const double wx[6] = {-a, c - gamma, c - gamma, c + gamma, c + gamma, b};
  std::copy(wx, wx + 6, Wx);

  // sample look-up tables on a uniform distance grid over [-a, b]. The grid spacing follows the
  // inner zone scale a where the curve is sharpest. The factor 500 gives Hermite interpolation
  // errors of ~4e-11/4e-10/7e-9 (in units of a) for H/H'/H'', two orders below a 1e-7 optimizer
  // tolerance. The cap bounds memory for extreme extensions, where the flat tail tolerates it
  n_tab = std::min(int(500.0 * b / a), 20000);
  H_tab.Resize(n_tab + 1);
  dH_tab.Resize(n_tab + 1);
  ddH_tab.Resize(n_tab + 1);
  dddH_tab.Resize(n_tab + 1);
  delta = (a + b) / n_tab;
  for(int i = 0; i <= n_tab; i++)
  {
    double t = SolveT(-a + i * delta);
    H_tab[i]    = EvalByT(Wy, t);
    dH_tab[i]   = BezierDyDx(Wx, Wy, t);
    ddH_tab[i]  = BezierD2yDx2(Wx, Wy, t);
    dddH_tab[i] = BezierD3yDx3(Wx, Wy, t);
  }
}

inline double FeatureMappingDesign::Bezier::Interpolate(const Vector<double>& tab, const Vector<double>& dtab, double dist) const
{
  double s = (dist + a) / delta;
  int i = std::min(int(s), n_tab - 1); // guard floating point rounding at dist -> b
  double w = s - i;
  // cubic Hermite basis on the unit interval, the slopes are scaled by the grid spacing
  double w2 = w * w;
  double w3 = w2 * w;
  return (2*w3 - 3*w2 + 1) * tab[i]   + (w3 - 2*w2 + w) * delta * dtab[i]
       + (3*w2 - 2*w3)     * tab[i+1] + (w3 - w2)       * delta * dtab[i+1];
}



double FeatureMappingDesign::Bezier::Bernstein(int n, const double* C, const double* W, double t)
{
  // the straightforward implementation is 5-7 times slower as std::pow is expensive
  // even for small integer exponents. The iterative powers below are bit-identical:
  //   double r = 0.0;
  //   for(int i = 0; i <= n; i++)
  //     r += C[i] * std::pow(t, i) * std::pow(1.0 - t, n - i) * W[i];
  //   return r;
  assert(n <= 5);
  double sp[6]; // powers of (1-t) up to n
  sp[0] = 1.0;
  for(int i = 1; i <= n; i++)
    sp[i] = sp[i-1] * (1.0 - t);

  double r = 0.0;
  double tp = 1.0; // t^i
  for(int i = 0; i <= n; i++)
  {
    r += C[i] * tp * sp[n-i] * W[i];
    tp *= t;
  }
  return r;
}

double FeatureMappingDesign::Bezier::EvalByT(const double* W, double t)
{
  const double C[6] = {1, 5, 10, 10, 5, 1};
  return Bernstein(5, C, W, t);
}

double FeatureMappingDesign::Bezier::FirstDerivT(const double* W, double t)
{
  const double C[5] = {1, 4, 6, 4, 1};
  const double dW[5] = {W[1]-W[0], W[2]-W[1], W[3]-W[2], W[4]-W[3], W[5]-W[4]};
  return 5.0 * Bernstein(4, C, dW, t);
}

double FeatureMappingDesign::Bezier::SecondDerivT(const double* W, double t)
{
  const double C[4] = {1, 3, 3, 1};
  const double d2W[4] = {W[2]-2*W[1]+W[0], W[3]-2*W[2]+W[1], W[4]-2*W[3]+W[2], W[5]-2*W[4]+W[3]};
  return 20.0 * Bernstein(3, C, d2W, t);
}

double FeatureMappingDesign::Bezier::ThirdDerivT(const double* W, double t)
{
  const double C[3] = {1, 2, 1};
  const double d3W[3] = {W[3]-3*W[2]+3*W[1]-W[0], W[4]-3*W[3]+3*W[2]-W[1], W[5]-3*W[4]+3*W[3]-W[2]};
  return 60.0 * Bernstein(2, C, d3W, t);
}

double FeatureMappingDesign::Bezier::BezierDyDx(const double* Wx, const double* Wy, double t)
{
  return FirstDerivT(Wy, t) / FirstDerivT(Wx, t);
}

double FeatureMappingDesign::Bezier::BezierD2yDx2(const double* Wx, const double* Wy, double t)
{
  double dx = FirstDerivT(Wx, t);
  return (SecondDerivT(Wy, t) * dx - FirstDerivT(Wy, t) * SecondDerivT(Wx, t)) / (dx*dx*dx);
}

double FeatureMappingDesign::Bezier::BezierD3yDx3(const double* Wx, const double* Wy, double t)
{
  double dx  = FirstDerivT(Wx, t);
  double dy  = FirstDerivT(Wy, t);
  double ddx = SecondDerivT(Wx, t);
  return (dx * (ThirdDerivT(Wy, t) * dx - dy * ThirdDerivT(Wx, t))
          - 3.0 * ddx * (SecondDerivT(Wy, t) * dx - dy * ddx)) / std::pow(dx, 5);
}



double FeatureMappingDesign::Bezier::FindGamma() const
{
  // the minimax objective is flat around gamma*, this resolution stays within 0.3% of the optimal
  // max |H''| while being 25 times cheaper than the 500/1000 of the original python implementation
  const int n_grid = 100; // gamma candidates
  const int n_t = 200;    // t samples for monotonicity check and max |H''|
  const double c = (a - b) / 30.0;
  const double g_max = 2.25 * a; // empirical search range, validity is checked via x'(t)

  double best_gamma = 0.0;
  double best_val = std::numeric_limits<double>::max();
  for(int k = 0; k < n_grid; k++)
  {
    double g = -g_max + k * 2.0 * g_max / (n_grid - 1);
    const double wx[6] = {-a, c - g, c - g, c + g, c + g, b};
    double val = 0.0;
    for(int i = 0; i <= n_t; i++)
    {
      double t = double(i) / n_t;
      if(FirstDerivT(wx, t) <= 1e-9) // skip gamma where x(t) is not strictly monotone
      {
        val = std::numeric_limits<double>::max();
        break;
      }
      val = std::max(val, std::abs(BezierD2yDx2(wx, Wy, t)));
    }
    if(val < best_val)
    {
      best_val = val;
      best_gamma = g;
    }
  }
  assert(best_val < std::numeric_limits<double>::max());
  return best_gamma;
}

double FeatureMappingDesign::Bezier::SolveT(double dist, double tol) const
{
  assert(dist >= -a && dist <= b);
  assert(tol > 0.0);
  double lo = 0.0;
  double hi = 1.0;
  // bisection, x(t) is strictly monotone by FindGamma(). The iteration cap guards against
  // tol below one ulp (~1.1e-16 near t=1), where the interval cannot shrink any further
  for(int i = 0; i < 60 && hi - lo > tol; i++)
  {
    double mid = 0.5 * (lo + hi);
    if(EvalByT(Wx, mid) < dist)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

double FeatureMappingDesign::Bezier::Eval(double dist, bool exact) const
{
  if(dist <= -a)
    return 1.0;
  if(dist >= b)
    return rhomin;

  double H = exact ? EvalByT(Wy, SolveT(dist)) : Interpolate(H_tab, dH_tab, dist);
  return rhomin + (1.0 - rhomin) * H;
}

double FeatureMappingDesign::Bezier::Grad(double dist, bool exact) const
{
  if(dist <= -a || dist >= b)
    return 0.0;

  double dH = exact ? BezierDyDx(Wx, Wy, SolveT(dist)) : Interpolate(dH_tab, ddH_tab, dist);
  return (1.0 - rhomin) * dH;
}

double FeatureMappingDesign::Bezier::Hessian(double dist, bool exact) const
{
  if(dist <= -a || dist >= b) // H'' is continuous here as the quintic curve has H''(-a) = H''(b) = 0
    return 0.0;

  double ddH = exact ? BezierD2yDx2(Wx, Wy, SolveT(dist)) : Interpolate(ddH_tab, dddH_tab, dist);
  return (1.0 - rhomin) * ddH;
}

inline double FeatureMappingDesign::Max::Eval(const Vector<double>& projected) const
{
  return projected.Max();
}

inline double FeatureMappingDesign::Max::Grad(const Vector<double>& projected, unsigned int num) const
{
  // Max combination is not really differentiable as only one feature can contribute, the others are 0
  // with more features 1, it is almost random. std::max_element() find another max than our Vector::Max() 
  // but as we are not differentiable anyway, we don't care
  assert(num < projected.GetSize());
  unsigned int idx = std::distance(projected.GetPointer(), std::max_element(projected.GetPointer(), projected.GetPointer()+projected.GetSize()));
  LOG_DBG3(fm) << "Max:G prj=" << projected.ToString() << " by " << num << " idx=" << idx << " -> " << (num == idx ? 1.0 : 0.0);
  return num == idx ? 1.0 : 0.0; // when we have 1 feature only, this is the plain rho
}

double FeatureMappingDesign::P_Norm::Eval(const Vector<double>& projected) const
{
   double sum = 0.0;
   for(double v : projected)
     sum += std::pow(v, p);
   double res = std::pow(sum, 1.0/p);
   LOG_DBG2(fm) << "PN:E " << projected.ToString() << " p=" << p << " sum=" << sum << " res=" << res;
   return res;
}

double FeatureMappingDesign::P_Norm::Grad(const Vector<double>& projected, unsigned int num)  const
{
  assert(projected.GetSize() > num);
  double sum = 0.0;
  for(double v : projected)
    sum += std::pow(v, p);
  if(sum == 0.0)
    return 0.0; // the p-norm is not differentiable at the origin (all rho zero, e.g. alpha=0 with rhomin=0), avoid inf*0=NaN
  double res = std::pow(sum, 1.0/p - 1.0) * std::pow(projected[num],p-1); // to be multiplied with d_rho_s/d_s
  LOG_DBG3(fm) << "PN:G by " << num << " v=" << projected.ToString() << " p=" << p << " -> " << res;
  return res;
}

double FeatureMappingDesign::P_Norm::Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const
{
  assert(f < projected.GetSize() && g < projected.GetSize());
  double sum = 0.0;
  for(double v : projected)
    sum += std::pow(v, p);
  // (p-1) * (delta_fg * v_f^(p-2) * mrho^(1-p) - (v_f v_g)^(p-1) * mrho^(1-2p)) with mrho = sum^(1/p)
  double res = -(p-1.0) * std::pow(projected[f] * projected[g], p-1) * std::pow(sum, 1.0/p - 2.0);
  if(f == g)
    res += (p-1.0) * std::pow(projected[f], p-2) * std::pow(sum, 1.0/p - 1.0);
  return res;
}

double FeatureMappingDesign::SoftMax::Eval(const Vector<double>& projected) const
{
  // (sum_f rho_f*exp(beta*rho_f)) / (sum_f exp(beta*rho_f))
  double num = 0.0;
  double denum = 0.0;

  for(double x : projected)
  {
      double exp = std::exp(beta * x);
      num += x * exp;
      denum += exp;
  }
  return num/denum;
}

double FeatureMappingDesign::SoftMax::Grad(const Vector<double>& projected, unsigned int num)  const
{
  double sum_x_exp = 0.0;
  double sum_exp = 0.0;

  for(double x : projected)
  {
      double exp = std::exp(beta * x);
      sum_x_exp += x * exp;
      sum_exp += exp;
  }
  // u = sum_x_exp (sum x_i*exp(beta_x_i))
  // v = sum_exp (sum exp(beta*x_i))
  // u' = exp(beta*x_j)+x_j*beta*exp(beta*x_j)=(1+beta*x_j)*exp(beta*x_j)
  // v' = beta*exp(beta*x_j)
  assert(projected.GetSize() > num);
  double x_j = projected[num];
  return std::exp(beta * x_j)*((1+beta*x_j)*sum_exp - beta*sum_x_exp)/(sum_exp*sum_exp);
}

double FeatureMappingDesign::KreisselmeierSteinhauser::Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const
{
  // with the softmax weights sigma_j = exp(beta*x_j)/sum_exp the gradient is sigma_f
  // and the Hessian beta * sigma_f * (delta_fg - sigma_g)
  double sigma_f = DerivSmoothMax(projected, beta, f);
  double sigma_g = f == g ? sigma_f : DerivSmoothMax(projected, beta, g);
  return beta * sigma_f * ((f == g ? 1.0 : 0.0) - sigma_g);
}

double FeatureMappingDesign::SoftMax::Hessian(const Vector<double>& projected, unsigned int f, unsigned int g) const
{
  assert(f < projected.GetSize() && g < projected.GetSize());
  double sum_x_exp = 0.0;
  double sum_exp = 0.0;

  for(double x : projected)
  {
      double exp = std::exp(beta * x);
      sum_x_exp += x * exp;
      sum_exp += exp;
  }
  // with the weights sigma_j = exp(beta*x_j)/sum_exp, the value m = sum_i x_i*sigma_i and
  // G_j = 1 + beta*(x_j - m), the gradient is sigma_f*G_f and the Hessian
  // beta * sigma_f * (delta_fg*(1 + G_f) - sigma_g*(G_f + G_g)), symmetric in f and g
  double m = sum_x_exp / sum_exp;
  double sigma_f = std::exp(beta * projected[f]) / sum_exp;
  double sigma_g = std::exp(beta * projected[g]) / sum_exp;
  double G_f = 1.0 + beta * (projected[f] - m);
  double G_g = 1.0 + beta * (projected[g] - m);
  return beta * sigma_f * ((f == g ? 1.0 + G_f : 0.0) - sigma_g * (G_f + G_g));
}

void FeatureMappingDesign::AnisoCurrentVariable::Reset()
{
  var_type = DE::NO_DERIVATIVE;
  var_idx = -1;
  feature_idx = -1;
  assert(dim_ == 2);
  RCdR.Resize(3,3);
  RCdR.Init(); 
}

} // end of namespace

