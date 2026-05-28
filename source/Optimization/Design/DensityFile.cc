#include <boost/algorithm/string.hpp>
#include <iostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/BaseFE.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/FeaturedDesign.hh"
#include "Optimization/Design/ShapeMapDesign.hh"
#ifdef USE_EMBEDDED_PYTHON // currently only the python version
  #include "Optimization/Design/SpaghettiDesign.hh"
#endif
#include "Optimization/Design/SplineBoxDesign.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "PDE/BasePDE.hh"
#include "Utils/Point.hh"
#include "Utils/StdVector.hh"

using namespace CoupledField;
using std::string;

DEFINE_LOG(density, "density")

DensityFile::DensityFile(DesignSpace* designSpace,
                            PtrParamNode export_pn,
                            ParamNodeList& des,
                            ParamNodeList& tfs,
                            PtrParamNode regulize_pn)
{
  this->space_ = designSpace;
  this->compress_ = export_pn->Get("compress")->As<bool>();

  name_ = export_pn->Get("file")->As<string>();
  if(name_ == "[problem]")
    name_ = progOpts->GetSimName() + ".density.xml";
  data = Create(des, tfs, regulize_pn, designSpace->DoNonDesignVicinity());
  all_iterations_ = export_pn->Get("save")->As<string>() == "all";
  finally_only_   = export_pn->Get("write")->As<string>() == "finally";
  write_density_  = export_pn->Get("density")->As<bool>();
  // append .gz if compress=true and not already file ends with it
  // note that gz is much faster than the better compressing bz2. ParamNode::ToFile() automatically compresses on .gz name
  if(compress_ && !boost::algorithm::ends_with(name_, ".gz"))
    name_ = name_ + ".gz";
}


DensityFile::~DensityFile()
{
  // in CommitIteration or here?
  if(finally_only_)
    data->ToFile(name_, true);

  data.reset();
}

bool DensityFile::NeedLoadErsatzMaterial()
{
  return domain->GetParamRoot()->Has("loadErsatzMaterial") || progOpts->GetErsatzMaterialStr() != "";
}

DesignSpace* DensityFile::CreateDesignSpace(bool force_region, const PtrParamNode& pn, const ParamNodeList& elems, const PtrParamNode& xml)
{
  Grid* grid = domain->GetGrid();
  PtrParamNode info = domain->GetInfoRoot();
  const unsigned int elsize = elems.GetSize();

  // only if the design space does not already exist (created by optimization)
  // the regions are normally implicitly defined by the element numbers. The exception
  // is force_region from <loadErsatzMaterial>
  StdVector<RegionIdType> regionIds;
  // check if we ignore the element numbers
  if(force_region)
  {
    regionIds.Push_back(grid->GetRegion().Parse(pn->Get("force_region")->As<string>()));
  }
  else
  {
    // find the regions by ourselves
    for (unsigned int e = 0; e < elsize; ++e)
    {
      unsigned int nr = elems[e]->Get("nr")->As<unsigned int>();
      if (!regionIds.Contains(grid->GetElem(nr)->regionId))
        regionIds.Push_back(grid->GetElem(nr)->regionId);
    }
  }
  // create the design space -> data has initial values!
  DesignSpace* space = new DesignSpace(regionIds, xml->Get("header"), ErsatzMaterial::SIMP_METHOD);
  space->PostInit(0, 0); // no objectives, no constraints
  // is cheap - for density filtering
  DesignStructure ds(space, space->GetRegionIds());
  PtrParamNode reg = xml->Get("header/filters/filter", ParamNode::PASS);
  if(reg)
  {
    ds.SetFilter(reg);
    ds.WriteFilterInfo(info->Get("ersatzMaterial"));
    space->SetFilterType(ds.GetCommonFilterType());
  }

  space->ToInfo(NULL);
  return space;
}

bool DensityFile::ReadDensity(PtrParamNode pn, const ParamNodeList& elems, bool force_region, DesignSpace* space,
    double& lower_violation, double& upper_violation)
{

  DesignElement::Type last_dt = DesignElement::NO_TYPE;
  bool enforce_bounds = false;
  double relative_bound = -1.0;
  const unsigned int elsize = elems.GetSize();

  // check the the dimensions! the number of design variables comes from the regions and designs
  if(space->data.GetSize() != elsize)
  {
    string msg = "the number of elements in the density file does not match the number of elements of the region!\n"\
                 "         check the results carefully!";
    domain->GetInfoRoot()->Get("ersatzMaterial")->SetWarning(msg);
  }

  string name = "design";
  if (pn != NULL && pn->Has("name"))
    name = pn->Get("name")->As<string>();

  for (unsigned int e = 0; e < elsize; ++e)
  {
    // the design set consists of entries like
    // <element nr="401" type="density" design="0.886466" physical="0.800454" />
    // only the combination nr and type is unique. E.g. in piezo we might have types density and polarization
    unsigned int nr = elems[e]->Get("nr")->As<unsigned int>();
    DesignElement::Type dt = (DesignElement::Type) DesignElement::type.Parse(elems[e]->Get("type")->As<string>());

    if (dt != last_dt)
    {
      // we don't want to have different enforce_bounds for the different designs. What is with the regions anyway??
      assert(!(last_dt != DesignElement::NO_TYPE && enforce_bounds != space->design[space->FindDesign(dt)].enforce_bounds));
      last_dt = dt;
      enforce_bounds = space->design[space->FindDesign(dt)].enforce_bounds;
      relative_bound = space->design[space->FindDesign(dt)].relative_bound;
    }

    double val;
    if (name != "design" && elems[e]->Has(name))
      val = elems[e]->Get(name)->As<double>();
    else
      val = elems[e]->Get("design")->As<double>();
    int idx = dt == DesignElement::MULTIMATERIAL ? elems[e]->Get("index")->As<int>() : -1;

    // replace the value of the DesignElement
    // we call Find(..,..,false) for meshes with two regions (e. g. cube and void)
    // where we want to ignore the "void"-region completely
    DesignElement* de = force_region ? &(space->data[e]) : space->Find(nr, dt, false, false, idx);
    assert(de == NULL || de->GetType() == dt);

    if (dt == DesignElement::MULTIMATERIAL)
    {
      de->multimaterial = &(space->GetMultiMaterials()[idx]);
      assert(de->multimaterial->index == idx);
    }

    // this is also for the void-region! mainly for computing high resolution inv hom problems
    if (de != NULL) // && regionIds.Find(de->elem->regionId) >= 0)
    {
      lower_violation = std::max(lower_violation, de->GetLowerBound() - val);
      upper_violation = std::max(upper_violation, val - de->GetUpperBound());

      if (enforce_bounds)
        de->SetDesign(std::min(de->GetUpperBound(), std::max(de->GetLowerBound(), val)));
      else
        de->SetDesign(val);

      // Get value of the relative bound for current design variable. If value not set, db = -1.
      if(relative_bound >= 0.0)
      {
        // if a relative_bound is set in the xml file, upper and lower bound are overwritten
        de->SetUpperBound(std::min(de->GetUpperBound(), val + relative_bound/2));
        de->SetLowerBound(std::max(de->GetLowerBound(), val - relative_bound/2));
      }
    }
  }
  return enforce_bounds;
}

DesignSpace* DensityFile::ReadErsatzMaterial(DesignSpace* space)
{
  PtrParamNode info = domain->GetInfoRoot();

  // perhaps Optimization has already called the SetEnums
  if(DesignElement::type.map.empty())
    DesignElement::SetEnums();
  if(Objective::type.map.empty())
    Optimization::SetEnums();

  // do we have a command line switch? we then use the filename and the last set
  bool cmd = progOpts->GetErsatzMaterialStr() != "";
  PtrParamNode pn; // some default auto pointer NULL stuff
  if(!cmd)
    pn = domain->GetParamRoot()->Get("loadErsatzMaterial");

  string file = cmd ? progOpts->GetErsatzMaterialStr() : pn->Get("file")->As<string>();

  // to be appended by the set name
  std::cout << "++ Load ersatz material file: '" << file << "'" << std::flush;

  PtrParamNode in = space ? space->GetInfo()->Get("ersatzMaterialFile") : info->Get("ersatzMaterialFile");
  in->Get("file")->SetValue(file);
  in->Get("source")->SetValue(cmd ? "command line" : "problem file");

  // we read something like <loadErsatzMaterial region="piezo" file="piezo_density.xml" set="last"/>
  // Initialize our xml parser to handle the external xml file
  // set the global ParamNode tree pointer
  PtrParamNode xml = XmlReader::ParseFile(file);

  // check this file
  if (xml->Count("set") == 0)
    throw Exception("There are no design sets in the ersatz material file");

  // find the proper design set. This is either 'first', 'last' or the * in <set id="*"> ...
  PtrParamNode set;
  string key = cmd ? "last" : pn->Get("set")->As<string>();
  if (key == "first")
    set = xml->GetList("set")[0];
  if (key == "last")
    set = xml->GetList("set").Last();
  if (set == NULL)
    set = xml->GetByVal("set", "id", key);

  // finish the output as we have now the set information
  std::cout << "/'" << set->Get("id")->As<string>() + "'" << std::endl;

  // read the set and replace the initial values for the optimization
  ParamNodeList elems = set->GetList("element"); // we be 0 for shape map

  // shall the bounds be enforced?
  bool force_region = pn != NULL && pn->Has("force_region");

  if(!space)
  {
    // only if the design space does not already exist (created by optimization)
    // the regions are normally implicitly defined by the element numbers. The exception
    // is force_region from <loadErsatzMaterial>
    space = CreateDesignSpace(force_region, pn, elems, xml);

    // In case where we read density file from external file the matrix based filtering is buggy.
    // Dirty fix of disabling it for now
    space->is_matrix_filt = false;

    if (domain->GetParamRoot()->Has("filters/use_mat_filt") && domain->GetParamRoot()->Get("filters/use_mat_filt")->As<bool>())
      EXCEPTION("Using Matrix based filtering is not tested when loading material from external file")
  }

  // check bound violations
  double lower_violation = 0.0;
  double upper_violation = 0.0;

  // shape map, spaghetti and spline box
  FeaturedDesign* fd = dynamic_cast<FeaturedDesign*>(space);
  if(fd)
    fd->ReadDensityXml(set, lower_violation, upper_violation);

  // parse slack variable
  PtrParamNode slack = set->Get("slack", ParamNode::ActionType::PASS);
  if (space->HasSlackVariable() && slack != nullptr) {
    AuxDesign* auxd = dynamic_cast<AuxDesign*>(space);
    assert(auxd);
    auxd->GetSlackDesign()->SetDesign(slack->Get("design")->As<double>());
    LOG_DBG(density) << "REM: slack=" <<  slack->Get("design")->As<double>() << " set=" << space->GetSlackVariable();
  }
  
  // enforce bounds
  bool enforce_bounds = fd ? false : ReadDensity(pn, elems, force_region, space, lower_violation, upper_violation);

  if(lower_violation > 1e-5) {
    std::string msg = "the external design violates lower design bounds up to " + std::to_string(lower_violation);
    if(enforce_bounds)
      in->Get("bound_violation")->SetValue(msg + " but has been cropped due to 'enforce_bounds'");
    else
      in->SetWarning(msg);
  }
  if(upper_violation > 1e-5) {
    std::string msg = "the external design violates upper design bounds up to " + std::to_string(upper_violation);
    if(enforce_bounds)
      in->Get("bound_violation", ParamNode::APPEND)->SetValue(msg + " but has been cropped due to 'enforce_bounds'");
    else
      in->SetWarning(msg);
  }

  return space;

}

PtrParamNode DensityFile::Create(ParamNodeList& des, ParamNodeList& tfs, PtrParamNode regularize, bool non_desing_vicinity)
{
   PtrParamNode in = PtrParamNode(new ParamNode(ParamNode::INSERT));
   in->SetName("cfsErsatzMaterial");

   // write header
   PtrParamNode in_ = in->Get("header");

   LOG_DBG(density) << "Create: regular=" << this->space_->IsRegular();

   // design space can be regular, but grid is probably not
   StdVector<unsigned int> grid = domain->GetGrid()->CalcRegulardGridDiscretization();
   if(!grid.IsEmpty())
   {
     PtrParamNode mesh = in_->Get("mesh");
     mesh->Get("x")->SetValue(grid[0]);
     mesh->Get("y")->SetValue(grid[1]);
     mesh->Get("z")->SetValue(grid[2]);
   }

   in_->Get("designSpace/non_design_vicinity")->SetValue(non_desing_vicinity);

   for(unsigned int i = 0, s = des.GetSize(); i < s; ++i)
     in_->Get("dummy", ParamNode::APPEND)->SetValue(des[i], true); // name is overwritten

   for(unsigned int i = 0, s = tfs.GetSize(); i < s; ++i)
     in_->Get("dummy", ParamNode::APPEND)->SetValue(tfs[i], true);

   if(regularize)
     in_->Get("dummy")->SetValue(regularize, true);

   domain->ToInfo(in_, true); // coordinateSystems and enforce min_x, max_x, ...
   const Matrix<double>& bounds = space_->domainBounds;
   PtrParamNode t = in_->Get("optimizationDomain");
   //  <designBounds min_x="0" max_x="1" min_y="0" max_y="1"/>
   t->Get("min_x")->SetValue(bounds[0][0]);
   t->Get("max_x")->SetValue(bounds[0][1]);
   t->Get("min_y")->SetValue(bounds[1][0]);
   t->Get("max_y")->SetValue(bounds[1][1]);
   if(domain->GetDim() > 2) {
     t->Get("min_z")->SetValue(bounds[2][0]);
     t->Get("max_z")->SetValue(bounds[2][1]);
   }


   // off the design space to add stuff to the header. Done by SpaghettiDesing
   space_->AddToDensityHeader(in_);

   return in;
}


void DensityFile::SetAndWriteCurrent(int current_iteration)
{
  //FIXME
  // check if design changed since last write to avoid writing same data again

  PtrParamNode in = data->Get("set", all_iterations_ ? ParamNode::APPEND : ParamNode::INSERT);
  // add the entry, note that the iteration counter was incremented in base implementation
  in->Get("id")->SetValue(current_iteration);

  // we use the fast (dirty) bulk block to be (measurable) faster
  StdVector<std::string>& block = in->GetFastBulkBlock();


  // for shape map we also want to export DesignSpace::data even if this are no design variables
  assert((space_->GetNumberOfFeatureMappingVariables() >  0 && (dynamic_cast<ShapeMapDesign*>(space_) != NULL || dynamic_cast<FeaturedDesign*>(space_) != NULL) ) ||
         (space_->GetNumberOfFeatureMappingVariables() == 0 &&  dynamic_cast<ShapeMapDesign*>(space_) == NULL && dynamic_cast<FeaturedDesign*>(space_) == NULL));

  unsigned int size = (write_density_ ? space_->data.GetSize() : 0) + space_->GetNumberOfFeatureMappingVariables();
  if(space_->HasSlackVariable())
    size++;
  if(space_->HasAlphaVariable())
    size++;

  block.Resize(size);
  unsigned int base = 0;

  // exporting densities can be switched off in "export"
  if(write_density_)
  {
    for(unsigned int i = 0, n = space_->data.GetSize(); i < n; ++i)
    {
      DesignElement* de = &space_->data[i];
      std::stringstream ss;
      ss << "<element nr=\"" << de->elem->elemNum;
      ss << "\" type=\"" << DesignElement::type.ToString(de->GetType());
      if(de->GetType() == DesignElement::MULTIMATERIAL)
        ss << "\" index=\"" << de->multimaterial->index;
      ss << "\" design=\"";
      ss.precision(11);
      ss << de->GetDesign(DesignElement::PLAIN) << "\"";
      if(de->HasPhysicalDesign())
        ss << " physical=\"" << de->GetPhysicalDesign(Optimization::context) << "\"";
      ss << "/>";
      block[i] = ss.str();
    }
    base += space_->data.GetSize();
  }

  if(space_->HasSlackVariable())
  {
    std::stringstream ss;
    ss.precision(11);
    ss << "<slack nr=\"0\" type=\"slack\" design=\"" << space_->GetSlackVariable() << "\"/>";
    block[base] = ss.str();
    base += 1;
  }

  if(space_->HasAlphaVariable())
  {
    assert(space_->HasSlackVariable());
    std::stringstream ss;
    ss << "<slack nr=\"1\" type=\"alpha\" design=\"" << space_->GetAlphaVariable() << "\"/>";
    block[base] = ss.str();
    base += 1;
  }

  // add shape map design if we have it. Can be visualized by shape_map.py.
  // Also in the SMD case the above design is of interest!
  if(space_->GetMethod() == ErsatzMaterial::Method::SHAPE_MAP)
  {
    ShapeMapDesign* smd = dynamic_cast<ShapeMapDesign*>(space_);
    // skip the aux variables slack and alpha -> they are written above the info.xml
    // this are all variables form shape_param_ and not only the ones from opt_shape_param_ - important for visualization!
    for(unsigned int i = 0, n = space_->GetNumberOfFeatureMappingVariables(); i < n; i++)
    {
      ShapeParamElement* spe = dynamic_cast<ShapeParamElement*>(smd->GetFeaturedDesignElement(i));
      assert(spe != NULL);

      ShapeMapDesign::ShapeParam* shape = smd->GetShape(spe);

      std::stringstream ss;
      ss << "<shapeParamElement nr=\"" << spe->GetIndex();
      ss << "\" type=\"" << DesignElement::type.ToString(spe->GetType());
      if(spe->GetType() == DesignElement::NODE)
        ss << "\" dof=\"" << spe->dof.ToString(spe->dof_);
      ss << "\" shape=\"" << shape->idx; // legacy density.xml files don't have this attribute
      ss << "\" ref=\"" << shape->GetReferenceId(); // legacy density.xml files don't have this attribute
      ss.precision(11);
      ss << "\" design=\"" << spe->GetDesign(BaseDesignElement::PLAIN);
      ss << "\"/>";
      block[base + i] = ss.str();
    }
  }

  if(space_->GetMethod() == ErsatzMaterial::Method::FEATURE_MAPPING ||
     space_->GetMethod() == ErsatzMaterial::Method::FEATURE_MAPPING_PARAM_MAT)
  {
    FeaturedDesign* fd = dynamic_cast<FeaturedDesign*>(space_);

    for(int i = 0; i < fd->GetNumberOfFeatureMappingVariables(); i++)
    {
      FeatureVariable* var = dynamic_cast<FeatureVariable*>(fd->GetFeaturedDesignElement((unsigned int) i));

      std::stringstream ss;
      ss << "<shapeParamElement nr=\"" << var->GetIndex();
      ss << "\" type=\"" << DesignElement::type.ToString(var->GetType());
      if(var->GetType() == DesignElement::NODE) {
        ss << "\" dof=\"" << var->dof.ToString(var->dof_);
        ss << "\" tip=\"" << FeatureVariable::tip_enum.ToString(var->tip);
      }
      ss << "\" shape=\"" << var->feature; // legacy density.xml files don't have this attribute
      ss << "\" design=\"" << var->GetDesign(BaseDesignElement::PLAIN);
      ss << "\"/>";
      block[base + i] = ss.str();
    }
  }

  // spaghetti.py can visualize the noodles
  if(ErsatzMaterial::IsSpaghetti(space_->GetMethod()))
  {
#ifdef USE_EMBEDDED_PYTHON // currently only the python version

    SpaghettiDesign* sd = dynamic_cast<SpaghettiDesign*>(space_);
    // skip the aux variables slack and alpha -> they are written to the info.xml
    for(unsigned int i = 0, n = space_->GetNumberOfFeatureMappingVariables(); i < n; i++)
    {
      FeatureVariable* spe = dynamic_cast<FeatureVariable*>(sd->GetFeaturedDesignElement(i));
      assert(spe != NULL);

      LOG_DBG2(density) << "SAWC: " << spe->GetLabel() << " " << spe->ToString();

      std::stringstream ss;
      ss << "<shapeParamElement nr=\"" << spe->GetIndex();
      ss << "\" type=\"" << DesignElement::type.ToString(spe->GetType());
      if(spe->GetType() == DesignElement::NODE)
      {
        ss << "\" dof=\"" << spe->dof.ToString(spe->dof_);
        ss << "\" tip=\"" << FeatureVariable::tip_enum.ToString(spe->tip);
        const SpaghettiDesign::Noodle& n = sd->spaghetti[spe->feature];
        ss << "\" num=\"" << (unsigned int) (spe->GetIndex() - n.points[0][0].GetIndex()) / domain->GetDim();
      }
      if(spe->GetType() == DesignElement::RADIUS)
      {
        const SpaghettiDesign::Noodle& n = sd->spaghetti[spe->feature];
        ss << "\" con=\"" << spe->GetIndex() - n.r[0].GetIndex();
      }

      ss << "\" shape=\"" << spe->feature; // legacy density.xml files don't have this attribute
      ss << "\" design=\"" << spe->GetDesign(BaseDesignElement::PLAIN);
      ss << "\"/>";
      block[base + i] = ss.str();
    }
#else
  EXCEPTION("compile with USE_EMBEDDED_PYTHON");
#endif
  }

  // add shape map design if we have it.
  if(space_->GetMethod() == ErsatzMaterial::Method::SPLINE_BOX)
  {
    SplineBoxDesign* sbd = dynamic_cast<SplineBoxDesign*>(space_);
    // skip the aux variables slack and alpha -> they are written to the info.xml
    for(unsigned int i = 0, n = space_->GetNumberOfFeatureMappingVariables(); i < n; i++)
    {
      BaseDesignElement* de = sbd->GetFeaturedDesignElement(i);
      assert(de != NULL);

      std::stringstream ss;
      ss << "<splineBoxParam nr=\"" << de->GetIndex();
      ss << "\" type=\"" << DesignElement::type.ToString(de->GetType());
      ss << "\" design=\"" << de->GetDesign();
      ss << "\"/>";
      block[base + i] = ss.str();
    }
  }

  // do we need to write?
  if(!finally_only_) // here or on destructor
    data->ToFile(name_);
}
