#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/DesignElement.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "Domain/domain.hh"

using namespace CoupledField;
using std::string;

DensityFile::DensityFile(DesignSpace* designSpace,
                            PtrParamNode export_pn,
                            ParamNodeList& des,
                            ParamNodeList& tfs,
                            PtrParamNode regulize_pn)
{
  this->ersatzMaterial_ = designSpace;

  name_ = export_pn->Get("file")->As<string>();
  if(name_ == "[problem]") name_ = progOpts->GetSimName() + ".density.xml";
  data = Create(des, tfs, regulize_pn);
  all_iterations_ = export_pn->Get("save")->As<string>() == "all";
  finally_only_   = export_pn->Get("write")->As<string>() == "finally";
}


DensityFile::~DensityFile()
{
  // in CommitIteration or here?
  if(finally_only_)
    data->ToFile(name_);

  data.reset();
}

bool DensityFile::NeedLoadErsatzMaterial()
{
  return param->Has("loadErsatzMaterial") || progOpts->GetErsatzMaterialStr() != "";
}

DesignSpace* DensityFile::ReadErsatzMaterial(DesignSpace* ersatzMaterial)
{
  Grid* grid = domain->GetGrid();

  // perhaps Optimization has already called the SetEnums
  if (DesignElement::type.map.empty())
    DesignElement::SetEnums();
  if (Objective::type.map.empty())
    Optimization::SetEnums();

  // do we have a command line switch? we then use the filename and the last set
  bool cmd = progOpts->GetErsatzMaterialStr() != "";
  PtrParamNode pn; // some default auto pointer NULL stuff
  if(!cmd) pn = param->Get("loadErsatzMaterial");

  string file = cmd ? progOpts->GetErsatzMaterialStr() : pn->Get("file")->As<string>();

  // to be appended by the set name
  std::cout << "++ Load ersatz material file: '" << file << "'";

  // we read something like <loadErsatzMaterial region="piezo" file="piezo_density.xml" set="last"/>
  // Initialize our xerces dom parser to handle the external xml file
  Xerces* xerces = new Xerces(file);
  // set the global ParamNode tree pointer
  PtrParamNode xml = xerces->CreateParamNodeInstance();
  // release the xerces ressources, param is not affected
  delete xerces;
  // check this file
  if (xml->Count("set") == 0)
    throw Exception("There are no design sets in the ersatz material file");

  // check if we ignore the element numbers
  bool ignore_numbers = pn != NULL ? pn->Get("ignore_element_numbers")->As<bool>() : false;

  // the region list is either from the loadErsatzMaterial element or if not given
  // and the command line argument is used there must be only a single region
  StdVector<RegionIdType> regionIds;

  if(pn != NULL)
  {
    ParamNodeList region_list = pn->GetList("region");
    for (unsigned int i = 0; i < region_list.GetSize(); i++)
    {
      string reg = region_list[i]->Get("name")->As<string>();
      if (!grid->GetRegion().IsValid(reg))
        throw Exception("region given in loadErsatzMaterial is invalid");
      regionIds.Push_back(grid->GetRegion().Parse(reg));
    }

    if (ignore_numbers && region_list.GetSize() != 1)
      EXCEPTION("'ignore_element_numbers' in 'loadErsatzMaterial' only allowed for a single region");
  }
  else
  {
    if(grid->GetNumRegions() != 1)
      throw Exception("can only load single region meshes with 'ersatz' command line option");
    regionIds.Push_back(0);
  }

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

  if (!ersatzMaterial)
  { // only if the design space does not already exist (created by optimization)
    // the header is like
    // <header>
    //   <design name="density" initial="1.0"/>
    //   ...
    //   <transferFunction type="simp" application="mech" design="density" param="1.0"/>
    //   ..
    //   <regularization type="filter">
    //    <filter contribution="linear" neighborhood="maxEdge" type="density" value="1.7"/>
    //   </regularization>
    // </header>

    // the design set consists of entries like
    // <element nr="401" type="density" design="0.886466" physical="0.800454" />
    // only the combination nr and type is unique. E.g. in piezo we might have types density and polarization
    ParamNodeList des = xml->Get("header")->GetList("design");
    ParamNodeList tfs = xml->Get("header")->GetList("transferFunction");
    PtrParamNode  reg = xml->Get("header/regularization/filter", ParamNode::PASS);
    ParamNodeList res(0); // empty

    // create the design space -> data has initial values!
    ersatzMaterial = new DesignSpace(regionIds, des, tfs, res, ErsatzMaterial::SIMP_METHOD);
    ersatzMaterial->PostInit(0, 0); // no objectives, no constraints
    // is cheap - for density filtering
    DesignStructure filter(ersatzMaterial, regionIds);
    if(reg) filter.SetFilters(reg, info->Get("ersatzMaterial"));

    ersatzMaterial->ToInfo(info->Get("ersatzMaterial")->Get(ParamNode::HEADER));
  }

  // read the set and replace the initial values
  ParamNodeList elems = set->GetList("element");

  // check the the dimensions! the number of design variables comes from the regions and designs
  if (ersatzMaterial->data.GetSize() != elems.GetSize())
  {
    if(grid->GetNumElems() == elems.GetSize())
      std::cout << "\033[01;31m" << "Warning:" << "\033[0m" << std::endl
                << "the number of elements in the region you are trying to read the densities into is not equal to"
                << " the number of elements in the density-file but matches the number of all elements!"
                << " We ignore this... I hope you know what you are doing!"
                << std::endl;
    else
      EXCEPTION("ErsatzMaterialFile '" << file << "' has " << elems.GetSize()
          << " entries, the mesh has "<< ersatzMaterial->data.GetSize() << " design elements");
  }

  for (unsigned int e = 0; e < elems.GetSize(); e++)
  {
    unsigned int nr = elems[e]->Get("nr")->As<Integer>();
    DesignElement::Type dt = (DesignElement::Type) DesignElement::type.Parse(
        elems[e]->Get("type")->As<string>());
    double val = elems[e]->Get("design")->As<Double>();

    // replace the value of the DesignElement
    DesignElement* de = ignore_numbers ? &(ersatzMaterial->data[e])
        : ersatzMaterial->Find(nr, dt, false);
    // it should be possible to specify less regions then specified during optimization and saving of results
    // if the element can not be found (e.g. lying in a not specified region) it is not set
    if (de != NULL)
    {
      // and the region can be set in optimization, thus exist, but not specified here, we should not set as well
      if (regionIds.Find(de->elem->regionId) >= 0)
      {
        de->SetDesign(val);
      }
    }
  }

  return ersatzMaterial;

}

PtrParamNode DensityFile::Create(ParamNodeList& des, ParamNodeList& tfs, PtrParamNode regularize)
{
   PtrParamNode in = PtrParamNode(new ParamNode(ParamNode::INSERT));
   in->SetName("cfsErsatzMaterial");

   // write header
   PtrParamNode in_ = in->Get("header");

   for(unsigned int i = 0; i < des.GetSize(); i++)
     in_->Get("dummy", ParamNode::APPEND)->SetValue(des[i], true); // name is overwritten

   for(unsigned int i = 0; i < tfs.GetSize(); i++)
     in_->Get("dummy", ParamNode::APPEND)->SetValue(tfs[i], true);

   if(regularize)
     in_->Get("dummy")->SetValue(regularize, true);

   return in;
}


void DensityFile::SetCurrent(int current_iteration)
{
  PtrParamNode in = data->Get("set", all_iterations_ ? ParamNode::APPEND : ParamNode::INSERT);
  // add the entry, note that the iteration counter was incremented in base implementation
  in->Get("id")->SetValue(current_iteration);

  // the loop can be quite long. When it is filled already, ParamNode::Get() is cached
  // but when it is filled it would search. Do it faster!

  // even for !exportDesignAllIterations in is empty in the first commit

  ParamNodeList& list = in->GetChildren();

  const unsigned int bias = 1; // there one "id" element as first element

  assert(list.GetSize() == bias || list.GetSize() == ersatzMaterial_->data.GetSize() + bias);
  assert(list[0]->GetName() == "id");

  bool append = list.GetSize() == bias;
  assert(!(all_iterations_ && ! append));

  if(append)
    list.Resize(ersatzMaterial_->data.GetSize() + 1);

  for(unsigned int i = 0, s = ersatzMaterial_->data.GetSize(); i < s; ++i)
  {
    DesignElement* de = &ersatzMaterial_->data[i];
    PtrParamNode el = append ? in->SetNewChild("element", i + bias) : list[i + bias];
    el->Get("nr")->SetValue(de->elem->elemNum);
    el->Get("type")->SetValue(DesignElement::type.ToString(de->GetType()));
    el->Get("design")->SetValue(de->GetDesign(DesignElement::PLAIN), 11);
    el->Get("physical")->SetValue(de->GetPhysicalDesign(), 11);
  }

  // do we need to write?
  if(!finally_only_) // here or on destructor
    data->ToFile(name_);
}
