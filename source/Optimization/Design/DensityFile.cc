#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/DesignElement.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Domain/domain.hh"

using namespace CoupledField;
using std::string;

DECLARE_LOG(density)
DEFINE_LOG(density, "density")

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
  ParamNodeList elems = set->GetList("element");
  const unsigned int elsize(elems.GetSize());
  bool force_region = pn != NULL && pn->Has("force_region");

  if(!ersatzMaterial)
  {
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
      for(unsigned int e = 0; e < elsize; ++e)
      {
        unsigned int nr = elems[e]->Get("nr")->As<unsigned int>();
        if(!regionIds.Contains(grid->GetElem(nr)->regionId))
          regionIds.Push_back(grid->GetElem(nr)->regionId);
      }
    }

    ParamNodeList des = xml->Get("header")->GetList("design");
    ParamNodeList tfs = xml->Get("header")->GetList("transferFunction");
    PtrParamNode  reg = xml->Get("header/regularization/filter", ParamNode::PASS);
    ParamNodeList res(0); // empty

    // create the design space -> data has initial values!
    ersatzMaterial = new DesignSpace(regionIds, des, tfs, res, ErsatzMaterial::SIMP_METHOD);
    ersatzMaterial->PostInit(0, 0); // no objectives, no constraints
    // is cheap - for density filtering
    DesignStructure filter(ersatzMaterial, ersatzMaterial->GetRegionIds());
    if(reg) filter.SetFilters(reg, info->Get("ersatzMaterial"));

    ersatzMaterial->ToInfo(info->Get("ersatzMaterial")->Get(ParamNode::HEADER));
  }


  // check the the dimensions! the number of design variables comes from the regions and designs
  if (ersatzMaterial->data.GetSize() != elsize)
  {
    string msg = "the number of elements in the density file does not match the number of elements of the region!\n"\
                 "         check the results carefully!";
    info->Get("ersatzMaterial")->Get(ParamNode::WARNING)->SetValue(msg);
  }

  for (unsigned int e = 0; e < elsize; ++e)
  {
    // the design set consists of entries like
    // <element nr="401" type="density" design="0.886466" physical="0.800454" />
    // only the combination nr and type is unique. E.g. in piezo we might have types density and polarization

    unsigned int nr = elems[e]->Get("nr")->As<unsigned int>();
    DesignElement::Type dt = (DesignElement::Type) DesignElement::type.Parse(elems[e]->Get("type")->As<string>());
    double val = elems[e]->Get("design")->As<double>();

    // replace the value of the DesignElement
    // we call Find(..,..,false) for meshes with two regions (e. g. cube and void)
    // where we want to ignore the "void"-region completely
    DesignElement* de = force_region ? &(ersatzMaterial->data[e]) : ersatzMaterial->Find(nr, dt, false);

    // this is also for the void-region! mainly for computing high resolution inv hom problems
    if(de != NULL) // && regionIds.Find(de->elem->regionId) >= 0)
      de->SetDesign(val);
  }

  return ersatzMaterial;

}

PtrParamNode DensityFile::Create(ParamNodeList& des, ParamNodeList& tfs, PtrParamNode regularize)
{
   PtrParamNode in = PtrParamNode(new ParamNode(ParamNode::INSERT));
   in->SetName("cfsErsatzMaterial");

   // write header
   PtrParamNode in_ = in->Get("header");

   LOG_TRACE(density) << "Create: regular=" << this->ersatzMaterial_->IsRegular();

   if(this->ersatzMaterial_->IsRegular())
   {
     DesignElement& de = ersatzMaterial_->data[0];
     Matrix<double> coords;
     domain->GetGrid()->GetElemNodesCoord(coords, de.elem->connect, false );
     StdVector<double> edges;
     de.elem->ptElem->GetEdgeLength(coords, edges);

     Point dim;
     domain->GetGrid()->CalcVolumeSpannedByNamedNodes(&dim);

     LOG_TRACE(density) << " dim=" << dim.ToString() << " edges=" << edges.ToString();
     PtrParamNode mesh = in_->Get("mesh");
     mesh->Get("x")->SetValue(dim.data[0] / edges[0]);
     mesh->Get("y")->SetValue(dim.data[1] / edges[1]);
     mesh->Get("z")->SetValue(domain->GetGrid()->GetDim() == 3 ? dim.data[2] / edges[2]: 1);
   }

   for(unsigned int i = 0, s = des.GetSize(); i < s; ++i)
     in_->Get("dummy", ParamNode::APPEND)->SetValue(des[i], true); // name is overwritten

   for(unsigned int i = 0, s = tfs.GetSize(); i < s; ++i)
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
    if(de->HasPhysicalDesign()){
      el->Get("physical")->SetValue(de->GetPhysicalDesign(), 11);
    }
  }

  // do we need to write?
  if(!finally_only_) // here or on destructor
    data->ToFile(name_);
}
