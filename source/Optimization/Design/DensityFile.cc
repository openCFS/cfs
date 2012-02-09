#include <boost/algorithm/string.hpp>
#include <iostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/programOptions.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Utils/Point.hh"
#include "Utils/StdVector.hh"

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
  data = Create(des, tfs, regulize_pn, designSpace->DoNonDesignVicinity());
  all_iterations_ = export_pn->Get("save")->As<string>() == "all";
  finally_only_   = export_pn->Get("write")->As<string>() == "finally";
  // append .bz2 if compress=true and not already file ends with it
  if(export_pn->Get("compress")->As<bool>() && !boost::algorithm::ends_with(name_, ".bz2")){
    name_ = name_ + ".bz2";
  }
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
  std::cout << "++ Load ersatz material file: '" << file << "'" << std::flush;

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

    // create the design space -> data has initial values!
    ersatzMaterial = new DesignSpace(regionIds, xml->Get("header"), ErsatzMaterial::SIMP_METHOD);
    ersatzMaterial->PostInit(0, 0); // no objectives, no constraints
    // is cheap - for density filtering
    DesignStructure filter(ersatzMaterial, ersatzMaterial->GetRegionIds());
    PtrParamNode  reg = xml->Get("header/regularization/filter", ParamNode::PASS);
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

PtrParamNode DensityFile::Create(ParamNodeList& des, ParamNodeList& tfs, PtrParamNode regularize, bool non_desing_vicinity)
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

   in_->Get("designSpace/non_design_vicinity")->SetValue(non_desing_vicinity);

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

  // we use the fast (dirty) bulk block to be (measurable) faster
  StdVector<std::string>& block = in->GetFastBulkBlock();
  block.Resize(ersatzMaterial_->data.GetSize());

  for(unsigned int i = 0, n = ersatzMaterial_->data.GetSize(); i < n; ++i)
  {
    DesignElement* de = &ersatzMaterial_->data[i];
    std::stringstream ss;
    ss << "<element nr=\"" << de->elem->elemNum
       << "\" type=\"" << DesignElement::type.ToString(de->GetType())
       << "\" design=\"";
    ss.precision(11);
    ss << de->GetDesign(DesignElement::PLAIN) << "\"";
    if(de->HasPhysicalDesign())
      ss << " physical=\"" << de->GetPhysicalDesign() << "\"";
    ss << "/>";
    block[i] = ss.str();
  }

  // do we need to write?
  if(!finally_only_) // here or on destructor
    data->ToFile(name_);
}
